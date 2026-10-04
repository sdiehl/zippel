//! Declarative YAML and Mathematica family interchange.
//!
//! Complete quadratic families support polynomial kinematics and rational
//! momentum coefficients. Unsupported definitions are
//! errors, never evaluated as Mathematica programs or silently approximated.

#![allow(clippy::too_many_lines)]

use crate::{family::PolynomialFamily, ibp::System, Plan};
use num_rational::BigRational;
use num_traits::Zero;
use polycore::{dense, Monomial, Order, Poly, Ring};
use serde_yaml_ng::Value;

use std::collections::BTreeMap;
use std::fmt::{self, Write};
use zippel_lift::Fraction;

type Q = BigRational;
type Result<T> = std::result::Result<T, FormatError>;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FormatError(pub String);
impl fmt::Display for FormatError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}
impl std::error::Error for FormatError {}
fn error(s: impl Into<String>) -> FormatError {
    FormatError(s.into())
}

#[derive(Clone, Debug)]
pub struct ImportedFamily {
    pub name: String,
    pub family: PolynomialFamily,
    /// Original propagator = sign * the normalized internal propagator.
    pub signs: Vec<i64>,
}

fn seq(v: &Value) -> Result<&[Value]> {
    v.as_sequence()
        .map(Vec::as_slice)
        .ok_or_else(|| error("expected YAML sequence"))
}
fn atom(v: &Value) -> Result<String> {
    match v {
        Value::String(s) => Ok(s.clone()),
        Value::Number(n) => Ok(n.to_string()),
        _ => Err(error("expected expression string or number")),
    }
}
fn names(v: &Value) -> Result<Vec<String>> {
    seq(v)?.iter().map(atom).collect()
}
fn parse(r: &Ring, s: &str) -> Result<Poly<Q>> {
    r.parse(&explicit_products(s))
        .map_err(|e| error(e.to_string()))
}
fn ring(names: Vec<String>) -> Result<Ring> {
    Ring::try_new(names, Order::Lex).map_err(|e| error(e.to_string()))
}
// Mathematica also writes products as `p q` or `2 (p+q)`.
fn explicit_products(s: &str) -> String {
    let chars: Vec<_> = s.chars().collect();
    let mut out = String::new();
    for (i, &c) in chars.iter().enumerate() {
        if c.is_whitespace() {
            let left = out.chars().last();
            let right = chars[i + 1..].iter().find(|c| !c.is_whitespace());
            if left.is_some_and(|c| c.is_ascii_alphanumeric() || c == '_' || c == ')')
                && right.is_some_and(|c| c.is_ascii_alphanumeric() || *c == '_' || *c == '(')
            {
                out.push('*');
            }
        } else {
            out.push(c);
        }
    }
    out
}

fn keys(v: &Value, allowed: &[&str]) -> Result<()> {
    for k in v
        .as_mapping()
        .ok_or_else(|| error("expected mapping"))?
        .keys()
    {
        let k = atom(k)?;
        if !allowed.contains(&k.as_str()) {
            return Err(error(format!("unsupported option {k}")));
        }
    }
    Ok(())
}

/// Read one named family from the YAML format’s `integralfamilies.yaml` and `kinematics.yaml`.
///
/// Supports sector unions, cuts, momentum conservation, polynomial scalar-product
/// rules involving momentum sums, and fixing one invariant to unity.
///
/// # Errors
/// Rejects malformed, unsupported, incomplete, or dependent family definitions.
pub fn read_family_yaml(families: &str, kinematics: &str, name: &str) -> Result<ImportedFamily> {
    let defs: Value = serde_yaml_ng::from_str(families).map_err(|e| error(e.to_string()))?;
    let kin: Value = serde_yaml_ng::from_str(kinematics).map_err(|e| error(e.to_string()))?;
    let k = &kin["kinematics"];
    keys(
        k,
        &[
            "incoming_momenta",
            "outgoing_momenta",
            "momentum_conservation",
            "kinematic_invariants",
            "scalarproduct_rules",
            "symbol_to_replace_by_one",
        ],
    )?;
    let f = seq(&defs["integralfamilies"])?
        .iter()
        .find(|f| f["name"].as_str() == Some(name))
        .ok_or_else(|| error("family not found"))?;
    keys(
        f,
        &[
            "name",
            "loop_momenta",
            "top_level_sectors",
            "propagators",
            "cut_propagators",
            "zero_sectors",
            "permutation_option",
        ],
    )?;
    let loops = names(&f["loop_momenta"])?;
    let mut external = Vec::new();
    for key in ["incoming_momenta", "outgoing_momenta"] {
        if !k[key].is_null() {
            external.extend(names(&k[key])?);
        }
    }
    let conservation = if k["momentum_conservation"].is_null() {
        None
    } else {
        let c = seq(&k["momentum_conservation"])?;
        if c.len() != 2 {
            return Err(error("momentum_conservation needs [momentum, expression]"));
        }
        Some((atom(&c[0])?, atom(&c[1])?))
    };
    if let Some((old, _)) = &conservation {
        if !external.contains(old) {
            return Err(error("eliminated momentum is not external"));
        }
        external.retain(|v| v != old);
    }
    let substitute = |s: String| {
        if let Some((old, new)) = &conservation {
            replace_ident(&s, old, &format!("({new})"))
        } else {
            s
        }
    };
    let vars = std::iter::once("d".to_owned())
        .chain(
            seq(&k["kinematic_invariants"])?
                .iter()
                .map(|v| {
                    seq(v)
                        .and_then(|v| v.first().ok_or_else(|| error("empty invariant")))
                        .and_then(atom)
                })
                .collect::<Result<Vec<_>>>()?,
        )
        .collect::<Vec<_>>();
    let vr = ring(vars.clone())?;
    let mr = ring(loops.iter().chain(&external).cloned().collect())?;
    let er = ring(external.clone())?;
    let rules: Vec<_> = seq(&k["scalarproduct_rules"])?
        .iter()
        .map(|v| {
            let v = seq(v)?;
            if v.len() != 2 {
                return Err(error("scalar rule requires two entries"));
            }
            let pair = seq(&v[0])?;
            if pair.len() != 2 {
                return Err(error("scalar rule requires a momentum pair"));
            }
            Ok((
                format!(
                    "({})*({})",
                    substitute(atom(&pair[0])?),
                    substitute(atom(&pair[1])?)
                ),
                atom(&v[1])?,
            ))
        })
        .collect::<Result<_>>()?;
    let products = scalar_products(&er, &vr, &rules)?;
    let full = ring(mr.names.iter().chain(&vars).cloned().collect())?;
    let propagators = seq(&f["propagators"])?
        .iter()
        .map(|v| {
            let v = seq(v)?;
            if v.len() != 2 {
                return Err(error(
                    "propagator requires [momentum or quadratic expression, mass_squared]",
                ));
            }
            let g = parse(&full, &substitute(atom(&v[0])?))?;
            let degree = g
                .terms
                .iter()
                .map(|(m, _)| m.exps()[..mr.nvars()].iter().sum::<u32>())
                .max()
                .unwrap_or(0);
            let g = if degree == 1 { g.pow(2) } else { g };
            Ok(&g - &parse(&full, &atom(&v[1])?)?)
        })
        .collect::<Result<Vec<_>>>()?;
    let sectors = seq(&f["top_level_sectors"])?
        .iter()
        .map(sector)
        .collect::<Result<Vec<_>>>()?;
    let zero_sectors = if f["zero_sectors"].is_null() {
        vec![]
    } else {
        seq(&f["zero_sectors"])?
            .iter()
            .map(sector)
            .collect::<Result<_>>()?
    };
    let cuts = if f["cut_propagators"].is_null() {
        vec![]
    } else {
        seq(&f["cut_propagators"])?
            .iter()
            .map(|v| {
                v.as_u64()
                    .and_then(|i| i.checked_sub(1))
                    .and_then(|i| usize::try_from(i).ok())
                    .ok_or_else(|| error("cut indices start at one"))
            })
            .collect::<Result<_>>()?
    };
    let mut family = PolynomialFamily {
        vars,
        loops: loops.len(),
        propagators,
        products,
        top_sectors: sectors,
        zero_sectors,
        cuts,
    };
    if !k["symbol_to_replace_by_one"].is_null() {
        let fixed = atom(&k["symbol_to_replace_by_one"])?;
        if fixed == "d" || !family.vars.contains(&fixed) {
            return Err(error("invalid fixed invariant"));
        }
        family = family.fix(&fixed, 1);
    }
    validate(&family)?;
    let signs = vec![1; family.propagators.len()];
    Ok(ImportedFamily {
        name: name.into(),
        family,
        signs,
    })
}

fn sector(v: &Value) -> Result<u32> {
    if let Some(n) = v.as_u64() {
        return u32::try_from(n).map_err(|_| error("sector exceeds 32 bits"));
    }
    let s = atom(v)?;
    if let Some(bits) = s.strip_prefix('b') {
        // In this notation the leftmost bit denotes the first propagator.
        return u32::from_str_radix(&bits.chars().rev().collect::<String>(), 2)
            .map_err(|_| error("invalid binary sector"));
    }
    Err(error("expected integer or binary sector"))
}

fn scalar_products(er: &Ring, vr: &Ring, rules: &[(String, String)]) -> Result<Vec<Vec<Poly<Q>>>> {
    let n = er.nvars();
    let pairs: Vec<_> = (0..n).flat_map(|i| (i..n).map(move |j| (i, j))).collect();
    let mut a = Vec::new();
    let mut rhs = Vec::new();
    for (lhs, rhs_expr) in rules {
        let mut row = vec![Q::zero(); pairs.len()];
        for (m, c) in parse(er, lhs)?.terms {
            let indices: Vec<_> = m
                .exps()
                .iter()
                .enumerate()
                .flat_map(|(i, &e)| std::iter::repeat_n(i, e as usize))
                .collect();
            if indices.len() != 2 {
                return Err(error("scalar-product rule must be quadratic in momenta"));
            }
            let pos = pairs
                .iter()
                .position(|&(i, j)| i == indices[0] && j == indices[1])
                .ok_or_else(|| error("unknown product"))?;
            row[pos] += c;
        }
        a.push(row);
        rhs.push(parse(vr, rhs_expr)?);
    }
    let mut support = std::collections::BTreeSet::new();
    for g in &rhs {
        for (m, _) in &g.terms {
            support.insert(m.exps().to_vec());
        }
    }
    support.insert(vec![0; vr.nvars()]);
    let mut values = vec![Vec::new(); pairs.len()];
    for exps in support {
        let b: Vec<_> = rhs
            .iter()
            .map(|g| {
                g.terms
                    .iter()
                    .find(|(m, _)| m.exps() == exps)
                    .map_or_else(Q::zero, |(_, c)| c.clone())
            })
            .collect();
        let x = dense::solve(&a, &b).map_err(|e| error(format!("scalar-product rules: {e}")))?;
        if x.len() != pairs.len() {
            return Err(error("incomplete scalar-product rules"));
        }
        for (v, c) in values.iter_mut().zip(x) {
            v.push((Monomial::new(exps.clone()), c));
        }
    }
    let zero = Poly::constant(Q::zero(), vr.nvars(), Order::Lex);
    let mut products = vec![vec![zero; n]; n];
    for ((i, j), terms) in pairs.into_iter().zip(values) {
        let g = Poly::new(terms, vr.nvars(), Order::Lex);
        products[i][j] = g.clone();
        products[j][i] = g;
    }
    Ok(products)
}

fn validate(f: &PolynomialFamily) -> Result<()> {
    if !f.valid() {
        return Err(error("invalid polynomial family dimensions or sectors"));
    }
    if !f.independent() {
        return Err(error("dependent propagators"));
    }
    Ok(())
}

fn replace_ident(s: &str, old: &str, new: &str) -> String {
    let mut out = String::new();
    let mut start = 0;
    for (i, c) in s.char_indices() {
        if !(c.is_ascii_alphanumeric() || c == '_') {
            let word = &s[start..i];
            out.push_str(if word == old { new } else { word });
            out.push(c);
            start = i + c.len_utf8();
        }
    }
    let word = &s[start..];
    out.push_str(if word == old { new } else { word });
    out
}

/// Split a Mathematica list or statement sequence at top-level delimiters.
fn split(s: &str, delimiter: char) -> Result<Vec<&str>> {
    let mut out = Vec::new();
    let mut depth = 0i32;
    let mut start = 0;
    for (i, c) in s.char_indices() {
        match c {
            '{' | '(' | '[' => depth += 1,
            '}' | ')' | ']' => {
                depth -= 1;
                if depth < 0 {
                    return Err(error("unbalanced expression"));
                }
            }
            _ => {}
        }
        if c == delimiter && depth == 0 {
            out.push(s[start..i].trim());
            start = i + 1;
        }
    }
    if depth != 0 {
        return Err(error("unbalanced expression"));
    }
    out.push(s[start..].trim());
    Ok(out)
}
fn list(s: &str) -> Result<Vec<&str>> {
    split(
        s.trim()
            .strip_prefix('{')
            .and_then(|s| s.strip_suffix('}'))
            .ok_or_else(|| error("expected Mathematica list"))?,
        ',',
    )
}

// Mathematica momentum components, e.g. p[1], are algebraic names, not calls.
fn indexed_names(s: &str) -> String {
    let chars: Vec<_> = s.chars().collect();
    let mut out = String::new();
    let mut i = 0;
    while i < chars.len() {
        if chars[i] == '[' && i > 0 && (chars[i - 1].is_ascii_alphanumeric() || chars[i - 1] == '_')
        {
            let mut j = i + 1;
            while j < chars.len() && chars[j].is_whitespace() {
                j += 1;
            }
            let start = j;
            while j < chars.len() && chars[j].is_ascii_digit() {
                j += 1;
            }
            let end = j;
            while j < chars.len() && chars[j].is_whitespace() {
                j += 1;
            }
            if end > start && chars.get(j) == Some(&']') {
                out.push('_');
                out.extend(chars[start..end].iter());
                i = j + 1;
                continue;
            }
        }
        out.push(chars[i]);
        i += 1;
    }
    out
}

/// Read Mathematica declarations: Internal, External, Propagators, Replacements.
///
/// `invariants` gives the independent scalar names (excluding d); `lines` excludes
/// trailing ISPs. Indexed momentum names, mapped squares, and `Thread` rules are
/// supported. General Mathematica evaluation and executable statements are rejected.
///
/// # Errors
/// Rejects executable input, unsupported propagators, and incomplete kinematics.
pub fn read_family_mathematica(
    src: &str,
    name: &str,
    invariants: &[&str],
    lines: usize,
) -> Result<ImportedFamily> {
    let mut clean = String::new();
    let mut rest = src;
    while let Some((before, after)) = rest.split_once("(*") {
        clean.push_str(before);
        let (_, tail) = after
            .split_once("*)")
            .ok_or_else(|| error("unterminated comment"))?;
        rest = tail;
    }
    clean.push_str(rest);
    let clean = indexed_names(&clean);
    let mut assignments = BTreeMap::new();
    for statement in split(&clean, ';')?.into_iter().filter(|s| !s.is_empty()) {
        let (key, value) = statement
            .split_once('=')
            .ok_or_else(|| error("only literal Mathematica assignments are supported"))?;
        let key = key.trim();
        if !["Internal", "External", "Propagators", "Replacements"].contains(&key)
            || assignments.insert(key, value.trim()).is_some()
        {
            return Err(error(format!("unsupported or duplicate assignment {key}")));
        }
    }
    let get = |key| {
        assignments
            .get(key)
            .copied()
            .ok_or_else(|| error(format!("missing {key}")))
    };
    let loops: Vec<_> = list(get("Internal")?)?
        .into_iter()
        .map(str::to_owned)
        .collect();
    let external: Vec<_> = list(get("External")?)?
        .into_iter()
        .filter(|s| !s.is_empty())
        .map(str::to_owned)
        .collect();
    let vars: Vec<_> = std::iter::once("d")
        .chain(invariants.iter().copied())
        .map(str::to_owned)
        .collect();
    let vr = ring(vars.clone())?;
    let er = ring(external.clone())?;
    let mr = ring(loops.iter().chain(&external).cloned().collect())?;
    let full = ring(mr.names.iter().chain(&vars).cloned().collect())?;
    let replacement = get("Replacements")?.trim();
    let rules = if let Some(inner) = replacement
        .strip_prefix("Thread[")
        .and_then(|s| s.strip_suffix(']'))
    {
        let (lhs, rhs) = inner
            .split_once("->")
            .ok_or_else(|| error("Thread requires a replacement arrow"))?;
        let lhs = list(lhs)?;
        let rhs = list(rhs)?;
        if lhs.len() != rhs.len() {
            return Err(error("Thread list lengths differ"));
        }
        lhs.into_iter()
            .zip(rhs)
            .map(|(a, b)| (a.to_owned(), b.to_owned()))
            .collect::<Vec<_>>()
    } else {
        list(replacement)?
            .into_iter()
            .filter(|s| !s.is_empty())
            .map(|s| {
                let (a, b) = s
                    .split_once("->")
                    .ok_or_else(|| error("expected replacement arrow"))?;
                Ok((a.trim().into(), b.trim().into()))
            })
            .collect::<Result<Vec<_>>>()?
    };
    let products = scalar_products(&er, &vr, &rules)?;
    let props = get("Propagators")?;
    let propagators = if let Some((function, values)) = props.split_once("/@") {
        let function = function
            .chars()
            .filter(|c| !c.is_whitespace())
            .collect::<String>();
        let sign = match function.as_str() {
            "#^2&" | "(#^2)&" | "(#)^2&" | "(-#)^2&" => 1,
            "-#^2&" | "(-#^2)&" => -1,
            _ => return Err(error("only mapped momentum squares are supported")),
        };
        list(values)?
            .iter()
            .map(|s| Ok(parse(&full, s)?.pow(2).scale(&Q::from_integer(sign.into()))))
            .collect::<Result<Vec<_>>>()?
    } else {
        list(props)?
            .iter()
            .map(|s| parse(&full, s))
            .collect::<Result<Vec<_>>>()?
    };
    if lines >= 32 || lines == 0 || lines > propagators.len() {
        return Err(error("invalid line count"));
    }
    let signs = vec![1; propagators.len()];
    let family = PolynomialFamily {
        vars,
        loops: loops.len(),
        propagators,
        products,
        top_sectors: vec![(1 << lines) - 1],
        zero_sectors: vec![],
        cuts: vec![],
    };
    validate(&family)?;
    Ok(ImportedFamily {
        name: name.into(),
        family,
        signs,
    })
}

impl ImportedFamily {
    /// Mathematica replacement rules, using the YAML format’s family head or the indexed format’s `G[id,{...}]`.
    ///
    /// # Errors
    /// Rejects incompatible coefficient dimensions or invalid variable names.
    pub fn rules(
        &self,
        system: &System,
        plan: &Plan,
        coefficients: &[Fraction],
        family_id: Option<u32>,
    ) -> Result<String> {
        if coefficients.len() != plan.targets.len() * plan.masters.len() {
            return Err(error("coefficient dimensions"));
        }
        let r = ring(system.vars.clone())?;
        let integral = |j: usize| {
            let indices: Vec<_> = system.integrals[j]
                .iter()
                .map(ToString::to_string)
                .collect();
            family_id.map_or_else(
                || format!("{}[{}]", self.name, indices.join(",")),
                |id| format!("G[{id},{{{}}}]", indices.join(",")),
            )
        };
        let sign = |j: usize| {
            system.integrals[j]
                .iter()
                .zip(&self.signs)
                .fold(1, |v, (&a, &s)| if a % 2 != 0 { v * s } else { v })
        };
        let mut out = String::from("{\n");
        for (i, &t) in plan.targets.iter().enumerate() {
            let mut terms = Vec::new();
            for (j, &m) in plan.masters.iter().enumerate() {
                let c = &coefficients[i * plan.masters.len() + j];
                if c.num.is_zero() {
                    continue;
                }
                let num = c.num.scale(&Q::from_integer((sign(t) * sign(m)).into()));
                terms.push(format!(
                    "({})/({})*{}",
                    r.show(&num),
                    r.show(&c.den),
                    integral(m)
                ));
            }
            writeln!(
                out,
                "{} -> {}{}",
                integral(t),
                if terms.is_empty() {
                    "0".into()
                } else {
                    terms.join(" + ")
                },
                if i + 1 == plan.targets.len() { "" } else { "," }
            )
            .unwrap();
        }
        out.push('}');
        Ok(out)
    }
}

impl ImportedFamily {
    /// Mathematica `.tables`: reduction entries plus an explicit integral-ID dictionary.
    ///
    /// # Errors
    /// Rejects incompatible coefficient dimensions or invalid variable names.
    pub fn reduction_tables(
        &self,
        system: &System,
        plan: &Plan,
        coefficients: &[Fraction],
        family_id: u32,
    ) -> Result<String> {
        if coefficients.len() != plan.targets.len() * plan.masters.len() {
            return Err(error("coefficient dimensions"));
        }
        let ring = ring(system.vars.clone())?;
        let sign = |j: usize| {
            system.integrals[j]
                .iter()
                .zip(&self.signs)
                .fold(1, |v, (&a, &s)| if a % 2 != 0 { v * s } else { v })
        };
        let mut ids: Vec<_> = plan.targets.iter().chain(&plan.masters).copied().collect();
        ids.sort_unstable();
        ids.dedup();
        let id = |j| ids.binary_search(&j).unwrap() + 1;
        let mut entries = Vec::new();
        for (i, &t) in plan.targets.iter().enumerate() {
            let terms: Vec<_> = plan
                .masters
                .iter()
                .enumerate()
                .filter_map(|(j, &m)| {
                    let c = &coefficients[i * plan.masters.len() + j];
                    if c.num.is_zero() {
                        return None;
                    }
                    let num = c.num.scale(&Q::from_integer((sign(t) * sign(m)).into()));
                    Some(format!(
                        "{{{},\"({})/({})\"}}",
                        id(m),
                        ring.show(&num),
                        ring.show(&c.den)
                    ))
                })
                .collect();
            entries.push(format!("{{{},{{{}}}}}", id(t), terms.join(",")));
        }
        for &m in &plan.masters {
            if !plan.targets.contains(&m) {
                entries.push(format!("{{{},{{{{{},\"1\"}}}}}}", id(m), id(m)));
            }
        }
        let dictionary: Vec<_> = ids
            .iter()
            .enumerate()
            .map(|(i, &j)| {
                let indices: Vec<_> = system.integrals[j]
                    .iter()
                    .map(ToString::to_string)
                    .collect();
                format!("{{{},{{{family_id},{{{}}}}}}}", i + 1, indices.join(","))
            })
            .collect();
        Ok(format!(
            "{{\n{{{}}},\n{{{}}}\n}}\n",
            entries.join(",\n"),
            dictionary.join(",\n")
        ))
    }
}

/// A Mathematica table with its integral dictionary retained for basis-independent checks.
#[derive(Clone, Debug)]
pub struct ReductionTable {
    pub integrals: BTreeMap<String, (u32, Vec<i32>)>,
    pub reductions: BTreeMap<String, Vec<(String, String)>>,
}

/// Read the two-list `.tables` interchange format without executing expressions.
///
/// # Errors
/// Rejects malformed lists, repeated IDs, and missing dictionary references.
pub fn read_reduction_tables(src: &str) -> Result<ReductionTable> {
    fn pair(s: &str) -> Result<(&str, &str)> {
        let v = list(s)?;
        if v.len() != 2 {
            return Err(error("expected pair"));
        }
        Ok((v[0], v[1]))
    }
    let id = |s: &str| {
        if !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit()) {
            Ok(s.to_owned())
        } else {
            Err(error("invalid integral ID"))
        }
    };
    let (entries, dictionary) = pair(src)?;
    let mut integrals = BTreeMap::new();
    for entry in list(dictionary)?.into_iter().filter(|s| !s.is_empty()) {
        let (i, value) = pair(entry)?;
        let (family, indices) = pair(value)?;
        let family = family
            .parse::<u32>()
            .map_err(|_| error("invalid family ID"))?;
        let indices = list(indices)?
            .iter()
            .map(|v| {
                v.parse::<i32>()
                    .map_err(|_| error("invalid integral index"))
            })
            .collect::<Result<Vec<_>>>()?;
        if integrals.insert(id(i)?, (family, indices)).is_some() {
            return Err(error("duplicate integral ID"));
        }
    }
    let mut reductions = BTreeMap::new();
    for entry in list(entries)?.into_iter().filter(|s| !s.is_empty()) {
        let (i, terms) = pair(entry)?;
        let i = id(i)?;
        if !integrals.contains_key(&i) {
            return Err(error("target missing from dictionary"));
        }
        let mut result = Vec::new();
        for term in list(terms)?.into_iter().filter(|s| !s.is_empty()) {
            let (j, c) = pair(term)?;
            let j = id(j)?;
            if !integrals.contains_key(&j) {
                return Err(error("master missing from dictionary"));
            }
            let c = c
                .strip_prefix('"')
                .and_then(|c| c.strip_suffix('"'))
                .ok_or_else(|| error("coefficient must be quoted"))?;
            result.push((j, c.into()));
        }
        if reductions.insert(i, result).is_some() {
            return Err(error("duplicate reduction"));
        }
    }
    Ok(ReductionTable {
        integrals,
        reductions,
    })
}

impl ReductionTable {
    /// Evaluate a table row over a prime, keyed by external integral indices.
    /// A missing row, unbound symbol, or singular coefficient rejects the probe.
    pub fn evaluate(
        &self,
        family: u32,
        indices: &[i32],
        vars: &[String],
        x: &[u64],
        p: u64,
    ) -> Option<BTreeMap<(u32, Vec<i32>), u64>> {
        use polycore::modp::add;
        use polycore::{Fp, Modular};
        if vars.len() != x.len() || x.iter().any(|&v| v >= p) {
            return None;
        }
        let (i, _) = self
            .integrals
            .iter()
            .find(|(_, v)| v.0 == family && v.1 == indices)?;
        let params = vars
            .iter()
            .zip(x)
            .map(|(s, &v)| (s.as_str(), Fp::new(v, p)))
            .collect::<Vec<_>>();
        let ring = Ring::new(Vec::<String>::new(), Order::Lex);
        let mut out = BTreeMap::new();
        for (j, c) in self.reductions.get(i)? {
            let c = ring
                .parse_with(
                    &explicit_products(c),
                    &|v| Fp::new(polycore::crt::reduce(&Q::from_integer(v), p).unwrap(), p),
                    &params,
                )
                .ok()?;
            let value = c.lc().map_or(0, |c| c.residue_mod(p));
            let entry = out.entry(self.integrals.get(j)?.clone()).or_insert(0);
            *entry = add(*entry, value, p);
        }
        out.retain(|_, c| *c != 0);
        Some(out)
    }
}
