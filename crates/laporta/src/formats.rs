//! Declarative Kira YAML and FIRE Mathematica family interchange.
//!
//! The IBP engine currently accepts complete quadratic families with affine
//! kinematics and integer momentum coefficients. Unsupported definitions are
//! errors, never evaluated as Mathematica programs or silently approximated.

#![allow(clippy::too_many_lines)]

use crate::{
    ibp::{Family, Lin, System},
    Plan,
};
use num_rational::BigRational;
use num_traits::{One, ToPrimitive, Zero};
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
    pub family: Family,
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

fn affine(g: &Poly<Q>) -> Result<Vec<Q>> {
    let mut out = vec![Q::zero(); g.nvars + 1];
    for (m, c) in &g.terms {
        let total: u32 = m.exps().iter().sum();
        if total > 1 {
            return Err(error(
                "non-affine kinematics: introduce squared-mass invariants",
            ));
        }
        let k = m.exps().iter().position(|&e| e == 1).map_or(0, |j| j + 1);
        out[k] += c;
    }
    Ok(out)
}
fn integers(cs: Vec<Q>) -> Result<Vec<i64>> {
    cs.into_iter()
        .map(|c| {
            if !c.is_integer() {
                return Err(error("fractional coefficient unsupported by IBP family"));
            }
            c.to_integer()
                .to_i64()
                .ok_or_else(|| error("coefficient exceeds i64"))
        })
        .collect()
}
fn momentum(r: &Ring, s: &str) -> Result<Vec<i64>> {
    let cs = integers(affine(&parse(r, s)?)?)?;
    if cs[0] != 0 {
        return Err(error("momentum contains scalar constant"));
    }
    Ok(cs[1..].to_vec())
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

/// Read one named family from Kira's `integralfamilies.yaml` and `kinematics.yaml`.
///
/// Supports one contiguous top sector, momentum conservation, scalar-product
/// rules involving momentum sums, and fixing one invariant to unity.
///
/// # Errors
/// Rejects malformed, unsupported, incomplete, or dependent family definitions.
pub fn read_kira(families: &str, kinematics: &str, name: &str) -> Result<ImportedFamily> {
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
        &["name", "loop_momenta", "top_level_sectors", "propagators"],
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
    let legs = scalar_products(&er, &vr, &rules)?;
    let props: Vec<_> = seq(&f["propagators"])?
        .iter()
        .map(|v| {
            let v = seq(v)?;
            if v.len() != 2 {
                return Err(error("propagator requires [momentum, mass_squared]"));
            }
            Ok((
                momentum(&mr, &substitute(atom(&v[0])?))?,
                integers(affine(&parse(&vr, &atom(&v[1])?)?)?)?,
            ))
        })
        .collect::<Result<_>>()?;
    let sectors = seq(&f["top_level_sectors"])?;
    if sectors.len() != 1 {
        return Err(error("exactly one top_level_sector is currently supported"));
    }
    let sector = sectors[0]
        .as_u64()
        .ok_or_else(|| error("sector must be an integer"))?;
    if sector == 0 || sector.checked_add(1).is_none_or(|s| !s.is_power_of_two()) {
        return Err(error("top sector must be contiguous: 2^lines-1"));
    }
    let lines = sector.count_ones() as usize;
    let signs = vec![1; props.len()];
    let mut family = Family {
        vars,
        loops: loops.len(),
        props,
        lines,
        legs,
        symmetries: vec![],
    };
    validate(&family)?;
    if !k["symbol_to_replace_by_one"].is_null() {
        let fixed = atom(&k["symbol_to_replace_by_one"])?;
        if fixed == "d" || !family.vars.contains(&fixed) {
            return Err(error("invalid fixed invariant"));
        }
        family = family.fix(&fixed, 1);
    }
    Ok(ImportedFamily {
        name: name.into(),
        family,
        signs,
    })
}

fn scalar_products(er: &Ring, vr: &Ring, rules: &[(String, String)]) -> Result<Vec<Vec<Lin>>> {
    let n = er.nvars();
    let pairs: Vec<_> = (0..n).flat_map(|i| (i..n).map(move |j| (i, j))).collect();
    let mut a = Vec::new();
    let mut rhs = Vec::new();
    for (lhs, rhs_expr) in rules {
        let g = parse(er, lhs)?;
        let mut row = vec![Q::zero(); pairs.len()];
        for (m, c) in g.terms {
            let indices: Vec<_> = m
                .exps()
                .iter()
                .enumerate()
                .flat_map(|(i, &e)| std::iter::repeat_n(i, e as usize))
                .collect();
            if indices.len() != 2 {
                return Err(error(
                    "scalar-product rule must be quadratic in external momenta",
                ));
            }
            let pos = pairs
                .iter()
                .position(|&(i, j)| i == indices[0] && j == indices[1])
                .ok_or_else(|| error("unknown scalar product"))?;
            row[pos] += c;
        }
        a.push(row);
        rhs.push(affine(&parse(vr, rhs_expr)?)?);
    }
    let mut values = vec![vec![Q::zero(); vr.nvars() + 1]; pairs.len()];
    for k in 0..=vr.nvars() {
        let b: Vec<_> = rhs.iter().map(|r| r[k].clone()).collect();
        let x = dense::solve(&a, &b).map_err(|e| error(format!("scalar-product rules: {e}")))?;
        if x.len() != pairs.len() {
            return Err(error("incomplete scalar-product rules"));
        }
        for (v, c) in values.iter_mut().zip(x) {
            v[k] = c * Q::from_integer(2.into());
        }
    }
    let mut legs = vec![vec![vec![]; n]; n];
    for ((i, j), v) in pairs.into_iter().zip(values) {
        let v = integers(v)?;
        legs[i][j].clone_from(&v);
        legs[j][i] = v;
    }
    Ok(legs)
}

fn validate(f: &Family) -> Result<()> {
    if f.loops == 0 || f.lines == 0 || f.lines >= 32 || f.lines > f.props.len() {
        return Err(error("invalid loop/line count"));
    }
    let n = f.loops + f.legs.len();
    let pairs: Vec<_> = (0..f.loops)
        .flat_map(|i| (i..n).map(move |j| (i, j)))
        .collect();
    if pairs.len() != f.props.len() {
        return Err(error(
            "propagators must span all loop scalar products, including ISPs",
        ));
    }
    let a: Vec<Vec<Q>> = f
        .props
        .iter()
        .map(|(q, _)| {
            pairs
                .iter()
                .map(|&(i, j)| Q::from_integer((q[i] * q[j] * if i == j { 1 } else { 2 }).into()))
                .collect()
        })
        .collect();
    if dense::invert(&a).is_none() {
        return Err(error("linearly dependent propagators"));
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

/// Read literal FIRE assignments: Internal, External, Propagators, Replacements.
///
/// `invariants` gives the independent scalar names (excluding d); `lines` excludes
/// trailing ISPs. General Mathematica evaluation and executable statements are rejected.
///
/// # Errors
/// Rejects executable input, unsupported propagators, and incomplete kinematics.
pub fn read_fire(
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
    let mut assignments = BTreeMap::new();
    for statement in split(&clean, ';')?.into_iter().filter(|s| !s.is_empty()) {
        let (key, value) = statement
            .split_once('=')
            .ok_or_else(|| error("only literal FIRE assignments are supported"))?;
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
    let rules: Vec<_> = list(get("Replacements")?)?
        .into_iter()
        .filter(|s| !s.is_empty())
        .map(|s| {
            let (a, b) = s
                .split_once("->")
                .ok_or_else(|| error("expected replacement arrow"))?;
            Ok((a.trim().into(), b.trim().into()))
        })
        .collect::<Result<_>>()?;
    let legs = scalar_products(&er, &vr, &rules)?;
    let mut props = Vec::new();
    let mut signs = Vec::new();
    for src in list(get("Propagators")?)? {
        let g = parse(&full, src)?;
        let coeff = |exps: &[u32]| {
            g.terms
                .iter()
                .find(|(m, _)| m.exps() == exps)
                .map_or_else(Q::zero, |(_, c)| c.clone())
        };
        let (base, sign) = (0..loops.len())
            .find_map(|i| {
                let mut e = vec![0; full.nvars()];
                e[i] = 2;
                let c = coeff(&e);
                if c == Q::one() {
                    Some((i, 1))
                } else if c == -Q::one() {
                    Some((i, -1))
                } else {
                    None
                }
            })
            .ok_or_else(|| error("propagator needs a loop-square coefficient +1 or -1"))?;
        let mut q = vec![0; mr.nvars()];
        q[base] = 1;
        for (j, v) in q.iter_mut().enumerate().filter(|(j, _)| *j != base) {
            let mut e = vec![0; full.nvars()];
            e[base] = 1;
            e[j] = 1;
            *v = integers(vec![coeff(&e) / Q::from_integer((2 * sign).into())])?[0];
        }
        let momentum = Poly::new(
            q.iter()
                .enumerate()
                .filter(|(_, v)| **v != 0)
                .map(|(i, &v)| {
                    let mut e = vec![0; full.nvars()];
                    e[i] = 1;
                    (Monomial::new(e), Q::from_integer(v.into()))
                })
                .collect(),
            full.nvars(),
            Order::Lex,
        );
        let mass = &momentum.pow(2) - &g.scale(&Q::from_integer(sign.into()));
        if mass
            .terms
            .iter()
            .any(|(m, _)| m.exps()[..mr.nvars()].iter().any(|&e| e != 0))
        {
            return Err(error(
                "propagator is not a signed momentum square minus affine mass",
            ));
        }
        let mass = Poly::new(
            mass.terms
                .iter()
                .map(|(m, c)| (Monomial::new(m.exps()[mr.nvars()..].to_vec()), c.clone()))
                .collect(),
            vars.len(),
            Order::Lex,
        );
        props.push((q, integers(affine(&mass)?)?));
        signs.push(sign);
    }
    let family = Family {
        vars,
        loops: loops.len(),
        props,
        lines,
        legs,
        symmetries: vec![],
    };
    validate(&family)?;
    Ok(ImportedFamily {
        name: name.into(),
        family,
        signs,
    })
}

impl ImportedFamily {
    /// Mathematica replacement rules, using Kira's family head or FIRE's `G[id,{...}]`.
    ///
    /// # Errors
    /// Rejects incompatible coefficient dimensions or invalid variable names.
    pub fn rules(
        &self,
        system: &System,
        plan: &Plan,
        coefficients: &[Fraction],
        fire_id: Option<u32>,
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
            fire_id.map_or_else(
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
    /// FIRE `.tables`: reduction entries plus an explicit integral-ID dictionary.
    ///
    /// # Errors
    /// Rejects incompatible coefficient dimensions or invalid variable names.
    pub fn fire_tables(
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
