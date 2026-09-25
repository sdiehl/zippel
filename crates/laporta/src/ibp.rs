//! Integration-by-parts identities of a family `I(a) = int prod_i d^d k_i / prod_m D_m^a_m`.
//!
//! Each `D_m` is quadratic in the momenta, so `d/dk_i . q_j` of it is linear in the scalar
//! products, and those involving a loop momentum are linear in the `D`s when the family is
//! complete. Every identity is then a relation between integrals whose indices differ by one,
//! with coefficients linear in `d` and the kinematics. Integrals are ordered as Laporta does:
//! more lines, then more dots, then more numerators is more complex.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write;

use groebner::{MonomialOrder, PolynomialRing};
use num_rational::{BigRational, Rational64};
use num_traits::{One, Zero};
use zippel_interp::modp::{add, mul, point};
use zippel_interp::Primes;
use zippel_lift::{lift, Fraction};

use crate::{Plan, Row};

/// `c_0 + c_1 x_1 + ...` in the variables `x = (d, kinematics..)`.
pub type Lin = Vec<i64>;

pub type Index = Vec<i32>;

type Form = Vec<Rational64>;

#[derive(Clone, Debug)]
pub struct Family {
    /// Variable names, `d` first.
    pub vars: Vec<&'static str>,
    pub loops: usize,
    /// Each `D_m` as its momentum on the loop then external momenta, and its squared mass. Only
    /// the first `lines` may be propagators; the rest are irreducible numerators.
    pub props: Vec<(Vec<i64>, Lin)>,
    pub lines: usize,
    /// `2 p_i . p_j` for the external momenta.
    pub legs: Vec<Vec<Lin>>,
    /// Permutations of the indices that map the family to itself.
    pub symmetries: Vec<Vec<usize>>,
}

/// A combination of the `D`s plus a constant.
#[derive(Clone, Debug)]
struct Expr {
    on: Vec<Rational64>,
    constant: Form,
}

#[derive(Debug)]
pub struct System {
    pub integrals: Vec<Index>,
    pub vars: Vec<&'static str>,
    seeds: BTreeSet<Index>,
    columns: BTreeMap<Index, usize>,
    eqs: Vec<Vec<(usize, Lin)>>,
}

impl Family {
    fn form(&self, l: &[i64], scale: Rational64) -> Form {
        let mut f = vec![Rational64::zero(); self.vars.len() + 1];
        for (a, &b) in f.iter_mut().zip(l) {
            *a = scale * b;
        }
        f
    }

    /// `d/dk_i . v_j` of each `D_m`, with `v` the loop then external momenta.
    fn derivatives(&self) -> Vec<Vec<Vec<Expr>>> {
        let (l, n) = (self.loops, self.props.len());
        let half = Rational64::new(1, 2);
        let mut sigma = BTreeMap::new();
        for u in 0..l {
            for v in u..l + self.legs.len() {
                let s = sigma.len();
                sigma.insert((u, v), s);
            }
        }
        assert_eq!(sigma.len(), n, "the family must span every scalar product");
        let product = |u: usize, v: usize| -> Result<usize, Form> {
            let (u, v) = (u.min(v), u.max(v));
            sigma
                .get(&(u, v))
                .copied()
                .ok_or_else(|| self.form(&self.legs[u - l][v - l], half))
        };
        let mut a = vec![vec![Rational64::zero(); n]; n];
        let mut c: Vec<Form> = Vec::new();
        for (m, (q, mass)) in self.props.iter().enumerate() {
            let mut cm = self.form(mass, -Rational64::one());
            for (u, &x) in q.iter().enumerate() {
                for (v, &y) in q.iter().enumerate() {
                    match product(u, v) {
                        Ok(s) => a[m][s] += x * y,
                        Err(f) => axpy(&mut cm, &f, Rational64::from(x * y)),
                    }
                }
            }
            c.push(cm);
        }
        let b = invert(a);
        let sp: Vec<Expr> = b
            .iter()
            .map(|row| {
                let mut constant = self.form(&[], Rational64::zero());
                for (bk, ck) in row.iter().zip(&c) {
                    axpy(&mut constant, ck, -*bk);
                }
                Expr {
                    on: row.clone(),
                    constant,
                }
            })
            .collect();
        let expr = |u, v| match product(u, v) {
            Ok(s) => sp[s].clone(),
            Err(constant) => Expr {
                on: vec![Rational64::zero(); n],
                constant,
            },
        };
        (0..l)
            .map(|i| {
                (0..l + self.legs.len())
                    .map(|j| {
                        self.props
                            .iter()
                            .map(|(q, _)| {
                                let mut t = Expr {
                                    on: vec![Rational64::zero(); n],
                                    constant: self.form(&[], Rational64::zero()),
                                };
                                for (u, &x) in q.iter().enumerate() {
                                    let e = expr(u, j);
                                    let k = Rational64::from(2 * q[i] * x);
                                    axpy(&mut t.on, &e.on, k);
                                    axpy(&mut t.constant, &e.constant, k);
                                }
                                t
                            })
                            .collect()
                    })
                    .collect()
            })
            .collect()
    }

    /// Whether every loop momentum flows through a line of the sector, without which the integral
    /// is scaleless.
    fn nonzero(&self, a: &[i32]) -> bool {
        (0..self.loops).all(|i| {
            a.iter()
                .zip(&self.props)
                .any(|(&x, (q, _))| x > 0 && q[i] != 0)
        })
    }

    /// The family with the kinematic variable `name` set to `value`, one fewer to reconstruct.
    /// Setting a scale to 1 loses nothing, since dimensional analysis puts it back.
    #[must_use]
    pub fn fix(mut self, name: &str, value: i64) -> Self {
        let k = 1 + self
            .vars
            .iter()
            .position(|&v| v == name)
            .expect("a variable");
        assert!(k > 1, "d stays symbolic");
        let sub = |l: &mut Lin| {
            if l.len() > k {
                l[0] += l[k] * value;
                l.remove(k);
            }
        };
        self.props.iter_mut().for_each(|(_, m)| sub(m));
        self.legs.iter_mut().flatten().for_each(sub);
        self.vars.remove(k - 1);
        self
    }

    /// Reduces `targets` to masters. Seeds carry fewer dots and numerators the fewer lines their
    /// sector has, up to `dots` and `numerators` in the top sector, and a sector is widened
    /// whenever one of its integrals comes out a master without being a seed, the mark of an
    /// identity missing there.
    pub fn reduce(
        &self,
        targets: &[Index],
        dots: i32,
        numerators: i32,
    ) -> Option<(System, Plan, Vec<Fraction>)> {
        let mut extra = BTreeMap::new();
        loop {
            let system = self.seeded(dots, numerators, &extra);
            let plan = system.learn(targets, 1);
            let mut widened = false;
            for &m in &plan.masters {
                let a = &system.integrals[m];
                let w = extra.entry(self.sector(a)).or_insert(0);
                if !system.seeds.contains(a) && *w < dots.max(numerators) {
                    *w += 1;
                    widened = true;
                }
            }
            if !widened {
                let coefficients = system.lift(&plan, 1)?;
                return Some((system, plan, coefficients));
            }
        }
    }

    fn sector(&self, a: &[i32]) -> u32 {
        (0..self.lines).filter(|&m| a[m] > 0).map(|m| 1 << m).sum()
    }

    /// Identities seeded at every index with up to `dots` extra powers on its lines and
    /// `numerators` powers of numerators in every sector, plus the symmetries.
    pub fn system(&self, dots: i32, numerators: i32) -> System {
        let full = (1..1u32 << self.lines)
            .map(|s| (s, dots.max(numerators)))
            .collect();
        self.seeded(dots, numerators, &full)
    }

    fn seeded(&self, dots: i32, numerators: i32, extra: &BTreeMap<u32, i32>) -> System {
        let n = self.props.len();
        let table = self.derivatives();
        let mut raw: Vec<BTreeMap<Index, Form>> = Vec::new();
        let mut all = BTreeSet::new();
        for sector in 1..1u32 << self.lines {
            let (on, off): (Vec<usize>, Vec<usize>) =
                (0..n).partition(|&m| m < self.lines && sector >> m & 1 == 1);
            let level = -i32::try_from(self.lines - on.len()).unwrap();
            let widen = extra.get(&sector).copied().unwrap_or(0);
            let trim = |max: i32| max.min(max.min(1).max(max + level) + widen);
            let mut seeds = Vec::new();
            for up in compositions(on.len(), trim(dots)) {
                for down in compositions(off.len(), trim(numerators)) {
                    let mut a = vec![0; n];
                    on.iter().zip(&up).for_each(|(&m, &x)| a[m] = 1 + x);
                    off.iter().zip(&down).for_each(|(&m, &x)| a[m] = -x);
                    seeds.push(a);
                }
            }
            seeds.retain(|a| self.nonzero(a));
            for a in &seeds {
                for (i, row) in table.iter().enumerate() {
                    for (j, ds) in row.iter().enumerate() {
                        let mut eq = BTreeMap::new();
                        if i == j {
                            let mut d = self.form(&[], Rational64::zero());
                            d[1] = Rational64::one();
                            eq.insert(a.clone(), d);
                        }
                        for (m, t) in ds.iter().enumerate().filter(|t| a[t.0] != 0) {
                            let k = -Rational64::from(i64::from(a[m]));
                            let mut b = a.clone();
                            b[m] += 1;
                            term(&mut eq, &b, &t.constant, k);
                            for (kk, &c) in t.on.iter().enumerate().filter(|t| !t.1.is_zero()) {
                                let mut e = b.clone();
                                e[kk] -= 1;
                                let f = self.form(&[1], c);
                                term(&mut eq, &e, &f, k);
                            }
                        }
                        eq.retain(|b, f| self.nonzero(b) && f.iter().any(|c| !c.is_zero()));
                        if !eq.is_empty() {
                            raw.push(eq);
                        }
                    }
                }
            }
            all.extend(seeds);
        }
        let mut seen: BTreeSet<Index> = raw.iter().flat_map(|eq| eq.keys().cloned()).collect();
        for a in seen.clone() {
            for s in &self.symmetries {
                let mut b = vec![0; n];
                for (m, &x) in a.iter().enumerate() {
                    b[s[m]] = x;
                }
                if b != a {
                    let one = self.form(&[1], Rational64::one());
                    let mut eq = BTreeMap::new();
                    term(&mut eq, &a, &one, Rational64::one());
                    term(&mut eq, &b, &one, -Rational64::one());
                    seen.insert(b);
                    raw.push(eq);
                }
            }
        }
        let mut integrals: Vec<Index> = seen.into_iter().collect();
        integrals.sort_by_key(|a| key(a, self.lines));
        let columns: BTreeMap<Index, usize> = integrals
            .iter()
            .enumerate()
            .map(|(i, a)| (a.clone(), i))
            .collect();
        let eqs = raw
            .into_iter()
            .map(|eq| {
                let l = eq
                    .values()
                    .flatten()
                    .fold(1i64, |l, c| num_integer::lcm(l, *c.denom()));
                eq.into_iter()
                    .map(|(a, f)| {
                        (
                            columns[&a],
                            f.iter().map(|c| (c * l).to_integer()).collect(),
                        )
                    })
                    .collect()
            })
            .collect();
        System {
            integrals,
            vars: self.vars.clone(),
            seeds: all,
            columns,
            eqs,
        }
    }
}

/// Every vector of `k` naturals summing to at most `max`.
fn compositions(k: usize, max: i32) -> Vec<Vec<i32>> {
    (0..k).fold(vec![vec![]], |acc, _| {
        acc.into_iter()
            .flat_map(|v: Vec<i32>| {
                let left = max - v.iter().sum::<i32>();
                (0..=left).map(move |x| {
                    let mut w = v.clone();
                    w.push(x);
                    w
                })
            })
            .collect()
    })
}

fn weight(a: &[i32]) -> (i32, i32) {
    a.iter().fold((0, 0), |(pos, neg), &x| {
        if x > 0 {
            (pos + x, neg)
        } else {
            (pos, neg - x)
        }
    })
}

/// Laporta's order: lines, then dots, then numerators, then the indices from the last.
fn key(a: &[i32], lines: usize) -> (usize, i32, i32, Vec<i32>) {
    let (pos, neg) = weight(a);
    let t = a[..lines].iter().filter(|&&x| x > 0).count();
    (t, pos, neg, a.iter().rev().copied().collect())
}

fn term(eq: &mut BTreeMap<Index, Form>, a: &Index, f: &Form, k: Rational64) {
    let e = eq
        .entry(a.clone())
        .or_insert_with(|| vec![Rational64::zero(); f.len()]);
    axpy(e, f, k);
}

fn axpy(a: &mut [Rational64], b: &[Rational64], k: Rational64) {
    for (x, y) in a.iter_mut().zip(b) {
        *x += *y * k;
    }
}

fn invert(mut a: Vec<Vec<Rational64>>) -> Vec<Vec<Rational64>> {
    let n = a.len();
    let mut b: Vec<Vec<Rational64>> = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| Rational64::from(i64::from(i == j)))
                .collect()
        })
        .collect();
    for c in 0..n {
        let r = (c..n)
            .find(|&r| !a[r][c].is_zero())
            .expect("the propagators must be independent");
        a.swap(c, r);
        b.swap(c, r);
        let l = a[c][c].recip();
        a[c].iter_mut().chain(&mut b[c]).for_each(|x| *x *= l);
        for r in (0..n).filter(|&r| r != c) {
            let k = a[r][c];
            let (ac, bc) = (a[c].clone(), b[c].clone());
            axpy(&mut a[r], &ac, -k);
            axpy(&mut b[r], &bc, -k);
        }
    }
    b
}

impl System {
    pub const fn len(&self) -> usize {
        self.eqs.len()
    }

    pub const fn is_empty(&self) -> bool {
        self.eqs.is_empty()
    }

    pub fn column(&self, a: &Index) -> usize {
        self.columns[a]
    }

    /// Equation `e` at the numeric point `x`.
    pub fn row(&self, e: usize, x: &[u64], p: u64) -> Row {
        let residue = |v: i64| u64::try_from(i128::from(v).rem_euclid(i128::from(p))).unwrap();
        self.eqs[e]
            .iter()
            .map(|(j, c)| {
                let v = c[1..].iter().zip(x).fold(residue(c[0]), |acc, (&ci, &xi)| {
                    add(acc, mul(residue(ci), xi, p), p)
                });
                (*j, v)
            })
            .collect()
    }

    /// Eliminates everything at a random point and keeps what `targets` need.
    pub fn learn(&self, targets: &[Index], seed: u64) -> Plan {
        let targets: Vec<usize> = targets.iter().map(|t| self.column(t)).collect();
        let p = Primes::new().next().unwrap();
        let x: Vec<u64> = (0..self.vars.len() as u64)
            .map(|i| point(&[seed, 7, i], p))
            .collect();
        let rows: Vec<Row> = (0..self.len()).map(|e| self.row(e, &x, p)).collect();
        Plan::learn(&rows, &targets, p)
    }

    /// Every coefficient of `plan`, lifted to Q from replays.
    pub fn lift(&self, plan: &Plan, seed: u64) -> Option<Vec<Fraction>> {
        lift(
            |x: &[u64], p| plan.replay(|e| self.row(e, x, p), p),
            self.vars.len(),
            seed,
        )
    }

    pub fn name(&self, j: usize) -> String {
        let a: Vec<String> = self.integrals[j].iter().map(ToString::to_string).collect();
        format!("I({})", a.join(","))
    }

    /// `I(a) = c * I(..) + ...` per target.
    pub fn render(&self, plan: &Plan, coefficients: &[Fraction]) -> String {
        let ring = PolynomialRing::<BigRational>::new(self.vars.clone(), MonomialOrder::Lex)
            .expect("valid variable names");
        let show = |f: &Fraction| {
            let num = ring.format(&f.num).unwrap();
            match ring.format(&f.den).unwrap().as_str() {
                "1" => format!("({num})"),
                den => format!("({num}) / ({den})"),
            }
        };
        let mut out = String::new();
        let width = plan.masters.len().max(1);
        for (t, row) in plan.targets.iter().zip(coefficients.chunks(width)) {
            let terms: Vec<String> = row
                .iter()
                .zip(&plan.masters)
                .filter(|(f, _)| !f.num.is_zero())
                .map(|(f, &m)| format!("{} * {}", show(f), self.name(m)))
                .collect();
            let rhs = if terms.is_empty() {
                "0".into()
            } else {
                terms.join("\n    + ")
            };
            writeln!(out, "{} = {rhs}", self.name(*t)).unwrap();
        }
        out
    }
}
