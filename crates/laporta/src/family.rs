//! General quadratic propagators with polynomial kinematics, evaluated numerically.
//!
//! Scalar-product changes of coordinates are inverted once per sample; IBP rows
//! retain only index shifts and derivative references. Cuts and sector unions are
//! applied while seeding and to every term of the resulting equations.

use crate::{
    ibp::{compositions, key, System},
    Row,
};
use num_rational::BigRational;
use num_traits::Zero;
use polycore::modp::{add, from_signed, mul, sub};
use polycore::{crt, modp_echelon, Fp, Monomial, Order, Poly};
use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::sync::{Arc, Mutex};
use zippel_interp::ModPoly;

mod analysis;
pub use analysis::FamilyAnalysis;

type Q = BigRational;
type SampleValue = Arc<Option<Vec<u64>>>;

/// A complete set of denominators, each linear in loop scalar products.
#[derive(Clone, Debug)]
pub struct PolynomialFamily {
    pub vars: Vec<String>,
    pub loops: usize,
    /// Propagator polynomials in `(loop momenta, external momenta, vars...)`.
    pub propagators: Vec<Poly<Q>>,
    /// External scalar products (not twice the products), in `vars`.
    pub products: Vec<Vec<Poly<Q>>>,
    /// Allowed top sectors, bit i corresponding to propagator i (including holes).
    pub top_sectors: Vec<u32>,
    /// Zero sectors and all their subsectors.
    pub zero_sectors: Vec<u32>,
    /// Zero-based propagators required to have positive indices.
    pub cuts: Vec<usize>,
}

#[derive(Clone, Debug)]
struct Term {
    column: usize,
    derivative: Option<usize>,
    factor: i32,
}
#[derive(Debug)]
pub(super) struct NumericRows {
    rows: Vec<Vec<Term>>,
    coordinates: Vec<Vec<Poly<Q>>>,
    derivatives: Vec<Vec<Poly<Q>>>,
    prime: Mutex<BTreeMap<u64, Arc<ModularFamily>>>,
    samples: Mutex<VecDeque<(Vec<u64>, SampleValue)>>,
}
#[derive(Debug)]
struct ModularFamily {
    coordinates: Vec<Vec<ModPoly>>,
    derivatives: Vec<Vec<ModPoly>>,
}

impl PolynomialFamily {
    const fn nprops(&self) -> usize {
        self.propagators.len()
    }
    fn allowed(&self, a: &[i32]) -> bool {
        let sector = a
            .iter()
            .enumerate()
            .filter(|(_, v)| **v > 0)
            .fold(0u32, |s, (i, _)| s | (1 << i));
        sector != 0
            && (0..self.loops).all(|i| {
                a.iter()
                    .zip(&self.propagators)
                    .any(|(&power, g)| power > 0 && g.terms.iter().any(|(m, _)| m.exps()[i] > 0))
            })
            && self.top_sectors.iter().any(|&top| sector & !top == 0)
            && !self.zero_sectors.iter().any(|&zero| sector & !zero == 0)
            && self.cuts.iter().all(|&i| a[i] > 0)
    }
    /// Check dimensions and polynomial momentum degrees without sampling.
    pub fn valid(&self) -> bool {
        let l = self.loops;
        let e = self.products.len();
        let n = l * (l + 1) / 2 + l * e;
        l > 0
            && n < 32
            && self.nprops() == n
            && self.vars.first().is_some_and(|v| v == "d")
            && self
                .products
                .iter()
                .all(|r| r.len() == e && r.iter().all(|g| g.nvars == self.vars.len()))
            && (0..e).all(|i| (0..e).all(|j| self.products[i][j] == self.products[j][i]))
            && self.propagators.iter().all(|g| {
                g.nvars == l + e + self.vars.len()
                    && g.terms
                        .iter()
                        .all(|(m, _)| matches!(m.exps()[..l + e].iter().sum::<u32>(), 0 | 2))
            })
            && !self.top_sectors.is_empty()
            && !self.top_sectors.contains(&0)
            && self
                .top_sectors
                .iter()
                .chain(&self.zero_sectors)
                .all(|&s| s < (1 << n))
            && self.cuts.iter().all(|&i| i < n)
    }
    /// Seed a union of sectors; supports polynomial kinematics with rational coefficients.
    #[allow(clippy::too_many_lines)]
    pub fn system(&self, dots: i32, numerators: i32) -> Option<System> {
        if !self.valid() || dots < 0 || numerators < 0 {
            return None;
        }
        let analysis = self.analyze().unwrap_or_default();
        let mut family = self.clone();
        family.zero_sectors.extend(&analysis.zero_sectors);
        let allowed = |a: &[i32]| family.allowed(a);
        let n = self.nprops();
        let nm = self.loops + self.products.len();
        let pairs: Vec<_> = (0..self.loops)
            .flat_map(|i| (i..nm).map(move |j| (i, j)))
            .collect();
        let coordinates = self
            .propagators
            .iter()
            .map(|g| self.scalar_form(g, &pairs))
            .collect::<Option<Vec<_>>>()?;
        let mut derivatives = Vec::new();
        for i in 0..self.loops {
            for j in 0..nm {
                for g in &self.propagators {
                    // On scalar products, d/dk_i . v_j replaces one momentum occurrence.
                    let dg = &g.derivative(i) * &Poly::var(j, g.nvars, Order::Lex);
                    derivatives.push(self.scalar_form(&dg, &pairs)?);
                }
            }
        }
        let mut sectors = BTreeSet::new();
        for &top in &self.top_sectors {
            let mut s = top;
            loop {
                if s != 0 {
                    sectors.insert(s);
                }
                if s == 0 {
                    break;
                }
                s = (s - 1) & top;
            }
        }
        let mut raw = Vec::new();
        let mut zeros = Vec::new();
        let mut seeds = BTreeSet::new();
        let mut seen = BTreeSet::new();
        for sector in sectors {
            let (on, off): (Vec<_>, Vec<_>) = (0..n).partition(|&i| sector >> i & 1 == 1);
            for up in compositions(on.len(), dots) {
                for down in compositions(off.len(), numerators) {
                    let mut a = vec![0; n];
                    for (&i, &v) in on.iter().zip(&up) {
                        a[i] = v + 1;
                    }
                    for (&i, &v) in off.iter().zip(&down) {
                        a[i] = -v;
                    }
                    if !allowed(&a) {
                        if self.allowed(&a) {
                            seen.insert(a.clone());
                            seeds.insert(a.clone());
                            zeros.push(a);
                        }
                        continue;
                    }
                    seeds.insert(a.clone());
                    seen.insert(a.clone());
                    for i in 0..self.loops {
                        for j in 0..nm {
                            let mut terms = Vec::new();
                            if i == j {
                                terms.push((a.clone(), None, 1));
                            }
                            for m in 0..n {
                                if a[m] == 0 {
                                    continue;
                                }
                                let mut b = a.clone();
                                b[m] += 1;
                                for k in 0..=n {
                                    let mut shifted = b.clone();
                                    if k < n {
                                        shifted[k] -= 1;
                                    }
                                    if allowed(&shifted) {
                                        seen.insert(shifted.clone());
                                        terms.push((
                                            shifted,
                                            Some(((i * nm + j) * n + m) * (n + 1) + k),
                                            -a[m],
                                        ));
                                    }
                                }
                            }
                            if !terms.is_empty() {
                                raw.push(terms);
                            }
                        }
                    }
                }
            }
        }
        let mut symmetry_rows = Vec::new();
        for a in seen.clone() {
            for b in analysis.images(&a) {
                if allowed(&b) {
                    seen.insert(b.clone());
                    symmetry_rows.push((a.clone(), b));
                }
            }
        }
        let mut integrals: Vec<_> = seen.into_iter().collect();
        integrals.sort_by_key(|a| key(a, n));
        let columns: BTreeMap<_, _> = integrals
            .iter()
            .cloned()
            .enumerate()
            .map(|(i, a)| (a, i))
            .collect();
        let rows = raw
            .into_iter()
            .map(|r| {
                r.into_iter()
                    .map(|(a, derivative, factor)| Term {
                        column: columns[&a],
                        derivative,
                        factor,
                    })
                    .collect()
            })
            .collect();
        let numeric = NumericRows {
            rows,
            coordinates,
            derivatives,
            prime: Mutex::default(),
            samples: Mutex::default(),
        };
        let mut eqs = symmetry_rows
            .into_iter()
            .map(|(a, b)| vec![(columns[&a], vec![1]), (columns[&b], vec![-1])])
            .collect::<Vec<_>>();
        eqs.extend(zeros.into_iter().map(|a| vec![(columns[&a], vec![1])]));
        Some(System {
            integrals,
            vars: self.vars.clone(),
            seeds,
            columns,
            eqs,
            numeric: Some(numeric),
        })
    }
    fn scalar_form(&self, g: &Poly<Q>, pairs: &[(usize, usize)]) -> Option<Vec<Poly<Q>>> {
        let nm = self.loops + self.products.len();
        let mut out = vec![Poly::constant(Q::zero(), self.vars.len(), Order::Lex); pairs.len() + 1];
        for (m, c) in &g.terms {
            let coefficient = Poly::new(
                vec![(Monomial::new(m.exps()[nm..].to_vec()), c.clone())],
                self.vars.len(),
                Order::Lex,
            );
            let factors: Vec<_> = m.exps()[..nm]
                .iter()
                .enumerate()
                .flat_map(|(j, &e)| std::iter::repeat_n(j, e as usize))
                .collect();
            let (index, value) = match factors.as_slice() {
                [] => (pairs.len(), coefficient),
                &[i, j] if i >= self.loops => (
                    pairs.len(),
                    &coefficient * &self.products[i - self.loops][j - self.loops],
                ),
                &[i, j] => (pairs.iter().position(|&v| v == (i, j))?, coefficient),
                _ => return None,
            };
            out[index] = &out[index] + &value;
        }
        Some(out)
    }
}

impl NumericRows {
    pub(super) const fn len(&self) -> usize {
        self.rows.len()
    }
    fn evaluate(&self, x: &[u64], p: u64) -> Arc<Option<Vec<u64>>> {
        let mut key = x.to_vec();
        key.push(p);
        if let Some((_, v)) = self.samples.lock().unwrap().iter().find(|(k, _)| *k == key) {
            return Arc::clone(v);
        }
        let value = (|| {
            let modular = {
                let mut cache = self.prime.lock().unwrap();
                if let Some(m) = cache.get(&p) {
                    Arc::clone(m)
                } else {
                    let convert = |rows: &Vec<Vec<Poly<Q>>>| {
                        rows.iter()
                            .map(|r| {
                                r.iter()
                                    .map(|g| {
                                        Some(ModPoly::from_poly(
                                            &g.try_map(|c| Some(Fp::new(crt::reduce(c, p)?, p)))?,
                                            p,
                                        ))
                                    })
                                    .collect::<Option<Vec<_>>>()
                            })
                            .collect::<Option<Vec<_>>>()
                    };
                    let m = Arc::new(ModularFamily {
                        coordinates: convert(&self.coordinates)?,
                        derivatives: convert(&self.derivatives)?,
                    });
                    cache.insert(p, Arc::clone(&m));
                    m
                }
            };
            let coordinates: Vec<Vec<_>> = modular
                .coordinates
                .iter()
                .map(|r| r.iter().map(|g| g.eval(x, p)).collect())
                .collect();
            let n = coordinates.len();
            let a = coordinates
                .iter()
                .map(|r| r[..n].to_vec())
                .collect::<Vec<_>>();
            let inverse = modp_echelon::invert(&a, p)?;
            let mut table = Vec::new();
            for row in &modular.derivatives {
                let v: Vec<_> = row.iter().map(|g| g.eval(x, p)).collect();
                let on: Vec<_> = (0..n)
                    .map(|j| (0..n).fold(0, |acc, k| add(acc, mul(v[k], inverse[k][j], p), p)))
                    .collect();
                let constant = on
                    .iter()
                    .zip(&coordinates)
                    .fold(v[n], |acc, (&c, r)| sub(acc, mul(c, r[n], p), p));
                table.extend(on);
                table.push(constant);
            }
            Some(table)
        })();
        let value = Arc::new(value);
        let mut cache = self.samples.lock().unwrap();
        if cache.len() == 32 {
            cache.pop_front();
        }
        cache.push_back((key, Arc::clone(&value)));
        value
    }
    pub(super) fn defined(&self, x: &[u64], p: u64) -> bool {
        self.evaluate(x, p).is_some()
    }
    pub(super) fn row(&self, e: usize, x: &[u64], p: u64) -> Option<Row> {
        let table = self.evaluate(x, p);
        let table = table.as_ref().as_ref()?;
        let mut row = BTreeMap::new();
        for t in &self.rows[e] {
            let c = t.derivative.map_or(x[0], |k| table[k]);
            let c = mul(c, from_signed(i64::from(t.factor), p), p);
            let value = row.entry(t.column).or_insert(0);
            *value = add(*value, c, p);
        }
        Some(row.into_iter().filter(|(_, c)| *c != 0).collect())
    }
}

impl PolynomialFamily {
    #[must_use]
    pub fn fix(mut self, name: &str, value: i64) -> Option<Self> {
        let k = self.vars.iter().position(|v| v == name)?;
        if k == 0 || !self.valid() {
            return None;
        }
        let remove = |g: &Poly<Q>, k| {
            let g = g.eval_var(k, &Q::from_integer(value.into()));
            Poly::new(
                g.terms
                    .iter()
                    .map(|(m, c)| {
                        let mut e = m.exps().to_vec();
                        e.remove(k);
                        (Monomial::new(e), c.clone())
                    })
                    .collect(),
                g.nvars - 1,
                Order::Lex,
            )
        };
        let nm = self.loops + self.products.len();
        self.propagators = self.propagators.iter().map(|g| remove(g, nm + k)).collect();
        self.products = self
            .products
            .iter()
            .map(|r| r.iter().map(|g| remove(g, k)).collect())
            .collect();
        self.vars.remove(k);
        Some(self)
    }
    /// Reject a dependent scalar-product coordinate system at several primes/points.
    pub fn independent(&self) -> bool {
        if !self.valid() {
            return false;
        }
        let nm = self.loops + self.products.len();
        let pairs: Vec<_> = (0..self.loops)
            .flat_map(|i| (i..nm).map(move |j| (i, j)))
            .collect();
        let Some(coordinates) = self
            .propagators
            .iter()
            .map(|g| self.scalar_form(g, &pairs))
            .collect::<Option<Vec<_>>>()
        else {
            return false;
        };
        let rows = NumericRows {
            rows: vec![],
            coordinates,
            derivatives: vec![],
            prime: Mutex::default(),
            samples: Mutex::default(),
        };
        polycore::modp::Primes::new().take(3).any(|p| {
            (0..3).any(|i| {
                rows.defined(
                    &(0..self.vars.len() as u64)
                        .map(|j| polycore::sample::point(&[113, i, j], p))
                        .collect::<Vec<_>>(),
                    p,
                )
            })
        })
    }
    pub fn reduce(
        &self,
        targets: &[crate::ibp::Index],
        dots: i32,
        numerators: i32,
    ) -> Option<(System, crate::Plan, Vec<zippel_lift::Fraction>)> {
        let system = self.system(dots, numerators)?;
        if targets.iter().any(|t| !system.integrals.contains(t)) {
            return None;
        }
        let plan = system.learn(targets, 1)?;
        let cs = system.lift(&plan, 1)?;
        Some((system, plan, cs))
    }
}

impl PolynomialFamily {
    /// A bounded, optional pool: dimension shifts, kinematic linear factors,
    /// scalar products, masses, and principal external Gram determinants.
    /// Guessing tests exact divisibility on independent slices; an incomplete
    /// pool is harmless because reconstruction retains a residual denominator.
    pub fn denominator_candidates(&self) -> Vec<Poly<Q>> {
        use num_traits::One;
        let n = self.vars.len();
        let zero = Poly::constant(Q::zero(), n, Order::Lex);
        let mut pool = Vec::new();
        let mut insert = |g: Poly<Q>| {
            if !g.is_zero()
                && g.terms
                    .iter()
                    .any(|(m, _)| m.exps().iter().any(|&e| e != 0))
            {
                let g = g.monic();
                if !pool.contains(&g) {
                    pool.push(g);
                }
            }
        };
        let d = Poly::var(0, n, Order::Lex);
        for denominator in 1..=2 {
            for numerator in -2 * denominator..=8 * denominator {
                insert(
                    &d - &Poly::constant(
                        Q::new(numerator.into(), denominator.into()),
                        n,
                        Order::Lex,
                    ),
                );
            }
        }
        for i in 1..n {
            let a = Poly::var(i, n, Order::Lex);
            insert(a.clone());
            insert(&a - &Poly::constant(Q::one(), n, Order::Lex));
            insert(&a + &Poly::constant(Q::one(), n, Order::Lex));
            for j in i + 1..n {
                let b = Poly::var(j, n, Order::Lex);
                insert(&a + &b);
                insert(&a - &b);
            }
        }
        for row in &self.products {
            for g in row {
                insert(g.clone());
            }
        }
        let nm = self.loops + self.products.len();
        for g in &self.propagators {
            insert(Poly::new(
                g.terms
                    .iter()
                    .filter(|(m, _)| m.exps()[..nm].iter().all(|&e| e == 0))
                    .map(|(m, c)| (Monomial::new(m.exps()[nm..].to_vec()), c.clone()))
                    .collect(),
                n,
                Order::Lex,
            ));
        }
        // Determinants by subset dynamic programming: O(r 2^r), no division.
        let e = self.products.len();
        if e <= 6 {
            for subset in 1usize..1 << e {
                let indices: Vec<_> = (0..e).filter(|&i| subset >> i & 1 == 1).collect();
                let r = indices.len();
                let mut det = vec![zero.clone(); 1 << r];
                det[0] = Poly::constant(Q::one(), n, Order::Lex);
                for mask in 1usize..1 << r {
                    let row = mask.count_ones() as usize - 1;
                    let mut position = 0;
                    for j in 0..r {
                        if mask >> j & 1 == 1 {
                            let term =
                                &det[mask ^ (1 << j)] * &self.products[indices[row]][indices[j]];
                            det[mask] = if (row + position).is_multiple_of(2) {
                                &det[mask] + &term
                            } else {
                                &det[mask] - &term
                            };
                            position += 1;
                        }
                    }
                }
                insert(det.pop().unwrap());
            }
        }
        pool
    }
}
