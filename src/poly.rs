//! Sparse polynomials in lex order (variable 0 most significant), terms sorted descending.

use crate::modp::{add, inv, mul, pow};
use crate::univariate::{self as uni, Dense};
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use std::collections::BTreeMap;

pub(crate) type Exps = Vec<u32>;

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct ModPoly {
    pub(crate) n: usize,
    pub(crate) terms: Vec<(Exps, u64)>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct IntPoly {
    pub(crate) n: usize,
    pub(crate) terms: Vec<(Exps, BigInt)>,
}

impl ModPoly {
    pub(crate) fn from_int(f: &IntPoly, p: u64) -> Self {
        let m = BigInt::from(p);
        let terms = f
            .terms
            .iter()
            .filter_map(|(e, c)| {
                let r = u64::try_from(c.mod_floor(&m)).unwrap_or(0);
                (r != 0).then(|| (e.clone(), r))
            })
            .collect();
        Self { n: f.n, terms }
    }

    pub(crate) fn lm(&self) -> &[u32] {
        &self.terms[0].0
    }

    pub(crate) fn scale(&mut self, c: u64, p: u64) {
        self.terms.iter_mut().for_each(|t| t.1 = mul(t.1, c, p));
    }

    pub(crate) fn support(&self) -> Vec<Exps> {
        self.terms.iter().map(|t| t.0.clone()).collect()
    }

    pub(crate) fn degree(&self, k: usize) -> usize {
        self.terms
            .iter()
            .map(|t| t.0[k] as usize)
            .max()
            .unwrap_or(0)
    }

    /// Coefficients in `x_k` of each monomial in `x_0..x_{k-1}`, assuming later variables are gone.
    pub(crate) fn groups(&self, k: usize) -> Vec<(Exps, Dense)> {
        let mut out: Vec<(Exps, Dense)> = Vec::new();
        for (e, c) in &self.terms {
            let mut key = e.clone();
            key[k] = 0;
            if out.last().is_none_or(|g| g.0 != key) {
                out.push((key, Vec::new()));
            }
            let d = &mut out.last_mut().expect("just pushed").1;
            let i = e[k] as usize;
            if d.len() <= i {
                d.resize(i + 1, 0);
            }
            d[i] = *c;
        }
        out
    }

    pub(crate) fn from_groups(n: usize, k: usize, groups: Vec<(Exps, Dense)>) -> Self {
        let mut terms = Vec::new();
        for (key, d) in groups {
            for (i, &c) in d.iter().enumerate().rev().filter(|t| *t.1 != 0) {
                let mut e = key.clone();
                e[k] = i as u32;
                terms.push((e, c));
            }
        }
        Self { n, terms }
    }

    /// Content in `GF(p)[x_k]` and the primitive part.
    pub(crate) fn primitive(&self, k: usize, p: u64) -> (Dense, Self) {
        let groups = self.groups(k);
        let c = groups
            .iter()
            .fold(Vec::new(), |acc, g| uni::gcd(&acc, &g.1, p));
        if c.len() == 1 {
            return (c, self.clone());
        }
        let groups = groups
            .into_iter()
            .map(|(e, d)| (e, uni::div_rem(&d, &c, p).0))
            .collect();
        (c, Self::from_groups(self.n, k, groups))
    }

    pub(crate) fn eval_var(&self, k: usize, a: u64, p: u64) -> Self {
        let terms = self
            .groups(k)
            .into_iter()
            .map(|(e, d)| (e, uni::eval(&d, a, p)))
            .filter(|t| t.1 != 0)
            .collect();
        Self { n: self.n, terms }
    }

    /// Substitute `point[i]` for every `x_i` with `i != k`, leaving a dense polynomial in `x_k`.
    pub(crate) fn eval_except(&self, k: usize, point: &[u64], p: u64) -> Dense {
        let mut d = vec![0; self.degree(k) + 1];
        for (e, c) in &self.terms {
            let v = e
                .iter()
                .enumerate()
                .filter(|&(i, &x)| i != k && x != 0)
                .fold(*c, |acc, (i, &x)| {
                    mul(acc, pow(point[i], u64::from(x), p), p)
                });
            d[e[k] as usize] = add(d[e[k] as usize], v, p);
        }
        uni::trim(&mut d);
        d
    }

    pub(crate) fn monic(mut self, p: u64) -> Self {
        let l = inv(self.terms[0].1, p);
        self.scale(l, p);
        self
    }
}

impl IntPoly {
    pub(crate) fn new(n: usize, terms: impl IntoIterator<Item = (Exps, BigInt)>) -> Self {
        let mut map: BTreeMap<Exps, BigInt> = BTreeMap::new();
        for (e, c) in terms {
            *map.entry(e).or_default() += c;
        }
        let terms = map.into_iter().rev().filter(|t| !t.1.is_zero()).collect();
        Self { n, terms }
    }

    pub(crate) fn constant(n: usize, c: BigInt) -> Self {
        Self::new(n, [(vec![0; n], c)])
    }

    pub(crate) const fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    pub(crate) fn is_constant(&self) -> bool {
        self.terms.len() == 1 && self.terms[0].0.iter().all(|&e| e == 0)
    }

    pub(crate) fn lc(&self) -> &BigInt {
        &self.terms[0].1
    }

    pub(crate) fn degree(&self, k: usize) -> u32 {
        self.terms.iter().map(|t| t.0[k]).max().unwrap_or(0)
    }

    pub(crate) fn content(&self) -> BigInt {
        self.terms.iter().fold(BigInt::zero(), |g, t| g.gcd(&t.1))
    }

    pub(crate) fn min_exps(&self) -> Exps {
        (0..self.n)
            .map(|k| self.terms.iter().map(|t| t.0[k]).min().unwrap_or(0))
            .collect()
    }

    pub(crate) fn map(&self, f: impl Fn(&Exps, &BigInt) -> (Exps, BigInt)) -> Self {
        Self::new(self.n, self.terms.iter().map(|(e, c)| f(e, c)))
    }

    /// Divide by the integer content and make the leading coefficient positive.
    pub(crate) fn primitive(&self) -> Self {
        let mut c = self.content();
        if self.lc().is_negative() {
            c = -c;
        }
        self.map(|e, x| (e.clone(), x / &c))
    }

    pub(crate) fn mul(&self, other: &Self) -> Self {
        let terms = self.terms.iter().flat_map(|(a, x)| {
            other
                .terms
                .iter()
                .map(move |(b, y)| (add_exps(a, b), x * y))
        });
        Self::new(self.n, terms)
    }

    /// Variable `i` of the result is variable `perm[i]` of `self`.
    pub(crate) fn permute(&self, perm: &[usize]) -> Self {
        self.map(|e, c| (perm.iter().map(|&i| e[i]).collect(), c.clone()))
    }

    pub(crate) fn div_exact(&self, d: &Self) -> Option<Self> {
        let (dm, dc) = &d.terms[0];
        let mut rem: BTreeMap<Exps, BigInt> = self.terms.iter().cloned().collect();
        let mut q = Vec::new();
        while let Some((m, c)) = rem.pop_last() {
            let e = m
                .iter()
                .zip(dm)
                .map(|(a, b)| a.checked_sub(*b))
                .collect::<Option<Exps>>()?;
            let (qc, r) = c.div_rem(dc);
            if !r.is_zero() {
                return None;
            }
            for (tm, tc) in &d.terms[1..] {
                let key = add_exps(&e, tm);
                let v = rem.entry(key.clone()).or_default();
                *v -= &qc * tc;
                if v.is_zero() {
                    rem.remove(&key);
                }
            }
            q.push((e, qc));
        }
        Some(Self {
            n: self.n,
            terms: q,
        })
    }

    /// Coefficients of `x_0^i` as polynomials in the remaining variables.
    pub(crate) fn coefficients_in_x0(&self) -> Vec<Self> {
        let mut groups: BTreeMap<u32, Vec<(Exps, BigInt)>> = BTreeMap::new();
        for (e, c) in &self.terms {
            let mut k = e.clone();
            k[0] = 0;
            groups.entry(e[0]).or_default().push((k, c.clone()));
        }
        groups
            .into_values()
            .map(|t| Self {
                n: self.n,
                terms: t,
            })
            .collect()
    }
}

pub(crate) fn add_exps(a: &[u32], b: &[u32]) -> Exps {
    a.iter().zip(b).map(|(x, y)| x + y).collect()
}

pub(crate) fn one(n: usize) -> IntPoly {
    IntPoly::constant(n, BigInt::one())
}
