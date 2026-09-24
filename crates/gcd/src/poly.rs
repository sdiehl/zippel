//! Sparse integer polynomials in the same lex layout as `ModPoly`.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use std::collections::BTreeMap;
pub(crate) use zippel_interp::poly::{Exps, ModPoly};

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct IntPoly {
    pub(crate) n: usize,
    pub(crate) terms: Vec<(Exps, BigInt)>,
}

impl IntPoly {
    pub(crate) fn reduce(&self, p: u64) -> ModPoly {
        let m = BigInt::from(p);
        let terms = self
            .terms
            .iter()
            .filter_map(|(e, c)| {
                let r = u64::try_from(c.mod_floor(&m)).unwrap_or(0);
                (r != 0).then(|| (e.clone(), r))
            })
            .collect();
        ModPoly { n: self.n, terms }
    }

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
