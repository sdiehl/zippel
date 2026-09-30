//! Sparse integer polynomials in the same lex layout as `ModPoly`.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use std::collections::binary_heap::PeekMut;
use std::collections::{BTreeMap, BinaryHeap};
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
        self.terms.iter().fold(BigInt::zero(), |g, t| gcd(&g, &t.1))
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

    pub(crate) fn degrees(&self) -> Exps {
        let mut d = vec![0; self.n];
        for (e, _) in &self.terms {
            d.iter_mut().zip(e).for_each(|(d, &x)| *d = (*d).max(x));
        }
        d
    }

    fn sum(&self) -> BigInt {
        self.terms.iter().map(|t| &t.1).sum()
    }

    /// `self / d` when `d` divides exactly. Cheap necessary conditions go first, then heap
    /// division on packed exponents.
    pub(crate) fn div_exact(&self, d: &Self) -> Option<Self> {
        if self.is_zero() {
            return Some(self.clone());
        }
        let (fd, dd) = (self.degrees(), d.degrees());
        let (fl, dl) = (
            &self.terms[self.terms.len() - 1],
            &d.terms[d.terms.len() - 1],
        );
        let ds = d.sum();
        if dd.iter().zip(&fd).any(|(a, b)| a > b)
            || fl.0.iter().zip(&dl.0).any(|(a, b)| a < b)
            || !self.lc().is_multiple_of(d.lc())
            || !fl.1.is_multiple_of(&dl.1)
            || (!ds.is_zero() && !self.sum().is_multiple_of(&ds))
        {
            return None;
        }
        let Some(pk) = Packing::new(&fd) else {
            return self.div_map(d);
        };
        let bound = pk.pack(&fd.iter().zip(&dd).map(|(a, b)| a - b).collect::<Vec<_>>());
        let dp: Vec<u64> = d.terms.iter().map(|t| pk.pack(&t.0)).collect();
        let mut heap: BinaryHeap<(u64, usize, usize)> = BinaryHeap::new();
        let mut q: Vec<(u64, BigInt)> = Vec::new();
        let mut f = self.terms.iter().map(|(e, c)| (pk.pack(e), c)).peekable();
        loop {
            let m = match (f.peek().map(|t| t.0), heap.peek().map(|t| t.0)) {
                (None, None) => break,
                (a, b) => a.max(b).expect("one is some"),
            };
            let mut c = f
                .next_if(|t| t.0 == m)
                .map_or_else(BigInt::zero, |t| t.1.clone());
            while let Some((_, i, j)) = heap.peek_mut().filter(|t| t.0 == m).map(PeekMut::pop) {
                c -= &q[i].1 * &d.terms[j].1;
                if j + 1 < dp.len() {
                    heap.push((q[i].0 + dp[j + 1], i, j + 1));
                }
            }
            if c.is_zero() {
                continue;
            }
            let e = pk.sub(m, dp[0]).filter(|&e| pk.sub(bound, e).is_some())?;
            let (qc, r) = c.div_rem(d.lc());
            if !r.is_zero() {
                return None;
            }
            if dp.len() > 1 {
                heap.push((e + dp[1], q.len(), 1));
            }
            q.push((e, qc));
        }
        let terms = q.into_iter().map(|(e, c)| (pk.unpack(e), c)).collect();
        Some(Self { n: self.n, terms })
    }

    fn div_map(&self, d: &Self) -> Option<Self> {
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

/// Nonnegative integer gcd by Lehmer's algorithm.
pub(crate) fn gcd(a: &BigInt, b: &BigInt) -> BigInt {
    polycore::lehmer::gcd(a.magnitude(), b.magnitude()).into()
}

/// Exponents in one word, variable 0 highest, each field topped by a guard bit that makes
/// componentwise subtraction checkable. Lex order is integer order.
struct Packing {
    shifts: Vec<u32>,
    guard: u64,
}

impl Packing {
    fn new(degrees: &[u32]) -> Option<Self> {
        let mut shifts = vec![0; degrees.len()];
        let (mut at, mut guard) = (0u32, 0u64);
        for (s, d) in shifts.iter_mut().zip(degrees).rev() {
            let w = 33 - d.leading_zeros();
            *s = at;
            at += w;
            if at > 64 {
                return None;
            }
            guard |= 1 << (at - 1);
        }
        Some(Self { shifts, guard })
    }

    fn pack(&self, e: &[u32]) -> u64 {
        e.iter()
            .zip(&self.shifts)
            .map(|(&x, s)| u64::from(x) << s)
            .sum()
    }

    fn unpack(&self, x: u64) -> Exps {
        let mut hi = 64;
        self.shifts
            .iter()
            .map(|&s| {
                let v = (x & (u64::MAX >> (64 - hi))) >> s;
                hi = s;
                v as u32
            })
            .collect()
    }

    /// `a - b` componentwise, or `None` if some component would go negative.
    fn sub(&self, a: u64, b: u64) -> Option<u64> {
        let r = (a | self.guard) - b;
        (r & self.guard == self.guard).then_some(r & !self.guard)
    }
}

pub(crate) fn add_exps(a: &[u32], b: &[u32]) -> Exps {
    a.iter().zip(b).map(|(x, y)| x + y).collect()
}

pub(crate) fn one(n: usize) -> IntPoly {
    IntPoly::constant(n, BigInt::one())
}
