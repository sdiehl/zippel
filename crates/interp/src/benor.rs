//! Ben-Or and Tiwari's interpolation from about `2t` evaluations for `t` terms, with no degree
//! bounds.
//!
//! At `x_i = s_i q_i^j`, with `q_i` the `i`-th prime, a term `c m` contributes `c m(s) m(q)^j`, so
//! the values form a sum of geometric sequences. Berlekamp-Massey finds the polynomial whose roots
//! are the ratios `m(q)`, and factoring each ratio over the `q_i` reads off the exponents. That
//! only works while `m(q) < p`, so it suits very sparse polynomials of modest total degree, where
//! Zippel would spend `t` points per degree of every variable.

use crate::modp::{add, hash, inv, mul, point, pow, sub};
use crate::poly::{Exps, ModPoly};
use crate::univariate::{self as uni, Dense};
use crate::vandermonde;
use crate::zippel::BlackBox;

const MAX_TERMS: usize = 1 << 12;

const BATCH: usize = 8;

/// Consecutive values the recurrence must predict before it is trusted.
const MARGIN: usize = 2;

/// The polynomial in `n` variables behind `f`, or `None` if its terms are too many or of too
/// high a degree. Monte Carlo, like [`crate::interpolate`].
pub fn interpolate(f: &impl BlackBox, n: usize, p: u64, seed: u64) -> Option<ModPoly> {
    (0..3).find_map(|attempt| attempt_once(f, n, p, hash(&[seed, attempt])))
}

fn attempt_once(f: &impl BlackBox, n: usize, p: u64, seed: u64) -> Option<ModPoly> {
    let q = small_primes(n);
    let s: Vec<u64> = (0..n as u64).map(|i| point(&[seed, 0, i], p)).collect();
    let mut bm = Massey::new();
    let mut a = Vec::new();
    while !bm.settled(a.len()) {
        if a.len() > 2 * MAX_TERMS {
            return None;
        }
        let xs: Vec<Vec<u64>> = (a.len()..a.len() + BATCH)
            .map(|j| {
                s.iter()
                    .zip(&q)
                    .map(|(&s, &q)| mul(s, pow(q, j as u64, p), p))
                    .collect()
            })
            .collect();
        for y in f.eval_many(&xs, p) {
            a.push(y?);
            bm.push(&a, p);
        }
    }
    let lambda = bm.generator();
    let z = uni::div_rem(&[0, 1], &lambda, p).1;
    if powmod(&z, p, &lambda, p) != z {
        return None;
    }
    let ratios = roots(&lambda, p, seed);
    let exps: Vec<Exps> = ratios
        .iter()
        .map(|&r| exponents(r, &q))
        .collect::<Option<_>>()?;
    let scaled = vandermonde::solve(&ratios, &lambda, &a[1..=ratios.len()], p);
    let mut terms: Vec<(Exps, u64)> = exps
        .into_iter()
        .zip(scaled)
        .map(|(e, c)| {
            let shift = e
                .iter()
                .zip(&s)
                .fold(1, |m, (&k, &s)| mul(m, pow(s, k.into(), p), p));
            (e, mul(c, inv(shift, p), p))
        })
        .filter(|t| t.1 != 0)
        .collect();
    terms.sort_by(|a, b| b.0.cmp(&a.0));
    let g = ModPoly { n, terms };
    (0..8)
        .find_map(|i| {
            let x: Vec<u64> = (0..n as u64).map(|k| point(&[seed, 1, i, k], p)).collect();
            f.eval(&x, p).map(|y| y == g.eval(&x, p))
        })?
        .then_some(g)
}

/// Berlekamp-Massey, one value at a time: the shortest recurrence `c` generating the sequence.
#[derive(Debug)]
struct Massey {
    c: Dense,
    b: Dense,
    l: usize,
    shift: usize,
    last: u64,
    zeros: usize,
}

impl Massey {
    fn new() -> Self {
        Self {
            c: vec![1],
            b: vec![1],
            l: 0,
            shift: 1,
            last: 1,
            zeros: 0,
        }
    }

    fn push(&mut self, a: &[u64], p: u64) {
        let k = a.len() - 1;
        let d = self
            .c
            .iter()
            .zip(a.iter().rev())
            .fold(0, |acc, (&c, &x)| add(acc, mul(c, x, p), p));
        if d == 0 {
            self.shift += 1;
            self.zeros += 1;
            return;
        }
        self.zeros = 0;
        let r = mul(d, inv(self.last, p), p);
        let mut c = self.c.clone();
        c.resize(c.len().max(self.b.len() + self.shift), 0);
        for (i, &x) in self.b.iter().enumerate() {
            c[i + self.shift] = sub(c[i + self.shift], mul(r, x, p), p);
        }
        if 2 * self.l <= k {
            self.l = k + 1 - self.l;
            self.b = std::mem::replace(&mut self.c, c);
            self.last = d;
            self.shift = 1;
        } else {
            self.c = c;
            self.shift += 1;
        }
    }

    const fn settled(&self, len: usize) -> bool {
        self.zeros >= MARGIN && len >= 2 * self.l + MARGIN
    }

    /// `prod (z - r)` over the ratios `r`: the recurrence reversed.
    fn generator(&self) -> Dense {
        let mut c = self.c.clone();
        c.resize(self.l + 1, 0);
        c.reverse();
        c
    }
}

fn powmod(a: &[u64], mut e: u64, f: &[u64], p: u64) -> Dense {
    let (mut r, mut a) = (uni::div_rem(&[1], f, p).1, a.to_vec());
    while e > 0 {
        if e & 1 == 1 {
            r = uni::div_rem(&uni::mul_poly(&r, &a, p), f, p).1;
        }
        a = uni::div_rem(&uni::mul_poly(&a, &a, p), f, p).1;
        e >>= 1;
    }
    r
}

/// The roots of a monic `f` with distinct roots, all in GF(p), by Cantor-Zassenhaus: a random
/// shift of `z^((p-1)/2) - 1` catches about half of them.
fn roots(f: &[u64], p: u64, seed: u64) -> Vec<u64> {
    if f.len() <= 2 {
        return f.get(1).map(|_| vec![sub(0, f[0], p)]).unwrap_or_default();
    }
    (0..)
        .find_map(|i| {
            let mut g = powmod(&[point(&[seed, 2, i], p), 1], (p - 1) / 2, f, p);
            g.resize(g.len().max(1), 0);
            g[0] = sub(g[0], 1, p);
            uni::trim(&mut g);
            let g = uni::gcd(&g, f, p);
            (g.len() > 1 && g.len() < f.len()).then(|| {
                let h = uni::div_rem(f, &g, p).0;
                [
                    roots(&g, p, hash(&[seed, i, 0])),
                    roots(&h, p, hash(&[seed, i, 1])),
                ]
                .concat()
            })
        })
        .unwrap()
}

fn exponents(mut r: u64, q: &[u64]) -> Option<Exps> {
    if r == 0 {
        return None;
    }
    let e = q
        .iter()
        .map(|&q| {
            let mut k = 0;
            while r.is_multiple_of(q) {
                r /= q;
                k += 1;
            }
            k
        })
        .collect();
    (r == 1).then_some(e)
}

fn small_primes(n: usize) -> Vec<u64> {
    (2..)
        .filter(|&k: &u64| (2..k).take_while(|d| d * d <= k).all(|d| k % d != 0))
        .take(n)
        .collect()
}
