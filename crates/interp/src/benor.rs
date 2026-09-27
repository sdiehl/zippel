//! Ben-Or and Tiwari's interpolation from about `2t` evaluations for `t` terms, with no degree
//! bounds.
//!
//! At `x_i = s_i q_i^j`, with `q_i` the `i`-th prime, a term `c m` contributes `c m(s) m(q)^j`, so
//! the values form a sum of geometric sequences. Berlekamp-Massey finds the polynomial whose roots
//! are the ratios `m(q)`, and factoring each ratio over the `q_i` reads off the exponents. That
//! only works while `m(q) < p`, so it suits very sparse polynomials of modest total degree, where
//! Zippel would spend `t` points per degree of every variable.

use polycore::interp::{solve, Massey};
use polycore::modp::{inv, mul, pow};
use polycore::sample::{hash, point, BlackBox};
use polycore::{Fp, Modular};

use crate::poly::{Exps, ModPoly};

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
    let mut bm = Massey::default();
    let mut a = Vec::new();
    while !bm.settled(MARGIN) {
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
            bm.push(Fp::new(a[a.len() - 1], p));
        }
    }
    let lambda = bm.generator();
    let ratios = lambda.roots();
    if ratios.len() != lambda.deg() {
        return None;
    }
    let exps: Vec<Exps> = ratios
        .iter()
        .map(|r| exponents(r.residue_mod(p), &q))
        .collect::<Option<_>>()?;
    let w: Vec<Fp> = a[1..=ratios.len()].iter().map(|&v| Fp::new(v, p)).collect();
    let scaled = solve(&ratios, &lambda, &w);
    let mut terms: Vec<(Exps, u64)> = exps
        .into_iter()
        .zip(scaled)
        .map(|(e, c)| {
            let shift = e
                .iter()
                .zip(&s)
                .fold(1, |m, (&k, &s)| mul(m, pow(s, k.into(), p), p));
            (e, mul(c.residue_mod(p), inv(shift, p), p))
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
