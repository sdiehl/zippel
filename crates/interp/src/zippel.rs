//! Zippel's sparse interpolation of a black box polynomial, one variable at a time.
//!
//! Stage `k` knows `f(x_0..x_{k-1}, a_k..)` at a random anchor `a` and lifts it to `x_k`. For each
//! fresh value `b` of `x_k` the known support is the skeleton, so evaluating `x_0..x_{k-1}` at the
//! powers `r^j` leaves a transposed Vandermonde system for the coefficients. Newton interpolation
//! in `x_k` stops at the first value it already predicts, so no degree bounds are needed.
//!
//! Every point is named by `(seed, stage, index)` rather than drawn from a stream, so two
//! interpolations with one seed ask the same questions for as long as their shapes agree.

use crate::modp::{add, hash, mul, point, pow};
use crate::poly::{Exps, ModPoly};
use crate::univariate::Newton;
use crate::vandermonde;

/// A polynomial known only through evaluation modulo `p`. `None` marks a point where the
/// evaluation broke down (a vanishing pivot, say) and should be avoided.
pub trait BlackBox {
    fn eval(&self, x: &[u64], p: u64) -> Option<u64>;
}

impl<F: Fn(&[u64], u64) -> Option<u64>> BlackBox for F {
    fn eval(&self, x: &[u64], p: u64) -> Option<u64> {
        self(x, p)
    }
}

/// Points tried per variable before giving up on a degree.
const MAX_POINTS: usize = 1 << 12;

/// The polynomial in `n` variables behind `f`, or `None` if it looks like no polynomial of
/// degree below `MAX_POINTS` in each variable. Monte Carlo, with error probability about
/// `deg / p` per step.
pub fn interpolate(f: &impl BlackBox, n: usize, p: u64, seed: u64) -> Option<ModPoly> {
    (0..3).find_map(|attempt| {
        let seed = hash(&[seed, attempt]);
        let anchor: Vec<u64> = (0..n as u64).map(|i| point(&[seed, 0, i], p)).collect();
        let mut h = first(f, &anchor, p, seed)?;
        for k in 1..n {
            h = lift(f, &h, k, &anchor, p, seed)?;
        }
        Some(h)
    })
}

/// Dense interpolation in `x_0` with every other variable at the anchor.
fn first(f: &impl BlackBox, anchor: &[u64], p: u64, seed: u64) -> Option<ModPoly> {
    let mut x = anchor.to_vec();
    if x.is_empty() {
        return Some(ModPoly {
            n: 0,
            terms: f
                .eval(&x, p)
                .filter(|&c| c != 0)
                .map(|c| (vec![], c))
                .into_iter()
                .collect(),
        });
    }
    let mut newton = Newton::default();
    for i in 0..MAX_POINTS as u64 {
        x[0] = point(&[seed, 1, i], p);
        let Some(y) = f.eval(&x, p).filter(|_| !newton.contains(x[0])) else {
            continue;
        };
        if !newton.add(x[0], y, p) && newton.len() > 1 {
            let d = newton.poly(p);
            let mut terms = Vec::new();
            for (i, &c) in d.iter().enumerate().rev().filter(|t| *t.1 != 0) {
                let mut e = vec![0; anchor.len()];
                e[0] = i as u32;
                terms.push((e, c));
            }
            return Some(ModPoly {
                n: anchor.len(),
                terms,
            });
        }
    }
    None
}

/// Lift `h = f(x_0..x_{k-1}, a_k..)` to `f(x_0..x_k, a_{k+1}..)`.
fn lift(
    f: &impl BlackBox,
    h: &ModPoly,
    k: usize,
    anchor: &[u64],
    p: u64,
    seed: u64,
) -> Option<ModPoly> {
    let skeleton = h.support();
    let t = skeleton.len();
    if t == 0 {
        return Some(h.clone());
    }
    let (r, vals) = distinct_values(&skeleton, k, p, seed)?;
    let master = vandermonde::master(&vals, p);
    let mut newton: Vec<Newton> = vec![Newton::default(); t];
    for (nw, (_, c)) in newton.iter_mut().zip(&h.terms) {
        nw.add(anchor[k], *c, p);
    }
    let mut x = anchor.to_vec();
    for i in 0..MAX_POINTS as u64 {
        x[k] = point(&[seed, 2, k as u64, i], p);
        if newton[0].contains(x[k]) {
            continue;
        }
        let mut w = Vec::with_capacity(t + 1);
        x[..k].fill(1);
        for _ in 0..=t {
            x[..k]
                .iter_mut()
                .zip(&r)
                .for_each(|(xi, &ri)| *xi = mul(*xi, ri, p));
            match f.eval(&x, p) {
                Some(y) => w.push(y),
                None => break,
            }
        }
        if w.len() <= t {
            continue;
        }
        let c = vandermonde::solve(&vals, &master, &w[..t], p);
        let e = t as u64 + 1;
        let check = c
            .iter()
            .zip(&vals)
            .fold(0, |acc, (&ci, &v)| add(acc, mul(ci, pow(v, e, p), p), p));
        if check != w[t] {
            return None;
        }
        let mut changed = false;
        for (nw, ci) in newton.iter_mut().zip(c) {
            changed |= nw.add(x[k], ci, p);
        }
        if !changed {
            return Some(assemble(h.n, k, &skeleton, &newton, p));
        }
    }
    None
}

/// Random `r` at which the skeleton's monomials take distinct values, so the system is solvable.
fn distinct_values(skeleton: &[Exps], k: usize, p: u64, seed: u64) -> Option<(Vec<u64>, Vec<u64>)> {
    (0..3).find_map(|a| {
        let r: Vec<u64> = (0..k as u64)
            .map(|i| point(&[seed, 3, k as u64, a, i], p))
            .collect();
        let vals: Vec<u64> = skeleton
            .iter()
            .map(|e| {
                r.iter()
                    .zip(e)
                    .fold(1, |acc, (&ri, &d)| mul(acc, pow(ri, u64::from(d), p), p))
            })
            .collect();
        let mut sorted = vals.clone();
        sorted.sort_unstable();
        sorted.dedup();
        (sorted.len() == vals.len()).then_some((r, vals))
    })
}

fn assemble(n: usize, k: usize, skeleton: &[Exps], newton: &[Newton], p: u64) -> ModPoly {
    let mut terms = Vec::new();
    for (e, nw) in skeleton.iter().zip(newton) {
        let d = nw.poly(p);
        for (i, &c) in d.iter().enumerate().rev().filter(|t| *t.1 != 0) {
            let mut e = e.clone();
            e[k] = i as u32;
            terms.push((e, c));
        }
    }
    ModPoly { n, terms }
}
