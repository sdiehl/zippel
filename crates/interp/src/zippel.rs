//! Zippel's sparse interpolation of a black box polynomial, one variable at a time.
//!
//! Stage `k` knows `f(x_0..x_{k-1}, a_k..)` at a random anchor `a` and lifts it to `x_k`. For each
//! fresh value `b` of `x_k` the known support is the skeleton, so evaluating `x_0..x_{k-1}` at the
//! powers `r^j` leaves a transposed Vandermonde system for the coefficients. Newton interpolation
//! in `x_k` stops at the first value it already predicts, so no degree bounds are needed.
//!
//! Every point is named by `(seed, stage, index)` rather than drawn from a stream, so two
//! interpolations with one seed ask the same questions for as long as their shapes agree.

use crate::poly::{power_table, Exps, ModPoly};
use crate::recovery::{RecoveryError, RecoveryOptions, RecoveryReport, Tracked};
use polycore::interp::{master, solve, Newton};
use polycore::modp::{add, mul, pow};
use polycore::sample::{hash, point, BlackBox};
use polycore::{Fp, Modular};

/// Points evaluated together when the next one depends on the last, at the risk of a few extra.
const BATCH: u64 = 8;

/// The polynomial in `n` variables behind `f`, or `None` if it looks like no polynomial of
/// degree below 4096 in each variable. Monte Carlo, with error probability about
/// `deg / p` per step.
pub fn interpolate(f: &impl BlackBox, n: usize, p: u64, seed: u64) -> Option<ModPoly> {
    interpolate_impl(f, n, p, seed, &RecoveryOptions::default(), &mut 0)
}

/// Reconstruct with configurable point, attempt and oracle budgets.
///
/// # Errors
/// Returns the failure category and work counters when reconstruction does not succeed.
pub fn interpolate_with_options(
    f: &impl BlackBox,
    n: usize,
    p: u64,
    seed: u64,
    options: &RecoveryOptions,
) -> Result<(ModPoly, RecoveryReport), RecoveryError> {
    let tracked = Tracked::new(f, options.max_probes);
    let mut attempts = 0;
    let value = interpolate_impl(&tracked, n, p, seed, options, &mut attempts);
    tracked.finish(value, attempts)
}

pub(crate) fn interpolate_impl(
    f: &impl BlackBox,
    n: usize,
    p: u64,
    seed: u64,
    options: &RecoveryOptions,
    attempts: &mut usize,
) -> Option<ModPoly> {
    (0..options.max_attempts).find_map(|attempt| {
        *attempts += 1;
        let seed = hash(&[seed, attempt as u64]);
        let anchor: Vec<u64> = (0..n as u64).map(|i| point(&[seed, 0, i], p)).collect();
        let mut h = first(f, &anchor, p, seed, options.max_points)?;
        for k in 1..n {
            h = lift(f, &h, k, &anchor, p, seed, options.max_points)?;
        }
        Some(h)
    })
}

/// Dense interpolation in `x_0` with every other variable at the anchor.
fn first(
    f: &impl BlackBox,
    anchor: &[u64],
    p: u64,
    seed: u64,
    max_points: usize,
) -> Option<ModPoly> {
    if anchor.is_empty() {
        let c = f.eval(&[], p)?;
        return Some(ModPoly {
            n: 0,
            terms: if c == 0 { vec![] } else { vec![(vec![], c)] },
        });
    }
    let mut newton = Newton::default();
    let points = (0..max_points as u64)
        .step_by(BATCH as usize)
        .flat_map(|start| {
            let xs: Vec<Vec<u64>> = (start..start.saturating_add(BATCH).min(max_points as u64))
                .map(|i| {
                    let mut x = anchor.to_vec();
                    x[0] = point(&[seed, 1, i], p);
                    x
                })
                .collect();
            let ys = f.eval_many(&xs, p);
            xs.into_iter().zip(ys)
        });
    for (x, y) in points {
        let x0 = Fp::new(x[0], p);
        let Some(y) = y.filter(|_| !newton.contains(&x0)) else {
            continue;
        };
        if !newton.add(x0, Fp::new(y, p)) && newton.len() > 1 {
            let n = anchor.len();
            return Some(assemble(n, 0, &[vec![0; n]], &[newton], p));
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
    max_points: usize,
) -> Option<ModPoly> {
    let skeleton = h.support();
    let t = skeleton.len();
    if t == 0 {
        return Some(h.clone());
    }
    let (r, vals) = distinct_values(&skeleton, &h.degrees()[..k], p, seed)?;
    let vals: Vec<Fp> = vals.iter().map(|&v| Fp::new(v, p)).collect();
    let master = master(&vals);
    let mut newton: Vec<Newton<Fp>> = vec![Newton::default(); t];
    for (nw, (_, c)) in newton.iter_mut().zip(&h.terms) {
        nw.add(Fp::new(anchor[k], p), Fp::new(*c, p));
    }
    let mut x = anchor.to_vec();
    for i in 0..max_points as u64 {
        x[k] = point(&[seed, 2, k as u64, i], p);
        let xk = Fp::new(x[k], p);
        if newton[0].contains(&xk) {
            continue;
        }
        x[..k].fill(1);
        let xs: Vec<Vec<u64>> = (0..=t)
            .map(|_| {
                x[..k]
                    .iter_mut()
                    .zip(&r)
                    .for_each(|(xi, &ri)| *xi = mul(*xi, ri, p));
                x.clone()
            })
            .collect();
        let Some(w) = f
            .eval_many(&xs, p)
            .into_iter()
            .map(|y| y.map(|y| Fp::new(y, p)))
            .collect::<Option<Vec<Fp>>>()
        else {
            continue;
        };
        let c = solve(&vals, &master, &w[..t]);
        let e = t as u64 + 1;
        let check = c.iter().zip(&vals).fold(0, |acc, (ci, v)| {
            add(
                acc,
                mul(ci.residue_mod(p), pow(v.residue_mod(p), e, p), p),
                p,
            )
        });
        if check != w[t].residue_mod(p) {
            return None;
        }
        let mut changed = false;
        for (nw, ci) in newton.iter_mut().zip(c) {
            changed |= nw.add(xk, ci);
        }
        if !changed {
            return Some(assemble(h.n, k, &skeleton, &newton, p));
        }
    }
    None
}

/// Random `r` at which the skeleton's monomials take distinct values, so the system is solvable.
fn distinct_values(
    skeleton: &[Exps],
    degrees: &[usize],
    p: u64,
    seed: u64,
) -> Option<(Vec<u64>, Vec<u64>)> {
    let k = degrees.len() as u64;
    (0..3).find_map(|a| {
        let r: Vec<u64> = (0..k).map(|i| point(&[seed, 3, k, a, i], p)).collect();
        let table = power_table(&r, degrees, p);
        let vals: Vec<u64> = skeleton
            .iter()
            .map(|e| {
                table
                    .iter()
                    .zip(e)
                    .fold(1, |acc, (ri, &d)| mul(acc, ri[d as usize], p))
            })
            .collect();
        let mut sorted = vals.clone();
        sorted.sort_unstable();
        sorted.dedup();
        (sorted.len() == vals.len()).then_some((r, vals))
    })
}

fn assemble(n: usize, k: usize, skeleton: &[Exps], newton: &[Newton<Fp>], p: u64) -> ModPoly {
    let groups = skeleton
        .iter()
        .cloned()
        .zip(newton.iter().map(Newton::poly));
    ModPoly::from_groups(n, k, groups.collect(), p)
}
