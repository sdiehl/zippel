//! Rational functions from a black box, after Klappert and Lange (2019).
//!
//! Along the ray `x = t*z + s`, with `z_0 = 1` and a fixed random shift `s`, the black box is a
//! univariate rational function of `t`. Thiele recovers it, scaled so the denominator is 1 at
//! `t = 0`, which is `D(s)` and therefore the same scale on every ray. The coefficient of `t^r` is
//! then the degree `r` homogeneous part of `N` (or `D`) at `z`, plus spill from the higher parts
//! through the shift. Recovering parts from the top down, that spill is already known, so each
//! part is a sparse polynomial for Zippel. Every part asks the same points, so each ray is
//! reconstructed once.

use std::collections::HashMap;
use std::sync::Mutex;

use crate::poly::{dense, Dense, ModPoly};
use crate::recovery::{RecoveryError, RecoveryOptions, RecoveryReport, Tracked};
use crate::zippel::interpolate_impl;
use polycore::interp::Thiele;
use polycore::modp::{add, inv, mul, sub};
use polycore::sample::{hash, point, BlackBox};
use polycore::{Fp, Modular};

const MAX_POINTS: u64 = 1 << 12;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RatFunc {
    pub num: ModPoly,
    pub den: ModPoly,
}

type Ray = Option<(Dense, Dense)>;

/// The rational function in `n` variables behind `f`, in lowest terms with the leading
/// coefficient of the denominator 1. Monte Carlo, like [`crate::interpolate`].
pub fn reconstruct(f: &(impl BlackBox + Sync), n: usize, p: u64, seed: u64) -> Option<RatFunc> {
    reconstruct_impl(f, n, p, seed, &RecoveryOptions::default(), &mut 0)
}

/// Reconstruct with configurable limits, counting actual oracle calls across all rays.
///
/// # Errors
/// Returns the failure category and work counters when reconstruction does not succeed.
pub fn reconstruct_with_options(
    f: &(impl BlackBox + Sync),
    n: usize,
    p: u64,
    seed: u64,
    options: &RecoveryOptions,
) -> Result<(RatFunc, RecoveryReport), RecoveryError> {
    let tracked = Tracked::new(f, options.max_probes);
    let mut attempts = 0;
    let value = reconstruct_impl(&tracked, n, p, seed, options, &mut attempts);
    tracked.finish(value, attempts)
}

fn reconstruct_impl(
    f: &(impl BlackBox + Sync),
    n: usize,
    p: u64,
    seed: u64,
    options: &RecoveryOptions,
    attempts: &mut usize,
) -> Option<RatFunc> {
    (0..options.max_attempts).find_map(|attempt| {
        *attempts += 1;
        let seed = hash(&[seed, attempt as u64]);
        let r = if n == 0 {
            let c = f.eval(&[], p)?;
            RatFunc {
                num: constant(0, c),
                den: constant(0, 1),
            }
        } else {
            Shifted::new(f, n, p, seed, *options).reconstruct()?
        };
        let x: Vec<u64> = (0..n as u64).map(|i| point(&[seed, 4, i], p)).collect();
        let y = f.eval(&x, p)?;
        (mul(y, r.den.eval(&x, p), p) == r.num.eval(&x, p)).then_some(r)
    })
}

struct Shifted<'a, F> {
    f: &'a F,
    n: usize,
    p: u64,
    seed: u64,
    s: Vec<u64>,
    rays: Mutex<HashMap<Vec<u64>, Ray>>,
    options: RecoveryOptions,
}

impl<F: BlackBox + Sync> Shifted<'_, F> {
    fn new(f: &F, n: usize, p: u64, seed: u64, options: RecoveryOptions) -> Shifted<'_, F> {
        let s = (0..n as u64).map(|i| point(&[seed, 0, i], p)).collect();
        Shifted {
            f,
            n,
            p,
            seed,
            s,
            rays: Mutex::default(),
            options,
        }
    }

    fn reconstruct(&self) -> Option<RatFunc> {
        let y: Vec<u64> = (1..self.n as u64)
            .map(|i| point(&[self.seed, 1, i], self.p))
            .collect();
        let first = self.ray(&y)?;
        let mut r = RatFunc {
            num: self.side(&first, |ray| &ray.0)?,
            den: self.side(&first, |ray| &ray.1)?,
        };
        let l = inv(r.den.terms[0].1, self.p);
        r.num.scale(l, self.p);
        r.den.scale(l, self.p);
        Some(r)
    }

    /// `f(t*z + s)` for `z = (1, y)` as `num / den` with `den(0) = 1`, once per ray. Rays
    /// asked together run on separate threads.
    fn ray(&self, y: &[u64]) -> Ray {
        if let Some(r) = self.rays.lock().unwrap().get(y) {
            return r.clone();
        }
        let z = ray_dir(y);
        let mut x = vec![0; self.n];
        let p = self.p;
        let ray = thiele_with_limit(
            |t| {
                for ((xi, &zi), &si) in x.iter_mut().zip(&z).zip(&self.s) {
                    *xi = add(mul(t, zi, p), si, p);
                }
                self.f.eval(&x, p)
            },
            p,
            hash(&[self.seed, 2]),
            self.options.max_points,
        )
        .and_then(|(num, den)| {
            let l = den.0.first()?.try_inverse()?;
            Some((num.scale(&l), den.scale(&l)))
        });
        self.rays.lock().unwrap().insert(y.to_vec(), ray.clone());
        ray
    }

    /// Numerator or denominator, one homogeneous part at a time from the top degree down.
    fn side(&self, first: &(Dense, Dense), pick: fn(&(Dense, Dense)) -> &Dense) -> Option<ModPoly> {
        let shape = (first.0 .0.len(), first.1 .0.len());
        let mut known = constant(self.n, 0);
        for r in (0..pick(first).0.len()).rev() {
            let part = |y: &[u64], p| {
                let ray = self
                    .ray(y)
                    .filter(|(a, b)| (a.0.len(), b.0.len()) == shape)?;
                let spill = self.spill(&known, &ray_dir(y));
                let at = |d: &Dense| d.0.get(r).map_or(0, |c| c.residue_mod(p));
                Some(sub(at(pick(&ray)), at(&spill), p))
            };
            let h = interpolate_impl(
                &part,
                self.n - 1,
                self.p,
                hash(&[self.seed, 3]),
                &self.options,
                &mut 0,
            )?;
            for (e, c) in h.terms {
                let d = r.checked_sub(e.iter().sum::<u32>() as usize)?;
                known
                    .terms
                    .push((std::iter::once(d as u32).chain(e).collect(), c));
            }
        }
        known.terms.sort_by(|a, b| b.0.cmp(&a.0));
        Some(known)
    }

    /// `g(t*z + s)` as a polynomial in `t`.
    fn spill(&self, g: &ModPoly, z: &[u64]) -> Dense {
        let p = self.p;
        let mut out = Dense::zero();
        for (e, c) in &g.terms {
            let mut term = dense([*c], p);
            for ((&d, &zi), &si) in e.iter().zip(z).zip(&self.s) {
                let line = dense([si, zi], p);
                for _ in 0..d {
                    term = &term * &line;
                }
            }
            out = &out + &term;
        }
        out
    }
}

/// The rational function behind `g` by Thiele interpolation, from points named by `seed`. `g`
/// may refuse points.
pub(crate) fn thiele(
    g: impl FnMut(u64) -> Option<u64>,
    p: u64,
    seed: u64,
) -> Option<(Dense, Dense)> {
    thiele_with_limit(g, p, seed, MAX_POINTS as usize)
}

fn thiele_with_limit(
    mut g: impl FnMut(u64) -> Option<u64>,
    p: u64,
    seed: u64,
    max_points: usize,
) -> Option<(Dense, Dense)> {
    let mut th = Thiele::default();
    for i in 0..max_points as u64 {
        let t = point(&[seed, i], p);
        let ft = Fp::new(t, p);
        if th.contains(&ft) {
            continue;
        }
        let Some(y) = g(t).map(|y| Fp::new(y, p)) else {
            continue;
        };
        if th.eval(&ft) == Some(y) {
            return Some(th.rational());
        }
        th.add(ft, y);
    }
    None
}

fn constant(n: usize, c: u64) -> ModPoly {
    ModPoly {
        n,
        terms: if c == 0 {
            vec![]
        } else {
            vec![(vec![0; n], c)]
        },
    }
}

fn ray_dir(y: &[u64]) -> Vec<u64> {
    std::iter::once(1).chain(y.iter().copied()).collect()
}
