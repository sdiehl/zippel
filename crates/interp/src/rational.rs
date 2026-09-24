//! Rational functions from a black box, after Klappert and Lange's `FireFly` (2019).
//!
//! Along the ray `x = t*z + s`, with `z_0 = 1` and a fixed random shift `s`, the black box is a
//! univariate rational function of `t`. Thiele recovers it, scaled so the denominator is 1 at
//! `t = 0`, which is `D(s)` and therefore the same scale on every ray. The coefficient of `t^r` is
//! then the degree `r` homogeneous part of `N` (or `D`) at `z`, plus spill from the higher parts
//! through the shift. Recovering parts from the top down, that spill is already known, so each
//! part is a sparse polynomial for Zippel. Every part asks the same points, so each ray is
//! reconstructed once.

use std::cell::RefCell;
use std::collections::HashMap;

use crate::modp::{add, hash, inv, mul, point, sub};
use crate::poly::ModPoly;
use crate::thiele;
use crate::univariate::{self as uni, Dense};
use crate::zippel::{interpolate, BlackBox};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RatFunc {
    pub num: ModPoly,
    pub den: ModPoly,
}

type Ray = Option<(Dense, Dense)>;

/// The rational function in `n` variables behind `f`, in lowest terms with the leading
/// coefficient of the denominator 1. Monte Carlo, like [`interpolate`].
pub fn reconstruct(f: &impl BlackBox, n: usize, p: u64, seed: u64) -> Option<RatFunc> {
    (0..3).find_map(|attempt| {
        let seed = hash(&[seed, attempt]);
        let r = if n == 0 {
            let c = f.eval(&[], p)?;
            RatFunc {
                num: constant(0, c),
                den: constant(0, 1),
            }
        } else {
            Shifted::new(f, n, p, seed).reconstruct()?
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
    rays: RefCell<HashMap<Vec<u64>, Ray>>,
}

impl<F: BlackBox> Shifted<'_, F> {
    fn new(f: &F, n: usize, p: u64, seed: u64) -> Shifted<'_, F> {
        let s = (0..n as u64).map(|i| point(&[seed, 0, i], p)).collect();
        Shifted {
            f,
            n,
            p,
            seed,
            s,
            rays: RefCell::default(),
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

    /// `f(t*z + s)` for `z = (1, y)` as `num / den` with `den(0) = 1`, once per ray.
    fn ray(&self, y: &[u64]) -> Ray {
        if let Some(r) = self.rays.borrow().get(y) {
            return r.clone();
        }
        let z = ray_dir(y);
        let mut x = vec![0; self.n];
        let p = self.p;
        let ray = thiele::reconstruct(
            |t| {
                for ((xi, &zi), &si) in x.iter_mut().zip(&z).zip(&self.s) {
                    *xi = add(mul(t, zi, p), si, p);
                }
                self.f.eval(&x, p)
            },
            p,
            hash(&[self.seed, 2]),
        )
        .and_then(|(mut num, mut den)| {
            let l = inv(*den.first().filter(|&&d| d != 0)?, p);
            num.iter_mut()
                .chain(&mut den)
                .for_each(|c| *c = mul(*c, l, p));
            Some((num, den))
        });
        self.rays.borrow_mut().insert(y.to_vec(), ray.clone());
        ray
    }

    /// Numerator or denominator, one homogeneous part at a time from the top degree down.
    fn side(&self, first: &(Dense, Dense), pick: fn(&(Dense, Dense)) -> &Dense) -> Option<ModPoly> {
        let shape = (first.0.len(), first.1.len());
        let mut known = constant(self.n, 0);
        for r in (0..pick(first).len()).rev() {
            let part = |y: &[u64], p| {
                let ray = self.ray(y).filter(|(a, b)| (a.len(), b.len()) == shape)?;
                let spill = self.spill(&known, &ray_dir(y));
                Some(sub(pick(&ray)[r], spill.get(r).copied().unwrap_or(0), p))
            };
            let h = interpolate(&part, self.n - 1, self.p, hash(&[self.seed, 3]))?;
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
        let mut out = Dense::new();
        for (e, c) in &g.terms {
            let mut term = vec![*c];
            for ((&d, &zi), &si) in e.iter().zip(z).zip(&self.s) {
                for _ in 0..d {
                    term = uni::mul_poly(&term, &[si, zi], p);
                }
            }
            out.resize(out.len().max(term.len()), 0);
            for (a, b) in out.iter_mut().zip(term) {
                *a = add(*a, b, p);
            }
        }
        out
    }
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
