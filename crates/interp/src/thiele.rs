//! Univariate rational interpolation by Thiele's continued fraction.
//!
//! The fraction `c_0 + (t - t_0) / (c_1 + (t - t_1) / (c_2 + ...))` grows one point at a time
//! until it predicts the next one, so no degree bounds are needed.

use crate::modp::{add, inv, mul, point, sub};
use crate::univariate::{self as uni, Dense};

const MAX_POINTS: u64 = 1 << 12;

#[derive(Clone, Debug, Default)]
pub struct Thiele {
    ts: Vec<u64>,
    cs: Vec<u64>,
}

impl Thiele {
    /// The fraction at `t`, or `None` at a pole (and before the first point).
    pub fn eval(&self, t: u64, p: u64) -> Option<u64> {
        let mut v = None;
        for (&ti, &ci) in self.ts.iter().zip(&self.cs).rev() {
            v = match v {
                None => Some(ci),
                Some(0) => None,
                Some(x) => Some(add(ci, mul(sub(t, ti, p), inv(x, p), p), p)),
            };
        }
        v
    }

    /// Adds `(t, y)` by inverted differences, returning `false` if they break down there.
    pub fn add(&mut self, t: u64, y: u64, p: u64) -> bool {
        let mut v = y;
        for (&ti, &ci) in self.ts.iter().zip(&self.cs) {
            if v == ci {
                return false;
            }
            v = mul(sub(t, ti, p), inv(sub(v, ci, p), p), p);
        }
        self.ts.push(t);
        self.cs.push(v);
        true
    }

    /// Collapses the fraction to `num / den` in lowest terms.
    pub fn rational(&self, p: u64) -> (Dense, Dense) {
        let mut num: Dense = vec![*self.cs.last().unwrap()];
        let mut den: Dense = vec![1];
        for (&ti, &ci) in self.ts.iter().zip(&self.cs).rev().skip(1) {
            let mut next = uni::mul_poly(&den, &[sub(0, ti, p), 1], p);
            next.resize(next.len().max(num.len()), 0);
            for (a, &b) in next.iter_mut().zip(&num) {
                *a = add(*a, mul(ci, b, p), p);
            }
            uni::trim(&mut next);
            den = std::mem::replace(&mut num, next);
        }
        let g = uni::gcd(&num, &den, p);
        (uni::div_rem(&num, &g, p).0, uni::div_rem(&den, &g, p).0)
    }
}

/// The rational function behind `g`, from points named by `seed`. `g` may refuse points.
pub fn reconstruct(
    mut g: impl FnMut(u64) -> Option<u64>,
    p: u64,
    seed: u64,
) -> Option<(Dense, Dense)> {
    let mut th = Thiele::default();
    for i in 0..MAX_POINTS {
        let t = point(&[seed, i], p);
        if th.ts.contains(&t) {
            continue;
        }
        let Some(y) = g(t) else { continue };
        if th.eval(t, p) == Some(y) {
            return Some(th.rational(p));
        }
        th.add(t, y, p);
    }
    None
}
