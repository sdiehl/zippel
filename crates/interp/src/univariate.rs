//! Dense polynomials over GF(p), coefficient `i` of `x^i`, no trailing zeros.

use crate::modp::{add, inv, mul, sub};

pub type Dense = Vec<u64>;

pub fn trim(a: &mut Dense) {
    while a.last() == Some(&0) {
        a.pop();
    }
}

pub const fn deg(a: &[u64]) -> usize {
    a.len().saturating_sub(1)
}

pub fn eval(a: &[u64], x: u64, p: u64) -> u64 {
    a.iter().rev().fold(0, |acc, &c| add(mul(acc, x, p), c, p))
}

pub fn mul_poly(a: &[u64], b: &[u64], p: u64) -> Dense {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut r = vec![0; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        for (j, &y) in b.iter().enumerate() {
            r[i + j] = add(r[i + j], mul(x, y, p), p);
        }
    }
    r
}

pub fn div_rem(a: &[u64], b: &[u64], p: u64) -> (Dense, Dense) {
    let mut r = a.to_vec();
    if r.len() < b.len() {
        return (Vec::new(), r);
    }
    let lead = inv(b[b.len() - 1], p);
    let mut q = vec![0; r.len() - b.len() + 1];
    for i in (0..q.len()).rev() {
        let c = mul(r[i + b.len() - 1], lead, p);
        q[i] = c;
        for (j, &y) in b.iter().enumerate() {
            r[i + j] = sub(r[i + j], mul(c, y, p), p);
        }
    }
    trim(&mut r);
    (q, r)
}

pub fn monic(mut a: Dense, p: u64) -> Dense {
    if let Some(&l) = a.last() {
        let l = inv(l, p);
        for c in &mut a {
            *c = mul(*c, l, p);
        }
    }
    a
}

pub fn gcd(a: &[u64], b: &[u64], p: u64) -> Dense {
    let (mut a, mut b) = (a.to_vec(), b.to_vec());
    while !b.is_empty() {
        let r = div_rem(&a, &b, p).1;
        a = b;
        b = r;
    }
    monic(a, p)
}

/// Newton interpolation through `(xs[i], ys[i])`.
pub fn interpolate(xs: &[u64], ys: &[u64], p: u64) -> Dense {
    let mut newton = Newton::default();
    xs.iter().zip(ys).for_each(|(&x, &y)| {
        newton.add(x, y, p);
    });
    newton.poly(p)
}

/// Newton's divided differences, one point at a time, at distinct abscissae.
#[derive(Clone, Debug, Default)]
pub struct Newton {
    xs: Vec<u64>,
    cs: Vec<u64>,
}

impl Newton {
    /// Adds `(x, y)`, returning `false` when the interpolant already passed through it.
    pub fn add(&mut self, x: u64, y: u64, p: u64) -> bool {
        let (mut v, mut w) = (0, 1);
        for (&xi, &c) in self.xs.iter().zip(&self.cs) {
            v = add(v, mul(c, w, p), p);
            w = mul(w, sub(x, xi, p), p);
        }
        let c = mul(sub(y, v, p), inv(w, p), p);
        self.xs.push(x);
        self.cs.push(c);
        c != 0
    }

    #[must_use]
    pub const fn len(&self) -> usize {
        self.xs.len()
    }

    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.xs.is_empty()
    }

    pub fn contains(&self, x: u64) -> bool {
        self.xs.contains(&x)
    }

    pub fn poly(&self, p: u64) -> Dense {
        let mut r: Dense = Vec::new();
        for (&xi, &c) in self.xs.iter().zip(&self.cs).rev() {
            r = mul_poly(&r, &[sub(0, xi, p), 1], p);
            match r.first_mut() {
                Some(r0) => *r0 = add(*r0, c, p),
                None => r.push(c),
            }
        }
        trim(&mut r);
        r
    }
}
