//! Dense polynomials over GF(p), coefficient `i` of `x^i`, no trailing zeros.

use crate::modp::{add, inv, mul, sub};

pub(crate) type Dense = Vec<u64>;

pub(crate) fn trim(a: &mut Dense) {
    while a.last() == Some(&0) {
        a.pop();
    }
}

pub(crate) const fn deg(a: &[u64]) -> usize {
    a.len().saturating_sub(1)
}

pub(crate) fn eval(a: &[u64], x: u64, p: u64) -> u64 {
    a.iter().rev().fold(0, |acc, &c| add(mul(acc, x, p), c, p))
}

pub(crate) fn mul_poly(a: &[u64], b: &[u64], p: u64) -> Dense {
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

pub(crate) fn div_rem(a: &[u64], b: &[u64], p: u64) -> (Dense, Dense) {
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

pub(crate) fn monic(mut a: Dense, p: u64) -> Dense {
    if let Some(&l) = a.last() {
        let l = inv(l, p);
        for c in &mut a {
            *c = mul(*c, l, p);
        }
    }
    a
}

pub(crate) fn gcd(a: &[u64], b: &[u64], p: u64) -> Dense {
    let (mut a, mut b) = (a.to_vec(), b.to_vec());
    while !b.is_empty() {
        let r = div_rem(&a, &b, p).1;
        a = b;
        b = r;
    }
    monic(a, p)
}

/// Newton interpolation through `(xs[i], ys[i])`.
pub(crate) fn interpolate(xs: &[u64], ys: &[u64], p: u64) -> Dense {
    let n = xs.len();
    let mut c = ys.to_vec();
    for j in 1..n {
        for i in (j..n).rev() {
            c[i] = mul(sub(c[i], c[i - 1], p), inv(sub(xs[i], xs[i - j], p), p), p);
        }
    }
    let mut r = vec![c[n - 1]];
    for i in (0..n - 1).rev() {
        r = mul_poly(&r, &[sub(0, xs[i], p), 1], p);
        r[0] = add(r[0], c[i], p);
    }
    trim(&mut r);
    r
}
