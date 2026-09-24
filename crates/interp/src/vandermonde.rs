//! Transposed Vandermonde systems `sum_l c_l v_l^(j+1) = w_j`, the linear algebra of sparse
//! interpolation: a monomial with value `v` at a point contributes `v^j` at that point's `j`-th
//! power.

use crate::modp::{add, inv, mul, sub};
use crate::univariate::{self as uni, Dense};

/// `prod_l (z - v_l)`.
pub fn master(vals: &[u64], p: u64) -> Dense {
    vals.iter().fold(vec![1], |acc, &v| {
        uni::mul_poly(&acc, &[sub(0, v, p), 1], p)
    })
}

/// Solves for `c` given distinct nonzero `vals` and `w_0..w_{n-1}`, in quadratic time: the
/// cofactor `q = master / (z - v_l)` annihilates every column but the `l`-th.
pub fn solve(vals: &[u64], master: &[u64], w: &[u64], p: u64) -> Vec<u64> {
    vals.iter()
        .map(|&v| {
            let q = uni::div_rem(master, &[sub(0, v, p), 1], p).0;
            let num = q
                .iter()
                .zip(w)
                .fold(0, |acc, (&a, &b)| add(acc, mul(a, b, p), p));
            mul(num, inv(mul(uni::eval(&q, v, p), v, p), p), p)
        })
        .collect()
}
