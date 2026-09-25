//! Recover a symbolic determinant from numeric Gaussian elimination mod p, the same shape as
//! reconstructing IBP reduction coefficients from numeric linear solves.
#![allow(clippy::many_single_char_names)]

use std::sync::atomic::{AtomicUsize, Ordering};
use zippel_interp::modp::{inv, mul, sub};
use zippel_interp::{interpolate, Primes};

const NAMES: [&str; 4] = ["x", "y", "z", "w"];

/// Plain elimination without pivoting: a vanishing pivot is a bad point, not a zero determinant.
fn determinant(v: &[u64], p: u64) -> Option<u64> {
    let [x, y, z, w] = [v[0], v[1], v[2], v[3]];
    let mut m = [[x, y, 0, 1], [z, x, y, 0], [0, z, x, y], [w, p - 1, z, x]];
    let mut det = 1;
    for i in 0..4 {
        let pivot = m[i][i];
        if pivot == 0 {
            return None;
        }
        det = mul(det, pivot, p);
        let pinv = inv(pivot, p);
        for r in i + 1..4 {
            let factor = mul(m[r][i], pinv, p);
            let pivot_row = m[i];
            for (a, b) in m[r][i..].iter_mut().zip(&pivot_row[i..]) {
                *a = sub(*a, mul(factor, *b, p), p);
            }
        }
    }
    Some(det)
}

fn main() {
    let p = Primes::new().next().unwrap();
    let calls = AtomicUsize::new(0);
    let black_box = |x: &[u64], p| {
        calls.fetch_add(1, Ordering::Relaxed);
        determinant(x, p)
    };
    let det = interpolate(&black_box, 4, p, 1).unwrap();
    println!("det = {}", det.show(&NAMES, p));
    println!(
        "{} terms from {} numeric eliminations mod {p}",
        det.terms.len(),
        calls.load(Ordering::Relaxed)
    );
}
