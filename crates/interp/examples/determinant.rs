//! Recover a symbolic determinant from numeric Gaussian elimination mod p, the same shape as
//! reconstructing IBP reduction coefficients from numeric linear solves.
#![allow(clippy::many_single_char_names)]

use std::cell::Cell;
use zippel_interp::modp::{inv, mul, sub};
use zippel_interp::{interpolate, ModPoly, Primes, Rng};

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

/// Coefficients in the symmetric range, which shows small integers as themselves.
fn show(f: &ModPoly, p: u64) -> String {
    let terms = f.terms.iter().map(|(e, c)| {
        let c = if *c > p / 2 {
            -i128::from(p - c)
        } else {
            i128::from(*c)
        };
        let vars: Vec<String> = e
            .iter()
            .zip(NAMES)
            .filter(|t| *t.0 > 0)
            .map(|(&d, v)| {
                if d == 1 {
                    v.to_string()
                } else {
                    format!("{v}^{d}")
                }
            })
            .collect();
        let vars = vars.join("*");
        let body = match (c.abs(), vars.is_empty()) {
            (1, false) => vars,
            (a, true) => a.to_string(),
            (a, false) => format!("{a}*{vars}"),
        };
        (c < 0, body)
    });
    terms
        .enumerate()
        .fold(String::new(), |out, (i, (neg, t))| match (i, neg) {
            (0, false) => t,
            (0, true) => format!("-{t}"),
            (_, true) => format!("{out} - {t}"),
            (_, false) => format!("{out} + {t}"),
        })
}

fn main() {
    let p = Primes::new().next().unwrap();
    let calls = Cell::new(0);
    let black_box = |x: &[u64], p| {
        calls.set(calls.get() + 1);
        determinant(x, p)
    };
    let det = interpolate(&black_box, 4, p, &mut Rng::new(1)).unwrap();
    println!("det = {}", show(&det, p));
    println!(
        "{} terms from {} numeric eliminations mod {p}",
        det.terms.len(),
        calls.get()
    );
}
