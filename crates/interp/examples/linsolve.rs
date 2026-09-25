//! Reconstruct a coefficient of a linear solve whose entries depend on the dimension `d`, the
//! invariant `s` and the mass `m`, the way IBP reduction coefficients are rational functions of
//! kinematics recovered from numeric solves.
#![allow(clippy::many_single_char_names)]

use std::sync::atomic::{AtomicUsize, Ordering};
use zippel_interp::modp::{add, inv, mul, sub};
use zippel_interp::{reconstruct, Primes};

const NAMES: [&str; 3] = ["d", "s", "m"];

/// `x_0` of `A x = e_0` by elimination without pivoting; a vanishing pivot is a bad point.
fn solve(v: &[u64], p: u64) -> Option<u64> {
    let [d, s, m] = [v[0], v[1], v[2]];
    let mut a = [
        [sub(d, 4, p), s, 1, 1],
        [m, sub(d, 3, p), s, 0],
        [add(s, m, p), m, sub(d, 5, p), 0],
    ];
    for i in 0..3 {
        let pinv = inv(Some(a[i][i]).filter(|&x| x != 0)?, p);
        for r in 0..3 {
            if r != i {
                let factor = mul(a[r][i], pinv, p);
                let pivot_row = a[i];
                for (x, y) in a[r].iter_mut().zip(pivot_row) {
                    *x = sub(*x, mul(factor, y, p), p);
                }
            }
        }
    }
    Some(mul(a[0][3], inv(a[0][0], p), p))
}

fn main() {
    let p = Primes::new().next().unwrap();
    let calls = AtomicUsize::new(0);
    let black_box = |x: &[u64], p| {
        calls.fetch_add(1, Ordering::Relaxed);
        solve(x, p)
    };
    let r = reconstruct(&black_box, 3, p, 1).unwrap();
    println!(
        "x0 = ({}) / ({})",
        r.num.show(&NAMES, p),
        r.den.show(&NAMES, p)
    );
    println!(
        "from {} numeric solves mod {p}",
        calls.load(Ordering::Relaxed)
    );
}
