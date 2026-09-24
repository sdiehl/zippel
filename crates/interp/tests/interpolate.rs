#![allow(clippy::cast_possible_truncation)]

use std::cell::Cell;
mod common;

use common::{random_poly, unlucky};
use zippel_interp::{interpolate, BlackBox, Primes, Rng};

#[test]
fn recovers_random_sparse_polynomials() {
    let p = Primes::new().next().unwrap();
    let mut rng = Rng::new(1);
    for n in 0..7 {
        for terms in [0, 1, 5, 40] {
            let f = random_poly(&mut rng, n, terms, 12, p);
            let bb = |x: &[u64], p| Some(f.eval(x, p));
            assert_eq!(interpolate(&bb, n, p, rng.nonzero(p)), Some(f));
        }
    }
}

#[test]
fn high_degree_and_many_variables() {
    let p = Primes::new().next().unwrap();
    let mut rng = Rng::new(2);
    let f = random_poly(&mut rng, 8, 40, 60, p);
    let calls = Cell::new(0);
    let bb = |x: &[u64], p| {
        calls.set(calls.get() + 1);
        Some(f.eval(x, p))
    };
    assert_eq!(interpolate(&bb, 8, p, 2).as_ref(), Some(&f));
    let zippel_bound = 8 * (60 + 2) * (f.terms.len() + 1);
    assert!(calls.get() <= zippel_bound, "{} evaluations", calls.get());
}

/// Points where the black box fails are skipped rather than trusted.
#[test]
fn tolerates_failing_evaluations() {
    let p = Primes::new().next().unwrap();
    let mut rng = Rng::new(3);
    let f = random_poly(&mut rng, 4, 20, 6, p);
    let bb = |x: &[u64], p| (!unlucky(x)).then(|| f.eval(x, p));
    assert_eq!(interpolate(&bb, 4, p, 3), Some(f));
}

/// A rational function is not a polynomial, so the budget runs out.
#[test]
fn rejects_non_polynomials() {
    let p = Primes::new().next().unwrap();
    let bb = |x: &[u64], p| Some(zippel_interp::modp::inv(x[0], p));
    assert!(bb.eval(&[2], p).is_some());
    assert_eq!(interpolate(&bb, 1, p, 4), None);
}
