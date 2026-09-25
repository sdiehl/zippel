#![allow(clippy::cast_possible_truncation)]

use std::sync::atomic::{AtomicUsize, Ordering};
mod common;

use common::{random_poly, unlucky};
use zippel_interp::modp::pow;
use zippel_interp::{benor, interpolate, BlackBox, Primes, Rng};

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
    let calls = AtomicUsize::new(0);
    let bb = |x: &[u64], p| {
        calls.fetch_add(1, Ordering::Relaxed);
        Some(f.eval(x, p))
    };
    assert_eq!(interpolate(&bb, 8, p, 2).as_ref(), Some(&f));
    let zippel_bound = 8 * (60 + 2) * (f.terms.len() + 1);
    let calls = calls.into_inner();
    assert!(calls <= zippel_bound, "{calls} evaluations");
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

/// Ben-Or/Tiwari needs about two points per term, where Zippel pays per degree of each variable.
#[test]
fn ben_or_tiwari_on_very_sparse_polynomials() {
    let p = Primes::new().next().unwrap();
    let mut rng = Rng::new(5);
    for (n, terms, deg) in [(0, 1, 0), (1, 0, 3), (3, 6, 8), (6, 30, 4), (10, 12, 2)] {
        let f = random_poly(&mut rng, n, terms, deg, p);
        let calls = [AtomicUsize::new(0), AtomicUsize::new(0)];
        let bb = |k: usize| {
            let (f, calls) = (&f, &calls);
            move |x: &[u64], p| {
                calls[k].fetch_add(1, Ordering::Relaxed);
                Some(f.eval(x, p))
            }
        };
        assert_eq!(benor::interpolate(&bb(0), n, p, 5).as_ref(), Some(&f));
        assert_eq!(interpolate(&bb(1), n, p, 5).as_ref(), Some(&f));
        let [bt, zippel] = calls.map(AtomicUsize::into_inner);
        assert!(
            bt <= 2 * f.terms.len() + 16 && (n < 3 || bt < zippel),
            "{bt} vs {zippel}"
        );
    }
}

/// Monomials worth more than `p` at the primes no longer factor, so the result is refused.
#[test]
fn ben_or_tiwari_refuses_high_degree() {
    let p = Primes::new().next().unwrap();
    let bb = |x: &[u64], p| Some(pow(x[0], 70, p));
    assert_eq!(benor::interpolate(&bb, 1, p, 6), None);
}
