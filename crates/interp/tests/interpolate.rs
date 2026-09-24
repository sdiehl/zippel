#![allow(clippy::cast_possible_truncation)]

use std::cell::Cell;
use zippel_interp::{interpolate, BlackBox, ModPoly, Primes, Rng};

fn random_poly(rng: &mut Rng, n: usize, terms: usize, deg: u64, p: u64) -> ModPoly {
    let mut t: Vec<_> = (0..terms)
        .map(|_| {
            (
                (0..n)
                    .map(|_| (rng.nonzero(deg + 2) - 1) as u32)
                    .collect::<Vec<_>>(),
                rng.nonzero(p),
            )
        })
        .collect();
    t.sort_by(|a, b| b.0.cmp(&a.0));
    t.dedup_by(|a, b| a.0 == b.0);
    ModPoly { n, terms: t }
}

#[test]
fn recovers_random_sparse_polynomials() {
    let p = Primes::new().next().unwrap();
    let mut rng = Rng::new(1);
    for n in 0..7 {
        for terms in [0, 1, 5, 40] {
            let f = random_poly(&mut rng, n, terms, 12, p);
            let bb = |x: &[u64], p| Some(f.eval(x, p));
            assert_eq!(interpolate(&bb, n, p, &mut rng), Some(f));
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
    assert_eq!(interpolate(&bb, 8, p, &mut rng).as_ref(), Some(&f));
    let zippel_bound = 8 * (60 + 2) * (f.terms.len() + 1);
    assert!(calls.get() <= zippel_bound, "{} evaluations", calls.get());
}

/// Points where the black box fails are skipped rather than trusted.
#[test]
fn tolerates_failing_evaluations() {
    let p = Primes::new().next().unwrap();
    let mut rng = Rng::new(3);
    let f = random_poly(&mut rng, 4, 20, 6, p);
    let hash = |x: &[u64]| {
        x.iter()
            .fold(17u64, |h, &v| (h ^ v).wrapping_mul(0x100_0000_01b3))
    };
    let bb = |x: &[u64], p| (hash(x) % 7 != 0).then(|| f.eval(x, p));
    assert_eq!(interpolate(&bb, 4, p, &mut rng), Some(f));
}

/// A rational function is not a polynomial, so the budget runs out.
#[test]
fn rejects_non_polynomials() {
    let p = Primes::new().next().unwrap();
    let mut rng = Rng::new(4);
    let bb = |x: &[u64], p| Some(zippel_interp::modp::inv(x[0], p));
    assert!(bb.eval(&[2], p).is_some());
    assert_eq!(interpolate(&bb, 1, p, &mut rng), None);
}
