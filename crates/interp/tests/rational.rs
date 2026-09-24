#![allow(clippy::many_single_char_names)]

mod common;

use common::{random_poly, unlucky};
use zippel_interp::modp::{inv, mul};
use zippel_interp::{reconstruct, ModPoly, Primes, RatFunc, Rng};

/// `num / den` as a black box that fails at poles, normalized as `reconstruct` returns it.
fn quotient(
    mut num: ModPoly,
    mut den: ModPoly,
    p: u64,
) -> (impl Fn(&[u64], u64) -> Option<u64>, RatFunc) {
    let l = inv(den.terms[0].1, p);
    num.scale(l, p);
    den.scale(l, p);
    let r = RatFunc { num, den };
    let f = r.clone();
    let bb = move |x: &[u64], p| {
        let d = f.den.eval(x, p);
        (d != 0).then(|| mul(f.num.eval(x, p), inv(d, p), p))
    };
    (bb, r)
}

/// A constant term keeps the random pair coprime: no monomial divides the denominator.
fn random_den(rng: &mut Rng, n: usize, terms: usize, deg: u64, p: u64) -> ModPoly {
    let mut d = random_poly(rng, n, terms, deg, p);
    if d.terms.last().is_none_or(|t| t.0.iter().any(|&e| e > 0)) {
        d.terms.push((vec![0; n], rng.nonzero(p)));
    }
    d
}

#[test]
fn recovers_random_rational_functions() {
    let p = Primes::new().next().unwrap();
    let mut rng = Rng::new(5);
    for n in 0..5 {
        for (nt, dt) in [(0, 0), (1, 0), (4, 3), (12, 6), (30, 10)] {
            let num = random_poly(&mut rng, n, nt, 7, p);
            let den = random_den(&mut rng, n, dt, 5, p);
            let (bb, r) = quotient(num, den, p);
            assert_eq!(reconstruct(&bb, n, p, rng.nonzero(p)), Some(r), "n = {n}");
        }
    }
}

/// `(x y - z^3) / (x^2 + y z)`: no constant term, so `s` is what keeps `t = 0` off the pole, and
/// a black box that also fails at scattered points.
#[test]
fn shift_avoids_the_origin() {
    let p = Primes::new().next().unwrap();
    let m = p - 1;
    let num = ModPoly {
        n: 3,
        terms: vec![(vec![1, 1, 0], 1), (vec![0, 0, 3], m)],
    };
    let den = ModPoly {
        n: 3,
        terms: vec![(vec![2, 0, 0], 1), (vec![0, 1, 1], 1)],
    };
    let (f, r) = quotient(num, den, p);
    let bb = |x: &[u64], p| f(x, p).filter(|_| !unlucky(x));
    assert_eq!(reconstruct(&bb, 3, p, 7), Some(r));
}
