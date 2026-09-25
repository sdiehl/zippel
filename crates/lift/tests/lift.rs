#![allow(clippy::many_single_char_names)]

use groebner::{MonomialOrder, PolynomialRing};
use num_integer::Integer;
use num_rational::BigRational;
use zippel_interp::modp::{add, inv, mul, pow};
use zippel_lift::{lift, Poly};

type Ring = PolynomialRing<BigRational>;

fn eval_mod(f: &Poly, x: &[u64], p: u64) -> Option<u64> {
    f.terms.iter().try_fold(0, |acc, t| {
        let c = &t.coefficient;
        let den = (c.denom() % p).try_into().ok().filter(|&d: &u64| d != 0)?;
        let num: u64 = (c.numer().mod_floor(&p.into())).try_into().ok()?;
        let v = t
            .monomial
            .exponents()
            .iter()
            .zip(x)
            .fold(mul(num, inv(den, p), p), |v, (&e, &xi)| {
                mul(v, pow(xi, u64::from(e), p), p)
            });
        Some(add(acc, v, p))
    })
}

#[test]
fn lifts_to_rationals() {
    let ring = Ring::new(["x", "y", "z"], MonomialOrder::Lex).unwrap();
    let cases = [
        ("x^2 - 3*y*z + 12345678901234567890123456789", "2*x*y - 5"),
        ("7/3*y - 1/2", "x + z^3"),
        ("0", "x - y"),
        ("x*y*z", "1/31415926535897932384626433 + x^2"),
    ];
    let fs: Vec<(Poly, Poly)> = cases
        .iter()
        .map(|(n, d)| (ring.parse(n).unwrap(), ring.parse(d).unwrap()))
        .collect();
    let bb = |x: &[u64], p| {
        fs.iter()
            .map(|(n, d)| {
                let d = eval_mod(d, x, p).filter(|&d| d != 0)?;
                Some(mul(eval_mod(n, x, p)?, inv(d, p), p))
            })
            .collect()
    };
    let out: Vec<String> = lift(bb, 3, 1)
        .unwrap()
        .iter()
        .map(|f| {
            format!(
                "({}) / ({})",
                ring.format(&f.num).unwrap(),
                ring.format(&f.den).unwrap()
            )
        })
        .collect();
    insta::assert_snapshot!(out.join("\n"));
}
