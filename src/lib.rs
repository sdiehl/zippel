//! Sparse multivariate polynomial GCD by Zippel's modular interpolation.
#![allow(
    clippy::multiple_crate_versions,
    clippy::redundant_pub_crate,
    clippy::cast_possible_truncation,
    clippy::missing_panics_doc,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::module_name_repetitions
)]

mod gcd;
mod linzip;
mod modp;
mod pgcd;
mod poly;
mod univariate;

use groebner::{Monomial, Polynomial, Term};
use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Signed};
use poly::IntPoly;

type Poly = Polynomial<BigRational>;

/// `f = F / d` with `F` integral.
fn integral(f: &Poly) -> (IntPoly, BigInt) {
    let d = f
        .terms
        .iter()
        .fold(BigInt::one(), |d, t| d.lcm(t.coefficient.denom()));
    let terms = f.terms.iter().map(|t| {
        let c = &t.coefficient;
        (
            t.monomial.exponents().to_vec(),
            c.numer() * (&d / c.denom()),
        )
    });
    (IntPoly::new(f.nvars, terms), d)
}

fn rational(h: &IntPoly, d: &BigInt, like: &Poly) -> Poly {
    let terms = h.terms.iter().map(|(e, c)| {
        Term::new(
            BigRational::new(c.clone(), d.clone()),
            Monomial::new(e.clone()),
        )
    });
    Polynomial::new(terms.collect(), like.nvars, like.order.clone())
}

/// Primitive over Z with a positive leading coefficient in `like`'s monomial order.
fn normalized(h: &IntPoly, like: &Poly) -> Poly {
    let r = rational(
        &if h.is_zero() {
            h.clone()
        } else {
            h.primitive()
        },
        &BigInt::one(),
        like,
    );
    match r.leading_coefficient() {
        Some(c) if c.is_negative() => r.multiply_scalar(&-BigRational::one()),
        _ => r,
    }
}

fn gcd_int(f: &Poly, g: &Poly) -> IntPoly {
    assert_eq!(f.nvars, g.nvars, "polynomials must share a ring");
    gcd::gcd_z(&integral(f).0, &integral(g).0)
}

/// The greatest common divisor, as a primitive integer polynomial with positive leading
/// coefficient. `gcd(0, 0) = 0`.
#[must_use]
pub fn gcd(f: &Poly, g: &Poly) -> Poly {
    normalized(&gcd_int(f, g), f)
}

/// `(h, f / h, g / h)` where `h = gcd(f, g)`.
#[must_use]
pub fn cofactors(f: &Poly, g: &Poly) -> (Poly, Poly, Poly) {
    let h = gcd_int(f, g);
    let zero = Polynomial::zero(f.nvars, f.order.clone());
    if h.is_zero() {
        return (zero.clone(), zero.clone(), zero);
    }
    let h = normalized(&h, f);
    let hi = integral(&h).0;
    let (hf, df) = integral(f);
    let (hg, dg) = integral(g);
    let qf = hf.div_exact(&hi).expect("gcd divides f");
    let qg = hg.div_exact(&hi).expect("gcd divides g");
    (h, rational(&qf, &df, f), rational(&qg, &dg, g))
}

/// The least common multiple, normalized like [`gcd`]. `lcm(f, 0) = 0`.
#[must_use]
pub fn lcm(f: &Poly, g: &Poly) -> Poly {
    let h = gcd_int(f, g);
    if f.is_zero() || g.is_zero() {
        return Polynomial::zero(f.nvars, f.order.clone());
    }
    let (hg, _) = integral(g);
    normalized(
        &integral(f).0.mul(&hg.div_exact(&h).expect("gcd divides g")),
        f,
    )
}
