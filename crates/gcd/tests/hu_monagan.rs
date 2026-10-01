use num_bigint::BigInt;
use num_traits::One;
use polycore::{Order, Ring};
use zippel_gcd::{hu_monagan, hu_monagan_bivariate};

#[test]
fn public_entry_points_preserve_content_and_rational_normalization() {
    let ring = Ring::new(["x", "y", "z", "w"], Order::GRevLex);
    let h = ring.parse("(y+z+1)^2*((y+w)*x^3+z*x+1)").unwrap();
    let a = &h * &ring.parse("-2/7*x^2*y*(x+z+2)").unwrap();
    let b = &h * &ring.parse("3/11*x*y^2*(x+w+3)").unwrap();
    let expected = (&h * &ring.parse("x*y").unwrap()).primitive();
    for gcd in [hu_monagan, hu_monagan_bivariate] {
        assert_eq!(gcd(&a, &b).unwrap(), expected);
        assert_eq!(gcd(&b, &a).unwrap(), expected);
    }
}

#[test]
fn reconstructs_large_coefficients_across_smooth_primes() {
    let ring = Ring::new(["x", "y", "z"], Order::Lex);
    let c = (BigInt::one() << 190) + BigInt::from(321);
    let h = ring.parse(&format!("{c}*x^3*y-17*y*z+31*z+1")).unwrap();
    let a = &h * &ring.parse("x+y+2").unwrap();
    let b = &h * &ring.parse("y+z+3").unwrap();
    for gcd in [hu_monagan, hu_monagan_bivariate] {
        assert_eq!(gcd(&a, &b).unwrap(), h.primitive());
    }
}

#[test]
fn retained_variables_cover_univariate_and_bivariate_rings() {
    let ring = Ring::new(["x"], Order::Lex);
    let a = ring.parse("(3*x+2)^3*(x+1)").unwrap();
    let b = ring.parse("(3*x+2)^2*(x+2)").unwrap();
    for gcd in [hu_monagan, hu_monagan_bivariate] {
        assert_eq!(gcd(&a, &b).unwrap(), ring.parse("(3*x+2)^2").unwrap());
    }
    let ring = Ring::new(["x", "y"], Order::Lex);
    let h = ring.parse("(y+1)*x^2+y*x+3").unwrap();
    let a = &h * &ring.parse("x+y").unwrap();
    let b = &h * &ring.parse("x+2*y+1").unwrap();
    assert_eq!(hu_monagan_bivariate(&a, &b).unwrap(), h);
}
