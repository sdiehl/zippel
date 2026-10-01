use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use polycore::{Order, Ring};

#[test]
fn reconstructs_large_coefficients_contents_and_monomials() {
    let ring = Ring::new(["x", "y", "z"], Order::GRevLex);
    let huge = (BigInt::one() << 190) + BigInt::from(123_457);
    let h = ring
        .parse(&format!("{huge}*x^2*y + 17*y*z - 31*z + 1"))
        .unwrap();
    let monomial = ring.parse("x^3*y^2").unwrap();
    let expected = &h * &monomial;
    let a = (&expected * &ring.parse("y + z + 1").unwrap())
        .scale(&BigRational::new(BigInt::from(-6), BigInt::from(35)));
    let b = (&expected * &ring.parse("x*z + 2").unwrap())
        .scale(&BigRational::new(BigInt::from(10), BigInt::from(77)));
    let (g, ca, cb) = zippel_gcd::cofactors(&a, &b);
    assert_eq!(g, expected.primitive());
    assert_eq!(&g * &ca, a);
    assert_eq!(&g * &cb, b);
}

#[test]
fn preserves_content_in_variables_outside_the_main_variable() {
    let ring = Ring::new(["x", "y", "z"], Order::Lex);
    let h = ring.parse("(y + z + 1)*(2*x*y - 3*z)").unwrap();
    let a = &h * &ring.parse("x^4 + y + 3").unwrap();
    let b = &h * &ring.parse("x^3 + z + 5").unwrap();
    assert_eq!(zippel_gcd::gcd(&a, &b), h.primitive());
}

#[test]
fn zero_and_constant_rings_keep_existing_normalization() {
    let ring = Ring::new(["x", "y"], Order::Lex);
    let zero = ring.parse("0").unwrap();
    let f = ring.parse("-6*x^2 + 12*y").unwrap();
    assert_eq!(zippel_gcd::gcd(&zero, &f), f.primitive());
    assert_eq!(zippel_gcd::gcd(&f, &zero), f.primitive());
    assert_eq!(zippel_gcd::gcd(&zero, &zero), zero);
    let constants = Ring::new(std::iter::empty::<&str>(), Order::Lex);
    assert_eq!(
        zippel_gcd::gcd(
            &constants.parse("6").unwrap(),
            &constants.parse("15").unwrap()
        ),
        constants.parse("1").unwrap()
    );
}
