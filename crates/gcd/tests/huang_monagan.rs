use num_bigint::BigInt;
use num_traits::One;
use polycore::{Order, Ring};
use zippel_gcd::{cofactors_with_algorithm, gcd_with_algorithm, lcm_with_algorithm, GcdAlgorithm};

const ALGORITHMS: [GcdAlgorithm; 5] = [
    GcdAlgorithm::HuangGao,
    GcdAlgorithm::HuangMonagan,
    GcdAlgorithm::Zippel,
    GcdAlgorithm::HuMonagan,
    GcdAlgorithm::HuMonaganBivariate,
];

#[test]
fn selectable_backends_agree_on_gcd_cofactors_and_lcm() {
    let ring = Ring::new(["x", "y", "z"], Order::GRevLex);
    let h = ring.parse("(x + y + 1)^2*(2*y - 3*z)").unwrap();
    let a = &h * &ring.parse("-2/7*x^2*y*(x + z + 2)").unwrap();
    let b = &h * &ring.parse("3/11*x*y^2*(y + z + 3)").unwrap();
    let expected = (&h * &ring.parse("x*y").unwrap()).primitive();
    let expected_lcm = zippel_gcd::lcm(&a, &b);
    for algorithm in ALGORITHMS {
        let (g, ca, cb) = cofactors_with_algorithm(&a, &b, algorithm);
        assert_eq!(g, expected, "{algorithm:?}");
        assert_eq!(&g * &ca, a);
        assert_eq!(&g * &cb, b);
        assert_eq!(lcm_with_algorithm(&a, &b, algorithm), expected_lcm);
    }
}

#[test]
fn reconstructs_large_integer_coefficients_with_prime_substitution() {
    let ring = Ring::new(["x", "y", "z"], Order::Lex);
    let c = (BigInt::one() << 190) + BigInt::from(321);
    let h = ring
        .parse(&format!("{c}*x^3*y - 17*y*z + 31*z + 1"))
        .unwrap();
    let a = &h * &ring.parse("x + y + 2").unwrap();
    let b = &h * &ring.parse("y + z + 3").unwrap();
    assert_eq!(
        gcd_with_algorithm(&a, &b, GcdAlgorithm::HuangMonagan),
        h.primitive()
    );
}

#[test]
fn selectors_preserve_zero_and_constant_conventions() {
    let ring = Ring::new(["x"], Order::Lex);
    let zero = ring.parse("0").unwrap();
    let f = ring.parse("-6*x + 12").unwrap();
    let constant = ring.parse("15").unwrap();
    for algorithm in ALGORITHMS {
        assert_eq!(gcd_with_algorithm(&zero, &f, algorithm), f.primitive());
        assert_eq!(gcd_with_algorithm(&f, &zero, algorithm), f.primitive());
        assert_eq!(gcd_with_algorithm(&zero, &zero, algorithm), zero);
        assert_eq!(
            gcd_with_algorithm(&constant, &f, algorithm),
            ring.parse("1").unwrap()
        );
        assert_eq!(
            cofactors_with_algorithm(&zero, &zero, algorithm),
            (zero.clone(), zero.clone(), zero.clone())
        );
        assert_eq!(lcm_with_algorithm(&zero, &f, algorithm), zero);
    }
    assert_eq!(GcdAlgorithm::default(), GcdAlgorithm::HuangGao);
}
