use polycore::{Order, Ring};
use zippel_gcd::{
    gcd, gcd_with_algorithm, hu_monagan, hu_monagan_bivariate, huang_gao, huang_monagan, zippel,
    GcdAlgorithm,
};

#[test]
fn named_algorithms_recover_the_same_normalized_gcd() {
    let ring = Ring::new(["x", "y", "z"], Order::GRevLex);
    let h = ring.parse("(2*x+y+1)^2").unwrap();
    let a = &h * &ring.parse("-2/7*x^2*y*(x+z+2)").unwrap();
    let b = &h * &ring.parse("3/11*x*y^2*(y+z+3)").unwrap();
    let expected = (&h * &ring.parse("x*y").unwrap()).primitive();
    assert_eq!(gcd(&a, &b), expected);
    assert_eq!(zippel(&a, &b), expected);
    for algorithm in [huang_gao, huang_monagan, hu_monagan, hu_monagan_bivariate] {
        assert_eq!(algorithm(&a, &b).unwrap(), expected);
        assert_eq!(algorithm(&b, &a).unwrap(), expected);
    }
}

#[test]
fn direct_recovery_declines_while_aggregate_apis_fall_back() {
    let ring = Ring::new(["x"], Order::Lex);
    let f = ring.parse("x^65537+1").unwrap();
    for (direct, selected) in [
        (huang_gao as fn(_, _) -> _, GcdAlgorithm::HuangGao),
        (huang_monagan, GcdAlgorithm::HuangMonagan),
        (hu_monagan, GcdAlgorithm::HuMonagan),
        (hu_monagan_bivariate, GcdAlgorithm::HuMonaganBivariate),
    ] {
        assert!(direct(&f, &f).is_none(), "{selected:?} must not fall back");
        assert_eq!(gcd_with_algorithm(&f, &f, selected), f);
    }
    assert_eq!(gcd(&f, &f), f);
    assert_eq!(zippel(&f, &f), f);
}

#[test]
fn named_algorithms_preserve_zero_and_constant_conventions() {
    let ring = Ring::new(["x", "y"], Order::Lex);
    let zero = ring.parse("0").unwrap();
    let f = ring.parse("-6*x+12*y").unwrap();
    let constant = ring.parse("15").unwrap();
    for algorithm in [huang_gao, huang_monagan, hu_monagan, hu_monagan_bivariate] {
        assert_eq!(algorithm(&zero, &zero), Some(zero.clone()));
        assert_eq!(algorithm(&zero, &f), Some(f.primitive()));
        assert_eq!(algorithm(&f, &zero), Some(f.primitive()));
        assert_eq!(algorithm(&constant, &f), Some(ring.parse("1").unwrap()));
    }
    assert_eq!(zippel(&zero, &zero), zero);
    assert_eq!(zippel(&zero, &f), f.primitive());
    assert_eq!(zippel(&f, &zero), f.primitive());
    assert_eq!(zippel(&constant, &f), ring.parse("1").unwrap());
}
