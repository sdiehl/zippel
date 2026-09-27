use num_rational::BigRational;
use num_traits::{One, Zero};
use polycore::{Order, Ring};

type Poly = polycore::Poly<BigRational>;
use zippel_gcd::Frac;

fn frac(ring: &Ring, num: &str, den: &str) -> Frac {
    Frac::new(ring.parse(num).unwrap(), ring.parse(den).unwrap()).unwrap()
}

fn q(n: i64, d: i64) -> BigRational {
    BigRational::new(n.into(), d.into())
}

#[test]
fn lowest_terms() {
    let r = Ring::new(["a", "b"], Order::GRevLex);
    let f = frac(&r, "a^2 - b^2", "2*a - 2*b");
    assert_eq!(r.show(f.num()), "1/2*a + 1/2*b");
    assert_eq!(r.show(f.den()), "1");
    assert_eq!(f, frac(&r, "a + b", "2"));
}

#[test]
fn field_operations() {
    let r = Ring::new(["a", "b"], Order::GRevLex);
    let x = frac(&r, "a", "a + b");
    let y = frac(&r, "b - 1", "a*b + 3");
    let at = [q(3, 2), q(-5, 7)];
    let v = |f: &Frac| f.eval(&at).unwrap();
    assert_eq!(v(&(x.clone() + y.clone())), v(&x) + v(&y));
    assert_eq!(v(&(x.clone() * y.clone())), v(&x) * v(&y));
    assert_eq!(v(&(x.clone() / y.clone())), v(&x) / v(&y));
    assert_eq!(x.clone() / x.clone(), Frac::one());
    assert_eq!(
        frac(&r, "1", "a + b") + frac(&r, "1", "a - b"),
        frac(&r, "2*a", "a^2 - b^2")
    );
    assert!((y.clone() - y).is_zero());
    assert_eq!(Frac::zero() + x.clone(), x);
    assert!(Frac::new(r.parse("1").unwrap(), Poly::zero(2, Order::GRevLex)).is_none());
}
