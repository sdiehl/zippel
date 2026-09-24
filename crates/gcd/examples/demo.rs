use groebner::{MonomialOrder, PolynomialRing};
use num_rational::BigRational;

fn main() {
    let ring = PolynomialRing::<BigRational>::new(["x", "y", "z"], MonomialOrder::Lex).unwrap();
    let common = ring.parse("x^3*y^2 - 7*x*y*z + 3*z^5 - 2").unwrap();
    let left = ring.parse("x^2 + y*z^4 - 1").unwrap();
    let right = ring.parse("y^3 - x*z + 5").unwrap();
    let f = common.multiply(&left);
    let g = common.multiply(&right);

    let (h, cf, cg) = zippel_gcd::cofactors(&f, &g);
    let show = |p| ring.format(p).unwrap();
    println!("f        = {}", show(&f));
    println!("g        = {}", show(&g));
    println!("gcd      = {}", show(&h));
    println!("f / gcd  = {}", show(&cf));
    println!("g / gcd  = {}", show(&cg));
    println!("lcm      = {}", show(&zippel_gcd::lcm(&f, &g)));
}
