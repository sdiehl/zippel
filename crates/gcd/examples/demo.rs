use polycore::{Order, Ring};

fn main() {
    let ring = Ring::new(["x", "y", "z"], Order::Lex);
    let common = ring.parse("x^3*y^2 - 7*x*y*z + 3*z^5 - 2").unwrap();
    let left = ring.parse("x^2 + y*z^4 - 1").unwrap();
    let right = ring.parse("y^3 - x*z + 5").unwrap();
    let f = &common * &left;
    let g = &common * &right;

    let (h, cf, cg) = zippel_gcd::cofactors(&f, &g);
    let show = |p| ring.show(p);
    println!("f        = {}", show(&f));
    println!("g        = {}", show(&g));
    println!("gcd      = {}", show(&h));
    println!("f / gcd  = {}", show(&cf));
    println!("g / gcd  = {}", show(&cg));
    println!("lcm      = {}", show(&zippel_gcd::lcm(&f, &g)));
}
