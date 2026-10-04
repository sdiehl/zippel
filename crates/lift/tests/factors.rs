use polycore::modp::{mul, try_inv};
use polycore::{crt, Fp, Modular, Order, Ring};
use std::sync::atomic::{AtomicUsize, Ordering};
use zippel_lift::{lift, lift_with_factors, Poly};

#[test]
fn factors_survive_lifting_and_reduce_total_probes() {
    let ring = Ring::new(["x", "y", "z"], Order::Lex);
    let num: Poly = ring.parse("x^3+y^2+z+1").unwrap();
    let den: Poly = ring.parse("(x+y+z+1)^4*(x+2*y+3*z+2)^3").unwrap();
    let candidates: Vec<Poly> = ["x+y+z+1", "x+2*y+3*z+2", "x+19"]
        .iter()
        .map(|s| ring.parse(s).unwrap())
        .collect();
    let calls = AtomicUsize::new(0);
    let f = |x: &[u64], p| {
        calls.fetch_add(1, Ordering::Relaxed);
        let x: Vec<_> = x.iter().map(|&v| Fp::new(v, p)).collect();
        let eval = |f: &Poly| {
            f.map(|c| Fp::new(crt::reduce(c, p).unwrap(), p))
                .eval(&x)
                .residue_mod(p)
        };
        let v = mul(eval(&num), try_inv(eval(&den), p)?, p);
        Some(vec![v, mul(v, 2, p), 0])
    };
    let want = lift(f, 3, 11).unwrap();
    let before = calls.swap(0, Ordering::Relaxed);
    assert_eq!(lift_with_factors(f, 3, &candidates, 11), Some(want));
    let after = calls.load(Ordering::Relaxed);
    eprintln!("all-prime probes: {before} -> {after}");
    assert!(after < before, "{before} -> {after}");
}
