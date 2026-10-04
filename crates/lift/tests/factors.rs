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

#[test]
fn discovers_dimension_factors_without_a_supplied_pool() {
    let ring = Ring::new(["d", "s"], Order::Lex);
    let num: Poly = ring.parse("d+s+1").unwrap();
    let den: Poly = ring.parse("(2*d-7)^2*(s^2+1)").unwrap();
    let f = |x: &[u64], p| {
        let x = x.iter().map(|&v| Fp::new(v, p)).collect::<Vec<_>>();
        let eval = |g: &Poly| {
            g.map(|c| Fp::new(crt::reduce(c, p).unwrap(), p))
                .eval(&x)
                .residue_mod(p)
        };
        Some(vec![mul(eval(&num), try_inv(eval(&den), p)?, p)])
    };
    assert_eq!(
        zippel_lift::lift_with_discovered_factors(f, 2, &[], 23),
        lift(f, 2, 23)
    );
}

#[test]
fn discovery_skips_an_unusable_first_prime() {
    let first = polycore::modp::Primes::new().next().unwrap();
    let f = |x: &[u64], p| (p != first).then(|| vec![polycore::modp::add(x[0], 1, p)]);
    assert_eq!(
        zippel_lift::lift_with_discovered_factors(f, 1, &[], 29),
        lift(f, 1, 29)
    );
}
