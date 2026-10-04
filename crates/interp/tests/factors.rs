use polycore::modp::{mul, try_inv, Primes};
use polycore::{Fp, Order, Ring};
use std::sync::atomic::{AtomicUsize, Ordering};
use zippel_interp::factors::{guess, reconstruct_with_factors, univariate_factors};
use zippel_interp::{reconstruct, ModPoly};

fn parse(s: &str, p: u64) -> ModPoly {
    let ring = Ring::new(["x", "y", "z"], Order::Lex);
    let q = ring.parse(s).unwrap();
    ModPoly::from_poly(
        &q.map(|c| Fp::new(polycore::crt::reduce(c, p).unwrap(), p)),
        p,
    )
}

#[test]
fn candidates_repeated_missing_false_and_duplicate() {
    let p = Primes::new().next().unwrap();
    let num = parse("x^2+y+1", p);
    let den = parse("(x+y+1)^3*(y-2)*(z+3)", p);
    let f = |x: &[u64], p| Some(mul(num.eval(x, p), try_inv(den.eval(x, p), p)?, p));
    let candidates: Vec<_> = ["x+y+1", "y-2", "z+3", "x+y+1", "x+17"]
        .iter()
        .map(|s| parse(s, p))
        .collect();
    let g = guess(&f, &candidates, 3, p, 7).unwrap();
    assert_eq!(g.powers, vec![3, 1, 1, 0, 0]);
    assert!(g.complete);
    let want = reconstruct(&f, 3, p, 9).unwrap();
    assert_eq!(
        reconstruct_with_factors(&f, &candidates, 3, p, 7),
        Some(want.clone())
    );
    assert_eq!(
        reconstruct_with_factors(&f, &candidates[..1], 3, p, 7),
        Some(want.clone())
    );
    assert_eq!(
        reconstruct_with_factors(&f, &[parse("0", p)], 3, p, 7),
        Some(want)
    );
    let auto = univariate_factors(&f, 3, p, 7).unwrap();
    assert!(auto.contains(&candidates[1]));
    assert!(auto.contains(&candidates[2]));
    assert_eq!(auto.len(), 2);
}

#[test]
fn complete_pool_saves_probes_including_discovery() {
    let p = Primes::new().next().unwrap();
    let num = parse("x^3+y^2+z+1", p);
    let den = parse("(x+y+z+1)^4*(x+2*y+3*z+2)^3", p);
    let count = AtomicUsize::new(0);
    let f = |x: &[u64], p| {
        count.fetch_add(1, Ordering::Relaxed);
        Some(mul(num.eval(x, p), try_inv(den.eval(x, p), p)?, p))
    };
    let baseline = reconstruct(&f, 3, p, 3).unwrap();
    let before = count.swap(0, Ordering::Relaxed);
    let factors: Vec<_> = ["x+y+z+1", "x+2*y+3*z+2"]
        .iter()
        .map(|s| parse(s, p))
        .collect();
    assert_eq!(
        reconstruct_with_factors(&f, &factors, 3, p, 3),
        Some(baseline)
    );
    let after = count.load(Ordering::Relaxed);
    eprintln!("reconstruction probes: {before} -> {after}");
    assert!(after < before, "{before} -> {after}");
}
