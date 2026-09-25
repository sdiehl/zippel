#[path = "../examples/ibp/mod.rs"]
mod ibp;

use zippel_interp::{Primes, Rng};
use zippel_laporta::{Plan, Row};

const TARGETS: [ibp::Index; 5] = [(2, 1), (2, 2), (3, 1), (3, 0), (1, -1)];

#[test]
fn reduces_the_bubble() {
    let bubble = ibp::Bubble::new(3);
    let (plan, coefficients) = bubble.reduce(&TARGETS, Primes::new().next().unwrap());
    insta::assert_snapshot!(bubble.render(&plan, &coefficients));
}

/// The trimmed replay agrees with eliminating everything from scratch at each new point.
#[test]
fn replay_matches_full_elimination() {
    let bubble = ibp::Bubble::new(3);
    let targets: Vec<usize> = TARGETS.iter().map(|&t| bubble.column(t)).collect();
    let mut primes = Primes::new();
    let mut rng = Rng::new(9);
    let rows = |x: &[u64], p| {
        (0..bubble.len())
            .map(|e| bubble.row(e, x, p))
            .collect::<Vec<Row>>()
    };
    let mut point = |p| (0..3).map(|_| rng.nonzero(p)).collect::<Vec<u64>>();
    let p = primes.next().unwrap();
    let plan = Plan::learn(&rows(&point(p), p), &targets, p);
    for p in primes.take(3) {
        let x = point(p);
        let fresh = Plan::learn(&rows(&x, p), &targets, p);
        assert_eq!(fresh.masters, plan.masters);
        let all = rows(&x, p);
        assert_eq!(
            plan.replay(|e| bubble.row(e, &x, p), p),
            fresh.replay(|e| all[e].clone(), p)
        );
    }
}
