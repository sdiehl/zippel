#[path = "../examples/families/mod.rs"]
mod families;

use zippel_interp::{Primes, Rng};
use zippel_laporta::{ibp::Index, Plan, Row};

fn reduce(system: &zippel_laporta::ibp::System, targets: &[Index]) -> String {
    let (plan, coefficients) = system.reduce(targets, 1).unwrap();
    system.render(&plan, &coefficients)
}

#[test]
fn reduces_the_bubble() {
    let targets = [[2, 1], [2, 2], [3, 1], [3, 0], [1, -1]].map(|a| a.to_vec());
    insta::assert_snapshot!(reduce(&families::bubble().system(2, 1), &targets));
}

#[test]
fn reduces_the_box() {
    let targets = [[2, 1, 1, 1], [2, 2, 1, 1], [1, 1, 2, 1], [2, 1, 1, 0]].map(|a| a.to_vec());
    insta::assert_snapshot!(reduce(&families::one_loop_box().system(2, 2), &targets));
}

#[test]
fn reduces_the_double_box() {
    let targets = [
        [1, 1, 1, 1, 1, 1, 1, -1, 0],
        [1, 1, 1, 1, 1, 1, 1, 0, -1],
        [2, 1, 1, 1, 1, 1, 1, 0, 0],
    ]
    .map(|a| a.to_vec());
    insta::assert_snapshot!(reduce(&families::double_box().system(1, 1), &targets));
}

/// The trimmed replay agrees with eliminating everything from scratch at each new point.
#[test]
fn replay_matches_full_elimination() {
    let system = families::one_loop_box().system(2, 2);
    let targets: Vec<usize> = [[2, 1, 1, 1], [2, 2, 1, 1], [1, 1, 1, 0]]
        .iter()
        .map(|a| system.column(&a.to_vec()))
        .collect();
    let rows = |x: &[u64], p| {
        (0..system.len())
            .map(|e| system.row(e, x, p))
            .collect::<Vec<Row>>()
    };
    let mut primes = Primes::new();
    let mut rng = Rng::new(9);
    let n = system.vars.len();
    let mut point = |p| (0..n).map(|_| rng.nonzero(p)).collect::<Vec<u64>>();
    let p = primes.next().unwrap();
    let plan = Plan::learn(&rows(&point(p), p), &targets, p);
    for p in primes.take(3) {
        let x = point(p);
        let all = rows(&x, p);
        let fresh = Plan::learn(&all, &targets, p);
        assert_eq!(fresh.masters, plan.masters);
        assert_eq!(
            plan.replay(|e| system.row(e, &x, p), p),
            fresh.replay(|e| all[e].clone(), p)
        );
    }
}
