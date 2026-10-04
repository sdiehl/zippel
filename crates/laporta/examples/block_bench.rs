//! Probe benchmarks with setup/refit costs and fresh-prime agreement checks.
//! Optional arguments: `double-box`, or `nonplanar-cut` (maximal cut of the
//! five-point nonplanar double pentagon, four invariants fixed to 1,3,5,7).
mod families;
use polycore::{modp::Primes, sample::Rng};
use std::{
    hint::black_box,
    sync::atomic::{AtomicUsize, Ordering},
    time::Instant,
};
use zippel_laporta::{
    block::{AdaptiveSearch, BlockPlan, Search},
    formats::{read_family_mathematica, read_family_yaml},
};

#[allow(clippy::too_many_lines)]
fn main() {
    let imported = read_family_yaml(
        include_str!("../tests/fixtures/box.yaml"),
        include_str!("../tests/fixtures/kinematics.yaml"),
        "box",
    )
    .unwrap();
    let family = imported.family.fix("s", 1);
    let mode = std::env::args().nth(1).unwrap_or_default();
    let double = mode == "double-box";
    let (system, targets) = if double {
        (
            families::double_box().fix("s", 1).system(1, 1),
            vec![
                vec![2, 1, 1, 1, 1, 1, 1, 0, 0],
                vec![1, 2, 1, 1, 1, 1, 1, 0, 0],
            ],
        )
    } else if mode == "nonplanar-cut" {
        let mut family = read_family_mathematica(
            include_str!("../tests/fixtures/nonplanar-double-pentagon.m"),
            "nonplanar",
            &["s12", "s23", "s34", "s45", "s51"],
            8,
        )
        .unwrap()
        .family;
        family.cuts = (0..8).collect();
        for (v, value) in [("s12", 1), ("s34", 3), ("s45", 5), ("s51", 7)] {
            family = family.fix(v, value);
        }
        (
            family.system(2, 2).unwrap(),
            vec![
                vec![2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0],
                vec![1, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0],
            ],
        )
    } else {
        (
            family.system(2, 2).unwrap(),
            vec![vec![2, 1, 1, 1], vec![1, 2, 1, 1], vec![2, 2, 1, 1]],
        )
    };
    let mut plan = system.learn(&targets, 7);
    eprintln!(
        "{} equations; {} retained; {} masters",
        system.len(),
        plan.len(),
        plan.masters.len()
    );
    plan.targets.sort_unstable();
    let calls = AtomicUsize::new(0);
    let oracle = |x: &[u64], p| {
        calls.fetch_add(1, Ordering::Relaxed);
        plan.replay(|e| system.row(e, x, p), p)
    };
    let p = Primes::new().next().unwrap();
    let start = Instant::now();
    let (blocks, form, report) = BlockPlan::learn_adaptive(
        oracle,
        2,
        plan.targets.len(),
        plan.masters.len(),
        p,
        7,
        &AdaptiveSearch {
            limits: Search {
                max_degree: 8,
                ..Search::default()
            },
            ..AdaptiveSearch::default()
        },
    )
    .expect("block search budget");
    let learn = start.elapsed();
    eprintln!("adaptive search: {report:?}");
    let learning_calls = calls.swap(0, Ordering::Relaxed);
    let prime = Primes::new().nth(1).unwrap();
    let start = Instant::now();
    let fitted = blocks.fit(oracle, prime, 7).unwrap();
    let fit = start.elapsed();
    let fitting_calls = calls.swap(0, Ordering::Relaxed);
    let mut rng = Rng::new(31);
    let points: Vec<_> = (0..1000)
        .map(|_| vec![rng.nonzero(prime), rng.nonzero(prime)])
        .collect();
    for x in points.iter().take(10) {
        assert_eq!(fitted.eval(x), oracle(x, prime));
    }
    let start = Instant::now();
    for x in &points {
        black_box(oracle(x, prime).unwrap());
    }
    let laporta = start.elapsed();
    let start = Instant::now();
    for x in &points {
        black_box(fitted.eval(x).unwrap());
    }
    let block = start.elapsed();
    println!(
        "{} equations, {} retained, {} masters; {} blocks, {} terms",
        system.len(),
        plan.len(),
        plan.masters.len(),
        form.blocks(),
        form.terms()
    );
    println!(
        "learn: {learning_calls} oracle probes, {learn:?}; refit: {fitting_calls} probes, {fit:?}"
    );
    println!(
        "1000 probes: Laporta {laporta:?}, blocks {block:?}, {:.2}x",
        laporta.as_secs_f64() / block.as_secs_f64()
    );
}
