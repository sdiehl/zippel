//! Reproducible probe benchmark on an imported Kira box; no external reducer needed.
use polycore::{modp::Primes, sample::Rng};
use std::{
    hint::black_box,
    sync::atomic::{AtomicUsize, Ordering},
    time::Instant,
};
use zippel_laporta::{
    block::{BlockPlan, Search},
    formats::read_kira,
};

fn main() {
    let imported = read_kira(
        include_str!("../tests/fixtures/box.yaml"),
        include_str!("../tests/fixtures/kinematics.yaml"),
        "box",
    )
    .unwrap();
    let family = imported.family.fix("s", 1);
    let system = family.system(2, 2);
    let mut plan = system.learn(&[vec![2, 1, 1, 1], vec![1, 2, 1, 1], vec![2, 2, 1, 1]], 7);
    plan.targets.sort_unstable();
    let calls = AtomicUsize::new(0);
    let oracle = |x: &[u64], p| {
        calls.fetch_add(1, Ordering::Relaxed);
        plan.replay(|e| system.row(e, x, p), p)
    };
    let p = Primes::new().next().unwrap();
    let start = Instant::now();
    let (blocks, form) = BlockPlan::learn(
        oracle,
        2,
        plan.targets.len(),
        plan.masters.len(),
        p,
        7,
        Search {
            max_degree: 5,
            ..Search::default()
        },
    )
    .expect("block search budget");
    let learn = start.elapsed();
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
