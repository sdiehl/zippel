use polycore::modp::{add, mul, try_inv, Primes};
use polycore::sample::Rng;
use zippel_laporta::block::{BlockPlan, Search};
use zippel_laporta::formats::read_family_yaml;

fn oracle(x: &[u64], p: u64) -> Option<Vec<u64>> {
    let a = mul(add(x[0], 1, p), try_inv(add(x[1], 2, p), p)?, p);
    let b = mul(add(a, x[1], p), try_inv(add(x[0], 3, p), p)?, p);
    Some(vec![a, b])
}

#[test]
fn learns_refits_and_lifts_triangular_relations() {
    let p = Primes::new().next().unwrap();
    for block_size in [1, 2] {
        let search = Search {
            block_size,
            max_degree: 2,
            ..Search::default()
        };
        let (plan, form) = BlockPlan::learn(oracle, 2, 2, 1, p, 1, search).unwrap();
        let mut rng = Rng::new(5);
        for _ in 0..20 {
            let x = [rng.nonzero(p), rng.nonzero(p)];
            assert_eq!(form.eval(&x), oracle(&x, p));
        }
        assert_eq!(form.eval(&[p, 0]), None);
        assert_eq!(form.eval(&[p - 3, 17]), None);
        for prime in Primes::new().skip(1).take(3) {
            let fitted = plan.fit(oracle, prime, 2).unwrap();
            for _ in 0..20 {
                let x = [rng.nonzero(prime), rng.nonzero(prime)];
                assert_eq!(fitted.eval(&x), oracle(&x, prime));
            }
        }
        assert_eq!(
            zippel_laporta::block::lift(oracle, 2, 2, 1, 1, search),
            zippel_lift::lift(oracle, 2, 1)
        );
    }
}

#[test]
fn bounded_failure_and_incompatible_prime() {
    let p = Primes::new().next().unwrap();
    assert!(BlockPlan::learn(
        oracle,
        2,
        2,
        1,
        p,
        1,
        Search {
            max_probes: 1,
            ..Search::default()
        }
    )
    .is_none());
    let (plan, _) = BlockPlan::learn(oracle, 2, 2, 1, p, 1, Search::default()).unwrap();
    assert!(plan.fit(|_, _| Some(vec![1]), p, 2).is_none());
    assert!(plan.fit(|x, _p| Some(vec![x[0], x[1]]), p, 2).is_none());
}

#[test]
fn imported_box_matches_laporta_on_fresh_primes() {
    let imported = read_family_yaml(
        include_str!("fixtures/box.yaml"),
        include_str!("fixtures/kinematics.yaml"),
        "box",
    )
    .unwrap();
    let family = imported.family.fix("s", 1);
    let system = family.system(1, 1).unwrap();
    let mut plan = system.learn(&[vec![2, 1, 1, 1], vec![1, 2, 1, 1]], 9);
    plan.targets.sort_unstable();
    let oracle = |x: &[u64], p| plan.replay(|e| system.row(e, x, p), p);
    let p = Primes::new().next().unwrap();
    let (blocks, first) = BlockPlan::learn(
        oracle,
        2,
        plan.targets.len(),
        plan.masters.len(),
        p,
        9,
        Search {
            max_degree: 5,
            ..Search::default()
        },
    )
    .unwrap();
    let mut rng = Rng::new(31);
    for prime in Primes::new().take(3) {
        let form = if prime == p {
            first.clone()
        } else {
            blocks.fit(oracle, prime, 9).unwrap()
        };
        for _ in 0..10 {
            let x = [rng.nonzero(prime), rng.nonzero(prime)];
            assert_eq!(form.eval(&x), oracle(&x, prime));
        }
    }
}

#[test]
fn adaptive_weights_groups_and_coupled_blocks() {
    use polycore::modp::{pow, sub};
    use zippel_laporta::block::AdaptiveSearch;
    let p = Primes::new().next().unwrap();
    // x*a+b=1, a+y*b=0. A degree-one joint ansatz exists, but
    // either single reduction has a degree-two denominator.
    let coupled = |x: &[u64], p| {
        let den = try_inv(sub(mul(x[0], x[1], p), 1, p), p)?;
        Some(vec![mul(x[1], den, p), sub(0, den, p)])
    };
    let search = AdaptiveSearch {
        limits: Search {
            max_degree: 1,
            ..Search::default()
        },
        variable_weights: vec![1, 1],
        max_block_size: 2,
        ..AdaptiveSearch::default()
    };
    let (plan, form, report) = BlockPlan::learn_adaptive(coupled, 2, 2, 1, p, 4, &search).unwrap();
    assert_eq!(report.block_sizes, vec![2]);
    assert_eq!(form.eval(&[13, 17]), coupled(&[13, 17], p));
    let q = Primes::new().nth(1).unwrap();
    assert_eq!(
        plan.fit(coupled, q, 6).unwrap().eval(&[19, 23]),
        coupled(&[19, 23], q)
    );
    let anisotropic = |x: &[u64], p| {
        Some(vec![mul(
            add(pow(x[0], 4, p), 1, p),
            try_inv(add(x[1], 2, p), p)?,
            p,
        )])
    };
    let search = AdaptiveSearch {
        limits: Search {
            max_degree: 5,
            ..Search::default()
        },
        groups: vec![(vec![1], 1)],
        ..AdaptiveSearch::default()
    };
    let (_, form, report) = BlockPlan::learn_adaptive(anisotropic, 2, 1, 1, p, 8, &search).unwrap();
    assert_eq!(report.variable_weights, vec![1, 4]);
    assert_eq!(form.eval(&[31, 37]), anisotropic(&[31, 37], p));
    assert!(report.oracle_probes < search.limits.max_probes);
    let bounded = AdaptiveSearch {
        limits: Search {
            max_probes: 2,
            ..Search::default()
        },
        ..search
    };
    assert!(BlockPlan::learn_adaptive(anisotropic, 2, 1, 1, p, 8, &bounded).is_none());
}

#[test]
fn automatic_intermediates_and_factored_lifting() {
    use zippel_laporta::block::AdaptiveSearch;
    let imported = read_family_yaml(
        include_str!("fixtures/box.yaml"),
        include_str!("fixtures/kinematics.yaml"),
        "box",
    )
    .unwrap();
    let family = imported.family.fix("s", 1);
    let system = family.system(1, 1).unwrap();
    let plan = system.learn(&[vec![2, 1, 1, 1]], 9);
    let search = AdaptiveSearch {
        max_intermediates: 3,
        ..AdaptiveSearch::default()
    };
    let reduction = system.learn_blocks(&plan, 9, &search).unwrap();
    assert!(reduction.intermediates() > 0);
    let ring = polycore::Ring::new(["d", "t"], polycore::Order::Lex);
    let candidates = [
        ring.parse("t").unwrap(),
        ring.parse("1+t").unwrap(),
        ring.parse("d-4").unwrap(),
    ];
    let coefficients = reduction.lift(&system, &candidates, 9).unwrap();
    let indices: Vec<_> = reduction
        .masters()
        .iter()
        .map(|&j| system.integrals[j].clone())
        .collect();
    let reference = system.learn(&[vec![2, 1, 1, 1]], 9);
    assert_eq!(
        indices,
        reference
            .masters
            .iter()
            .map(|&j| system.integrals[j].clone())
            .collect::<Vec<_>>()
    );
    assert_eq!(coefficients, system.lift(&reference, 9).unwrap());
}

#[test]
fn repeated_targets_preserve_output_order() {
    use zippel_laporta::block::AdaptiveSearch;
    let imported = read_family_yaml(
        include_str!("fixtures/box.yaml"),
        include_str!("fixtures/kinematics.yaml"),
        "box",
    )
    .unwrap();
    let system = imported.family.fix("s", 1).system(1, 1).unwrap();
    let plan = system.learn(&[vec![2, 1, 1, 1], vec![2, 1, 1, 1]], 17);
    let reduction = system
        .learn_blocks(
            &plan,
            17,
            &AdaptiveSearch {
                max_intermediates: 0,
                ..AdaptiveSearch::default()
            },
        )
        .unwrap();
    assert_eq!(reduction.intermediates(), 0);
    let p = Primes::new().next().unwrap();
    let output = reduction.output_plan();
    assert_eq!(
        reduction.eval(&[13, 17]),
        output.replay(|e| system.row(e, &[13, 17], p), p)
    );
    assert_eq!(output.targets.len(), 2);
}
