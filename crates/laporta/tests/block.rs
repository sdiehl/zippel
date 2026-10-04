use polycore::modp::{add, mul, try_inv, Primes};
use polycore::sample::Rng;
use zippel_laporta::block::{BlockPlan, Search};
use zippel_laporta::formats::read_kira;

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
    let imported = read_kira(
        include_str!("fixtures/box.yaml"),
        include_str!("fixtures/kinematics.yaml"),
        "box",
    )
    .unwrap();
    let family = imported.family.fix("s", 1);
    let system = family.system(1, 1);
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
