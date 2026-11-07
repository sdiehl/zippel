use std::sync::atomic::{AtomicUsize, Ordering};

use polycore::modp::{add, mul, Primes};
use polycore::{Order, Ring};
use zippel_lift::{lift, lift_with_options, LiftFailureReason, LiftOptions};

#[test]
fn lifting_reports_refits_and_obeys_prime_and_probe_limits() {
    let calls = AtomicUsize::new(0);
    let oracle = |x: &[u64], p| {
        calls.fetch_add(1, Ordering::Relaxed);
        Some(vec![add(x[0], 1, p)])
    };
    let options = LiftOptions {
        max_primes: 1,
        ..LiftOptions::default()
    };
    let error = lift_with_options(oracle, 1, 7, &options).unwrap_err();
    assert_eq!(error.reason, LiftFailureReason::PrimeBudgetExceeded);
    assert_eq!(error.report.primes, 1);
    assert_eq!(calls.swap(0, Ordering::Relaxed), error.report.probes);

    let (value, report) = lift_with_options(oracle, 1, 7, &LiftOptions::default()).unwrap();
    assert_eq!(calls.swap(0, Ordering::Relaxed), report.probes);
    assert_eq!(report.primes, 2);
    assert_eq!(report.refits, 1);
    assert_eq!(report.failed_refits, 0);
    assert_eq!(report.reconstructions, 1);
    assert!(report.reconstruction_attempts > 0);
    assert_eq!(Some(value), lift(oracle, 1, 7));

    calls.store(0, Ordering::Relaxed);
    let options = LiftOptions {
        max_probes: 7,
        max_primes: 2,
        ..options
    };
    let error = lift_with_options(oracle, 1, 7, &options).unwrap_err();
    assert_eq!(error.reason, LiftFailureReason::ProbeBudgetExceeded);
    assert_eq!(error.report.probes, 7);
    assert_eq!(calls.load(Ordering::Relaxed), 7);
}

#[test]
fn lifting_distinguishes_refusals_from_unsuccessful_reconstruction() {
    let mut options = LiftOptions {
        max_primes: 2,
        ..LiftOptions::default()
    };
    let error = lift_with_options(|_, _| None, 0, 7, &options).unwrap_err();
    assert_eq!(error.reason, LiftFailureReason::NoUsableSamples);
    assert_eq!(error.report.probes, error.report.rejected_probes);
    assert_eq!(error.report.skipped_primes, 2);

    options.reconstruction.max_points = 1;
    options.reconstruction.max_attempts = 1;
    let error = lift_with_options(|x, _| Some(vec![x[0]]), 1, 7, &options).unwrap_err();
    assert_eq!(error.reason, LiftFailureReason::ReconstructionFailed);
    assert_eq!(error.report.reconstruction_attempts, 2);
    options.max_primes = 0;
    let error = lift_with_options(|_, _| panic!("no calls"), 1, 7, &options).unwrap_err();
    assert_eq!(error.reason, LiftFailureReason::PrimeBudgetExceeded);
    assert_eq!(error.report.probes, 0);
}

#[test]
fn an_unlucky_first_prime_records_the_fallback_to_full_reconstruction() {
    let coefficient = Primes::new().next().unwrap();
    let oracle = |x: &[u64], p| Some(vec![add(mul(coefficient % p, x[0], p), 1, p)]);
    let (value, report) = lift_with_options(oracle, 1, 7, &LiftOptions::default()).unwrap();
    let ring = Ring::new(["x"], Order::Lex);
    assert_eq!(
        value[0].num,
        ring.parse(&format!("{coefficient}*x+1")).unwrap()
    );
    assert_eq!(value[0].den, ring.parse("1").unwrap());
    assert!(report.failed_refits > 0);
    assert!(report.reconstructions > 1);
}
