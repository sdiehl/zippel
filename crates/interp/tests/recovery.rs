use std::sync::atomic::{AtomicUsize, Ordering};

use polycore::modp::{add, inv, mul, Primes};
use polycore::sample::BlackBox;
use zippel_interp::{
    interpolate, interpolate_with_options, reconstruct, reconstruct_with_options, FailureReason,
    RecoveryOptions,
};

#[test]
fn point_limits_can_be_raised_and_reports_match_oracle_calls() {
    let calls = AtomicUsize::new(0);
    let p = Primes::new().next().unwrap();
    let oracle = |x: &[u64], p| {
        calls.fetch_add(1, Ordering::Relaxed);
        Some(add(mul(x[0], x[0], p), 7, p))
    };
    let options = RecoveryOptions {
        max_points: 1,
        max_attempts: 2,
        ..RecoveryOptions::default()
    };
    let error = interpolate_with_options(&oracle, 1, p, 7, &options).unwrap_err();
    assert_eq!(error.reason, FailureReason::RecoveryFailed);
    assert_eq!(error.report.probes, 2);
    assert_eq!(error.report.attempts, 2);
    assert_eq!(calls.swap(0, Ordering::Relaxed), error.report.probes);

    let options = RecoveryOptions {
        max_points: 16,
        ..options
    };
    let (value, report) = interpolate_with_options(&oracle, 1, p, 7, &options).unwrap();
    assert_eq!(calls.load(Ordering::Relaxed), report.probes);
    assert_eq!(report.rejected_probes, 0);
    assert_eq!(report.attempts, 1);
    assert_eq!(Some(value), interpolate(&oracle, 1, p, 7));
}

struct Batched {
    calls: AtomicUsize,
    batches: AtomicUsize,
}

impl BlackBox for Batched {
    fn eval(&self, _: &[u64], _: u64) -> Option<u64> {
        panic!("the custom batch evaluator must be preserved");
    }

    fn eval_many(&self, xs: &[Vec<u64>], _: u64) -> Vec<Option<u64>> {
        self.calls.fetch_add(xs.len(), Ordering::Relaxed);
        self.batches.fetch_add(1, Ordering::Relaxed);
        vec![None; xs.len()]
    }
}

#[test]
fn a_partial_batch_never_exceeds_the_probe_budget() {
    let oracle = Batched {
        calls: AtomicUsize::new(0),
        batches: AtomicUsize::new(0),
    };
    let options = RecoveryOptions {
        max_probes: 3,
        max_points: 8,
        max_attempts: 2,
    };
    let error = interpolate_with_options(&oracle, 1, 101, 3, &options).unwrap_err();
    assert_eq!(error.reason, FailureReason::ProbeBudgetExceeded);
    assert_eq!(error.report.probes, 3);
    assert_eq!(error.report.rejected_probes, 3);
    assert_eq!(oracle.calls.load(Ordering::Relaxed), 3);
    assert_eq!(oracle.batches.load(Ordering::Relaxed), 1);
}

#[test]
fn wholly_refused_samples_and_disabled_search_are_distinct() {
    let options = RecoveryOptions::default();
    let error = interpolate_with_options(&|_: &[u64], _| None, 0, 101, 1, &options).unwrap_err();
    assert_eq!(error.reason, FailureReason::NoUsableSamples);
    assert_eq!(error.report.probes, options.max_attempts);
    let disabled = RecoveryOptions {
        max_attempts: 0,
        ..options
    };
    let error = interpolate_with_options(&|_: &[u64], _| panic!("no calls"), 1, 101, 1, &disabled)
        .unwrap_err();
    assert_eq!(error.reason, FailureReason::RecoveryFailed);
    assert_eq!(error.report.probes, 0);
}

#[test]
fn retries_are_counted_and_success_at_the_exact_budget_is_allowed() {
    let calls = AtomicUsize::new(0);
    let oracle = |_: &[u64], _| (calls.fetch_add(1, Ordering::Relaxed) > 0).then_some(7);
    let options = RecoveryOptions {
        max_probes: 2,
        ..RecoveryOptions::default()
    };
    let (value, report) = interpolate_with_options(&oracle, 0, 101, 1, &options).unwrap();
    assert_eq!(value.eval(&[], 101), 7);
    assert_eq!(report.probes, 2);
    assert_eq!(report.attempts, 2);
    assert_eq!(report.rejected_probes, 1);
}

#[test]
fn rational_budget_counts_actual_calls_across_nested_rays() {
    let calls = AtomicUsize::new(0);
    let p = Primes::new().next().unwrap();
    let oracle = |x: &[u64], p| {
        calls.fetch_add(1, Ordering::Relaxed);
        let den = add(x[1], 1, p);
        (den != 0).then(|| mul(add(mul(x[0], x[0], p), 7, p), inv(den, p), p))
    };
    let options = RecoveryOptions {
        max_probes: 5,
        max_points: 16,
        max_attempts: 2,
    };
    let error = reconstruct_with_options(&oracle, 2, p, 7, &options).unwrap_err();
    assert_eq!(error.reason, FailureReason::ProbeBudgetExceeded);
    assert_eq!(error.report.probes, 5);
    assert_eq!(calls.swap(0, Ordering::Relaxed), 5);
    let (value, report) =
        reconstruct_with_options(&oracle, 2, p, 7, &RecoveryOptions::default()).unwrap();
    assert_eq!(calls.load(Ordering::Relaxed), report.probes);
    assert_eq!(Some(value), reconstruct(&oracle, 2, p, 7));
}
