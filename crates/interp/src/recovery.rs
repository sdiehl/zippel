//! Optional work limits and diagnostics for Zippel and rational reconstruction.

use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

use polycore::sample::BlackBox;

/// Limits for one reconstruction. Defaults preserve the original search sizes.
#[derive(Clone, Copy, Debug)]
pub struct RecoveryOptions {
    /// Candidate points per variable or Thiele ray, including refused and repeated points.
    pub max_points: usize,
    /// Independent attempts, including the first. Also applies to each nested interpolant.
    pub max_attempts: usize,
    /// Total oracle evaluations, including refused samples and batched evaluations.
    pub max_probes: usize,
}

impl Default for RecoveryOptions {
    fn default() -> Self {
        Self {
            max_points: 4096,
            max_attempts: 3,
            max_probes: usize::MAX,
        }
    }
}

/// Work performed, available on both success and failure.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct RecoveryReport {
    /// Calls made to the supplied oracle; cached values do not count again.
    pub probes: usize,
    /// Oracle calls returning `None`; budget-denied calls are excluded.
    pub rejected_probes: usize,
    /// Top-level attempts, including the first; retries are `attempts.saturating_sub(1)`.
    pub attempts: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FailureReason {
    /// Recovery needed another oracle evaluation after the probe allowance was used.
    ProbeBudgetExceeded,
    /// Every oracle evaluation was refused (at least one evaluation was made).
    NoUsableSamples,
    /// The configured search did not recover a validated result. This does not prove
    /// that the oracle is non-polynomial or that a larger budget would succeed.
    RecoveryFailed,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RecoveryError {
    pub reason: FailureReason,
    pub report: RecoveryReport,
}

impl std::fmt::Display for RecoveryError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{:?} after {} probes and {} attempts",
            self.reason, self.report.probes, self.report.attempts
        )
    }
}

impl std::error::Error for RecoveryError {}

/// Preserve a custom oracle's batch implementation, while reserving its allowance atomically.
pub(crate) struct Tracked<'a, F> {
    f: &'a F,
    limit: usize,
    probes: AtomicUsize,
    rejected: AtomicUsize,
    denied: AtomicBool,
}

impl<'a, F> Tracked<'a, F> {
    pub(crate) const fn new(f: &'a F, limit: usize) -> Self {
        Self {
            f,
            limit,
            probes: AtomicUsize::new(0),
            rejected: AtomicUsize::new(0),
            denied: AtomicBool::new(false),
        }
    }

    fn reserve(&self, requested: usize) -> usize {
        // `try_update` is the new name on Rust 1.99, but is unavailable on our 1.88 MSRV.
        #[allow(deprecated)]
        let before = self
            .probes
            .fetch_update(Ordering::Relaxed, Ordering::Relaxed, |used| {
                Some(used + requested.min(self.limit - used))
            })
            .expect("reservation always succeeds");
        let granted = requested.min(self.limit - before);
        if granted < requested {
            self.denied.store(true, Ordering::Relaxed);
        }
        granted
    }

    pub(crate) fn finish<T>(
        &self,
        value: Option<T>,
        attempts: usize,
    ) -> Result<(T, RecoveryReport), RecoveryError> {
        let report = RecoveryReport {
            probes: self.probes.load(Ordering::Relaxed),
            rejected_probes: self.rejected.load(Ordering::Relaxed),
            attempts,
        };
        value.map(|v| (v, report)).ok_or_else(|| RecoveryError {
            reason: if self.denied.load(Ordering::Relaxed) {
                FailureReason::ProbeBudgetExceeded
            } else if report.probes > 0 && report.probes == report.rejected_probes {
                FailureReason::NoUsableSamples
            } else {
                FailureReason::RecoveryFailed
            },
            report,
        })
    }
}

impl<F: BlackBox> BlackBox for Tracked<'_, F> {
    fn eval(&self, x: &[u64], p: u64) -> Option<u64> {
        if self.reserve(1) == 0 {
            return None;
        }
        let value = self.f.eval(x, p);
        if value.is_none() {
            self.rejected.fetch_add(1, Ordering::Relaxed);
        }
        value
    }

    fn eval_many(&self, xs: &[Vec<u64>], p: u64) -> Vec<Option<u64>> {
        let granted = self.reserve(xs.len());
        let mut values = if granted == 0 {
            Vec::new()
        } else {
            self.f.eval_many(&xs[..granted], p)
        };
        self.rejected.fetch_add(
            values.iter().filter(|v| v.is_none()).count(),
            Ordering::Relaxed,
        );
        values.resize(xs.len(), None);
        values
    }
}
