//! Laporta elimination of sparse linear systems over GF(p), learned once and replayed.
//!
//! Unknowns are columns ordered by complexity, and elimination always removes the most complex
//! one first, so the columns left without a pivot are the simplest: the master integrals.
//! Learning runs on one numeric sample and records which equations the targets actually depend
//! on and where each pivots. Replay builds only those equations at every later sample and checks
//! that the pivots agree, which turns a large, redundant system into a lean black box.
#![allow(
    clippy::must_use_candidate,
    clippy::missing_panics_doc,
    clippy::multiple_crate_versions,
    clippy::many_single_char_names
)]

pub mod block;
pub mod family;
pub mod formats;
pub mod ibp;

#[cfg(doctest)]
#[doc = include_str!("../../../README.md")]
struct ReadmeDoctests;

use polycore::{modp_echelon::Echelon, Lead};
use std::collections::{BTreeMap, BTreeSet};

/// A sparse equation `sum c_j * u_j = 0` as `(j, c_j)` pairs.
pub type Row = Vec<(usize, u64)>;

#[derive(Clone, Debug)]
pub struct Plan {
    eqs: Vec<usize>,
    pivots: Vec<usize>,
    pub targets: Vec<usize>,
    pub masters: Vec<usize>,
}

impl Plan {
    /// Eliminates the whole system at one sample and keeps what `targets` need.
    pub fn learn(rows: &[Row], targets: &[usize], p: u64) -> Self {
        // Simplest equations first, so pivots are found before they are needed to reduce others.
        let mut order: Vec<usize> = (0..rows.len()).collect();
        order.sort_by_cached_key(|&i| {
            let mut cols: Vec<usize> = rows[i].iter().map(|t| t.0).collect();
            cols.sort_unstable_by(|a, b| b.cmp(a));
            (cols.first().copied(), cols.len(), cols)
        });
        let mut ech = Echelon::recording(Lead::High, p);
        let mut origin = BTreeMap::new();
        for (k, &i) in order.iter().enumerate() {
            if let Some(c) = ech.insert(&rows[i]) {
                origin.insert(c, (k, i));
            }
        }
        let mut masters = BTreeSet::new();
        let mut stack = Vec::new();
        for &t in targets {
            let (solution, used) = ech.solve(t);
            masters.extend(solution.iter().map(|m| m.0));
            stack.extend(used);
        }
        let mut needed = BTreeSet::new();
        while let Some(c) = stack.pop() {
            if needed.insert(c) {
                stack.extend(ech.uses(c));
            }
        }
        let kept: BTreeMap<(usize, usize), usize> =
            needed.iter().map(|c| (origin[c], *c)).collect();
        Self {
            eqs: kept.keys().map(|k| k.1).collect(),
            pivots: kept.values().copied().collect(),
            targets: targets.to_vec(),
            masters: masters.into_iter().collect(),
        }
    }

    /// Equations kept out of the learned system.
    pub const fn len(&self) -> usize {
        self.eqs.len()
    }

    pub const fn is_empty(&self) -> bool {
        self.eqs.is_empty()
    }

    /// The coefficient of every master in every target, targets major, from only the kept
    /// equations `row(i)`. `None` if a pivot vanished at this sample.
    pub fn replay(&self, row: impl Fn(usize) -> Row, p: u64) -> Option<Vec<u64>> {
        let mut ech = Echelon::new(Lead::High, p);
        for (&i, &c) in self.eqs.iter().zip(&self.pivots) {
            if ech.insert(&row(i))? != c {
                return None;
            }
        }
        let mut out = vec![0; self.targets.len() * self.masters.len()];
        for (&t, chunk) in self
            .targets
            .iter()
            .zip(out.chunks_mut(self.masters.len().max(1)))
        {
            for (m, v) in ech.solve(t).0 {
                chunk[self.masters.binary_search(&m).ok()?] = v;
            }
        }
        Some(out)
    }
}
