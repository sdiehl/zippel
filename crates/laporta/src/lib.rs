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

pub mod ibp;

use std::collections::{BTreeMap, BTreeSet, BinaryHeap};
use zippel_interp::modp::{add, inv, mul, sub};

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
        let mut ech = Echelon::default();
        let mut origin = BTreeMap::new();
        for (k, &i) in order.iter().enumerate() {
            if let Some(c) = ech.insert(&rows[i], p) {
                origin.insert(c, (k, i));
            }
        }
        let mut masters = BTreeSet::new();
        let mut stack = Vec::new();
        for &t in targets {
            let (solution, used) = ech.solve(t, p);
            masters.extend(solution.iter().map(|m| m.0));
            stack.extend(used);
        }
        let mut needed = BTreeSet::new();
        while let Some(c) = stack.pop() {
            if needed.insert(c) {
                stack.extend(&ech.uses[c]);
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
        let mut ech = Echelon::default();
        for (&i, &c) in self.eqs.iter().zip(&self.pivots) {
            if ech.insert(&row(i), p)? != c {
                return None;
            }
        }
        let mut out = vec![0; self.targets.len() * self.masters.len()];
        for (&t, chunk) in self
            .targets
            .iter()
            .zip(out.chunks_mut(self.masters.len().max(1)))
        {
            for (m, v) in ech.solve(t, p).0 {
                chunk[self.masters.binary_search(&m).ok()?] = v;
            }
        }
        Some(out)
    }
}

/// Monic rows by leading column, the pivots each was reduced by, and a dense scratch row.
#[derive(Debug, Default)]
struct Echelon {
    rows: Vec<Option<Row>>,
    uses: Vec<Vec<usize>>,
    acc: Vec<u64>,
    live: Vec<bool>,
}

impl Echelon {
    fn grow(&mut self, c: usize) {
        if c >= self.acc.len() {
            self.rows.resize(c + 1, None);
            self.uses.resize(c + 1, Vec::new());
            self.acc.resize(c + 1, 0);
            self.live.resize(c + 1, false);
        }
    }

    /// Reduces `row` column by column from the top: by every pivot if `full`, else only until
    /// the leading column has none. The remaining terms, descending, and the pivots used.
    fn reduce(&mut self, row: &[(usize, u64)], full: bool, p: u64) -> (Row, Vec<usize>) {
        let mut heap = BinaryHeap::new();
        for &(c, v) in row {
            self.grow(c);
            self.acc[c] = add(self.acc[c], v, p);
            if !self.live[c] {
                self.live[c] = true;
                heap.push(c);
            }
        }
        let Self {
            rows, acc, live, ..
        } = self;
        let (mut out, mut used) = (Row::new(), Vec::new());
        while let Some(c) = heap.pop() {
            live[c] = false;
            let v = std::mem::take(&mut acc[c]);
            match &rows[c] {
                _ if v == 0 => {}
                Some(pivot) if full || out.is_empty() => {
                    used.push(c);
                    for &(k, x) in &pivot[1..] {
                        acc[k] = sub(acc[k], mul(v, x, p), p);
                        if !live[k] {
                            live[k] = true;
                            heap.push(k);
                        }
                    }
                }
                _ => out.push((c, v)),
            }
        }
        (out, used)
    }

    /// Reduces `row` until its leading column has no pivot and makes it one. `None` if the row
    /// vanished, being a consequence of the rows before it.
    fn insert(&mut self, row: &[(usize, u64)], p: u64) -> Option<usize> {
        let (mut r, used) = self.reduce(row, false, p);
        let (c, v) = *r.first()?;
        let l = inv(v, p);
        for t in &mut r {
            t.1 = mul(t.1, l, p);
        }
        self.rows[c] = Some(r);
        self.uses[c] = used;
        Some(c)
    }

    /// Column `t` as a combination of columns without a pivot, the masters, and the pivots used.
    fn solve(&mut self, t: usize, p: u64) -> (Row, Vec<usize>) {
        self.grow(t);
        let Some(row) = self.rows[t].clone() else {
            return (vec![(t, 1)], Vec::new());
        };
        let (mut r, mut used) = self.reduce(&row[1..], true, p);
        for m in &mut r {
            m.1 = sub(0, m.1, p);
        }
        used.push(t);
        (r, used)
    }
}
