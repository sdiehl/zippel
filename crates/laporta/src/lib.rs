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
    clippy::many_single_char_names
)]

use std::collections::{BTreeMap, BTreeSet};
use zippel_interp::modp::{add, inv, mul, sub};

/// A sparse equation `sum c_j * u_j = 0` as `(j, c_j)` pairs.
pub type Row = Vec<(usize, u64)>;

type Sparse = BTreeMap<usize, u64>;

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
        let mut ech = Echelon::default();
        let mut origin = BTreeMap::new();
        for (i, row) in rows.iter().enumerate() {
            if let Some(c) = ech.insert(row, p) {
                origin.insert(c, i);
            }
        }
        ech.back_substitute(p);
        let mut needed = BTreeSet::new();
        let mut stack: Vec<usize> = targets.to_vec();
        while let Some(c) = stack.pop() {
            if let Some(uses) = ech.uses.get(&c).filter(|_| needed.insert(c)) {
                stack.extend(uses);
            }
        }
        let kept: BTreeMap<usize, usize> = needed.iter().map(|c| (origin[c], *c)).collect();
        let masters = targets
            .iter()
            .flat_map(|t| {
                ech.rows.get(t).map_or_else(
                    || vec![*t],
                    |r| r.keys().copied().filter(|k| k != t).collect(),
                )
            })
            .collect::<BTreeSet<_>>();
        Self {
            eqs: kept.keys().copied().collect(),
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
        ech.back_substitute(p);
        let coefficient = |t: &usize, m: &usize| {
            ech.rows.get(t).map_or_else(
                || u64::from(t == m),
                |r| sub(0, r.get(m).copied().unwrap_or(0), p),
            )
        };
        Some(
            self.targets
                .iter()
                .flat_map(|t| self.masters.iter().map(move |m| coefficient(t, m)))
                .collect(),
        )
    }
}

/// Monic rows keyed by leading column, and the pivots each row was reduced by.
#[derive(Debug, Default)]
struct Echelon {
    rows: BTreeMap<usize, Sparse>,
    uses: BTreeMap<usize, Vec<usize>>,
}

impl Echelon {
    /// Reduces `row` until its leading column has no pivot and makes it one. `None` if the row
    /// vanished, being a consequence of the rows before it.
    fn insert(&mut self, row: &Row, p: u64) -> Option<usize> {
        let mut r = Sparse::new();
        for &(c, v) in row {
            let e = r.entry(c).or_insert(0);
            *e = add(*e, v, p);
        }
        r.retain(|_, v| *v != 0);
        let mut used = Vec::new();
        let (c, v) = loop {
            let (&c, &v) = r.last_key_value()?;
            let Some(pivot) = self.rows.get(&c) else {
                break (c, v);
            };
            axpy(&mut r, pivot, v, p);
            used.push(c);
        };
        let l = inv(v, p);
        for x in r.values_mut() {
            *x = mul(*x, l, p);
        }
        self.rows.insert(c, r);
        self.uses.insert(c, used);
        Some(c)
    }

    /// Clears every pivot below the leading one, lowest rows first, leaving only masters.
    fn back_substitute(&mut self, p: u64) {
        let cols: Vec<usize> = self.rows.keys().copied().collect();
        for c in cols {
            let mut r = self.rows.remove(&c).unwrap();
            let lower: Vec<usize> = r
                .range(..c)
                .map(|(&k, _)| k)
                .filter(|k| self.rows.contains_key(k))
                .collect();
            for &k in &lower {
                let v = r[&k];
                axpy(&mut r, &self.rows[&k], v, p);
            }
            self.uses.get_mut(&c).unwrap().extend(lower);
            self.rows.insert(c, r);
        }
    }
}

/// `r -= v * s`.
fn axpy(r: &mut Sparse, s: &Sparse, v: u64, p: u64) {
    for (&c, &x) in s {
        let e = r.entry(c).or_insert(0);
        *e = sub(*e, mul(v, x, p), p);
        if *e == 0 {
            r.remove(&c);
        }
    }
}
