//! Budgeted exact sparse division over a prime field.

use crate::poly::ModPoly;
use polycore::modp::{inv, mul, sub};
use std::collections::BTreeMap;

pub(crate) fn quotient(f: &ModPoly, g: &ModPoly, p: u64) -> Option<ModPoly> {
    if f.n != g.n {
        return None;
    }
    let (lead, lc) = g.terms.first()?;
    let li = inv(*lc, p);
    let mut rem: BTreeMap<_, _> = f.terms.iter().cloned().collect();
    let mut terms = Vec::new();
    let mut budget = 100_000usize;
    while let Some((e, c)) = rem.pop_last() {
        budget = budget.checked_sub(g.terms.len())?;
        let q = e
            .iter()
            .zip(lead)
            .map(|(a, b)| a.checked_sub(*b))
            .collect::<Option<Vec<_>>>()?;
        let c = mul(c, li, p);
        for (e, v) in &g.terms[1..] {
            let e = e
                .iter()
                .zip(&q)
                .map(|(e, q)| e.checked_add(*q))
                .collect::<Option<Vec<_>>>()?;
            let value = sub(rem.get(&e).copied().unwrap_or(0), mul(c, *v, p), p);
            if value == 0 {
                rem.remove(&e);
            } else {
                rem.insert(e, value);
            }
        }
        terms.push((q, c));
    }
    Some(ModPoly { n: f.n, terms })
}
