//! GCD verification's work limit around polycore's exact sparse division.

use crate::poly::ModPoly;

pub(crate) fn quotient(f: &ModPoly, g: &ModPoly, p: u64) -> Option<ModPoly> {
    f.to_poly(p)
        .exact_with_budget(&g.to_poly(p), 100_000)
        .map(|q| ModPoly::from_poly(&q, p))
}
