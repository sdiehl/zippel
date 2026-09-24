#![allow(clippy::cast_possible_truncation, dead_code, unreachable_pub)]

use zippel_interp::{ModPoly, Rng};

pub fn random_poly(rng: &mut Rng, n: usize, terms: usize, deg: u64, p: u64) -> ModPoly {
    let mut t: Vec<_> = (0..terms)
        .map(|_| {
            (
                (0..n)
                    .map(|_| (rng.nonzero(deg + 2) - 1) as u32)
                    .collect::<Vec<_>>(),
                rng.nonzero(p),
            )
        })
        .collect();
    t.sort_by(|a, b| b.0.cmp(&a.0));
    t.dedup_by(|a, b| a.0 == b.0);
    ModPoly { n, terms: t }
}

/// A point-dependent coin, for black boxes that fail at scattered points.
pub fn unlucky(x: &[u64]) -> bool {
    x.iter()
        .fold(17u64, |h, &v| (h ^ v).wrapping_mul(0x100_0000_01b3))
        % 7
        == 0
}
