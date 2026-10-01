//! Degree certificates shared by the sparse GCD algorithms.

use crate::fast::Arithmetic;
use crate::poly::ModPoly;
use polycore::modp::{add, inv, mul};
use polycore::sample::Rng;
use zippel_interp::poly::power_table;

/// Degree-preserving univariate images give rigorous upper bounds in EVERY
/// variable. Weighted spans alone are insufficient: x-y has span zero when
/// the two weights coincide, so a constant image would miss that factor.
pub(crate) fn degree_bounds(f: &ModPoly, g: &ModPoly, p: u64, rng: &mut Rng) -> Option<Vec<usize>> {
    let field = Arithmetic { p };
    let point: Vec<_> = (0..f.n).map(|_| rng.nonzero(p)).collect();
    let aa = univariate_images(f, &point, p);
    let bb = univariate_images(g, &point, p);
    aa.iter()
        .zip(&bb)
        .map(|(a, b)| {
            if a.last() == Some(&0) || b.last() == Some(&0) {
                return None;
            }
            Some(field.gcd(a, b).len() - 1)
        })
        .collect()
}

/// All leave-one-variable-free images in O(n(T+D)) field operations. Compute
/// each term at the full point once, then undo one variable's contribution.
pub(crate) fn univariate_images(f: &ModPoly, point: &[u64], p: u64) -> Vec<Vec<u64>> {
    let degrees = f.degrees();
    let powers = power_table(point, &degrees, p);
    let inverse_point: Vec<_> = point.iter().map(|&b| inv(b, p)).collect();
    let inverses = power_table(&inverse_point, &degrees, p);
    let mut images: Vec<_> = degrees.iter().map(|&d| vec![0; d + 1]).collect();
    for (e, c) in &f.terms {
        let value = e
            .iter()
            .zip(&powers)
            .fold(*c, |v, (&e, pw)| mul(v, pw[e as usize], p));
        for (k, &e) in e.iter().enumerate() {
            images[k][e as usize] = add(
                images[k][e as usize],
                mul(value, inverses[k][e as usize], p),
                p,
            );
        }
    }
    images
}
