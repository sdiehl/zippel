//! Degree certificates shared by the sparse GCD algorithms.

use crate::poly::ModPoly;
use polycore::fast::PrimeField;
use polycore::sample::Rng;
use polycore::{Fp, Modular};

/// Degree-preserving univariate images give rigorous upper bounds in EVERY
/// variable. Weighted spans alone are insufficient: x-y has span zero when
/// the two weights coincide, so a constant image would miss that factor.
pub(crate) fn degree_bounds(f: &ModPoly, g: &ModPoly, p: u64, rng: &mut Rng) -> Option<Vec<usize>> {
    let field = PrimeField::new(p).expect("prime modulus");
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

/// Shared evaluation preserves coefficient slots so degree loss is observable.
pub(crate) fn univariate_images(f: &ModPoly, point: &[u64], p: u64) -> Vec<Vec<u64>> {
    let point: Vec<_> = point.iter().map(|&x| Fp::new(x, p)).collect();
    f.to_poly(p)
        .univariate_images(&point)
        .into_iter()
        .zip(f.degrees())
        .map(|(image, degree)| {
            let mut cs: Vec<_> = image.0.iter().map(|c| c.residue_mod(p)).collect();
            cs.resize(degree + 1, 0);
            cs
        })
        .collect()
}
