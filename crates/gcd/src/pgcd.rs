//! Zippel's recursive gcd over GF(p): dense interpolation in the last variable, with every image
//! after the first found by sparse interpolation against the first image's skeleton.

use crate::linzip::linzip;
use crate::poly::{Exps, ModPoly};
use std::cmp::Ordering;
use std::collections::BTreeMap;
use zippel_interp::modp::Rng;
use zippel_interp::univariate::{self as uni, Dense};

/// `gcd(f, g)` up to a scalar, where only `x_0..=x_k` occur.
pub(crate) fn pgcd(f: &ModPoly, g: &ModPoly, k: usize, p: u64, rng: &mut Rng) -> Option<ModPoly> {
    if k == 0 {
        let h = uni::gcd(&f.eval_except(0, &[], p), &g.eval_except(0, &[], p), p);
        return Some(ModPoly::from_groups(f.n, 0, vec![(vec![0; f.n], h)]));
    }
    let (cf, f) = f.primitive(k, p);
    let (cg, g) = g.primitive(k, p);
    let content = uni::gcd(&cf, &cg, p);
    let lf = f.groups(k).swap_remove(0).1;
    let lg = g.groups(k).swap_remove(0).1;
    let gamma = uni::gcd(&lf, &lg, p);
    let need = degree_bound(&f, &g, k, p, rng)? + uni::deg(&gamma) + 1;

    let mut images = Accumulator::default();
    for _ in 0..4 * need + 16 {
        let a = rng.nonzero(p);
        if uni::eval(&lf, a, p) == 0 || uni::eval(&lg, a, p) == 0 || images.points.contains(&a) {
            continue;
        }
        let (fa, ga) = (f.eval_var(k, a, p), g.eval_var(k, a, p));
        let sparse = images
            .skeleton
            .as_ref()
            .filter(|_| k >= 2)
            .and_then(|s| linzip(&fa, &ga, s, k - 1, p, rng));
        let Some(h) = sparse.or_else(|| pgcd(&fa, &ga, k - 1, p, rng)) else {
            continue;
        };
        images.offer(a, h, uni::eval(&gamma, a, p), p);
        if images.points.len() == need {
            let groups = images.interpolate(p);
            let h = ModPoly::from_groups(f.n, k, groups).primitive(k, p).1;
            let groups = h
                .groups(k)
                .into_iter()
                .map(|(e, d)| (e, uni::mul_poly(&d, &content, p)))
                .collect();
            return Some(ModPoly::from_groups(f.n, k, groups));
        }
    }
    None
}

/// `deg_{x_k} gcd(f, g)` is at most the degree of a gcd of univariate images in `x_k`.
fn degree_bound(f: &ModPoly, g: &ModPoly, k: usize, p: u64, rng: &mut Rng) -> Option<usize> {
    (0..8).find_map(|_| {
        let point: Vec<u64> = (0..f.n).map(|_| rng.nonzero(p)).collect();
        let (uf, ug) = (f.eval_except(k, &point, p), g.eval_except(k, &point, p));
        (uni::deg(&uf) == f.degree(k) && uni::deg(&ug) == g.degree(k))
            .then(|| uni::deg(&uni::gcd(&uf, &ug, p)))
    })
}

/// Images scaled to a common normalization, discarding those of the wrong shape.
#[derive(Default)]
pub(crate) struct Accumulator {
    pub(crate) points: Vec<u64>,
    pub(crate) images: Vec<ModPoly>,
    pub(crate) skeleton: Option<Vec<Exps>>,
}

impl Accumulator {
    pub(crate) fn offer(&mut self, point: u64, h: ModPoly, lead: u64, p: u64) {
        let mut h = h.monic(p);
        h.scale(lead, p);
        let Some(reset) = reshape(self.skeleton.as_deref(), &h) else {
            return;
        };
        if reset {
            self.points.clear();
            self.images.clear();
            self.skeleton = Some(h.support());
        }
        self.points.push(point);
        self.images.push(h);
    }

    fn interpolate(&self, p: u64) -> Vec<(Exps, Dense)> {
        let mut values: BTreeMap<&Exps, Vec<u64>> = BTreeMap::new();
        for (j, h) in self.images.iter().enumerate() {
            for (e, c) in &h.terms {
                values
                    .entry(e)
                    .or_insert_with(|| vec![0; self.images.len()])[j] = *c;
            }
        }
        values
            .into_iter()
            .rev()
            .map(|(e, ys)| (e.clone(), uni::interpolate(&self.points, &ys, p)))
            .collect()
    }
}

/// Compare an image against the current skeleton: `None` discards it, `Some(true)` starts over.
/// An unlucky image has a larger leading monomial than the true gcd; one with terms outside
/// the skeleton reveals that the skeleton came from an unlucky choice.
pub(crate) fn reshape(skeleton: Option<&[Exps]>, h: &ModPoly) -> Option<bool> {
    let Some(s) = skeleton else {
        return Some(true);
    };
    match h.lm().cmp(s[0].as_slice()) {
        Ordering::Greater => None,
        Ordering::Less => Some(true),
        Ordering::Equal => Some(
            !h.terms
                .iter()
                .all(|(e, _)| s.binary_search_by(|x| e.cmp(x)).is_ok()),
        ),
    }
}
