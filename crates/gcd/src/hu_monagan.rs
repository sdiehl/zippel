//! Hu–Monagan geometric sparse GCD: univariate or bivariate images, with
//! simultaneous recovery of the GCD and a cofactor (Monagan's 2022 update).
//!
//! References: Hu & Monagan, JSC 105 (2021), 28–63, Section 5.3;
//! Monagan, “Speeding up polynomial GCD, a crucial operation in Maple” (2022),
//! Section 3. The smooth subgroup and exact certificates are implementation
//! choices; these bounded paths do not claim the papers' parallel complexity.

use crate::fast::Arithmetic;
use crate::geometric::{Encoding, Orbit, Stream};
use crate::images::degree_bounds;
use crate::modular_division::quotient;
use crate::pgcd::pgcd;
use crate::poly::{Exps, ModPoly};
use polycore::modp::{add, mul};
use polycore::sample::Rng;
use std::collections::BTreeMap;

const MAX_SAMPLES: usize = 256;
const MAX_CELLS: usize = 2_000_000;
const MAX_IMAGE_DEGREE: usize = 4096;

#[derive(Default)]
struct Interpolant {
    streams: BTreeMap<Exps, Stream>,
    samples: usize,
}

impl Interpolant {
    fn push(&mut self, image: &ModPoly, p: u64) -> Option<()> {
        let new_streams = image
            .terms
            .iter()
            .filter(|(e, _)| !self.streams.contains_key(e))
            .count();
        if (self.streams.len() + new_streams).checked_mul(self.samples + 1)? > MAX_CELLS {
            return None;
        }
        for (e, _) in &image.terms {
            self.streams.entry(e.clone()).or_insert_with(|| {
                let mut stream = Stream::default();
                for _ in 0..self.samples {
                    stream.push(0, p);
                }
                stream
            });
        }
        self.samples += 1;
        for (e, stream) in &mut self.streams {
            let value = image
                .terms
                .binary_search_by(|t| e.cmp(&t.0))
                .map_or(0, |i| image.terms[i].1);
            stream.push(value, p);
        }
        Some(())
    }

    fn recover(&self, orbit: &Orbit<'_>, n: usize) -> Option<ModPoly> {
        let mut terms = Vec::new();
        for (prefix, stream) in &self.streams {
            for (suffix, c) in stream.recover(orbit)? {
                let e: Vec<_> = prefix.iter().chain(&suffix).copied().collect();
                if c != 0 {
                    terms.push((e, c));
                }
            }
        }
        terms.sort_unstable_by(|a, b| b.0.cmp(&a.0));
        (!terms.is_empty()).then_some(ModPoly { n, terms })
    }
}

/// Evaluate consecutive geometric samples by updating each term once.
/// Prefix groups stay fixed even when a coefficient vanishes at one sample.
struct ImageSequence {
    keep: usize,
    groups: BTreeMap<Exps, Vec<(u64, u64)>>,
}

impl ImageSequence {
    fn new(f: &ModPoly, keep: usize, orbit: &Orbit<'_>, p: u64) -> Self {
        let mut groups: BTreeMap<Exps, Vec<(u64, u64)>> = BTreeMap::new();
        for (e, c) in &f.terms {
            let (value, ratio) = orbit.monomial(&e[keep..]);
            groups
                .entry(e[..keep].to_vec())
                .or_default()
                .push((mul(*c, value, p), ratio));
        }
        Self { keep, groups }
    }

    fn advance(&mut self, p: u64) -> ModPoly {
        let terms = self
            .groups
            .iter_mut()
            .rev()
            .filter_map(|(e, terms)| {
                let c = terms.iter_mut().fold(0, |sum, (value, ratio)| {
                    *value = mul(*value, *ratio, p);
                    add(sum, *value, p)
                });
                (c != 0).then(|| (e.clone(), c))
            })
            .collect();
        ModPoly {
            n: self.keep,
            terms,
        }
    }
}

fn image_gcd(a: &ModPoly, b: &ModPoly, p: u64, rng: &mut Rng) -> Option<ModPoly> {
    if a.n == 2 {
        return pgcd(a, b, 1, p, rng).map(|h| h.monic(p));
    }
    let dense = |f: &ModPoly| {
        let mut cs = vec![0; f.degree(0) + 1];
        for (e, c) in &f.terms {
            cs[e[0] as usize] = *c;
        }
        cs
    };
    let h = Arithmetic { p }.gcd(&dense(a), &dense(b));
    Some(ModPoly {
        n: 1,
        terms: h
            .into_iter()
            .enumerate()
            .rev()
            .filter(|(_, c)| *c != 0)
            .map(|(e, c)| (vec![e as u32], c))
            .collect(),
    })
}

/// Remove LC(Abar)*G's or LC(G)*Abar's content in the retained variables.
/// Inputs to this backend are primitive in x0; hence the desired factor is
/// primitive with respect to either retained block as well.
fn primitive(f: &ModPoly, keep: usize, p: u64, rng: &mut Rng) -> Option<ModPoly> {
    let mut groups: BTreeMap<Exps, Vec<(Exps, u64)>> = BTreeMap::new();
    for (e, c) in &f.terms {
        let mut suffix = e.clone();
        suffix[..keep].fill(0);
        groups
            .entry(e[..keep].to_vec())
            .or_default()
            .push((suffix, *c));
    }
    let mut coefficients = groups.into_values().map(|terms| ModPoly { n: f.n, terms });
    let mut content = coefficients.next()?;
    for c in coefficients {
        if content.terms.len() == 1 && content.lm().iter().all(|&e| e == 0) {
            break;
        }
        content = pgcd(&content, &c, f.n - 1, p, rng)?;
    }
    quotient(f, &content, p).map(|h| h.monic(p))
}

#[derive(Debug, PartialEq, Eq)]
enum Recovered {
    Gcd,
    Cofactor,
}

struct Recovery<'a> {
    a: &'a ModPoly,
    b: &'a ModPoly,
    keep: usize,
    p: u64,
    bounds: &'a [usize],
}

impl Recovery<'_> {
    fn certify(&self, candidate: &ModPoly) -> bool {
        candidate.degrees() == self.bounds
            && quotient(self.a, candidate, self.p).is_some()
            && quotient(self.b, candidate, self.p).is_some()
    }

    fn attempt(&self, encoding: &Encoding, rng: &mut Rng) -> Option<(ModPoly, Recovered, usize)> {
        let Self { a, b, keep, p, .. } = *self;
        let orbit = Orbit::new(encoding, p, rng)?;
        let mut aa = ImageSequence::new(a, keep, &orbit, p);
        let mut bb = ImageSequence::new(b, keep, &orbit, p);
        let (mut h, mut c) = (Interpolant::default(), Interpolant::default());
        let mut leading = None;
        let (mut checkpoint, mut previous) = (4, 2);
        for j in 1..=MAX_SAMPLES {
            // Never skip a failed image: gaps invalidate the recurrence.
            let (aj, bj) = (aa.advance(p), bb.advance(p));
            if aj.terms.is_empty()
                || bj.terms.is_empty()
                || aj.lm() != &a.lm()[..keep]
                || bj.lm() != &b.lm()[..keep]
            {
                return None;
            }
            let mut gj = image_gcd(&aj, &bj, p, rng)?;
            if leading.as_ref().is_some_and(|e| *e != gj.lm()) {
                return None;
            }
            leading = Some(gj.lm().to_vec());
            c.push(&quotient(&aj, &gj, p)?, p)?;
            gj.scale(aj.terms[0].1, p);
            h.push(&gj, p)?;
            if j != checkpoint && j != MAX_SAMPLES {
                continue;
            }
            (checkpoint, previous) = (checkpoint + previous, checkpoint);
            for (interpolant, source) in [(&h, Recovered::Gcd), (&c, Recovered::Cofactor)] {
                let Some(recovered) = interpolant.recover(&orbit, a.n) else {
                    continue;
                };
                let Some(recovered) = primitive(&recovered, keep, p, rng) else {
                    continue;
                };
                let candidate = match source {
                    Recovered::Gcd => recovered,
                    Recovered::Cofactor => {
                        let Some(g) = quotient(a, &recovered, p) else {
                            continue;
                        };
                        g.monic(p)
                    }
                };
                if self.certify(&candidate) {
                    return Some((candidate, source, j));
                }
            }
        }
        None
    }
}

/// Exactly certified finite-field recovery. Both inputs must be nonzero and
/// primitive in x0. A declined attempt is handled by the integer driver's fallback.
pub(crate) fn gcd(a: &ModPoly, b: &ModPoly, keep: usize, p: u64, rng: &mut Rng) -> Option<ModPoly> {
    let encoding = Encoding::new(
        &a.degrees()[keep..]
            .iter()
            .map(|&d| d as u32)
            .collect::<Vec<_>>(),
    )?;
    let cells = a
        .degrees()
        .into_iter()
        .chain(b.degrees())
        .try_fold(0usize, |s, d| s.checked_add(d + 1))?;
    if cells > MAX_CELLS || (0..keep).any(|k| a.degree(k).max(b.degree(k)) > MAX_IMAGE_DEGREE) {
        return None;
    }
    let bounds = degree_bounds(a, b, p, rng)?;
    let recovery = Recovery {
        a,
        b,
        keep,
        p,
        bounds: &bounds,
    };
    (0..3).find_map(|_| recovery.attempt(&encoding, rng).map(|(h, _, _)| h))
}

#[cfg(test)]
mod tests {
    use super::*;
    use polycore::{Order, Ring};

    fn recover(h: &str, u: &str, v: &str, keep: usize) -> (Recovered, usize) {
        let ring = Ring::new(["x", "y", "z", "w"], Order::Lex);
        let h = ring.parse(h).unwrap();
        let a = crate::integral(&(&h * &ring.parse(u).unwrap())).0;
        let b = crate::integral(&(&h * &ring.parse(v).unwrap())).0;
        let encoding = Encoding::new(&a.degrees()[keep..]).unwrap();
        let p = encoding.primes().next().unwrap();
        let (a, b) = (a.reduce(p), b.reduce(p));
        let mut rng = Rng::new(31);
        let bounds = degree_bounds(&a, &b, p, &mut rng).unwrap();
        let recovery = Recovery {
            a: &a,
            b: &b,
            keep,
            p,
            bounds: &bounds,
        };
        let (g, source, samples) = recovery
            .attempt(&encoding, &mut rng)
            .expect("direct recovery, no Zippel fallback");
        assert_eq!(g.terms, crate::integral(&h).0.reduce(p).monic(p).terms);
        assert_eq!(gcd(&a, &b, keep, p, &mut rng).unwrap().terms, g.terms);
        (source, samples)
    }

    #[test]
    fn unlucky_leading_coefficient_restarts_the_entire_orbit() {
        let encoding = Encoding::new(&[1]).unwrap();
        let p = encoding.primes().next().unwrap();
        let orbit = Orbit::new(&encoding, p, &mut Rng::new(43)).unwrap();
        let (value, ratio) = orbit.monomial(&[1]);
        let root = mul(value, ratio, p);
        let a = ModPoly {
            n: 2,
            terms: vec![(vec![1, 1], 1), (vec![1, 0], p - root), (vec![0, 0], 1)],
        };
        let b = ModPoly {
            n: 2,
            terms: vec![(vec![1, 1], 1), (vec![1, 0], p - root), (vec![0, 0], 2)],
        };
        let recovery = Recovery {
            a: &a,
            b: &b,
            keep: 1,
            p,
            bounds: &[0, 0],
        };
        let mut rng = Rng::new(43);
        assert!(recovery.attempt(&encoding, &mut rng).is_none());
        assert_eq!(
            recovery.attempt(&encoding, &mut rng).unwrap().0.terms,
            vec![(vec![0, 0], 1)]
        );
    }

    #[test]
    fn geometric_images_match_evaluation_and_pad_vanishing_coefficients() {
        let encoding = Encoding::new(&[1]).unwrap();
        let p = encoding.primes().next().unwrap();
        let orbit = Orbit::new(&encoding, p, &mut Rng::new(43)).unwrap();
        let (value, ratio) = orbit.monomial(&[1]);
        let root = mul(value, ratio, p);
        let f = ModPoly {
            n: 2,
            terms: vec![(vec![1, 0], 1), (vec![0, 1], 1), (vec![0, 0], p - root)],
        };
        let mut images = ImageSequence::new(&f, 1, &orbit, p);
        let mut samples = Interpolant::default();
        let mut point = value;
        for _ in 0..6 {
            point = mul(point, ratio, p);
            let image = images.advance(p);
            let expected = f.eval_var(1, point, p);
            assert_eq!(
                image.terms,
                expected
                    .terms
                    .iter()
                    .map(|(e, c)| (e[..1].to_vec(), *c))
                    .collect::<Vec<_>>()
            );
            samples.push(&image, p).unwrap();
        }
        assert_eq!(samples.recover(&orbit, 2).unwrap().terms, f.terms);
    }

    #[test]
    fn both_image_dimensions_recover_nonmonic_repeated_factors() {
        for keep in [1, 2] {
            recover(
                "((y+z)*x^2+(y^2*w+3)*x+z+1)^2",
                "(z+w)*x+y+2",
                "(y+w)*x+z+3",
                keep,
            );
        }
    }

    #[test]
    fn cofactor_recovery_finishes_before_a_dense_gcd_coefficient() {
        let tail = (1..=18)
            .map(|e| format!("z^{e}"))
            .collect::<Vec<_>>()
            .join("+");
        for keep in [1, 2] {
            let (source, samples) = recover(&format!("x+{tail}"), "x+y+w", "x+2*y+3*w+1", keep);
            assert_eq!(source, Recovered::Cofactor);
            assert!(samples <= 6);
        }
    }

    #[test]
    fn bivariate_images_reduce_the_required_coefficient_samples() {
        let tail = |variable| {
            (1..=10)
                .map(|e| format!("y^{e}*{variable}^{e}"))
                .collect::<Vec<_>>()
                .join("+")
        };
        let h = format!("x+{}", tail("z"));
        let u = format!("x+{}", tail("w"));
        let (_, univariate) = recover(&h, &u, "x+1", 1);
        let (_, bivariate) = recover(&h, &u, "x+1", 2);
        assert!(univariate >= 22);
        assert_eq!(bivariate, 4);
    }
}
