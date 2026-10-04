//! Denominator candidates learned on independent slices. Guesses are Monte Carlo;
//! reconstruction is checked at fresh points and falls back to ordinary reconstruction.

use crate::rational::thiele;
use crate::{interpolate, reconstruct, ModPoly, RatFunc};
use polycore::modp::{add, inv, mul};
use polycore::sample::{hash, point, BlackBox};
use polycore::{Fp, Modular, Uni};

/// Multiplicities in the caller's candidate order. Duplicate/overlapping candidates
/// consume the slice denominator in that order. `complete` is a polynomial hint.
#[derive(Clone, Debug)]
pub struct FactorGuess {
    pub powers: Vec<u32>,
    pub complete: bool,
}

/// Estimate denominator multiplicities on three independent shifted lines.
///
/// Candidates must be nonconstant, nonzero polynomials in exactly `n` variables.
/// A failed/degree-dropping slice is retried. No black-box identity is certified here.
pub fn guess(
    f: &(impl BlackBox + Sync),
    candidates: &[ModPoly],
    n: usize,
    p: u64,
    seed: u64,
) -> Option<FactorGuess> {
    if candidates.iter().any(|g| {
        g.n != n
            || g.terms.is_empty()
            || g.terms.iter().any(|(e, c)| e.len() != n || *c >= p)
            || g.terms.iter().all(|(e, _)| e.iter().all(|&d| d == 0))
    }) {
        return None;
    }
    let polys: Vec<_> = candidates.iter().map(|g| g.to_poly(p)).collect();
    if polys.iter().any(|g| g.is_zero() || g.is_constant()) {
        return None;
    }
    let mut powers = vec![u32::MAX; candidates.len()];
    let mut complete = true;
    let mut accepted = 0;
    for attempt in 0..12 {
        let s: Vec<_> = (0..n as u64)
            .map(|j| point(&[seed, 80, attempt, j], p))
            .collect();
        let z: Vec<_> = (0..n as u64)
            .map(|j| point(&[seed, 81, attempt, j], p))
            .collect();
        let fp = |x: &u64| Fp::new(*x, p);
        let lines: Vec<_> = polys
            .iter()
            .map(|g| {
                g.on_line(
                    &z.iter().map(fp).collect::<Vec<_>>(),
                    &s.iter().map(fp).collect::<Vec<_>>(),
                )
            })
            .collect();
        if lines
            .iter()
            .zip(&polys)
            .any(|(l, g)| l.0.len() != g.total_degree() as usize + 1)
        {
            continue;
        }
        let Some((_, mut den)) = thiele(
            |t| {
                let x: Vec<_> = s
                    .iter()
                    .zip(&z)
                    .map(|(&s, &z)| add(s, mul(t, z, p), p))
                    .collect();
                f.eval(&x, p)
            },
            p,
            hash(&[seed, 82, attempt]),
        ) else {
            continue;
        };
        for ((power, line), g) in powers.iter_mut().zip(&lines).zip(&polys) {
            let mut count = 0;
            while den.deg() >= g.total_degree() as usize {
                let (q, r) = den.divrem(line);
                if !r.0.is_empty() {
                    break;
                }
                den = q;
                count += 1;
            }
            *power = (*power).min(count);
        }
        complete &= den.deg() == 0;
        accepted += 1;
        if accepted == 3 {
            return Some(FactorGuess { powers, complete });
        }
    }
    None
}

/// Discover factors depending on a single variable. Factors must recur unchanged
/// on three independent coordinate slices; mixed factors belong in a supplied pool.
pub fn univariate_factors(
    f: &(impl BlackBox + Sync),
    n: usize,
    p: u64,
    seed: u64,
) -> Option<Vec<ModPoly>> {
    let mut out = Vec::new();
    for k in 0..n {
        let mut common: Option<Vec<Uni<Fp>>> = None;
        for slice in 0..3 {
            let mut x: Vec<_> = (0..n as u64)
                .map(|j| point(&[seed, 83, k as u64, slice, j], p))
                .collect();
            let (_, den) = thiele(
                |t| {
                    x[k] = t;
                    f.eval(&x, p)
                },
                p,
                hash(&[seed, 84, k as u64, slice]),
            )?;
            let factors: Vec<_> = polyfactor::factor_mod(&den)
                .1
                .into_iter()
                .map(|(g, _)| g)
                .collect();
            common = Some(common.map_or_else(
                || factors.clone(),
                |old| old.into_iter().filter(|g| factors.contains(g)).collect(),
            ));
        }
        for g in common.unwrap_or_default() {
            let terms =
                g.0.iter()
                    .enumerate()
                    .rev()
                    .filter_map(|(i, c)| {
                        let c = c.residue_mod(p);
                        if c == 0 {
                            return None;
                        }
                        let mut e = vec![0; n];
                        e[k] = i as u32;
                        Some((e, c))
                    })
                    .collect();
            out.push(ModPoly { n, terms });
        }
    }
    Some(out)
}

/// Reconstruct using a candidate pool, returning the usual expanded, monic fraction.
/// Repeated factors are retained; erroneous overestimates are cancelled exactly.
pub fn reconstruct_with_factors(
    f: &(impl BlackBox + Sync),
    candidates: &[ModPoly],
    n: usize,
    p: u64,
    seed: u64,
) -> Option<RatFunc> {
    factored(f, candidates, n, p, seed).or_else(|| reconstruct(f, n, p, seed))
}

fn factored(
    f: &(impl BlackBox + Sync),
    candidates: &[ModPoly],
    n: usize,
    p: u64,
    seed: u64,
) -> Option<RatFunc> {
    let guessed = guess(f, candidates, n, p, seed)?;
    let mut known = ModPoly {
        n,
        terms: vec![(vec![0; n], 1)],
    }
    .to_poly(p);
    for (g, &m) in candidates.iter().zip(&guessed.powers) {
        known = &known * &g.to_poly(p).pow(m);
    }
    let known_mod = ModPoly::from_poly(&known, p);
    let scaled = |x: &[u64], p| {
        let d = known_mod.eval(x, p);
        if d == 0 {
            None
        } else {
            Some(mul(f.eval(x, p)?, d, p))
        }
    };
    let r = if guessed.complete {
        interpolate(&scaled, n, p, hash(&[seed, 85])).map(|num| RatFunc {
            num,
            den: ModPoly {
                n,
                terms: vec![(vec![0; n], 1)],
            },
        })
    } else {
        None
    }
    .or_else(|| reconstruct(&scaled, n, p, hash(&[seed, 86])))?;
    let mut num = r.num.to_poly(p);
    let mut den = &r.den.to_poly(p) * &known;
    for g in candidates {
        let g = g.to_poly(p);
        loop {
            if num.is_zero() {
                break;
            }
            let Some(a) = num.exact(&g) else { break };
            let Some(b) = den.exact(&g) else { break };
            num = a;
            den = b;
        }
    }
    if num.is_zero() {
        den = ModPoly {
            n,
            terms: vec![(vec![0; n], 1)],
        }
        .to_poly(p);
    }
    let mut result = RatFunc {
        num: ModPoly::from_poly(&num, p),
        den: ModPoly::from_poly(&den, p),
    };
    let l = inv(result.den.terms.first()?.1, p);
    result.num.scale(l, p);
    result.den.scale(l, p);
    let mut checked = 0;
    for i in 0..24 {
        let x: Vec<_> = (0..n as u64).map(|j| point(&[seed, 87, i, j], p)).collect();
        let Some(y) = f.eval(&x, p) else { continue };
        let d = result.den.eval(&x, p);
        if d == 0 {
            continue;
        }
        if mul(y, d, p) != result.num.eval(&x, p) {
            return None;
        }
        checked += 1;
        if checked == 3 {
            return Some(result);
        }
    }
    None
}
