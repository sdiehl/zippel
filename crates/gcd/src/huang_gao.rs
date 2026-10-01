//! Derivative-aided separated lifting (Huang–Gao, arXiv:2609.08074v1,
//! Algorithms 1, 3–6). Signed exponents represent the monic factor's Laurent
//! coefficients. Euler derivatives `x_k d/dx_k` avoid divisions by the base point.
//!
//! This is a bounded, exactly certified implementation of the recovery method,
//! using fast NTT/Newton/half-GCD arithmetic. Verification uses exact sparse
//! division, and small fields are declined rather than
//! extended. The integer driver falls back to Zippel when a budget is exhausted.

use crate::images::degree_bounds;
use crate::poly::{Exps, ModPoly};
use polycore::fast::PrimeField;
use polycore::modp::{add, from_signed, inv, mul, pow, sub, symmetric};
use polycore::sample::Rng;
use std::collections::BTreeMap;

// Bound memory and recovery work before returning to the Zippel fallback.
const MAX_DEGREE: usize = 65_536;
const MAX_TERMS: usize = 256;
const MAX_JET_CELLS: usize = 4_000_000;

type Laurent = BTreeMap<Vec<i64>, u64>;
type Series = BTreeMap<i64, Jet>;
type Term = (Vec<i64>, u64, i64);

#[derive(Clone, Debug)]
struct Jet {
    // Evaluations at b, b^2, b^3; Euler derivatives only at b.
    values: [u64; 3],
    derivatives: Vec<u64>,
}

impl Jet {
    fn zero(n: usize) -> Self {
        Self {
            values: [0; 3],
            derivatives: vec![0; n],
        }
    }
}

fn weight(e: &[i64], s: &[i64]) -> i64 {
    e.iter().zip(s).map(|(e, s)| e * s).sum()
}

fn monomial(e: &[i64], b: &[u64], p: u64) -> u64 {
    e.iter().zip(b).fold(1, |v, (&e, &b)| {
        mul(
            v,
            pow(if e < 0 { inv(b, p) } else { b }, e.unsigned_abs(), p),
            p,
        )
    })
}

fn term_jet(e: &[i64], c: u64, b: &[u64], p: u64) -> Jet {
    let m = monomial(e, b, p);
    let values = [
        mul(c, m, p),
        mul(c, pow(m, 2, p), p),
        mul(c, pow(m, 3, p), p),
    ];
    let derivatives = e
        .iter()
        .map(|&e| mul(from_signed(e, p), values[0], p))
        .collect();
    Jet {
        values,
        derivatives,
    }
}

fn accumulate(target: &mut Jet, source: &Jet, subtract: bool, p: u64) {
    let op = if subtract { sub } else { add };
    for (a, &b) in target.values.iter_mut().zip(&source.values) {
        *a = op(*a, b, p);
    }
    for (a, &b) in target.derivatives.iter_mut().zip(&source.derivatives) {
        *a = op(*a, b, p);
    }
}

/// `Phi_s(f)`, with its y valuation removed before any dense allocation.
fn project(f: &ModPoly, s: &[i64], b: &[u64], p: u64) -> Option<Series> {
    let terms: Vec<_> = f
        .terms
        .iter()
        .map(|(e, c)| {
            let e: Vec<_> = e.iter().map(|&x| i64::from(x)).collect();
            (weight(&e, s), e, *c)
        })
        .collect();
    let low = terms.iter().map(|t| t.0).min()?;
    let high = terms.iter().map(|t| t.0).max()?;
    let degree = usize::try_from(high - low).ok()?;
    if degree > MAX_DEGREE || (degree + 1).checked_mul(f.n + 3)? > MAX_JET_CELLS {
        return None;
    }
    let mut result = Series::new();
    for (d, e, c) in terms {
        accumulate(
            result.entry(d - low).or_insert_with(|| Jet::zero(f.n)),
            &term_jet(&e, c, b, p),
            false,
            p,
        );
    }
    // Preserve the projection's support endpoints at all three points.
    if [0, high - low]
        .iter()
        .any(|d| result[d].values.contains(&0))
    {
        return None;
    }
    Some(result)
}

fn index(d: i64) -> usize {
    usize::try_from(d).expect("nonnegative bounded image degree")
}

fn evaluation(series: &Series, j: usize) -> Vec<u64> {
    let degree = index(*series.last_key_value().expect("nonempty projection").0);
    let mut cs = vec![0; degree + 1];
    for (&d, jet) in series {
        cs[index(d)] = jet.values[j];
    }
    while cs.last() == Some(&0) {
        cs.pop();
    }
    cs
}

/// Algorithms 3 and 5, steps 5–17: the derivative of the monic GCD factor
/// is rem(H^{-1} * derivative(A + c B), G), including repeated GCD factors.
fn lift(a: &Series, b: &Series, n: usize, c: u64, p: u64) -> Option<Series> {
    let field = PrimeField::new(p).expect("prime modulus");
    let aa: Vec<_> = (0..3).map(|j| evaluation(a, j)).collect();
    let bb: Vec<_> = (0..3).map(|j| evaluation(b, j)).collect();
    let gg: Vec<_> = aa.iter().zip(&bb).map(|(a, b)| field.gcd(a, b)).collect();
    let degree = gg[0].len().checked_sub(1)?;
    if gg
        .iter()
        .any(|g| g.is_empty() || g.len() != degree + 1 || g[0] == 0)
    {
        return None;
    }
    let mut result: Series = (0..=degree)
        .map(|d| {
            (
                i64::try_from(d).expect("bounded degree"),
                Jet {
                    values: std::array::from_fn(|j| gg[j][d]),
                    derivatives: vec![0; n],
                },
            )
        })
        .collect();
    if degree == 0 {
        return Some(result);
    }
    let f = field.add(&aa[0], &field.scale(&bb[0], c));
    let (h, rem) = field.divrem(&f, &gg[0])?;
    if !rem.is_empty() {
        return None;
    }
    let v = field.inverse_mod(&h, &gg[0])?;
    let fdegree = aa[0].len().max(bb[0].len()) - 1;
    for k in 0..n {
        let mut cs = vec![0; fdegree + 1];
        for (&d, jet) in a {
            cs[index(d)] = jet.derivatives[k];
        }
        for (&d, jet) in b {
            cs[index(d)] = add(cs[index(d)], mul(c, jet.derivatives[k], p), p);
        }
        let derivative = field.divrem(&field.mul(&v, &cs), &gg[0])?.1;
        for (&d, jet) in &mut result {
            jet.derivatives[k] = derivative.get(index(d)).copied().unwrap_or(0);
        }
    }
    Some(result)
}

/// Algorithms 4 and 1: reject colliding buckets, then recover signed exponents
/// from Euler derivative / coefficient. Negative exponents arise from monicization.
fn recover(series: &Series, b: &[u64], bounds: &[i64], p: u64) -> Vec<Term> {
    series
        .iter()
        .filter_map(|(&d, jet)| {
            let [a, c, z] = jet.values;
            if a == 0 || mul(a, z, p) != mul(c, c, p) {
                return None;
            }
            let ai = inv(a, p);
            let e = jet
                .derivatives
                .iter()
                .zip(bounds)
                .map(|(&v, &bound)| {
                    let e = i64::try_from(symmetric(mul(v, ai, p), p)).ok()?;
                    (e.abs() <= bound).then_some(e)
                })
                .collect::<Option<Vec<_>>>()?;
            let coeff = mul(a, inv(monomial(&e, b, p), p), p);
            // Also check that the recovered term explains all three evaluations.
            (term_jet(&e, coeff, b, p).values == jet.values).then_some((e, coeff, d))
        })
        .collect()
}

/// Align by the lexicographically lowest recovered term, which is the common
/// anchor of successive approximations. Subtract known terms before recovering
/// the residual, so terms colliding with known terms can still be recovered.
fn improve(
    known: &Laurent,
    g: &Series,
    s: &[i64],
    b: &[u64],
    bounds: &[i64],
    p: u64,
) -> Option<Laurent> {
    let recovered = recover(g, b, bounds, p);
    let (anchor, ac, ad) = recovered.iter().min_by(|a, b| a.0.cmp(&b.0))?;
    let shift: Vec<_> = anchor.iter().map(|e| -e).collect();
    let gamma = term_jet(&shift, inv(*ac, p), b, p);
    let mut residual = Series::new();
    for (&d, jet) in g {
        let mut aligned = Jet::zero(b.len());
        for j in 0..3 {
            aligned.values[j] = mul(gamma.values[j], jet.values[j], p);
        }
        for (k, &e) in shift.iter().enumerate() {
            aligned.derivatives[k] = mul(
                gamma.values[0],
                add(
                    jet.derivatives[k],
                    mul(from_signed(e, p), jet.values[0], p),
                    p,
                ),
                p,
            );
        }
        residual.insert(d - ad, aligned);
    }
    for (e, &c) in known {
        accumulate(
            residual
                .entry(weight(e, s))
                .or_insert_with(|| Jet::zero(b.len())),
            &term_jet(e, c, b, p),
            true,
            p,
        );
    }
    let mut next = known.clone();
    for (e, c, d) in recover(&residual, b, bounds, p) {
        if weight(&e, s) != d {
            return None;
        }
        let v = next.entry(e).or_default();
        *v = add(*v, c, p);
    }
    next.retain(|_, c| *c != 0);
    Some(next)
}

fn polynomial(known: &Laurent, n: usize, p: u64) -> Option<ModPoly> {
    let mins: Vec<_> = (0..n)
        .map(|k| known.keys().map(|e| e[k]).min())
        .collect::<Option<_>>()?;
    let terms = known
        .iter()
        .rev()
        .map(|(e, &c)| {
            let e = e
                .iter()
                .zip(&mins)
                .map(|(e, m)| u32::try_from(e - m).ok())
                .collect::<Option<Exps>>()?;
            Some((e, c))
        })
        .collect::<Option<Vec<_>>>()?;
    Some(ModPoly { n, terms }.monic(p))
}

fn divides(f: &ModPoly, g: &ModPoly, p: u64) -> bool {
    crate::modular_division::quotient(f, g, p).is_some()
}

/// Inputs are nonzero and have no common monomial factor (the integer driver
/// strips these first). Matching the degree upper bounds in every variable
/// proves maximality once a candidate divides both inputs. No probabilistic
/// candidate escapes unchecked.
pub(crate) fn gcd(f: &ModPoly, g: &ModPoly, p: u64, rng: &mut Rng) -> Option<ModPoly> {
    let bounds: Vec<_> = (0..f.n)
        .map(|k| i64::try_from(f.degree(k).min(g.degree(k))).expect("u32 exponent"))
        .collect();
    if bounds.iter().any(|&d| p <= 2 * d.unsigned_abs()) {
        return None;
    }
    // Adaptive term guesses with fresh approximations after unsuccessful rounds.
    // The smaller sampling range than the paper's error-budget constants is safe
    // because certification is exact; it can only increase retries/fallbacks.
    if f.degrees()
        .into_iter()
        .chain(g.degrees())
        .any(|d| d > MAX_DEGREE)
    {
        return None;
    }
    let degrees = degree_bounds(f, g, p, rng)?;
    if degrees.iter().all(|&d| d == 0) {
        return Some(ModPoly {
            n: f.n,
            terms: vec![(vec![0; f.n], 1)],
        });
    }
    for guess in [2, 4, 8, 16, 32, 64, 128, MAX_TERMS] {
        let mut known = Laurent::new();
        for _ in 0..(guess.ilog2() + 2) {
            let s: Vec<_> = (0..f.n)
                .map(|_| {
                    i64::try_from(1 + rng.next_u64() % (4 * guess as u64)).expect("bounded weight")
                })
                .collect();
            let b: Vec<_> = (0..f.n).map(|_| rng.nonzero(p)).collect();
            let (Some(a), Some(bb)) = (project(f, &s, &b, p), project(g, &s, &b, p)) else {
                continue;
            };
            let Some(lifted) = lift(&a, &bb, f.n, rng.nonzero(p), p) else {
                continue;
            };
            let Some(next) = improve(&known, &lifted, &s, &b, &bounds, p) else {
                break;
            };
            known = next;
            let Some(candidate) = polynomial(&known, f.n, p) else {
                continue;
            };
            if candidate.degrees() == degrees
                && divides(f, &candidate, p)
                && divides(g, &candidate, p)
            {
                return Some(candidate);
            }
            if known.len() > guess {
                break;
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::integral;
    use polycore::{Order, Ring};

    const P: u64 = 1_000_000_007;

    fn input(text: &str) -> ModPoly {
        let ring = Ring::new(["x", "y", "z"], Order::Lex);
        integral(&ring.parse(text).unwrap()).0.reduce(P)
    }

    #[test]
    fn recovers_nonmonic_repeated_and_laurent_factors() {
        let ring = Ring::new(["x", "y", "z"], Order::Lex);
        for src in [
            "2*x^3 + 3*y^2 + 5*z + 7",
            "x*y + x*z + y*z",
            "(x + y + z + 1)^3",
            "x^9 + y^8 + z^7 + 1",
            "x + y",
            "1",
        ] {
            let h = ring.parse(src).unwrap();
            // One cofactor shares a factor with G, so lifting A alone fails.
            let a = &h * &(&h * &ring.parse("x + 2*y + 3").unwrap());
            let b = &h * &ring.parse("y + z + 5").unwrap();
            let found = gcd(
                &integral(&a).0.reduce(P),
                &integral(&b).0.reduce(P),
                P,
                &mut Rng::new(17),
            )
            .expect(src);
            assert_eq!(found, integral(&h).0.reduce(P).monic(P), "{src}");
        }
    }

    #[test]
    fn collisions_are_deferred_and_residuals_recovered() {
        let h = input("x^2 + x*y + y^2 + 1");
        let b = [2, 3, 5];
        let s = [1, 1, 1];
        let series = project(&h, &s, &b, P).unwrap();
        let terms = recover(&series, &b, &[2, 2, 0], P);
        assert_eq!(terms, vec![(vec![0, 0, 0], 1, 0)]);
        // Subtracting known colliding terms exposes the remaining x*y term.
        let mut residual = series;
        for e in [vec![2, 0, 0], vec![0, 2, 0]] {
            accumulate(
                residual.get_mut(&2).unwrap(),
                &term_jet(&e, 1, &b, P),
                true,
                P,
            );
        }
        let terms = recover(&residual, &b, &[2, 2, 0], P);
        assert!(terms.contains(&(vec![1, 1, 0], 1, 2)));
    }

    #[test]
    fn rejects_small_characteristic_and_oversized_images() {
        let f = input("x^9 + y + 1");
        assert!(gcd(&f, &f, 7, &mut Rng::new(1)).is_none());
        assert!(project(&f, &[i64::from(u32::MAX); 3], &[2, 3, 5], P).is_none());
    }

    #[test]
    fn modular_division_rejects_false_factors() {
        assert!(divides(&input("x^2 - y^2"), &input("x + y"), P));
        assert!(!divides(&input("x^2 + y^2"), &input("x + y"), P));
    }

    #[test]
    fn successive_lifts_recover_residual_colliding_with_known_terms() {
        let h = input("1 + x + y + x^2 + z^4");
        let b = [2, 3, 5];
        let first = project(&h, &[1, 1, 1], &b, P).unwrap();
        let lifted = lift(&first, &first, 3, 7, P).unwrap();
        let known = improve(&Laurent::new(), &lifted, &[1, 1, 1], &b, &[2, 1, 4], P).unwrap();
        assert_eq!(polynomial(&known, 3, P).unwrap(), input("1 + x^2 + z^4"));
        // Here y collides with the known x^2 term. Subtraction before recovery
        // is essential; taking the union of noncolliding buckets cannot find y.
        let second = project(&h, &[1, 2, 1], &b, P).unwrap();
        let lifted = lift(&second, &second, 3, 7, P).unwrap();
        let known = improve(&known, &lifted, &[1, 2, 1], &b, &[2, 1, 4], P).unwrap();
        assert_eq!(polynomial(&known, 3, P).unwrap(), h);
    }

    #[test]
    fn rejects_bad_points_and_noncoprime_lifts() {
        assert!(project(&input("x - y"), &[1, 1, 1], &[2, 2, 5], P).is_none());
        let h = project(&input("x + y + 1"), &[1, 2, 1], &[2, 3, 5], P).unwrap();
        assert!(lift(&h, &h, 3, P - 1, P).is_none());
    }

    #[test]
    fn degree_certificate_detects_a_factor_hidden_by_equal_weights() {
        let h = input("x - y");
        let image = project(&h, &[1, 1, 1], &[2, 3, 5], P).unwrap();
        assert_eq!(evaluation(&image, 0).len(), 1);
        assert_eq!(
            degree_bounds(&h, &h, P, &mut Rng::new(2)).unwrap(),
            vec![1, 1, 0]
        );
        for seed in 0..20 {
            assert_eq!(gcd(&h, &h, P, &mut Rng::new(seed)).unwrap(), h);
        }
    }

    #[test]
    fn all_univariate_images_match_direct_evaluation() {
        let h = input("2*x^3*y^4 + 3*x*y*z + 5*z^6 + 7");
        let b = [2, 3, 5];
        let images = crate::images::univariate_images(&h, &b, P);
        for (k, image) in images.iter().enumerate() {
            assert_eq!(
                zippel_interp::poly::dense(image.clone(), P),
                h.eval_except(k, &b, P)
            );
        }
    }

    #[test]
    fn seeded_sparse_products_use_recovery_without_fallback() {
        use crate::poly::IntPoly;
        use num_bigint::BigInt;

        let mut rng = Rng::new(2026);
        for case in 0..24 {
            let n = 2 + case % 4;
            let mut terms = vec![(vec![0; n], BigInt::from(1))];
            for _ in 0..(2 + case % 12) {
                let exps = (0..n).map(|_| (rng.next_u64() % 5) as u32).collect();
                terms.push((exps, BigInt::from(rng.nonzero(101))));
            }
            let h = IntPoly::new(n, terms);
            let linear = |k, c| {
                let mut e = vec![0; n];
                e[k] = 1;
                IntPoly::new(n, [(e, BigInt::from(1)), (vec![0; n], BigInt::from(c))])
            };
            let a = h.mul(&linear(0, 2)).reduce(P);
            let b = h.mul(&linear(1, 3)).reduce(P);
            let found = gcd(&a, &b, P, &mut rng).expect("separated recovery succeeds");
            assert_eq!(found, h.reduce(P).monic(P), "case {case}");
        }
    }

    #[test]
    fn high_degree_sparse_factor_uses_fast_lifting() {
        let ring = Ring::new(["x", "y", "z"], Order::Lex);
        let h = ring.parse("x^120 + 2*y^95 + 3*z^77 + 5").unwrap();
        let a = &h * &ring.parse("x + y + 1").unwrap();
        let b = &h * &ring.parse("y + z + 2").unwrap();
        let found = gcd(
            &integral(&a).0.reduce(P),
            &integral(&b).0.reduce(P),
            P,
            &mut Rng::new(99),
        )
        .unwrap();
        assert_eq!(found, integral(&h).0.reduce(P).monic(P));
    }
}
