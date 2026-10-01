//! Huang–Monagan, arXiv:2609.10626v1, Algorithms 1 and 6.
//!
//! Substitute `x_i = p_i y^s_i` for distinct small primes, compute a primitive
//! integer univariate GCD, and recover exponents by prime valuations. As in the
//! paper's Section 6 implementation, use small-prime modular reconstruction
//! rather than the alternative Hensel lifting subroutine of Algorithm 3.
//! Exact divisibility and degree certificates replace Monte Carlo verification.

use crate::images::degree_bounds;
use crate::poly::{one, Exps, IntPoly};
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::Zero;
use polycore::modp::{is_prime, Primes};
use polycore::sample::Rng;
use std::collections::BTreeSet;

const MAX_DEGREE: u32 = 65_536;
const MAX_COEFFICIENT_BITS: u64 = 65_536;
const MAX_IMAGE_BITS: u64 = 64_000_000;
const MAX_DEGREE_CELLS: usize = 4_000_000;

fn weight(e: &[u32], s: &[u64]) -> Option<u64> {
    e.iter().zip(s).try_fold(0u64, |v, (&e, &s)| {
        v.checked_add(u64::from(e).checked_mul(s)?)
    })
}

/// Remove the symbolic y valuation, then substitute the distinct integer primes.
/// Refuse degree loss and coefficient growth before doing large allocations.
fn substitute(f: &IntPoly, s: &[u64], primes: &[u64]) -> Option<IntPoly> {
    let weights = f
        .terms
        .iter()
        .map(|(e, _)| weight(e, s))
        .collect::<Option<Vec<_>>>()?;
    let low = *weights.iter().min()?;
    let high = *weights.iter().max()?;
    let degree = u32::try_from(high - low).ok()?;
    if degree > MAX_DEGREE {
        return None;
    }
    let prime_bits: Vec<_> = primes
        .iter()
        .map(|p| u64::from(u64::BITS - p.leading_zeros()))
        .collect();
    let mut total_bits = 0u64;
    for (e, c) in &f.terms {
        let bits = e
            .iter()
            .zip(&prime_bits)
            .try_fold(c.bits(), |v, (&e, &bits)| {
                v.checked_add(u64::from(e).checked_mul(bits)?)
            })?;
        if bits > MAX_COEFFICIENT_BITS {
            return None;
        }
        total_bits = total_bits.checked_add(bits)?;
        if total_bits > MAX_IMAGE_BITS {
            return None;
        }
    }
    let terms = f.terms.iter().zip(weights).map(|((e, c), w)| {
        let coefficient = e
            .iter()
            .zip(primes)
            .fold(c.clone(), |c, (&e, &p)| c * BigInt::from(p).pow(e));
        (
            vec![u32::try_from(w - low).expect("bounded weighted degree")],
            coefficient,
        )
    });
    let result = IntPoly::new(1, terms);
    if result.is_zero() || result.degree(0) != degree || result.terms.last()?.0[0] != 0 {
        return None;
    }
    Some(result)
}

/// Algorithm 1, steps 17–25. Coefficient factors equal to a substitution prime
/// are indistinguishable from exponents here, so the weighted-support check and
/// subsequent exact certificate are essential. No general integer factoring.
fn decode(u: &IntPoly, s: &[u64], primes: &[u64], bounds: &[u32]) -> Option<IntPoly> {
    let n = primes.len();
    let mut terms = Vec::with_capacity(u.terms.len());
    let mut shift = None;
    for (y, a) in &u.terms {
        let mut c = a.clone();
        let mut e = vec![0; n];
        for (k, &prime) in primes.iter().enumerate() {
            let prime = BigInt::from(prime);
            loop {
                let (q, r) = c.div_rem(&prime);
                if !r.is_zero() {
                    break;
                }
                e[k] += 1;
                if e[k] > bounds[k] {
                    return None;
                }
                c = q;
            }
        }
        // All decoded weights must differ from the y exponents by the same
        // nonnegative valuation. This also rules out duplicated decoded terms.
        let offset = weight(&e, s)?.checked_sub(u64::from(y[0]))?;
        if shift.is_some_and(|shift| shift != offset) {
            return None;
        }
        shift = Some(offset);
        terms.push((e, c));
    }
    let result = IntPoly::new(n, terms);
    (!result.is_zero()).then(|| result.primitive())
}

/// The practical O(n T log C) prime-size choice from Section 6, with independent
/// resampling on retries. Exact certification makes unlucky choices harmless.
fn sample_primes(n: usize, guess: u64, bits: u64, rng: &mut Rng) -> Option<Vec<u64>> {
    let size = u64::try_from(n)
        .ok()?
        .checked_mul(guess)?
        .checked_mul(bits.max(2))?
        .max(32);
    let mut seen = BTreeSet::new();
    let mut result = Vec::with_capacity(n);
    for _ in 0..n.checked_mul(256)? {
        let prime = size.checked_add(rng.next_u64() % size)?.checked_add(1)? | 1;
        if is_prime(prime) && seen.insert(prime) {
            result.push(prime);
            if result.len() == n {
                return Some(result);
            }
        }
    }
    None
}

/// Degree-preserving modular images give upper bounds on the INTEGER GCD's
/// degree in every variable. A prime or specialization may enlarge these
/// bounds, but cannot make them smaller; exact division then proves maximality.
fn integer_degrees(f: &IntPoly, g: &IntPoly, rng: &mut Rng) -> Option<Exps> {
    let (fd, gd) = (f.degrees(), g.degrees());
    if fd.iter().chain(&gd).any(|&d| d > MAX_DEGREE) {
        return None;
    }
    let cells = fd
        .iter()
        .chain(&gd)
        .try_fold(0usize, |a, &d| a.checked_add(d as usize + 1))?;
    if cells > MAX_DEGREE_CELLS {
        return None;
    }
    let mut best: Option<Exps> = None;
    let mut samples = 0;
    for p in Primes::new().take(8) {
        let (a, b) = (f.reduce(p), g.reduce(p));
        if a.terms.is_empty()
            || b.terms.is_empty()
            || a.degrees().iter().zip(&fd).any(|(&a, &b)| a != b as usize)
            || b.degrees().iter().zip(&gd).any(|(&a, &b)| a != b as usize)
        {
            continue;
        }
        if let Some(degrees) = degree_bounds(&a, &b, p, rng) {
            let degrees = degrees
                .into_iter()
                .map(|d| u32::try_from(d).ok())
                .collect::<Option<Exps>>()?;
            if let Some(best) = &mut best {
                best.iter_mut()
                    .zip(degrees)
                    .for_each(|(b, d)| *b = (*b).min(d));
            } else {
                best = Some(degrees);
            }
            samples += 1;
            if samples == 3 {
                break;
            }
        }
    }
    best
}

fn certify(f: &IntPoly, g: &IntPoly, h: &IntPoly, degrees: &[u32]) -> bool {
    h.degrees() == degrees && f.div_exact(h).is_some() && g.div_exact(h).is_some()
}

/// Nonzero, monomial-primitive inputs; the caller restores common contents.
/// Bounded term guesses and image sizes let the caller fall back to Zippel.
pub(crate) fn gcd(f: &IntPoly, g: &IntPoly) -> Option<IntPoly> {
    let mut rng = Rng::new(0x4855_414e_474d);
    let degrees = integer_degrees(f, g, &mut rng)?;
    if degrees.iter().all(|&d| d == 0) {
        return Some(one(f.n));
    }
    let bits = f
        .terms
        .iter()
        .chain(&g.terms)
        .map(|(_, c)| c.bits())
        .max()?;
    for guess in [2, 4, 8, 16, 32, 64] {
        // Theorem 2.9: fully separate a T-term GCD with weights in
        // [0, 9 T(T-1)/2]. Guesses double as in Algorithm 6.
        let range = 9 * guess * (guess - 1) / 2 + 1;
        for _ in 0..2 {
            let s: Vec<_> = (0..f.n).map(|_| rng.next_u64() % range).collect();
            let primes = sample_primes(f.n, guess, bits, &mut rng)?;
            let (Some(a), Some(b)) = (substitute(f, &s, &primes), substitute(g, &s, &primes))
            else {
                continue;
            };
            let Some(u) = crate::gcd::univariate(&a, &b) else {
                continue;
            };
            let Some(h) = decode(&u, &s, &primes, &degrees) else {
                continue;
            };
            if certify(f, g, &h, &degrees) {
                return Some(h);
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

    fn input(src: &str) -> IntPoly {
        integral(&Ring::new(["x", "y", "z"], Order::Lex).parse(src).unwrap()).0
    }

    #[test]
    fn prime_valuations_recover_signed_coefficients_and_shifted_support() {
        let h = input("4*x^2*y^3 - 5*x^5*y^6 + 10*x^4*y^4");
        let s = [2, 3, 1];
        let primes = [7, 11, 13];
        let u = substitute(&h, &s, &primes).unwrap();
        assert_eq!(
            u.terms.iter().map(|t| t.0[0]).collect::<Vec<_>>(),
            vec![15, 7, 0]
        );
        assert_eq!(decode(&u, &s, &primes, &[5, 6, 0]).unwrap(), h.primitive());
    }

    #[test]
    fn rejects_coefficients_containing_substitution_primes() {
        let h = input("7*x + y + 1");
        let s = [2, 3, 1];
        let u = substitute(&h, &s, &[7, 11, 13]).unwrap();
        assert!(decode(&u, &s, &[7, 11, 13], &[3, 3, 0]).is_none());
    }

    #[test]
    fn rejects_lost_endpoints_and_oversized_images() {
        assert!(substitute(&input("5*x - 3*y + 1"), &[1, 1, 1], &[3, 5, 7]).is_none());
        assert!(substitute(&input("x^2 + 1"), &[u64::MAX, 1, 1], &[3, 5, 7]).is_none());
        assert!(substitute(&input("x^65537 + 1"), &[1, 1, 1], &[3, 5, 7]).is_none());
    }

    #[test]
    fn rejects_proper_divisor_from_paper_example_3_4() {
        let f = input("(5*x - 3*y + 1)*(2*x*y + 7)");
        let h = input("2*x*y + 7");
        let degrees = integer_degrees(&f, &f, &mut Rng::new(12)).unwrap();
        assert!(f.div_exact(&h).is_some());
        assert!(!certify(&f, &f, &h, &degrees));
    }

    #[test]
    fn computes_repeated_nonmonic_and_sparse_factors_without_fallback() {
        for src in [
            "2*x^3 + 3*y^2 + 5*z + 7",
            "x*y + x*z + y*z",
            "(x + y + 1)^2",
            "x^12 + y^9 + z^7 + 1",
            "x - y",
            "1",
        ] {
            let h = input(src);
            let a = h.mul(&input("x + 2*y + 3"));
            let b = h.mul(&input("y + z + 5"));
            assert_eq!(gcd(&a, &b).expect(src), h.primitive(), "{src}");
        }
    }

    #[test]
    fn sampled_primes_are_distinct() {
        let primes = sample_primes(16, 2, 10, &mut Rng::new(4)).unwrap();
        assert!(primes.iter().all(|&p| is_prime(p)));
        assert_eq!(primes.iter().collect::<BTreeSet<_>>().len(), 16);
    }

    #[test]
    fn seeded_products_recover_without_fallback() {
        let mut rng = Rng::new(2026);
        for case in 0..20 {
            let mut terms = vec![(vec![0; 3], BigInt::from(1))];
            for _ in 0..(2 + case % 8) {
                let e = (0..3).map(|_| (rng.next_u64() % 4) as u32).collect();
                terms.push((e, BigInt::from(rng.nonzero(101))));
            }
            let h = IntPoly::new(3, terms).primitive();
            let a = h.mul(&input("x + 2"));
            let b = h.mul(&input("y + 3"));
            assert_eq!(
                gcd(&a, &b).expect("substitution recovery"),
                h,
                "case {case}"
            );
        }
    }

    #[test]
    fn univariate_crt_recovers_missing_terms_after_unlucky_primes() {
        let p = Primes::new().next().unwrap();
        let ring = Ring::new(["y"], Order::Lex);
        let h = integral(&ring.parse(&format!("y^2 + {p}*y + 1")).unwrap()).0;
        let a = h.mul(&integral(&ring.parse("y + 2").unwrap()).0);
        let b = h.mul(&integral(&ring.parse(&format!("y + {}", p + 2)).unwrap()).0);
        // Modulo the first prime the cofactors coincide and h loses a term.
        assert_eq!(crate::gcd::univariate(&a, &b).unwrap(), h);
    }

    #[test]
    fn degree_certificates_discard_unlucky_modular_factors() {
        let p = Primes::new().next().unwrap();
        let a = input("(x + y + 1)*(x + 2)");
        let b = input(&format!("(x + y + 1)*(x + {})", p + 2));
        assert_eq!(
            integer_degrees(&a, &b, &mut Rng::new(3)).unwrap(),
            vec![1, 1, 0]
        );
        assert_eq!(gcd(&a, &b).unwrap(), input("x + y + 1"));
    }
}
