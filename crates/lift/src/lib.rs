//! Rational functions over Q from black boxes modulo word-sized primes.
//!
//! Each prime gives every component as a rational function over GF(p) with a monic denominator.
//! Images with one support are combined coefficientwise by CRT, and Wang's rational number
//! reconstruction turns the residues into fractions. A candidate is accepted once a further
//! prime agrees with it. An image with fewer terms than the others came from an unlucky prime
//! and is skipped. Only the first prime pays for reconstruction: later ones fit coefficients to
//! the support already known, from as many samples as the largest component has terms.
#![allow(
    clippy::multiple_crate_versions,
    clippy::redundant_pub_crate,
    clippy::many_single_char_names,
    clippy::missing_panics_doc,
    clippy::must_use_candidate
)]

use std::collections::HashMap;
use std::sync::{Arc, Mutex};

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Zero};
use polycore::evaluation::power_table;
use polycore::modp::{add, mul, sub, Primes};
use polycore::sample::point;
use polycore::{crt, dense, Fp, Modular};
use polycore::{Monomial, Order};
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use zippel_interp::{reconstruct, Exps, RatFunc};

pub type Poly = polycore::Poly<BigRational>;

/// `num / den` in lowest terms over Z: coprime integer coefficients, `den` leading positive in lex.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Fraction {
    pub num: Poly,
    pub den: Poly,
}

type Shape = Vec<(Vec<Exps>, Vec<Exps>)>;

const MAX_PRIMES: usize = 64;

/// Every component of the vector black box `f` in `n` variables as a fraction over Q. `f` must
/// work modulo any prime it is given and may refuse points. Monte Carlo, like
/// [`zippel_interp::reconstruct`].
pub fn lift(
    f: impl Fn(&[u64], u64) -> Option<Vec<u64>> + Sync,
    n: usize,
    seed: u64,
) -> Option<Vec<Fraction>> {
    lift_impl(&f, n, seed, &[])
}

fn lift_impl(
    f: &(impl Fn(&[u64], u64) -> Option<Vec<u64>> + Sync),
    n: usize,
    seed: u64,
    polynomial: &[bool],
) -> Option<Vec<Fraction>> {
    let mut acc: Option<Residues> = None;
    let mut candidate: Option<Vec<BigRational>> = None;
    for p in Primes::new().take(MAX_PRIMES) {
        let fitted = acc.as_ref().and_then(|a| {
            let values = fit(f, &a.shape, n, p, seed)?;
            Some((a.shape.clone(), values))
        });
        let Some((shape, values)) = fitted.or_else(|| image(f, n, p, seed, polynomial)) else {
            continue;
        };
        if let (Some(a), Some(q)) = (&acc, &candidate) {
            if a.shape == shape
                && q.iter()
                    .zip(&values)
                    .all(|(r, &v)| crt::reduce(r, p) == Some(v))
            {
                return Some(fractions(n, &shape, q));
            }
        }
        match &mut acc {
            Some(a) if a.shape == shape => crt::garner(&mut a.xs, &mut a.m, &values, p),
            Some(a) if terms(&a.shape) > terms(&shape) => continue,
            _ => acc = Some(Residues::new(shape, &values, p)),
        }
        candidate = acc.as_ref().and_then(Residues::rationals);
    }
    None
}

/// All components modulo `p`, sharing one memo so that each point is sampled once.
fn image(
    f: &(impl Fn(&[u64], u64) -> Option<Vec<u64>> + Sync),
    n: usize,
    p: u64,
    seed: u64,
    polynomial: &[bool],
) -> Option<(Shape, Vec<u64>)> {
    let memo = Mutex::new(HashMap::new());
    let sample = |x: &[u64]| -> Option<Vec<u64>> {
        if let Some(y) = memo.lock().unwrap().get(x) {
            return Option::clone(y);
        }
        let y = f(x, p);
        memo.lock().unwrap().insert(x.to_vec(), y.clone());
        y
    };
    let len = (0..16)
        .find_map(|i| {
            sample(
                &(0..n as u64)
                    .map(|j| point(&[seed, 5, i, j], p))
                    .collect::<Vec<_>>(),
            )
        })?
        .len();
    let images: Vec<RatFunc> = (0..len)
        .map(|i| {
            let component = |x: &[u64], _| sample(x)?.get(i).copied();
            let poly = polynomial
                .get(i)
                .copied()
                .unwrap_or(false)
                .then(|| {
                    zippel_interp::interpolate(&component, n, p, seed).map(|num| RatFunc {
                        num,
                        den: zippel_interp::ModPoly {
                            n,
                            terms: vec![(vec![0; n], 1)],
                        },
                    })
                })
                .flatten();
            poly.or_else(|| reconstruct(&component, n, p, seed))
        })
        .collect::<Option<_>>()?;
    let shape = images
        .iter()
        .map(|r| (r.num.support(), r.den.support()))
        .collect();
    let values = images
        .iter()
        .flat_map(|r| r.num.terms.iter().chain(&r.den.terms).map(|t| t.1))
        .collect();
    Some((shape, values))
}

/// The coefficients on a known support modulo `p`, one linear solve per component through
/// samples taken in parallel. `None` if a solve is singular or misses a check point.
fn fit(
    f: &(impl Fn(&[u64], u64) -> Option<Vec<u64>> + Sync),
    shape: &Shape,
    n: usize,
    p: u64,
    seed: u64,
) -> Option<Vec<u64>> {
    let need = shape.iter().map(|(a, b)| a.len() + b.len()).max()?;
    let xs: Vec<Vec<u64>> = (0..need as u64 + need as u64 / 8 + 2)
        .map(|i| (0..n as u64).map(|j| point(&[seed, 6, i, j], p)).collect())
        .collect();
    #[cfg(feature = "parallel")]
    let xs = xs.par_iter();
    #[cfg(not(feature = "parallel"))]
    let xs = xs.iter();
    let mut degrees = vec![0; n];
    for e in shape.iter().flat_map(|(a, b)| a.iter().chain(b)) {
        for (d, &x) in degrees.iter_mut().zip(e) {
            *d = (*d).max(x as usize);
        }
    }
    let samples: Vec<(Vec<Vec<u64>>, Vec<u64>)> = xs
        .filter_map(|x| Some((power_table(x, &degrees, p), f(x, p)?)))
        .collect();
    let monomial = |e: &Exps, x: &[Vec<u64>]| {
        e.iter()
            .zip(x)
            .fold(1, |v, (&d, xi)| mul(v, xi[d as usize], p))
    };
    let mut values = Vec::new();
    for (i, (num, den)) in shape.iter().enumerate() {
        let unknowns = num.len() + den.len() - 1;
        let mut rows: Vec<Vec<u64>> = samples
            .iter()
            .take(unknowns + 1)
            .map(|(x, y)| {
                let y = *y.get(i)?;
                let mut row: Vec<u64> = num.iter().map(|e| monomial(e, x)).collect();
                row.extend(
                    den[1..]
                        .iter()
                        .map(|e| mul(sub(0, y, p), monomial(e, x), p)),
                );
                row.push(mul(y, monomial(&den[0], x), p));
                Some(row)
            })
            .collect::<Option<_>>()?;
        let check = rows.pop().filter(|_| rows.len() == unknowns)?;
        let fp = |v: &u64| Fp::new(*v, p);
        let a: Vec<Vec<Fp>> = rows
            .iter()
            .map(|r| r[..unknowns].iter().map(fp).collect())
            .collect();
        let b: Vec<Fp> = rows.iter().map(|r| fp(&r[unknowns])).collect();
        let c: Vec<u64> = dense::solve(&a, &b)
            .ok()?
            .iter()
            .map(|x| x.residue_mod(p))
            .collect();
        let lhs = c
            .iter()
            .zip(&check)
            .fold(0, |acc, (&cj, &a)| add(acc, mul(cj, a, p), p));
        if lhs != check[unknowns] {
            return None;
        }
        values.extend(&c[..num.len()]);
        values.push(1);
        values.extend(&c[num.len()..]);
    }
    Some(values)
}

fn terms(shape: &Shape) -> usize {
    shape.iter().map(|(a, b)| a.len() + b.len()).sum()
}

#[derive(Debug)]
struct Residues {
    shape: Shape,
    xs: Vec<BigInt>,
    m: BigInt,
}

impl Residues {
    fn new(shape: Shape, values: &[u64], p: u64) -> Self {
        Self {
            shape,
            xs: values.iter().map(|&v| BigInt::from(v)).collect(),
            m: BigInt::from(p),
        }
    }

    fn rationals(&self) -> Option<Vec<BigRational>> {
        self.xs.iter().map(|x| crt::wang(x, &self.m)).collect()
    }
}

fn fractions(n: usize, shape: &Shape, q: &[BigRational]) -> Vec<Fraction> {
    let mut q = q.iter();
    let mut poly = |es: &[Exps]| {
        let terms = es
            .iter()
            .zip(q.by_ref())
            .map(|(e, c)| (Monomial::new(e.clone()), c.clone()))
            .collect();
        Poly::new(terms, n, Order::Lex)
    };
    shape
        .iter()
        .map(|(ne, de)| integral(&poly(ne), &poly(de)))
        .collect()
}

/// Scales `num / den` to coprime integer coefficients.
fn integral(num: &Poly, den: &Poly) -> Fraction {
    let cs: Vec<&BigRational> = num.terms.iter().chain(&den.terms).map(|t| &t.1).collect();
    let l = cs.iter().fold(BigInt::one(), |l, c| l.lcm(c.denom()));
    let g = cs.iter().fold(BigInt::zero(), |g, c| {
        g.gcd(&(c.numer() * (&l / c.denom())))
    });
    let k = BigRational::new(l, g);
    Fraction {
        num: num.scale(&k),
        den: den.scale(&k),
    }
}

/// Lift with denominator candidates over Q, shared by all output components.
///
/// Three slices select multiplicities independently per component. Subsequent
/// primes fit the residual fraction, keeping known factors out of the linear solve.
/// Invalid candidates or unsuccessful learning fall back to ordinary lifting.
pub fn lift_with_factors(
    f: impl Fn(&[u64], u64) -> Option<Vec<u64>> + Sync,
    n: usize,
    candidates: &[Poly],
    seed: u64,
) -> Option<Vec<Fraction>> {
    lift_factored(&f, n, candidates, seed).or_else(|| lift(&f, n, seed))
}

fn lift_factored(
    f: &(impl Fn(&[u64], u64) -> Option<Vec<u64>> + Sync),
    n: usize,
    candidates: &[Poly],
    seed: u64,
) -> Option<Vec<Fraction>> {
    use zippel_interp::{factors::guess, ModPoly};
    if candidates.iter().any(|g| g.nvars != n || g.is_constant()) {
        return None;
    }
    let p = Primes::new().next()?;
    let modular: Vec<_> = candidates
        .iter()
        .map(|g| {
            let g = g.try_map(|c| Some(Fp::new(crt::reduce(c, p)?, p)))?;
            Some(ModPoly::from_poly(&g, p))
        })
        .collect::<Option<_>>()?;
    // Slice probes are shared by all coefficients, just as in ordinary lifting.
    let memo = Mutex::new(HashMap::new());
    let sample = |x: &[u64]| {
        if let Some(y) = memo.lock().unwrap().get(x) {
            return Option::clone(y);
        }
        let y = f(x, p);
        memo.lock().unwrap().insert(x.to_vec(), y.clone());
        y
    };
    let len = (0..16)
        .find_map(|i| {
            sample(
                &(0..n as u64)
                    .map(|j| point(&[seed, 90, i, j], p))
                    .collect::<Vec<_>>(),
            )
        })?
        .len();
    let guesses: Vec<_> = (0..len)
        .map(|i| {
            guess(
                &|x: &[u64], _| sample(x)?.get(i).copied(),
                &modular,
                n,
                p,
                seed,
            )
        })
        .collect::<Option<_>>()?;
    let known: Vec<_> = guesses
        .iter()
        .map(|g| {
            candidates.iter().zip(&g.powers).fold(
                Poly::constant(BigRational::one(), n, Order::Lex),
                |acc, (f, &m)| &acc * &f.pow(m),
            )
        })
        .collect();
    let polynomial: Vec<_> = guesses.iter().map(|g| g.complete).collect();
    // Convert the candidate pool once per prime; evaluate the factored product
    // without expanding it or repeatedly reducing rational coefficients.
    let pools = Mutex::new(HashMap::new());
    let scaled = |x: &[u64], prime| {
        let pool = {
            let mut cache = pools.lock().unwrap();
            cache
                .entry(prime)
                .or_insert_with(|| {
                    candidates
                        .iter()
                        .map(|g| {
                            let g = g.try_map(|c| Some(Fp::new(crt::reduce(c, prime)?, prime)))?;
                            Some(ModPoly::from_poly(&g, prime))
                        })
                        .collect::<Option<Vec<_>>>()
                        .map(Arc::new)
                })
                .clone()?
        };
        let factors: Vec<_> = pool.iter().map(|g| g.eval(x, prime)).collect();
        let y = if prime == p { sample(x)? } else { f(x, prime)? };
        if y.len() != len {
            return None;
        }
        y.iter()
            .zip(&guesses)
            .map(|(&v, g)| {
                let d = factors.iter().zip(&g.powers).fold(1, |acc, (&c, &m)| {
                    mul(acc, polycore::modp::pow(c, u64::from(m), prime), prime)
                });
                (d != 0).then(|| mul(v, d, prime))
            })
            .collect::<Option<Vec<_>>>()
    };
    let residual = lift_impl(&scaled, n, seed, &polynomial)?;
    Some(restore_factors(residual, &known, candidates, n))
}

fn restore_factors(
    residual: Vec<Fraction>,
    known: &[Poly],
    candidates: &[Poly],
    n: usize,
) -> Vec<Fraction> {
    let out: Vec<_> = residual
        .into_iter()
        .zip(known)
        .map(|(r, d)| {
            let mut num = r.num;
            let mut den = &r.den * d;
            for g in candidates {
                while !num.is_zero() {
                    let Some(a) = num.exact(g) else { break };
                    let Some(b) = den.exact(g) else { break };
                    num = a;
                    den = b;
                }
            }
            if num.is_zero() {
                den = Poly::constant(BigRational::one(), n, Order::Lex);
            }
            if den.lc().is_some_and(|c| c < &BigRational::zero()) {
                num = num.scale(&-BigRational::one());
                den = den.scale(&-BigRational::one());
            }
            integral(&num, &den)
        })
        .collect();
    out
}

/// Discover factors depending on one variable and combine them with supplied
/// multivariate candidates before lifting.
///
/// Discovery uses `polyfactor` on three
/// independent slices per coordinate and shares vector probes across components.
/// Small rational factor coefficients are guessed with Wang reconstruction;
/// unsuccessful guesses simply leave a residual denominator to reconstruct.
/// This optional preprocessing is most useful for expensive, factor-rich oracles.
pub fn lift_with_discovered_factors(
    f: impl Fn(&[u64], u64) -> Option<Vec<u64>> + Sync,
    n: usize,
    candidates: &[Poly],
    seed: u64,
) -> Option<Vec<Fraction>> {
    let p = Primes::new().next()?;
    let memo = Mutex::new(HashMap::new());
    let sample = |x: &[u64], prime| {
        if prime != p {
            return f(x, prime);
        }
        if let Some(y) = memo.lock().unwrap().get(x) {
            return Option::clone(y);
        }
        let y = f(x, p);
        memo.lock().unwrap().insert(x.to_vec(), y.clone());
        y
    };
    let initial = (0..16).find_map(|i| {
        sample(
            &(0..n as u64)
                .map(|j| point(&[seed, 85, i, j], p))
                .collect::<Vec<_>>(),
            p,
        )
    });
    let Some(initial) = initial else {
        return lift_with_factors(&f, n, candidates, seed);
    };
    let len = initial.len();
    let mut pool = candidates.to_vec();
    let modulus = BigInt::from(p);
    for component in 0..len {
        let oracle = |x: &[u64], p| sample(x, p)?.get(component).copied();
        for factor in
            zippel_interp::factors::univariate_factors(&oracle, n, p, seed).unwrap_or_default()
        {
            let terms = factor
                .terms
                .iter()
                .map(|(e, c)| {
                    Some((
                        Monomial::new(e.clone()),
                        crt::wang(&BigInt::from(*c), &modulus)?,
                    ))
                })
                .collect::<Option<Vec<_>>>();
            if let Some(terms) = terms {
                let g = Poly::new(terms, n, Order::Lex);
                if !pool.contains(&g) {
                    pool.push(g);
                }
            }
        }
    }
    lift_with_factors(sample, n, &pool, seed)
}
