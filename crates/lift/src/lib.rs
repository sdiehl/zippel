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

mod crt;

use std::collections::HashMap;
use std::sync::Mutex;

use groebner::{Monomial, MonomialOrder, Polynomial, Term};
use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Zero};
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use zippel_interp::modp::{add, inv, mul, point, pow, sub};
use zippel_interp::{reconstruct, Exps, Primes, RatFunc};

pub type Poly = Polynomial<BigRational>;

/// `num / den` in lowest terms over Z: coprime integer coefficients, `den` leading positive in lex.
#[derive(Clone, Debug, PartialEq)]
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
    let mut acc: Option<Residues> = None;
    let mut candidate: Option<Vec<BigRational>> = None;
    for p in Primes::new().take(MAX_PRIMES) {
        let fitted = acc.as_ref().and_then(|a| {
            let values = fit(&f, &a.shape, n, p, seed)?;
            Some((a.shape.clone(), values))
        });
        let Some((shape, values)) = fitted.or_else(|| image(&f, n, p, seed)) else {
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
        .map(|i| reconstruct(&|x: &[u64], _| sample(x)?.get(i).copied(), n, p, seed))
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
    let samples: Vec<(&Vec<u64>, Vec<u64>)> = xs.filter_map(|x| Some((x, f(x, p)?))).collect();
    let monomial = |e: &Exps, x: &[u64]| {
        e.iter()
            .zip(x)
            .fold(1, |v, (&d, &xi)| mul(v, pow(xi, u64::from(d), p), p))
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
        let c = solve(rows, p)?;
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

/// Gauss-Jordan on the augmented square system `a`, or `None` if it is singular.
fn solve(mut a: Vec<Vec<u64>>, p: u64) -> Option<Vec<u64>> {
    let m = a.len();
    for c in 0..m {
        let r = (c..m).find(|&r| a[r][c] != 0)?;
        a.swap(c, r);
        let l = inv(a[c][c], p);
        a[c].iter_mut().for_each(|x| *x = mul(*x, l, p));
        let pivot = a[c].clone();
        for row in a.iter_mut().enumerate().filter(|t| t.0 != c).map(|t| t.1) {
            let k = row[c];
            if k != 0 {
                row.iter_mut()
                    .zip(&pivot)
                    .for_each(|(x, &y)| *x = sub(*x, mul(k, y, p), p));
            }
        }
    }
    Some(a.into_iter().map(|r| r[m]).collect())
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
            .map(|(e, c)| Term::new(c.clone(), Monomial::new(e.clone())))
            .collect();
        Polynomial::new(terms, n, MonomialOrder::Lex)
    };
    shape
        .iter()
        .map(|(ne, de)| integral(&poly(ne), &poly(de)))
        .collect()
}

/// Scales `num / den` to coprime integer coefficients.
fn integral(num: &Poly, den: &Poly) -> Fraction {
    let cs: Vec<&BigRational> = num
        .terms
        .iter()
        .chain(&den.terms)
        .map(|t| &t.coefficient)
        .collect();
    let l = cs.iter().fold(BigInt::one(), |l, c| l.lcm(c.denom()));
    let g = cs.iter().fold(BigInt::zero(), |g, c| {
        g.gcd(&(c.numer() * (&l / c.denom())))
    });
    let k = BigRational::new(l, g);
    Fraction {
        num: num.multiply_scalar(&k),
        den: den.multiply_scalar(&k),
    }
}
