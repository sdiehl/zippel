//! Rational functions over Q from black boxes modulo word-sized primes.
//!
//! Each prime gives every component as a rational function over GF(p) with a monic denominator.
//! Images with one support are combined coefficientwise by CRT, and Wang's rational number
//! reconstruction turns the residues into fractions. A candidate is accepted once a further
//! prime agrees with it. An image with fewer terms than the others came from an unlucky prime
//! and is skipped.
#![allow(
    clippy::multiple_crate_versions,
    clippy::redundant_pub_crate,
    clippy::many_single_char_names,
    clippy::missing_panics_doc,
    clippy::must_use_candidate
)]

mod crt;

use std::cell::RefCell;
use std::collections::HashMap;

use groebner::{Monomial, MonomialOrder, Polynomial, Term};
use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Zero};
use zippel_interp::modp::point;
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
    f: impl Fn(&[u64], u64) -> Option<Vec<u64>>,
    n: usize,
    seed: u64,
) -> Option<Vec<Fraction>> {
    let mut acc: Option<Residues> = None;
    let mut candidate: Option<Vec<BigRational>> = None;
    for p in Primes::new().take(MAX_PRIMES) {
        let Some((shape, values)) = image(&f, n, p, seed) else {
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
    f: &impl Fn(&[u64], u64) -> Option<Vec<u64>>,
    n: usize,
    p: u64,
    seed: u64,
) -> Option<(Shape, Vec<u64>)> {
    let memo = RefCell::new(HashMap::new());
    let sample = |x: &[u64]| -> Option<Vec<u64>> {
        memo.borrow_mut()
            .entry(x.to_vec())
            .or_insert_with(|| f(x, p))
            .clone()
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
