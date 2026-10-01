//! Sparse multivariate polynomial GCD using Huang–Gao derivative-aided separated
//! Hensel lifting, Huang–Monagan prime substitution and Hu–Monagan geometric
//! interpolation, with a Zippel/LINZIP fallback.
//!
//! The default path uses weighted projections, one derivative lift per variable,
//! collision detection and iterative sparse recovery. Fast NTT convolution,
//! Newton division and half-GCD support the univariate computations. CRT and
//! rational reconstruction lift modular images to integers; exact division
//! certifies the result. [`GcdAlgorithm`] selects the preferred backend.
//!
//! The recovery backends use bounded retries and exact verification, rather than
//! the papers' Monte Carlo verification, and do not claim their full asymptotic
//! bit-complexity bounds. See [`GcdAlgorithm::HuangMonagan`] for substitution limits.
#![allow(
    clippy::multiple_crate_versions,
    clippy::redundant_pub_crate,
    clippy::cast_possible_truncation,
    clippy::missing_panics_doc,
    clippy::many_single_char_names,
    clippy::similar_names,
    clippy::module_name_repetitions
)]

mod fast;
mod frac;
mod gcd;
mod geometric;
mod hu_monagan;
mod huang_gao;
mod huang_monagan;
mod images;
mod linzip;
mod modular_division;
mod pgcd;
mod poly;

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::One;
use poly::IntPoly;
use polycore::Monomial;

pub use frac::Frac;

/// Preferred sparse GCD algorithm. Recovery backends use exact certification
/// and fall back to Zippel if their recovery budgets are exhausted.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum GcdAlgorithm {
    /// Derivative-aided separated Hensel lifting, then Huang–Monagan, then Zippel.
    #[default]
    HuangGao,
    /// Integer prime substitution and valuation decoding, then Zippel.
    ///
    /// Follows Algorithms 1 and 6 of Huang–Monagan,
    /// <https://arxiv.org/abs/2609.10626v1>, using the practical modular univariate
    /// GCD and small-prime sampling of Section 6. Trial division by the chosen
    /// primes recovers exponents; exact division and degree certificates reject
    /// unlucky substitutions.
    ///
    /// Tries two substitutions per term-bound guess from 2 through 64. Images
    /// are limited to degree 65,536, estimated coefficient size 65,536 bits and
    /// 64 million coefficient bits in total. Degree certificates use at most
    /// four million coefficient slots; each univariate CRT uses at most 256
    /// primes. Exhausting a budget invokes the Zippel fallback.
    HuangMonagan,
    /// Hu–Monagan geometric interpolation of a GCD and cofactor together.
    /// Uses univariate images, smooth-subgroup discrete logarithms and
    /// Ben-Or/Tiwari recovery, then exact certification and a Zippel fallback.
    ///
    /// Both Hu–Monagan variants allow three orbits of at most 256 samples per
    /// prime, 256 CRT primes, a 48-bit exponent box, retained degrees up to
    /// 4,096 and two million coefficient slots per interpolant. Exact modular
    /// division permits 100,000 term reductions. Exceeding a budget falls back.
    HuMonagan,
    /// Hu–Monagan simultaneous GCD/cofactor recovery using bivariate images.
    /// Keeps two variables symbolic; otherwise uses the same recovery and
    /// budgets as [`Self::HuMonagan`]. In one variable, uses univariate images.
    HuMonaganBivariate,
    /// Recursive modular interpolation with LINZIP support reuse.
    Zippel,
}

type Poly = polycore::Poly<BigRational>;

/// `f = F / d` with `F` integral.
fn integral(f: &Poly) -> (IntPoly, BigInt) {
    let d = f
        .terms
        .iter()
        .fold(BigInt::one(), |d, (_, c)| d.lcm(c.denom()));
    let terms = f
        .terms
        .iter()
        .map(|(m, c)| (m.exps().to_vec(), c.numer() * (&d / c.denom())));
    (IntPoly::new(f.nvars, terms), d)
}

fn rational(h: &IntPoly, d: &BigInt, like: &Poly) -> Poly {
    let terms = h.terms.iter().map(|(e, c)| {
        (
            Monomial::new(e.clone()),
            BigRational::new(c.clone(), d.clone()),
        )
    });
    Poly::new(terms.collect(), like.nvars, like.order.clone())
}

/// Primitive over Z with a positive leading coefficient in `like`'s monomial order.
fn normalized(h: &IntPoly, like: &Poly) -> Poly {
    rational(h, &BigInt::one(), like).primitive()
}

fn gcd_int(f: &Poly, g: &Poly, algorithm: GcdAlgorithm) -> IntPoly {
    assert_eq!(f.nvars, g.nvars, "polynomials must share a ring");
    gcd::gcd_z(&integral(f).0, &integral(g).0, algorithm)
}

/// The greatest common divisor, as a primitive integer polynomial with positive leading
/// coefficient. `gcd(0, 0) = 0`.
#[must_use]
pub fn gcd(f: &Poly, g: &Poly) -> Poly {
    gcd_with_algorithm(f, g, GcdAlgorithm::default())
}

/// Hu–Monagan GCD using univariate images and simultaneous cofactor recovery.
///
/// Has the normalization and zero conventions of [`gcd`], with a certified
/// Zippel fallback when sparse recovery exceeds its budgets.
#[must_use]
pub fn gcd_hu_monagan(f: &Poly, g: &Poly) -> Poly {
    gcd_with_algorithm(f, g, GcdAlgorithm::HuMonagan)
}

/// Hu–Monagan GCD retaining two symbolic variables in each image.
///
/// Uses univariate images for a one-variable ring. See [`gcd_hu_monagan`]
/// for normalization and fallback behavior.
#[must_use]
pub fn gcd_hu_monagan_bivariate(f: &Poly, g: &Poly) -> Poly {
    gcd_with_algorithm(f, g, GcdAlgorithm::HuMonaganBivariate)
}

/// `(h, f / h, g / h)` where `h = gcd(f, g)`.
#[must_use]
pub fn cofactors(f: &Poly, g: &Poly) -> (Poly, Poly, Poly) {
    cofactors_with_algorithm(f, g, GcdAlgorithm::default())
}

/// Like [`gcd`], using the selected algorithm and its certified fallback.
///
/// ```
/// use polycore::{Order, Ring};
/// use zippel_gcd::{gcd_with_algorithm, GcdAlgorithm};
/// let ring = Ring::new(["x", "y"], Order::Lex);
/// let a = ring.parse("(x + y)*(x + 1)").unwrap();
/// let b = ring.parse("(x + y)*(y + 2)").unwrap();
/// let h = gcd_with_algorithm(&a, &b, GcdAlgorithm::HuangMonagan);
/// assert_eq!(h, ring.parse("x + y").unwrap());
/// ```
#[must_use]
pub fn gcd_with_algorithm(f: &Poly, g: &Poly, algorithm: GcdAlgorithm) -> Poly {
    normalized(&gcd_int(f, g, algorithm), f)
}

/// Like [`cofactors`], using the selected algorithm and its certified fallback.
#[must_use]
pub fn cofactors_with_algorithm(f: &Poly, g: &Poly, algorithm: GcdAlgorithm) -> (Poly, Poly, Poly) {
    let h = gcd_int(f, g, algorithm);
    let zero = Poly::zero(f.nvars, f.order.clone());
    if h.is_zero() {
        return (zero.clone(), zero.clone(), zero);
    }
    let h = normalized(&h, f);
    let hi = integral(&h).0;
    let (hf, df) = integral(f);
    let (hg, dg) = integral(g);
    let qf = hf.div_exact(&hi).expect("gcd divides f");
    let qg = hg.div_exact(&hi).expect("gcd divides g");
    (h, rational(&qf, &df, f), rational(&qg, &dg, g))
}

/// The least common multiple, normalized like [`gcd`]. `lcm(f, 0) = 0`.
#[must_use]
pub fn lcm(f: &Poly, g: &Poly) -> Poly {
    lcm_with_algorithm(f, g, GcdAlgorithm::default())
}

/// Like [`lcm`], using the selected algorithm and its certified fallback.
#[must_use]
pub fn lcm_with_algorithm(f: &Poly, g: &Poly, algorithm: GcdAlgorithm) -> Poly {
    let h = gcd_int(f, g, algorithm);
    if f.is_zero() || g.is_zero() {
        return Poly::zero(f.nvars, f.order.clone());
    }
    let (hg, _) = integral(g);
    normalized(
        &integral(f).0.mul(&hg.div_exact(&h).expect("gcd divides g")),
        f,
    )
}
