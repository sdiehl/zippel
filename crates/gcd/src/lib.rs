//! Sparse multivariate polynomial GCD, cofactors and LCM over the rationals.
//!
//! [`gcd()`] automatically tries Huang–Gao, Huang–Monagan, then Zippel.
//! [`gcd_with_algorithm`] chooses the first algorithm in a documented fallback
//! chain. [`cofactors`] and [`lcm`] provide the same choice through their
//! `_with_algorithm` variants.
//!
//! Call an algorithm directly to control recovery:
//!
//! | Function | Algorithm | Result |
//! | --- | --- | --- |
//! | [`zippel()`] | Zippel modular GCD with LINZIP | `Poly` |
//! | [`huang_gao()`] | Huang–Gao separated Hensel lifting; asymptotic SOTA (2026) | `Option<Poly>` |
//! | [`huang_monagan()`] | Huang–Monagan prime substitution | `Option<Poly>` |
//! | [`hu_monagan()`] | Hu–Monagan GCD/cofactor interpolation | `Option<Poly>` |
//! | [`hu_monagan_bivariate()`] | Hu–Monagan with bivariate images | `Option<Poly>` |
//!
//! Named sparse methods return `None` when bounded recovery fails; they never
//! switch to another recovery backend. They share integer normalization, CRT,
//! exact certification and, for Hu–Monagan, Zippel content extraction.
//! Every successful result has the normalization and zero conventions of [`gcd()`].
//!
//! Huang–Gao's [2026 paper](https://arxiv.org/abs/2609.08074v1) is the
//! state of the art in asymptotic sparse integer GCD complexity. This bounded,
//! exactly certified implementation does not claim the paper's full complexity
//! bound or the fastest practical runtime for every input.
//!
//! ```
//! use polycore::{Order, Ring};
//! use zippel_gcd::{gcd, huang_gao, huang_monagan, hu_monagan,
//!                  hu_monagan_bivariate, zippel};
//!
//! let ring = Ring::new(["x", "y"], Order::Lex);
//! let a = ring.parse("(x + y)*(x + 1)").unwrap();
//! let b = ring.parse("(x + y)*(y + 2)").unwrap();
//! let expected = ring.parse("x + y").unwrap();
//! assert_eq!(gcd(&a, &b), expected);       // Automatic recovery and fallback.
//! assert_eq!(zippel(&a, &b), expected);    // Explicit Zippel.
//! for algorithm in [huang_gao, huang_monagan, hu_monagan, hu_monagan_bivariate] {
//!     assert_eq!(algorithm(&a, &b).unwrap(), expected);
//! }
//! ```
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

/// Preferred algorithm for [`gcd_with_algorithm`], [`cofactors_with_algorithm`]
/// and [`lcm_with_algorithm`], including each variant's documented fallback.
///
/// To run a single recovery backend without fallback, use its named function:
/// [`zippel()`], [`huang_gao()`], [`huang_monagan()`], [`hu_monagan()`] or
/// [`hu_monagan_bivariate()`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum GcdAlgorithm {
    /// Huang–Gao separated Hensel lifting (asymptotic SOTA, 2026), then
    /// Huang–Monagan, then Zippel. See [`huang_gao()`] for the paper and scope.
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

/// Automatic GCD: Huang–Gao, then Huang–Monagan, then Zippel.
///
/// Returns a primitive integer polynomial with positive leading coefficient
/// in the input monomial order. `gcd(0, 0) = 0`. Use the named algorithm
/// functions to run a specific recovery backend without fallback.
#[must_use]
pub fn gcd(f: &Poly, g: &Poly) -> Poly {
    gcd_with_algorithm(f, g, GcdAlgorithm::default())
}

fn recover(f: &Poly, g: &Poly, algorithm: GcdAlgorithm) -> Option<Poly> {
    assert_eq!(f.nvars, g.nvars, "polynomials must share a ring");
    gcd::try_gcd_z(&integral(f).0, &integral(g).0, algorithm).map(|h| normalized(&h, f))
}

/// Zippel GCD by recursive modular interpolation with LINZIP support reuse.
///
/// Uses only the Zippel recovery backend, with retries until exact reconstruction
/// succeeds. Has the normalization and zero conventions of [`gcd()`].
#[must_use]
pub fn zippel(f: &Poly, g: &Poly) -> Poly {
    recover(f, g, GcdAlgorithm::Zippel).expect("Zippel retries until reconstruction succeeds")
}

/// Huang–Gao GCD by derivative-aided separated Hensel lifting (asymptotic SOTA, 2026).
///
/// Runs Huang–Gao recovery only; returns `None` if bounded recovery fails.
/// Every `Some` result is exactly certified and normalized like [`gcd()`].
/// Use [`gcd()`] for automatic fallback.
///
/// The [2026 paper](https://arxiv.org/abs/2609.08074v1) gives expected integer
/// bit complexity `O~(n T D log(Hi) log(Ho))`, linear in each fundamental parameter.
/// The SOTA designation concerns that asymptotic result; this bounded implementation
/// does not claim the full bound or the best practical runtime on every input.
#[must_use]
pub fn huang_gao(f: &Poly, g: &Poly) -> Option<Poly> {
    recover(f, g, GcdAlgorithm::HuangGao)
}

/// Huang–Monagan GCD by integer prime substitution and valuation decoding.
///
/// Runs Huang–Monagan recovery only; returns `None` if bounded recovery fails.
/// Every `Some` result is exactly certified and normalized like [`gcd()`].
/// See [`GcdAlgorithm::HuangMonagan`] for the paper and recovery budgets.
#[must_use]
pub fn huang_monagan(f: &Poly, g: &Poly) -> Option<Poly> {
    recover(f, g, GcdAlgorithm::HuangMonagan)
}

/// Hu–Monagan GCD/cofactor interpolation using univariate images.
///
/// Runs Hu–Monagan recovery only; returns `None` if bounded recovery fails.
/// Shared preprocessing uses Zippel to extract polynomial contents. Every `Some`
/// result is exactly certified and normalized like [`gcd()`]. See
/// [`GcdAlgorithm::HuMonagan`] for recovery budgets.
#[must_use]
pub fn hu_monagan(f: &Poly, g: &Poly) -> Option<Poly> {
    recover(f, g, GcdAlgorithm::HuMonagan)
}

/// Hu–Monagan GCD/cofactor interpolation using bivariate images.
///
/// Retains two symbolic variables, or one for a univariate ring. Has the same
/// normalization, content preprocessing and `None` behavior as [`hu_monagan()`],
/// without switching recovery backends.
#[must_use]
pub fn hu_monagan_bivariate(f: &Poly, g: &Poly) -> Option<Poly> {
    recover(f, g, GcdAlgorithm::HuMonaganBivariate)
}

/// GCD and cofactors `(h, f / h, g / h)` using the automatic Huang–Gao-first chain.
///
/// `h = gcd(f, g)`. Both cofactors are zero when both inputs are zero.
#[must_use]
pub fn cofactors(f: &Poly, g: &Poly) -> (Poly, Poly, Poly) {
    cofactors_with_algorithm(f, g, GcdAlgorithm::default())
}

/// GCD using a chosen first algorithm and its documented fallback chain.
///
/// [`GcdAlgorithm::HuangGao`] tries Huang–Gao, Huang–Monagan, then Zippel.
/// Other sparse variants try the selected method, then Zippel; the Zippel
/// variant uses only Zippel. Named functions such as [`huang_gao()`] give
/// direct access without fallback. Normalizes the result like [`gcd()`].
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

/// GCD and cofactors using a chosen first algorithm and its documented fallback chain.
///
/// See [`gcd_with_algorithm`] for the exact chain and [`cofactors`] for the result.
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

/// LCM using the automatic Huang–Gao-first GCD chain.
///
/// Normalized like [`gcd()`]. `lcm(f, 0) = 0`.
#[must_use]
pub fn lcm(f: &Poly, g: &Poly) -> Poly {
    lcm_with_algorithm(f, g, GcdAlgorithm::default())
}

/// LCM using a chosen first algorithm and its documented fallback chain.
///
/// See [`gcd_with_algorithm`] for the exact chain. Normalizes the result like [`gcd()`].
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
