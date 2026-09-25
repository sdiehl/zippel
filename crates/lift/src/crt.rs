//! Chinese remaindering and Wang's rational number reconstruction.

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Signed, ToPrimitive, Zero};
use zippel_interp::modp::{inv, mul, sub};

fn residue(x: &BigInt, p: u64) -> u64 {
    x.mod_floor(&BigInt::from(p)).to_u64().unwrap()
}

/// Garner's step: residues modulo `m` and `values` modulo `p` to residues modulo `m * p`.
pub(crate) fn garner(xs: &mut [BigInt], m: &mut BigInt, values: &[u64], p: u64) {
    let minv = inv(residue(m, p), p);
    for (x, &v) in xs.iter_mut().zip(values) {
        *x += &*m * mul(sub(v, residue(x, p), p), minv, p);
    }
    *m *= p;
}

/// `q` modulo `p`, or `None` when `p` divides its denominator.
pub(crate) fn reduce(q: &BigRational, p: u64) -> Option<u64> {
    let s = residue(q.denom(), p);
    (s != 0).then(|| mul(residue(q.numer(), p), inv(s, p), p))
}

/// The fraction `r / s` congruent to `x` modulo `m` with `|r|, |s| <= sqrt(m / 2)`, unique when it
/// exists, found by stopping the extended Euclidean algorithm halfway.
pub(crate) fn wang(x: &BigInt, m: &BigInt) -> Option<BigRational> {
    let bound = (m / 2u32).sqrt();
    let (mut r0, mut r1) = (m.clone(), x.clone());
    let (mut s0, mut s1) = (BigInt::zero(), BigInt::one());
    while r1 > bound {
        let q = &r0 / &r1;
        let r = &r0 - &q * &r1;
        r0 = std::mem::replace(&mut r1, r);
        let s = &s0 - &q * &s1;
        s0 = std::mem::replace(&mut s1, s);
    }
    (!s1.is_zero() && s1.abs() <= bound && r1.gcd(&s1).is_one()).then(|| BigRational::new(r1, s1))
}
