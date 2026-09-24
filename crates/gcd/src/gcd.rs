//! Gcd over Z: strip contents, then join Zippel images modulo word-sized primes by CRT.

use crate::linzip::linzip;
use crate::pgcd::{pgcd, reshape};
use crate::poly::{add_exps, one, Exps, IntPoly, ModPoly};
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Zero};
use std::iter::once;
use zippel_interp::modp::{inv, mul, sub, Primes, Rng};
use zippel_interp::univariate as uni;

/// `gcd(f, g)` up to sign.
pub(crate) fn gcd_z(f: &IntPoly, g: &IntPoly) -> IntPoly {
    if f.is_zero() || g.is_zero() {
        let h = if f.is_zero() { g } else { f };
        return if h.is_zero() {
            h.clone()
        } else {
            h.primitive()
        };
    }
    let c = f.content().gcd(&g.content());
    let (mf, mg) = (f.min_exps(), g.min_exps());
    let m: Exps = mf.iter().zip(&mg).map(|(a, b)| *a.min(b)).collect();
    let shift = |h: &IntPoly, s: &Exps| {
        h.primitive()
            .map(|e, x| (e.iter().zip(s).map(|(a, b)| a - b).collect(), x.clone()))
    };
    let h = gcd_shifted(&shift(f, &mf), &shift(g, &mg));
    h.map(|e, x| (add_exps(e, &m), x * &c))
}

/// Main variable is the shared one of largest degree; its content comes from recursion on the
/// coefficients, leaving the modular stage a gcd that is primitive in `x_0`.
fn gcd_shifted(f: &IntPoly, g: &IntPoly) -> IntPoly {
    let n = f.n;
    let shared = (0..n).filter(|&k| f.degree(k) > 0 && g.degree(k) > 0);
    let Some(x0) = shared.max_by_key(|&k| f.degree(k).min(g.degree(k))) else {
        return one(n);
    };
    let perm: Vec<usize> = once(x0).chain((0..n).filter(|&k| k != x0)).collect();
    let mut back = vec![0; n];
    perm.iter().enumerate().for_each(|(i, &k)| back[k] = i);
    let (f, g) = (f.permute(&perm), g.permute(&perm));
    let (cf, cg) = (content_x0(&f), content_x0(&g));
    let cont = gcd_z(&cf, &cg);
    let pf = f.div_exact(&cf).expect("content divides");
    let pg = g.div_exact(&cg).expect("content divides");
    modular(&pf, &pg).mul(&cont).permute(&back)
}

fn content_x0(f: &IntPoly) -> IntPoly {
    let mut coeffs = f.coefficients_in_x0().into_iter();
    let first = coeffs.next().expect("nonzero");
    coeffs
        .try_fold(first, |acc, c| {
            let h = gcd_z(&acc, &c);
            if h.is_constant() {
                Err(())
            } else {
                Ok(h)
            }
        })
        .unwrap_or_else(|()| one(f.n))
}

fn residue(a: &BigInt, p: u64) -> u64 {
    u64::try_from(a.mod_floor(&BigInt::from(p))).expect("reduced")
}

/// Images are scaled so their leading coefficient is `gcd(lc f, lc g)`, which makes them agree
/// across primes. A candidate that divides both inputs is the gcd, since its leading monomial is
/// no smaller than the true one.
fn modular(f: &IntPoly, g: &IntPoly) -> IntPoly {
    let n = f.n;
    let mut rng = Rng::new(0x5eed);
    let gamma = f.lc().gcd(g.lc());
    let mut skeleton: Option<Vec<Exps>> = None;
    let mut acc: Vec<BigInt> = Vec::new();
    let mut modulus = BigInt::one();
    let mut last: Option<IntPoly> = None;
    for p in Primes::new() {
        let gp_ = residue(&gamma, p);
        let (fp, gp) = (f.reduce(p), g.reduce(p));
        if gp_ == 0 || fp.lm() != f.terms[0].0 || gp.lm() != g.terms[0].0 {
            continue;
        }
        if skeleton.is_none() && coprime(&fp, &gp, p, &mut rng) {
            return one(n);
        }
        let sparse = skeleton
            .as_deref()
            .filter(|_| n >= 2)
            .and_then(|s| linzip(&fp, &gp, s, n - 1, p, &mut rng));
        let Some(h) = sparse.or_else(|| pgcd(&fp, &gp, n - 1, p, &mut rng)) else {
            continue;
        };
        let mut h = h.monic(p);
        h.scale(gp_, p);
        match reshape(skeleton.as_deref(), &h) {
            None => continue,
            Some(true) => {
                skeleton = Some(h.support());
                acc = vec![BigInt::zero(); h.terms.len()];
                modulus = BigInt::one();
                last = None;
            }
            Some(false) => {}
        }
        let s = skeleton.as_ref().expect("set above");
        let m_inv = inv(residue(&modulus, p), p);
        for (e, a) in s.iter().zip(&mut acc) {
            let r = h
                .terms
                .binary_search_by(|t| e.cmp(&t.0))
                .map_or(0, |i| h.terms[i].1);
            let t = mul(sub(r, residue(a, p), p), m_inv, p);
            *a += &modulus * t;
        }
        modulus *= p;
        let half = &modulus >> 1;
        let candidate = IntPoly::new(
            n,
            s.iter()
                .zip(&acc)
                .map(|(e, a)| (e.clone(), if *a > half { a - &modulus } else { a.clone() })),
        );
        if last.as_ref() == Some(&candidate) {
            let h = candidate.primitive();
            if f.div_exact(&h).is_some() && g.div_exact(&h).is_some() {
                return h;
            }
        }
        last = Some(candidate);
    }
    unreachable!("ran out of primes")
}

/// With `x_0`-degrees preserved, a constant univariate image proves `gcd = 1`, because the gcd
/// is primitive in `x_0`.
fn coprime(f: &ModPoly, g: &ModPoly, p: u64, rng: &mut Rng) -> bool {
    let point: Vec<u64> = (0..f.n).map(|_| rng.nonzero(p)).collect();
    let (uf, ug) = (f.eval_except(0, &point, p), g.eval_except(0, &point, p));
    uni::deg(&uf) == f.degree(0) && uni::deg(&ug) == g.degree(0) && uni::gcd(&uf, &ug, p).len() == 1
}
