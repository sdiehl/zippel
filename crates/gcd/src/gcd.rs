//! GCD over Z: explicit recovery backends and an automatic fallback chain.

use crate::linzip::linzip;
use crate::pgcd::{pgcd, reshape};
use crate::poly::{add_exps, gcd, one, Exps, IntPoly, ModPoly};
use crate::GcdAlgorithm;
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Zero};
use polycore::crt::WangContext;
use polycore::fast::PrimeField;
use polycore::modp::{inv, mul, sub, Primes};
use polycore::sample::Rng;
use std::iter::once;

/// `gcd(f, g)` up to sign.
pub(crate) fn gcd_z(f: &IntPoly, g: &IntPoly, algorithm: GcdAlgorithm) -> IntPoly {
    try_gcd_z(f, g, algorithm)
        .or_else(|| {
            (algorithm == GcdAlgorithm::HuangGao)
                .then(|| try_gcd_z(f, g, GcdAlgorithm::HuangMonagan))
                .flatten()
        })
        .unwrap_or_else(|| {
            try_gcd_z(f, g, GcdAlgorithm::Zippel)
                .expect("Zippel retries until reconstruction succeeds")
        })
}

/// Run only the requested recovery backend, retaining shared normalization
/// and content extraction. None means recovery declined, not that the GCD is 1.
pub(crate) fn try_gcd_z(f: &IntPoly, g: &IntPoly, algorithm: GcdAlgorithm) -> Option<IntPoly> {
    if f.is_zero() || g.is_zero() {
        let h = if f.is_zero() { g } else { f };
        return Some(if h.is_zero() {
            h.clone()
        } else {
            h.primitive()
        });
    }
    let c = gcd(&f.content(), &g.content());
    let (mf, mg) = (f.min_exps(), g.min_exps());
    let m: Exps = mf.iter().zip(&mg).map(|(a, b)| *a.min(b)).collect();
    let shift = |h: &IntPoly, s: &Exps| {
        h.primitive()
            .map(|e, x| (e.iter().zip(s).map(|(a, b)| a - b).collect(), x.clone()))
    };
    let h = gcd_shifted(&shift(f, &mf), &shift(g, &mg), algorithm)?;
    Some(h.map(|e, x| (add_exps(e, &m), x * &c)))
}

/// Separated lifting handles all variables and their polynomial contents at once.
/// A declined attempt returns to the caller without switching recovery methods.
fn gcd_shifted(f: &IntPoly, g: &IntPoly, algorithm: GcdAlgorithm) -> Option<IntPoly> {
    let n = f.n;
    let shared = (0..n).filter(|&k| f.degree(k) > 0 && g.degree(k) > 0);
    let Some(x0) = shared.max_by_key(|&k| f.degree(k).min(g.degree(k))) else {
        return Some(one(n));
    };
    if algorithm == GcdAlgorithm::HuangGao {
        return modular(f, g, Backend::Separated);
    }
    if algorithm == GcdAlgorithm::HuangMonagan {
        return crate::huang_monagan::gcd(f, g);
    }
    let perm: Vec<usize> = once(x0).chain((0..n).filter(|&k| k != x0)).collect();
    let mut back = vec![0; n];
    perm.iter().enumerate().for_each(|(i, &k)| back[k] = i);
    let (f, g) = (f.permute(&perm), g.permute(&perm));
    let (cf, cg) = (content_x0(&f), content_x0(&g));
    let cont = gcd_z(&cf, &cg, GcdAlgorithm::Zippel);
    let pf = f.div_exact(&cf).expect("content divides");
    let pg = g.div_exact(&cg).expect("content divides");
    let keep = match algorithm {
        GcdAlgorithm::HuMonagan => Some(1),
        GcdAlgorithm::HuMonaganBivariate => Some(n.min(2)),
        _ => None,
    };
    let h = if let Some(keep) = keep {
        let (a, b) = if pf.terms.len() <= pg.terms.len() {
            (&pf, &pg)
        } else {
            (&pg, &pf)
        };
        modular(a, b, Backend::HuMonagan(keep))?
    } else {
        modular(&pf, &pg, Backend::Zippel)?
    };
    Some(h.mul(&cont).permute(&back))
}

fn content_x0(f: &IntPoly) -> IntPoly {
    let mut coeffs = f.coefficients_in_x0().into_iter();
    let first = coeffs.next().expect("nonzero");
    coeffs
        .try_fold(first, |acc, c| {
            let h = gcd_z(&acc, &c, GcdAlgorithm::Zippel);
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

#[derive(Clone, Copy, PartialEq, Eq)]
enum Backend {
    Separated,
    Univariate,
    HuMonagan(usize),
    Zippel,
}

/// Primitive univariate integer GCD for prime-substituted inputs. This calls
/// the modular CRT engine directly, never the multivariate dispatch recursively.
pub(crate) fn univariate(f: &IntPoly, g: &IntPoly) -> Option<IntPoly> {
    debug_assert!(f.n == 1 && g.n == 1 && !f.is_zero() && !g.is_zero());
    modular(f, g, Backend::Univariate)
}

/// Images are scaled so their leading coefficient is `gcd(lc f, lc g)`, which makes them agree
/// across primes. A candidate that divides both inputs is the gcd, since its leading monomial is
/// no smaller than the true one.
fn modular(f: &IntPoly, g: &IntPoly, backend: Backend) -> Option<IntPoly> {
    let n = f.n;
    let mut rng = Rng::new(0x5eed);
    let gamma = gcd(f.lc(), g.lc());
    let mut skeleton: Option<Vec<Exps>> = None;
    let mut acc: Vec<BigInt> = Vec::new();
    let mut modulus = BigInt::one();
    let mut images = 0usize;
    let mut tried: Option<IntPoly> = None;
    let primes: Box<dyn Iterator<Item = u64>> = if let Backend::HuMonagan(keep) = backend {
        Box::new(zippel_interp::geometric::Encoding::new(&f.degrees()[keep..])?.primes())
    } else {
        Box::new(Primes::new())
    };
    for (attempt, p) in primes.enumerate() {
        if matches!(backend, Backend::Univariate | Backend::HuMonagan(_)) && attempt >= 256 {
            return None;
        }
        let gp_ = residue(&gamma, p);
        let (fp, gp) = (f.reduce(p), g.reduce(p));
        if gp_ == 0 || fp.lm() != f.terms[0].0 || gp.lm() != g.terms[0].0 {
            continue;
        }
        if backend == Backend::Zippel && skeleton.is_none() && coprime(&fp, &gp, p, &mut rng) {
            return Some(one(n));
        }
        let h = if backend == Backend::Separated {
            // A declined field attempt restarts with the established Zippel path;
            // never combine a partial factor with the CRT accumulator.
            crate::huang_gao::gcd(&fp, &gp, p, &mut rng)?
        } else if let Backend::HuMonagan(keep) = backend {
            crate::hu_monagan::gcd(&fp, &gp, keep, p, &mut rng)?
        } else if backend == Backend::Univariate {
            let dense = |f: &ModPoly| {
                let mut cs = vec![0; f.degree(0) + 1];
                for (e, c) in &f.terms {
                    cs[e[0] as usize] = *c;
                }
                cs
            };
            let h = PrimeField::new(p)
                .expect("prime modulus")
                .gcd(&dense(&fp), &dense(&gp));
            ModPoly {
                n: 1,
                terms: h
                    .into_iter()
                    .enumerate()
                    .rev()
                    .filter(|(_, c)| *c != 0)
                    .map(|(e, c)| (vec![e as u32], c))
                    .collect(),
            }
        } else {
            let sparse = skeleton
                .as_deref()
                .filter(|_| n >= 2)
                .and_then(|s| linzip(&fp, &gp, s, n - 1, p, &mut rng));
            let Some(h) = sparse.or_else(|| pgcd(&fp, &gp, n - 1, p, &mut rng)) else {
                continue;
            };
            h
        };
        let mut h = h.monic(p);
        h.scale(gp_, p);
        match reshape(skeleton.as_deref(), &h) {
            None => continue,
            Some(true) => {
                skeleton = Some(h.support());
                acc = vec![BigInt::zero(); h.terms.len()];
                modulus = BigInt::one();
                images = 0;
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
        images += 1;
        let Some(h) = candidate(n, s, &acc, &modulus, images.is_power_of_two())
            .filter(|h| tried.as_ref() != Some(h))
        else {
            continue;
        };
        if f.div_exact(&h).is_some() && g.div_exact(&h).is_some() {
            return Some(h);
        }
        tried = Some(h);
    }
    None
}

/// Bits of headroom below the modulus that make a reconstruction worth a trial division.
const SLACK: u64 = 16;

/// A candidate once the residues settle well inside the modulus: as integers, or else as the
/// fractions `h / lc(h)` of the images divided by `gamma`, which settle first when `gamma` is
/// much larger than `lc(h)`. Fractions cost a reconstruction, so only some rounds try them.
fn candidate(n: usize, s: &[Exps], acc: &[BigInt], m: &BigInt, fractions: bool) -> Option<IntPoly> {
    let room = m.bits().saturating_sub(SLACK);
    let half = m >> 1;
    let lift: Vec<BigInt> = acc
        .iter()
        .map(|a| if *a > half { a - m } else { a.clone() })
        .collect();
    let poly = |cs: Vec<BigInt>| Some(IntPoly::new(n, s.iter().cloned().zip(cs)).primitive());
    if lift.iter().all(|a| a.bits() <= room) {
        return poly(lift);
    }
    if !fractions {
        return None;
    }
    let w = WangContext::new(m)?;
    let gi = acc[0].modinv(m)?;
    let fraction = |a: &BigInt| {
        w.reconstruct(&(a * &gi).mod_floor(m))
            .filter(|q| q.numer().bits() + q.denom().bits() <= room)
    };
    fraction(acc.last()?)?;
    let qs = acc.iter().map(fraction).collect::<Option<Vec<_>>>()?;
    let l = qs.iter().fold(BigInt::one(), |l, q| l.lcm(q.denom()));
    poly(qs.iter().map(|q| q.numer() * (&l / q.denom())).collect())
}

/// With `x_0`-degrees preserved, a constant univariate image proves `gcd = 1`, because the gcd
/// is primitive in `x_0`.
fn coprime(f: &ModPoly, g: &ModPoly, p: u64, rng: &mut Rng) -> bool {
    let point: Vec<u64> = (0..f.n).map(|_| rng.nonzero(p)).collect();
    let (uf, ug) = (f.eval_except(0, &point, p), g.eval_except(0, &point, p));
    uf.deg() == f.degree(0) && ug.deg() == g.degree(0) && uf.gcd(&ug).0.len() == 1
}

#[cfg(test)]
mod tests {
    use super::*;
    use polycore::{Order, Ring};

    #[test]
    fn hu_monagan_crt_recovers_large_and_initially_missing_coefficients() {
        let ring = Ring::new(["x", "y", "z"], Order::Lex);
        let prime = zippel_interp::geometric::Encoding::new(&[4, 4])
            .unwrap()
            .primes()
            .next()
            .unwrap();
        let large = (BigInt::one() << 190) + BigInt::from(321);
        let h = ring.parse(&format!("x^3+{large}*x*y+{prime}*z+1")).unwrap();
        let a = crate::integral(&(&h * &ring.parse("x+y+2").unwrap())).0;
        let b = crate::integral(&(&h * &ring.parse("x+z+3").unwrap())).0;
        for keep in [1, 2] {
            let g =
                modular(&a, &b, Backend::HuMonagan(keep)).expect("Hu–Monagan CRT without fallback");
            assert_eq!(g, crate::integral(&h).0);
        }
    }
}
