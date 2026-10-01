//! Ben-Or/Tiwari recovery in a smooth subgroup, with mixed-radix exponents.

use polycore::interp::{solve, Massey};
use polycore::modp::{add, inv, is_prime, mul, pow};
use polycore::sample::Rng;
use polycore::{Fp, Modular};

/// Injective Kronecker encoding within a power-of-two subgroup. The bound is
/// on the entire exponent box, not just the number of terms to recover.
pub(crate) struct Encoding {
    radices: Vec<u64>,
    strides: Vec<u64>,
    size: u64,
    bits: u32,
}

impl Encoding {
    pub(crate) fn new(degrees: &[u32]) -> Option<Self> {
        let mut size = 1u64;
        let mut strides = Vec::new();
        let radices: Vec<_> = degrees.iter().map(|&d| u64::from(d) + 1).collect();
        for &radix in &radices {
            strides.push(size);
            size = size.checked_mul(radix)?;
        }
        let bits = (64 - (size - 1).leading_zeros()).max(32);
        (bits <= 48).then_some(Self {
            radices,
            strides,
            size,
            bits,
        })
    }

    fn decode(&self, mut e: u64) -> Option<Vec<u32>> {
        if e >= self.size {
            return None;
        }
        Some(
            self.radices
                .iter()
                .map(|&r| {
                    let digit = (e % r) as u32;
                    e /= r;
                    digit
                })
                .collect(),
        )
    }

    pub(crate) const fn primes(&self) -> SmoothPrimes {
        let step = 1u64 << self.bits;
        let odd = (((1u64 << 62) - 1) / step) | 1;
        SmoothPrimes {
            next: odd * step + 1,
            step: 2 * step,
        }
    }
}

/// Distinct word-sized primes whose multiplicative groups contain the needed
/// smooth subgroup. Its exact order avoids factoring the odd part of p - 1.
pub(crate) struct SmoothPrimes {
    next: u64,
    step: u64,
}

impl Iterator for SmoothPrimes {
    type Item = u64;
    fn next(&mut self) -> Option<u64> {
        while self.next > self.step {
            let p = self.next;
            self.next -= self.step;
            if is_prime(p) {
                return Some(p);
            }
        }
        None
    }
}

pub(crate) struct Orbit<'a> {
    encoding: &'a Encoding,
    root: u64,
    p: u64,
    scale: Vec<u64>,
    ratios: Vec<u64>,
}

impl<'a> Orbit<'a> {
    pub(crate) fn new(encoding: &'a Encoding, p: u64, rng: &mut Rng) -> Option<Self> {
        let order = 1u64 << encoding.bits;
        if !(p - 1).is_multiple_of(order) {
            return None;
        }
        let root = (0..32).find_map(|_| {
            let r = pow(rng.nonzero(p), (p - 1) / order, p);
            (pow(r, order / 2, p) != 1).then_some(r)
        })?;
        let ratios = encoding.strides.iter().map(|&s| pow(root, s, p)).collect();
        let scale: Vec<_> = encoding.radices.iter().map(|_| rng.nonzero(p)).collect();
        Some(Self {
            encoding,
            root,
            p,
            scale,
            ratios,
        })
    }

    /// Initial monomial value and geometric ratio. Subsequent evaluations
    /// need one multiplication per input term (Hu–Monagan, Section 5.1).
    pub(crate) fn monomial(&self, exponents: &[u32]) -> (u64, u64) {
        let at = |point: &[u64]| {
            exponents.iter().zip(point).fold(1, |v, (&e, &x)| {
                mul(v, pow(x, u64::from(e), self.p), self.p)
            })
        };
        (at(&self.scale), at(&self.ratios))
    }

    /// Binary Pohlig-Hellman in the subgroup of order 2^bits.
    fn log(&self, mut value: u64) -> Option<u64> {
        let mut inverse = inv(self.root, self.p);
        let mut exponent = 0;
        for bit in 0..self.encoding.bits {
            match pow(value, 1u64 << (self.encoding.bits - bit - 1), self.p) {
                1 => {}
                minus_one if minus_one == self.p - 1 => {
                    exponent |= 1u64 << bit;
                    value = mul(value, inverse, self.p);
                }
                _ => return None,
            }
            inverse = mul(inverse, inverse, self.p);
        }
        (value == 1).then_some(exponent)
    }
}

#[derive(Default)]
pub(crate) struct Stream {
    recurrence: Massey<Fp>,
    values: Vec<u64>,
}

impl Stream {
    pub(crate) fn push(&mut self, value: u64, p: u64) {
        self.recurrence.push(Fp::new(value, p));
        self.values.push(value);
    }

    pub(crate) fn recover(&self, orbit: &Orbit<'_>) -> Option<Vec<(Vec<u32>, u64)>> {
        if !self.recurrence.settled(2) {
            return None;
        }
        let t = self.recurrence.complexity();
        if t == 0 {
            return Some(Vec::new());
        }
        let master = self.recurrence.generator();
        let roots = master.roots();
        if roots.len() != t {
            return None;
        }
        let p = orbit.p;
        let values: Vec<_> = self.values[..t].iter().map(|&v| Fp::new(v, p)).collect();
        let coefficients = solve(&roots, &master, &values);
        let residues: Vec<_> = roots.iter().map(|v| v.residue_mod(p)).collect();
        let coefficients: Vec<_> = coefficients.iter().map(|v| v.residue_mod(p)).collect();
        let mut powers = residues.clone();
        for &value in &self.values {
            let predicted = powers
                .iter()
                .zip(&coefficients)
                .fold(0, |a, (&v, &c)| add(a, mul(v, c, p), p));
            if predicted != value {
                return None;
            }
            for (v, &r) in powers.iter_mut().zip(&residues) {
                *v = mul(*v, r, p);
            }
        }
        residues
            .iter()
            .zip(coefficients)
            .map(|(&r, c)| {
                let e = orbit.encoding.decode(orbit.log(r)?)?;
                let scale = e
                    .iter()
                    .zip(&orbit.scale)
                    .fold(1, |s, (&e, &v)| mul(s, pow(v, u64::from(e), p), p));
                Some((e, mul(c, inv(scale, p), p)))
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn smooth_primes_and_binary_logs_cover_the_subgroup() {
        let encoding = Encoding::new(&[100_000, 100_000]).unwrap();
        let primes: Vec<_> = encoding.primes().take(3).collect();
        assert!(primes.windows(2).all(|w| w[0] > w[1]));
        for p in primes {
            assert!(is_prime(p));
            let orbit = Orbit::new(&encoding, p, &mut Rng::new(17)).unwrap();
            for e in [0, 1, 53, encoding.size - 1, (1u64 << encoding.bits) - 1] {
                assert_eq!(orbit.log(pow(orbit.root, e, p)), Some(e));
            }
            assert_eq!(orbit.log(0), None);
            assert!((2..100).any(|v| orbit.log(v).is_none()));
        }
    }

    #[test]
    fn recovers_multivariate_exponents_and_unscales_coefficients() {
        let encoding = Encoding::new(&[100_000, 100_000]).unwrap();
        let p = encoding.primes().next().unwrap();
        let orbit = Orbit::new(&encoding, p, &mut Rng::new(71)).unwrap();
        let mut point = orbit.scale.clone();
        let terms = [
            (vec![0u32, 0], 7),
            (vec![100_000, 100_000], p - 31),
            (vec![51, 13], 5),
        ];
        let mut stream = Stream::default();
        for _ in 0..10 {
            for (x, &r) in point.iter_mut().zip(&orbit.ratios) {
                *x = mul(*x, r, p);
            }
            let value = terms.iter().fold(0, |sum, (e, c)| {
                let term = e
                    .iter()
                    .zip(&point)
                    .fold(*c, |v, (&e, &x)| mul(v, pow(x, u64::from(e), p), p));
                add(sum, term, p)
            });
            stream.push(value, p);
        }
        let mut recovered = stream.recover(&orbit).unwrap();
        recovered.sort_unstable();
        let mut expected = terms.to_vec();
        expected.sort_unstable();
        assert_eq!(recovered, expected);
    }

    #[test]
    fn exponent_boxes_never_wrap_or_alias() {
        assert!(Encoding::new(&[u32::MAX, u32::MAX]).is_none());
        assert!(Encoding::new(&[1 << 24, 1 << 24]).is_none());
        let encoding = Encoding::new(&[3, 5]).unwrap();
        assert_eq!(encoding.decode(23), Some(vec![3, 5]));
        assert_eq!(encoding.decode(24), None);
    }
}
