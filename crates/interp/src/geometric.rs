//! Ben-Or/Tiwari recovery in a smooth subgroup, with mixed-radix exponents.

use polycore::interp::{solve, Massey};
use polycore::modp::{add, inv, mul, pow, PowerOfTwoSubgroup, SmoothPrimes};
use polycore::sample::Rng;
use polycore::{Fp, Modular};

/// Injective Kronecker encoding within a power-of-two subgroup. The bound is
/// on the entire exponent box, not just the number of terms to recover.
#[derive(Clone, Debug)]
pub struct Encoding {
    radices: Vec<u64>,
    strides: Vec<u64>,
    size: u64,
    bits: u32,
}

impl Encoding {
    /// A collision-free exponent box of at most 48 bits, or None if too large.
    pub fn new(degrees: &[u32]) -> Option<Self> {
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

    pub fn primes(&self) -> SmoothPrimes {
        SmoothPrimes::new(self.bits).expect("encoding uses 32..=48 bits")
    }
}

/// Randomly scaled geometric sampling in the encoding's smooth subgroup.
#[derive(Clone, Debug)]
pub struct Orbit<'a> {
    encoding: &'a Encoding,
    group: PowerOfTwoSubgroup,
    p: u64,
    scale: Vec<u64>,
    ratios: Vec<u64>,
}

impl<'a> Orbit<'a> {
    pub fn new(encoding: &'a Encoding, p: u64, rng: &mut Rng) -> Option<Self> {
        let group = PowerOfTwoSubgroup::new(p, encoding.bits, rng)?;
        let root = group.generator();
        let ratios = encoding.strides.iter().map(|&s| pow(root, s, p)).collect();
        let scale: Vec<_> = encoding.radices.iter().map(|_| rng.nonzero(p)).collect();
        Some(Self {
            encoding,
            group,
            p,
            scale,
            ratios,
        })
    }

    /// Initial monomial value and geometric ratio. Subsequent evaluations
    /// need one multiplication per input term (Hu–Monagan, Section 5.1).
    pub fn monomial(&self, exponents: &[u32]) -> (u64, u64) {
        assert_eq!(
            exponents.len(),
            self.scale.len(),
            "monomial dimension differs"
        );
        let at = |point: &[u64]| {
            exponents.iter().zip(point).fold(1, |v, (&e, &x)| {
                mul(v, pow(x, u64::from(e), self.p), self.p)
            })
        };
        (at(&self.scale), at(&self.ratios))
    }

    /// Scale and ratio vectors for successive geometric evaluations.
    pub fn scale(&self) -> &[u64] {
        &self.scale
    }
    pub fn ratios(&self) -> &[u64] {
        &self.ratios
    }

    fn log(&self, value: u64) -> Option<u64> {
        self.group.log(value)
    }
}

/// Incremental Ben-Or/Tiwari samples at consecutive powers starting at one.
#[derive(Clone, Debug, Default)]
pub struct Stream {
    recurrence: Massey<Fp>,
    values: Vec<u64>,
    modulus: Option<u64>,
}

impl Stream {
    /// Add one sample. All samples must belong to the same prime field.
    pub fn push(&mut self, value: u64, p: u64) {
        if let Some(previous) = self.modulus {
            assert_eq!(p, previous, "incompatible sample fields");
        } else {
            assert!(polycore::modp::is_prime(p), "modulus must be prime");
            self.modulus = Some(p);
        }
        let value = value % p;
        self.recurrence.push(Fp::new(value, p));
        self.values.push(value);
    }

    pub fn recover(&self, orbit: &Orbit<'_>) -> Option<Vec<(Vec<u32>, u64)>> {
        if self.modulus != Some(orbit.p) || !self.recurrence.settled(2) {
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
    #[should_panic(expected = "incompatible sample fields")]
    fn rejects_mixed_sample_fields() {
        let mut stream = Stream::default();
        stream.push(1, 101);
        stream.push(1, 103);
    }

    #[test]
    fn normalizes_samples_and_checks_the_recovery_field() {
        let encoding = Encoding::new(&[1]).unwrap();
        let mut primes = encoding.primes();
        let p = primes.next().unwrap();
        let q = primes.next().unwrap();
        let orbit = Orbit::new(&encoding, p, &mut Rng::new(9)).unwrap();
        let other = Orbit::new(&encoding, q, &mut Rng::new(9)).unwrap();
        let mut stream = Stream::default();
        for _ in 0..4 {
            stream.push(p + 7, p);
        }
        assert_eq!(stream.recover(&orbit), Some(vec![(vec![0], 7)]));
        assert!(stream.recover(&other).is_none());
    }

    #[test]
    fn smooth_primes_and_binary_logs_cover_the_subgroup() {
        let encoding = Encoding::new(&[100_000, 100_000]).unwrap();
        let primes: Vec<_> = encoding.primes().take(3).collect();
        assert!(primes.windows(2).all(|w| w[0] > w[1]));
        for p in primes {
            assert!(polycore::modp::is_prime(p));
            let orbit = Orbit::new(&encoding, p, &mut Rng::new(17)).unwrap();
            for e in [0, 1, 53, encoding.size - 1, (1u64 << encoding.bits) - 1] {
                assert_eq!(orbit.log(orbit.group.pow(e)), Some(e));
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
