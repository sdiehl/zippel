//! Fast univariate arithmetic over word-sized prime fields. Convolution uses
//! three auxiliary NTT primes and exact CRT (not a floating-point FFT). Newton
//! inversion supplies fast division; half-GCD batches Euclidean steps.

use polycore::modp::{add, inv, is_prime, mul, pow, sub};
use std::sync::OnceLock;

type Poly = Vec<u64>;
type Matrix = [[Poly; 2]; 2];
const CROSSOVER: usize = 64;
const TRANSFORM_BITS: u32 = 24;

#[derive(Clone, Copy)]
struct NttPrime {
    p: u64,
    root: u64,
}

fn primes() -> &'static [NttPrime; 3] {
    static PRIMES: OnceLock<[NttPrime; 3]> = OnceLock::new();
    PRIMES.get_or_init(|| {
        let step = 1 << TRANSFORM_BITS;
        let mut candidate = (1 << 61) + 1;
        std::array::from_fn(|_| {
            loop {
                candidate -= step;
                if is_prime(candidate) {
                    break;
                }
            }
            let p = candidate;
            assert!(
                p > 1 << 60,
                "CRT convolution requires 60-bit auxiliary primes"
            );
            let root = (2..)
                .find_map(|a| {
                    let r = pow(a, (p - 1) / step, p);
                    (pow(r, step / 2, p) != 1).then_some(r)
                })
                .expect("prime field has roots of unity");
            NttPrime { p, root }
        })
    })
}

fn ntt(a: &mut [u64], prime: NttPrime, inverse: bool) {
    let n = a.len();
    let mut j = 0;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            a.swap(i, j);
        }
    }
    let root = if inverse {
        inv(prime.root, prime.p)
    } else {
        prime.root
    };
    let mut width = 2;
    while width <= n {
        let step = pow(root, (1 << TRANSFORM_BITS) / width as u64, prime.p);
        for block in a.chunks_exact_mut(width) {
            let mut w = 1;
            let (left, right) = block.split_at_mut(width / 2);
            for (a, b) in left.iter_mut().zip(right) {
                let v = mul(*b, w, prime.p);
                (*a, *b) = (add(*a, v, prime.p), sub(*a, v, prime.p));
                w = mul(w, step, prime.p);
            }
        }
        width *= 2;
    }
    if inverse {
        let scale = inv(n as u64, prime.p);
        for a in a {
            *a = mul(*a, scale, prime.p);
        }
    }
}

fn trim(mut a: Poly) -> Poly {
    while a.last() == Some(&0) {
        a.pop();
    }
    a
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct Arithmetic {
    pub(crate) p: u64,
}

impl Arithmetic {
    pub(crate) fn scale(self, a: &[u64], c: u64) -> Poly {
        trim(a.iter().map(|&a| mul(a, c, self.p)).collect())
    }

    pub(crate) fn add(self, a: &[u64], b: &[u64]) -> Poly {
        trim(
            (0..a.len().max(b.len()))
                .map(|i| {
                    add(
                        a.get(i).copied().unwrap_or(0),
                        b.get(i).copied().unwrap_or(0),
                        self.p,
                    )
                })
                .collect(),
        )
    }

    fn sub(self, a: &[u64], b: &[u64]) -> Poly {
        trim(
            (0..a.len().max(b.len()))
                .map(|i| {
                    sub(
                        a.get(i).copied().unwrap_or(0),
                        b.get(i).copied().unwrap_or(0),
                        self.p,
                    )
                })
                .collect(),
        )
    }

    pub(crate) fn mul(self, a: &[u64], b: &[u64]) -> Poly {
        if a.is_empty() || b.is_empty() {
            return vec![];
        }
        let len = a.len() + b.len() - 1;
        if a.len().min(b.len()) < CROSSOVER {
            let mut result = vec![0; len];
            for (i, &a) in a.iter().enumerate() {
                for (j, &b) in b.iter().enumerate() {
                    result[i + j] = add(result[i + j], mul(a, b, self.p), self.p);
                }
            }
            return trim(result);
        }
        let n = len.next_power_of_two();
        assert!(n <= 1 << TRANSFORM_BITS, "NTT capacity exceeded");
        let qs = primes();
        let residues: Vec<_> = qs
            .iter()
            .map(|&q| {
                let mut x = vec![0; n];
                let mut y = vec![0; n];
                for (x, a) in x.iter_mut().zip(a) {
                    *x = a % q.p;
                }
                for (y, b) in y.iter_mut().zip(b) {
                    *y = b % q.p;
                }
                ntt(&mut x, q, false);
                ntt(&mut y, q, false);
                for (x, y) in x.iter_mut().zip(y) {
                    *x = mul(*x, y, q.p);
                }
                ntt(&mut x, q, true);
                x.truncate(len);
                x
            })
            .collect();
        let [q0, q1, q2] = qs.map(|q| q.p);
        let i01 = inv(q0 % q1, q1);
        let i012 = inv(mul(q0 % q2, q1 % q2, q2), q2);
        let q01p = mul(q0 % self.p, q1 % self.p, self.p);
        // Each exact coefficient is < 2^24 * (2^64)^2 = 2^152;
        // the product of these three primes is > 2^180. Garner's mixed
        // radix digits therefore recover the INTEGER convolution uniquely.
        trim(
            (0..len)
                .map(|i| {
                    let a = residues[0][i];
                    let b = mul(sub(residues[1][i], a % q1, q1), i01, q1);
                    let ab2 = add(a % q2, mul(q0 % q2, b % q2, q2), q2);
                    let c = mul(sub(residues[2][i], ab2, q2), i012, q2);
                    add(
                        add(a % self.p, mul(q0 % self.p, b % self.p, self.p), self.p),
                        mul(q01p, c % self.p, self.p),
                        self.p,
                    )
                })
                .collect(),
        )
    }

    fn inverse_series(self, a: &[u64], len: usize) -> Poly {
        let mut result = vec![inv(a[0], self.p)];
        let mut size = 1;
        while size < len {
            size = (2 * size).min(len);
            let mut correction = self.mul(&a[..a.len().min(size)], &result);
            correction.resize(size, 0);
            correction.truncate(size);
            for c in &mut correction {
                *c = sub(0, *c, self.p);
            }
            correction[0] = add(correction[0], 2 % self.p, self.p);
            result = self.mul(&result, &correction);
            result.resize(size, 0);
            result.truncate(size);
        }
        result
    }

    pub(crate) fn divrem(self, a: &[u64], b: &[u64]) -> (Poly, Poly) {
        assert!(!b.is_empty(), "division by zero polynomial");
        if a.len() < b.len() {
            return (vec![], a.to_vec());
        }
        let count = a.len() - b.len() + 1;
        if count.min(b.len()) < CROSSOVER {
            let mut r = a.to_vec();
            let mut q = vec![0; count];
            let li = inv(*b.last().unwrap(), self.p);
            for i in (0..count).rev() {
                q[i] = mul(r[i + b.len() - 1], li, self.p);
                for (j, &b) in b.iter().enumerate() {
                    r[i + j] = sub(r[i + j], mul(q[i], b, self.p), self.p);
                }
            }
            return (trim(q), trim(r));
        }
        let ar: Vec<_> = a.iter().rev().take(count).copied().collect();
        let br: Vec<_> = b.iter().rev().take(count).copied().collect();
        let mut q = self.mul(&ar, &self.inverse_series(&br, count));
        q.resize(count, 0);
        q.truncate(count);
        q.reverse();
        let r = self.sub(a, &self.mul(&q, b));
        debug_assert!(r.len() < b.len());
        (trim(q), r)
    }

    fn matrix_mul(self, a: &Matrix, b: &Matrix) -> Matrix {
        std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                self.add(&self.mul(&a[i][0], &b[0][j]), &self.mul(&a[i][1], &b[1][j]))
            })
        })
    }

    fn apply(self, m: &Matrix, a: &[u64], b: &[u64]) -> (Poly, Poly) {
        (
            self.add(&self.mul(&m[0][0], a), &self.mul(&m[0][1], b)),
            self.add(&self.mul(&m[1][0], a), &self.mul(&m[1][1], b)),
        )
    }

    /// Transformation reducing the second degree below ceil(deg(a)/2).
    fn half_gcd(self, a: &[u64], b: &[u64]) -> Matrix {
        let mid = (a.len() - 1).div_ceil(2);
        let identity = || [[vec![1], vec![]], [vec![], vec![1]]];
        if b.len() <= mid {
            return identity();
        }
        if a.len() < CROSSOVER {
            let (mut a, mut b) = (a.to_vec(), b.to_vec());
            let mut result = identity();
            while b.len() > mid {
                let (q, r) = self.divrem(&a, &b);
                let step = [[vec![], vec![1]], [vec![1], self.scale(&q, self.p - 1)]];
                result = self.matrix_mul(&step, &result);
                (a, b) = (b, r);
            }
            return result;
        }
        let r = self.half_gcd(&a[mid..], &b[mid..]);
        let (c, d) = self.apply(&r, a, b);
        if d.len() <= mid {
            return r;
        }
        let (q, e) = self.divrem(&c, &d);
        let step = [[vec![], vec![1]], [vec![1], self.scale(&q, self.p - 1)]];
        let rr = self.matrix_mul(&step, &r);
        if e.len() <= mid {
            return rr;
        }
        let shift = 2 * mid - (d.len() - 1);
        let next = self.half_gcd(&d[shift..], &e[shift..]);
        self.matrix_mul(&next, &rr)
    }

    /// Monic gcd, and optionally the coefficient of b in its Bezout identity.
    fn euclid(self, a: &[u64], b: &[u64], bezout: bool) -> (Poly, Poly) {
        let (mut a, mut b) = (a.to_vec(), b.to_vec());
        let (mut u, mut v) = (vec![], vec![1]);
        if a.len() < b.len() {
            std::mem::swap(&mut a, &mut b);
            std::mem::swap(&mut u, &mut v);
        }
        while !b.is_empty() {
            if a.len() >= CROSSOVER && b.len() > a.len().div_ceil(2) {
                let m = self.half_gcd(&a, &b);
                (a, b) = self.apply(&m, &a, &b);
                if bezout {
                    (u, v) = self.apply(&m, &u, &v);
                }
                if b.is_empty() {
                    break;
                }
            }
            let (q, r) = self.divrem(&a, &b);
            (a, b) = (b, r);
            if bezout {
                let next = self.sub(&u, &self.mul(&q, &v));
                (u, v) = (v, next);
            }
        }
        if let Some(&lead) = a.last() {
            let li = inv(lead, self.p);
            (
                self.scale(&a, li),
                if bezout { self.scale(&u, li) } else { vec![] },
            )
        } else {
            (vec![], vec![])
        }
    }

    pub(crate) fn gcd(self, a: &[u64], b: &[u64]) -> Poly {
        self.euclid(a, b, false).0
    }

    pub(crate) fn inverse_mod(self, h: &[u64], g: &[u64]) -> Option<Poly> {
        if g.len() <= 1 {
            return None;
        }
        let (d, v) = self.euclid(g, &self.divrem(h, g).1, true);
        (d == [1]).then_some(v)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use polycore::sample::Rng;
    use polycore::Modular;
    use zippel_interp::poly::dense;

    #[test]
    fn fast_arithmetic_matches_classical() {
        let mut rng = Rng::new(77);
        for p in [
            101,
            1_000_000_007,
            4_611_686_018_427_387_847,
            18_446_744_073_709_551_557,
        ] {
            assert!(is_prime(p));
            let field = Arithmetic { p };
            for (n, m) in [
                (1, 1),
                (63, 40),
                (65, 64),
                (130, 67),
                (257, 129),
                (450, 210),
            ] {
                let a: Vec<_> = (0..n).map(|_| rng.nonzero(p)).collect();
                let b: Vec<_> = (0..m).map(|_| rng.nonzero(p)).collect();
                let (aa, bb) = (dense(a.clone(), p), dense(b.clone(), p));
                let raw = |x: &zippel_interp::Dense| {
                    x.0.iter().map(|c| c.residue_mod(p)).collect::<Vec<_>>()
                };
                assert_eq!(field.mul(&a, &b), raw(&(&aa * &bb)));
                let (q, r) = aa.divrem(&bb);
                assert_eq!(field.divrem(&a, &b), (raw(&q), raw(&r)));
                assert_eq!(field.gcd(&a, &b), raw(&aa.gcd(&bb)));
                if let Some(v) = field.inverse_mod(&b, &a) {
                    assert_eq!(field.divrem(&field.mul(&v, &b), &a).1, vec![1]);
                } else {
                    assert!(aa.gcd(&bb).deg() > 0 || aa.deg() == 0);
                }
                let h = [1, 2, 3, 4, 1];
                let ah = field.mul(&a, &h);
                let bh = field.mul(&b, &h);
                assert_eq!(field.gcd(&ah, &bh), field.mul(&field.gcd(&a, &b), &h));
            }
        }
    }

    #[test]
    fn half_gcd_handles_sparse_degree_drops_and_zero_operands() {
        let field = Arithmetic { p: 1_000_000_007 };
        assert!(field.gcd(&[], &[]).is_empty());
        assert_eq!(field.gcd(&[], &[2, 2]), vec![1, 1]);
        assert_eq!(field.gcd(&[2, 2], &[]), vec![1, 1]);
        for (n, m) in [(512, 256), (513, 257), (700, 3), (1024, 511)] {
            let mut a = vec![0; n + 1];
            let mut b = vec![0; m + 1];
            a[0] = field.p - 1;
            a[n] = 1;
            b[0] = field.p - 1;
            b[m] = 1;
            let expected = dense(a.clone(), field.p).gcd(&dense(b.clone(), field.p));
            assert_eq!(dense(field.gcd(&a, &b), field.p), expected);
            assert!(field.inverse_mod(&b, &a).is_none());
        }
    }
}
