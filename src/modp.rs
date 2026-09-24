pub(crate) const fn add(a: u64, b: u64, p: u64) -> u64 {
    let s = a + b;
    if s >= p {
        s - p
    } else {
        s
    }
}

pub(crate) const fn sub(a: u64, b: u64, p: u64) -> u64 {
    if a >= b {
        a - b
    } else {
        a + p - b
    }
}

pub(crate) fn mul(a: u64, b: u64, p: u64) -> u64 {
    (u128::from(a) * u128::from(b) % u128::from(p)) as u64
}

pub(crate) fn pow(mut a: u64, mut e: u64, p: u64) -> u64 {
    let mut r = 1;
    while e > 0 {
        if e & 1 == 1 {
            r = mul(r, a, p);
        }
        a = mul(a, a, p);
        e >>= 1;
    }
    r
}

pub(crate) fn inv(a: u64, p: u64) -> u64 {
    pow(a, p - 2, p)
}

fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    let bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
    if let Some(&b) = bases.iter().find(|&&b| n.is_multiple_of(b)) {
        return n == b;
    }
    let s = (n - 1).trailing_zeros();
    let d = (n - 1) >> s;
    bases.iter().all(|&a| {
        let mut x = pow(a, d, n);
        if x == 1 || x == n - 1 {
            return true;
        }
        (1..s).any(|_| {
            x = mul(x, x, n);
            x == n - 1
        })
    })
}

/// Primes descending from 2^62, so sums of two residues never overflow a `u64`.
#[derive(Debug)]
pub(crate) struct Primes(u64);

impl Primes {
    pub(crate) const fn new() -> Self {
        Self(1 << 62)
    }
}

impl Iterator for Primes {
    type Item = u64;
    fn next(&mut self) -> Option<u64> {
        while self.0 > 3 {
            self.0 -= 1;
            if is_prime(self.0) {
                return Some(self.0);
            }
        }
        None
    }
}

/// `SplitMix64`: deterministic, so every run and snapshot is reproducible.
#[derive(Debug)]
pub(crate) struct Rng(u64);

impl Rng {
    pub(crate) const fn new(seed: u64) -> Self {
        Self(seed)
    }

    pub(crate) const fn nonzero(&mut self, p: u64) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        1 + (z ^ (z >> 31)) % (p - 1)
    }
}
