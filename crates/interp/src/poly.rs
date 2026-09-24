//! Sparse polynomials over GF(p) in lex order (variable 0 most significant), terms descending.

use crate::modp::{add, inv, mul, pow};
use crate::univariate::{self as uni, Dense};

pub type Exps = Vec<u32>;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ModPoly {
    pub n: usize,
    pub terms: Vec<(Exps, u64)>,
}

impl ModPoly {
    pub fn lm(&self) -> &[u32] {
        &self.terms[0].0
    }

    pub fn scale(&mut self, c: u64, p: u64) {
        self.terms.iter_mut().for_each(|t| t.1 = mul(t.1, c, p));
    }

    pub fn support(&self) -> Vec<Exps> {
        self.terms.iter().map(|t| t.0.clone()).collect()
    }

    pub fn degree(&self, k: usize) -> usize {
        self.terms
            .iter()
            .map(|t| t.0[k] as usize)
            .max()
            .unwrap_or(0)
    }

    /// Coefficients in `x_k` of each monomial in `x_0..x_{k-1}`, assuming later variables are gone.
    pub fn groups(&self, k: usize) -> Vec<(Exps, Dense)> {
        let mut out: Vec<(Exps, Dense)> = Vec::new();
        for (e, c) in &self.terms {
            let mut key = e.clone();
            key[k] = 0;
            if out.last().is_none_or(|g| g.0 != key) {
                out.push((key, Vec::new()));
            }
            let d = &mut out.last_mut().expect("just pushed").1;
            let i = e[k] as usize;
            if d.len() <= i {
                d.resize(i + 1, 0);
            }
            d[i] = *c;
        }
        out
    }

    pub fn from_groups(n: usize, k: usize, groups: Vec<(Exps, Dense)>) -> Self {
        let mut terms = Vec::new();
        for (key, d) in groups {
            for (i, &c) in d.iter().enumerate().rev().filter(|t| *t.1 != 0) {
                let mut e = key.clone();
                e[k] = i as u32;
                terms.push((e, c));
            }
        }
        Self { n, terms }
    }

    /// Content in `GF(p)[x_k]` and the primitive part.
    pub fn primitive(&self, k: usize, p: u64) -> (Dense, Self) {
        let groups = self.groups(k);
        let c = groups
            .iter()
            .fold(Vec::new(), |acc, g| uni::gcd(&acc, &g.1, p));
        if c.len() == 1 {
            return (c, self.clone());
        }
        let groups = groups
            .into_iter()
            .map(|(e, d)| (e, uni::div_rem(&d, &c, p).0))
            .collect();
        (c, Self::from_groups(self.n, k, groups))
    }

    #[must_use]
    pub fn eval_var(&self, k: usize, a: u64, p: u64) -> Self {
        let terms = self
            .groups(k)
            .into_iter()
            .map(|(e, d)| (e, uni::eval(&d, a, p)))
            .filter(|t| t.1 != 0)
            .collect();
        Self { n: self.n, terms }
    }

    /// Substitute `point[i]` for every `x_i` with `i != k`, leaving a dense polynomial in `x_k`.
    pub fn eval_except(&self, k: usize, point: &[u64], p: u64) -> Dense {
        let mut d = vec![0; self.degree(k) + 1];
        for (e, c) in &self.terms {
            let v = e
                .iter()
                .enumerate()
                .filter(|&(i, &x)| i != k && x != 0)
                .fold(*c, |acc, (i, &x)| {
                    mul(acc, pow(point[i], u64::from(x), p), p)
                });
            d[e[k] as usize] = add(d[e[k] as usize], v, p);
        }
        uni::trim(&mut d);
        d
    }

    pub fn eval(&self, x: &[u64], p: u64) -> u64 {
        self.terms.iter().fold(0, |acc, (e, c)| {
            let v = e
                .iter()
                .zip(x)
                .fold(*c, |v, (&d, &xi)| mul(v, pow(xi, u64::from(d), p), p));
            add(acc, v, p)
        })
    }

    #[must_use]
    pub fn monic(mut self, p: u64) -> Self {
        let l = inv(self.terms[0].1, p);
        self.scale(l, p);
        self
    }

    /// Coefficients in the symmetric range, which shows small integers as themselves.
    pub fn show(&self, names: &[&str], p: u64) -> String {
        let terms = self.terms.iter().map(|(e, c)| {
            let c = if *c > p / 2 {
                -i128::from(p - c)
            } else {
                i128::from(*c)
            };
            let vars: Vec<String> = e
                .iter()
                .zip(names)
                .filter(|t| *t.0 > 0)
                .map(|(&d, v)| {
                    if d == 1 {
                        (*v).to_string()
                    } else {
                        format!("{v}^{d}")
                    }
                })
                .collect();
            let vars = vars.join("*");
            let body = match (c.abs(), vars.is_empty()) {
                (1, false) => vars,
                (a, true) => a.to_string(),
                (a, false) => format!("{a}*{vars}"),
            };
            (c < 0, body)
        });
        let out = terms
            .enumerate()
            .fold(String::new(), |out, (i, (neg, t))| match (i, neg) {
                (0, false) => t,
                (0, true) => format!("-{t}"),
                (_, true) => format!("{out} - {t}"),
                (_, false) => format!("{out} + {t}"),
            });
        if out.is_empty() {
            "0".into()
        } else {
            out
        }
    }
}
