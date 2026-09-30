//! Sparse polynomials over GF(p) in lex order (variable 0 most significant), terms descending.

use polycore::modp::{add, inv, mul};
use polycore::{Fp, Modular, Uni};
use std::iter::successors;

pub type Exps = Vec<u32>;

/// A dense univariate polynomial over GF(p).
pub type Dense = Uni<Fp>;

pub fn dense(cs: impl IntoIterator<Item = u64>, p: u64) -> Dense {
    Uni::new(cs.into_iter().map(|c| Fp::new(c, p)).collect())
}

/// `x[i]^j` for `j <= d[i]`, so evaluating a monomial takes one multiply per variable.
pub fn power_table(x: &[u64], d: &[usize], p: u64) -> Vec<Vec<u64>> {
    x.iter()
        .zip(d)
        .map(|(&xi, &di)| {
            successors(Some(1 % p), |&v| Some(mul(v, xi, p)))
                .take(di + 1)
                .collect()
        })
        .collect()
}

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

    /// The degree in every variable.
    pub fn degrees(&self) -> Vec<usize> {
        let mut d = vec![0; self.n];
        for (e, _) in &self.terms {
            d.iter_mut()
                .zip(e)
                .for_each(|(d, &x)| *d = (*d).max(x as usize));
        }
        d
    }

    /// Coefficients in `x_k` of each monomial in `x_0..x_{k-1}`, assuming later variables are gone.
    pub fn groups(&self, k: usize, p: u64) -> Vec<(Exps, Dense)> {
        let mut out: Vec<(Exps, Vec<u64>)> = Vec::new();
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
        out.into_iter().map(|(e, d)| (e, dense(d, p))).collect()
    }

    pub fn from_groups(n: usize, k: usize, groups: Vec<(Exps, Dense)>, p: u64) -> Self {
        let mut terms = Vec::new();
        for (key, d) in groups {
            for (i, c) in d.0.iter().enumerate().rev() {
                let c = c.residue_mod(p);
                if c != 0 {
                    let mut e = key.clone();
                    e[k] = i as u32;
                    terms.push((e, c));
                }
            }
        }
        Self { n, terms }
    }

    /// Content in `GF(p)[x_k]` and the primitive part.
    pub fn primitive(&self, k: usize, p: u64) -> (Dense, Self) {
        let groups = self.groups(k, p);
        let c = groups
            .iter()
            .try_fold(Uni::zero(), |acc, g| {
                let c = acc.gcd(&g.1);
                if c.deg() == 0 {
                    Err(c)
                } else {
                    Ok(c)
                }
            })
            .unwrap_or_else(|c| c);
        if c.deg() == 0 {
            return (c, self.clone());
        }
        let groups = groups.into_iter().map(|(e, d)| (e, &d / &c)).collect();
        (c, Self::from_groups(self.n, k, groups, p))
    }

    #[must_use]
    pub fn eval_var(&self, k: usize, a: u64, p: u64) -> Self {
        let pw = &power_table(&[a], &[self.degree(k)], p)[0];
        let mut terms: Vec<(Exps, u64)> = Vec::new();
        for (e, c) in &self.terms {
            let v = mul(*c, pw[e[k] as usize], p);
            match terms.last_mut() {
                Some((l, s)) if l[..k] == e[..k] => *s = add(*s, v, p),
                _ => {
                    let mut l = e.clone();
                    l[k] = 0;
                    terms.push((l, v));
                }
            }
        }
        terms.retain(|t| t.1 != 0);
        Self { n: self.n, terms }
    }

    /// The leading coefficient in `GF(p)[x_k]`, assuming later variables are gone.
    pub fn lead(&self, k: usize, p: u64) -> Dense {
        let key = &self.terms[0].0[..k];
        let mut d = vec![0; self.terms[0].0[k] as usize + 1];
        for (e, c) in self.terms.iter().take_while(|t| &t.0[..k] == key) {
            d[e[k] as usize] = *c;
        }
        dense(d, p)
    }

    /// Substitute `point[i]` for every `x_i` with `i != k`, leaving a dense polynomial in `x_k`.
    pub fn eval_except(&self, k: usize, point: &[u64], p: u64) -> Dense {
        let mut degs = self.degrees();
        let mut d = vec![0; degs[k] + 1];
        degs[k] = 0;
        let pw = power_table(point, &degs, p);
        for (e, c) in &self.terms {
            let v = e
                .iter()
                .enumerate()
                .filter(|&(i, &x)| i != k && x != 0)
                .fold(*c, |acc, (i, &x)| mul(acc, pw[i][x as usize], p));
            d[e[k] as usize] = add(d[e[k] as usize], v, p);
        }
        dense(d, p)
    }

    pub fn eval(&self, x: &[u64], p: u64) -> u64 {
        let pw = power_table(x, &self.degrees(), p);
        self.terms.iter().fold(0, |acc, (e, c)| {
            let v = e
                .iter()
                .zip(&pw)
                .fold(*c, |v, (&d, xi)| mul(v, xi[d as usize], p));
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
