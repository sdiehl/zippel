//! Sparse interpolation of `gcd(f, g)` in `x_0..x_m` from a known skeleton.
//!
//! Evaluating `x_1..x_m` at the powers `b, b^2, ..` turns every monomial `M` into a geometric
//! sequence `M(b)^j`, so each `x_0` coefficient of the gcd solves a transposed Vandermonde system.
//! Univariate gcd images are only known up to a scalar `m_j`; those unknowns are recovered
//! alongside the coefficients (de Kleine, Monagan, Wittkopf 2005).

use crate::poly::{Exps, ModPoly};
use zippel_interp::modp::{inv, mul, pow, sub, Rng};
use zippel_interp::univariate::{self as uni, Dense};
use zippel_interp::vandermonde;

struct Block {
    e0: u32,
    monos: Vec<Exps>,
    vals: Vec<u64>,
    master: Dense,
}

enum Failure {
    Retry,
    Abort,
}

pub(crate) fn linzip(
    f: &ModPoly,
    g: &ModPoly,
    skeleton: &[Exps],
    m: usize,
    p: u64,
    rng: &mut Rng,
) -> Option<ModPoly> {
    for _ in 0..3 {
        match attempt(f, g, skeleton, m, p, rng) {
            Ok(h) => return Some(h),
            Err(Failure::Abort) => return None,
            Err(Failure::Retry) => {}
        }
    }
    None
}

fn attempt(
    f: &ModPoly,
    g: &ModPoly,
    skeleton: &[Exps],
    m: usize,
    p: u64,
    rng: &mut Rng,
) -> Result<ModPoly, Failure> {
    let b: Vec<u64> = (0..=m).map(|_| rng.nonzero(p)).collect();
    let value = |e: &Exps| (1..=m).fold(1, |acc, i| mul(acc, pow(b[i], u64::from(e[i]), p), p));
    let blocks = blocks(skeleton, &value, p)?;
    let n_max = blocks.iter().map(|bl| bl.monos.len()).max().unwrap_or(0);
    let total: usize = blocks.iter().map(|bl| bl.monos.len()).sum();
    let anchor = blocks.iter().position(|bl| bl.monos.len() == 1);
    let r = blocks.len();
    let mut s = match anchor {
        Some(_) => n_max,
        None if r > 1 => n_max.max((total - 1).div_ceil(r - 1)),
        None => return Err(Failure::Abort),
    } + 1;

    let mut images = Images::new(f, g, &value, p);
    let deg0 = blocks[0].e0 as usize;
    let mut coeffs: Vec<Vec<u64>> = Vec::new();
    loop {
        while coeffs.len() < s {
            let h = images.next().ok_or(Failure::Retry)?;
            if uni::deg(&h) != deg0 {
                return Err(if uni::deg(&h) > deg0 {
                    Failure::Retry
                } else {
                    Failure::Abort
                });
            }
            let expected = h.iter().filter(|&&c| c != 0).count();
            let row: Vec<u64> = blocks.iter().map(|bl| h[bl.e0 as usize]).collect();
            if row.iter().filter(|&&c| c != 0).count() != expected {
                return Err(Failure::Abort);
            }
            coeffs.push(row);
        }
        let scales = match anchor {
            Some(a) => anchored_scales(&blocks[a], a, &coeffs, p).ok_or(Failure::Abort)?,
            None => match solve_scales(&blocks, &coeffs, p) {
                Scales::Solved(v) => v,
                Scales::Deficient if s < 2 * total + 2 => {
                    s += 1;
                    continue;
                }
                Scales::Inconsistent | Scales::Deficient => return Err(Failure::Abort),
            },
        };
        return assemble(&blocks, &coeffs, &scales, p);
    }
}

fn blocks(skeleton: &[Exps], value: &impl Fn(&Exps) -> u64, p: u64) -> Result<Vec<Block>, Failure> {
    let mut out: Vec<Block> = Vec::new();
    for e in skeleton {
        if out.last().is_none_or(|bl| bl.e0 != e[0]) {
            out.push(Block {
                e0: e[0],
                monos: Vec::new(),
                vals: Vec::new(),
                master: vec![1],
            });
        }
        let bl = out.last_mut().expect("just pushed");
        bl.monos.push(e.clone());
        bl.vals.push(value(e));
    }
    for bl in &mut out {
        let mut sorted = bl.vals.clone();
        sorted.sort_unstable();
        sorted.dedup();
        if sorted.len() != bl.vals.len() {
            return Err(Failure::Retry);
        }
        bl.master = bl.vals.iter().fold(vec![1], |acc, &v| {
            uni::mul_poly(&acc, &[sub(0, v, p), 1], p)
        });
    }
    Ok(out)
}

/// Univariate images `gcd(f, g)(x_0, b^j)` for `j = 1, 2, ..`, made monic.
struct Images {
    f: Vec<(usize, u64, u64)>,
    g: Vec<(usize, u64, u64)>,
    df: usize,
    dg: usize,
    p: u64,
}

impl Images {
    fn new(f: &ModPoly, g: &ModPoly, value: &impl Fn(&Exps) -> u64, p: u64) -> Self {
        let prep = |h: &ModPoly| {
            h.terms
                .iter()
                .map(|(e, c)| (e[0] as usize, *c, value(e)))
                .collect()
        };
        Self {
            f: prep(f),
            g: prep(g),
            df: f.degree(0),
            dg: g.degree(0),
            p,
        }
    }

    fn step(terms: &mut [(usize, u64, u64)], d: usize, p: u64) -> Option<Dense> {
        let mut u = vec![0; d + 1];
        for (e, c, v) in terms.iter_mut() {
            *c = mul(*c, *v, p);
            u[*e] = zippel_interp::modp::add(u[*e], *c, p);
        }
        uni::trim(&mut u);
        (uni::deg(&u) == d && !u.is_empty()).then_some(u)
    }

    fn next(&mut self) -> Option<Dense> {
        let uf = Self::step(&mut self.f, self.df, self.p)?;
        let ug = Self::step(&mut self.g, self.dg, self.p)?;
        Some(uni::gcd(&uf, &ug, self.p))
    }
}

/// A block with a single monomial `c M` fixes every scale: `m_j = M(b)^j / h_j`.
fn anchored_scales(bl: &Block, a: usize, coeffs: &[Vec<u64>], p: u64) -> Option<Vec<u64>> {
    let v = bl.vals[0];
    let mut vj = 1;
    coeffs
        .iter()
        .map(|row| {
            vj = mul(vj, v, p);
            (row[a] != 0).then(|| mul(vj, inv(row[a], p), p))
        })
        .collect()
}

enum Scales {
    Solved(Vec<u64>),
    Inconsistent,
    Deficient,
}

/// Each block's master polynomial annihilates its Vandermonde columns, which leaves a linear
/// system in the scales alone. `m_1 = 1` fixes the overall normalization.
fn solve_scales(blocks: &[Block], coeffs: &[Vec<u64>], p: u64) -> Scales {
    let s = coeffs.len();
    let mut rows: Vec<Vec<u64>> = Vec::new();
    for (i, bl) in blocks.iter().enumerate() {
        let n = bl.monos.len();
        for shift in 0..s - n {
            let mut row = vec![0; s];
            for (t, &pt) in bl.master.iter().enumerate() {
                row[shift + t] = mul(pt, coeffs[shift + t][i], p);
            }
            rows.push(row);
        }
    }
    let mut pivot_row = 0;
    for col in 1..s {
        let Some(r) = (pivot_row..rows.len()).find(|&r| rows[r][col] != 0) else {
            return Scales::Deficient;
        };
        rows.swap(pivot_row, r);
        let l = inv(rows[pivot_row][col], p);
        rows[pivot_row].iter_mut().for_each(|x| *x = mul(*x, l, p));
        let pr = rows[pivot_row].clone();
        for (k, row) in rows.iter_mut().enumerate() {
            if k != pivot_row && row[col] != 0 {
                let c = row[col];
                row.iter_mut()
                    .zip(&pr)
                    .for_each(|(x, y)| *x = sub(*x, mul(c, *y, p), p));
            }
        }
        pivot_row += 1;
    }
    if rows[pivot_row..].iter().any(|r| r[0] != 0) {
        return Scales::Inconsistent;
    }
    let mut m = vec![1; s];
    for (col, row) in (1..s).zip(&rows) {
        m[col] = sub(0, row[0], p);
    }
    Scales::Solved(m)
}

fn assemble(
    blocks: &[Block],
    coeffs: &[Vec<u64>],
    scales: &[u64],
    p: u64,
) -> Result<ModPoly, Failure> {
    let mut terms = Vec::new();
    for (i, bl) in blocks.iter().enumerate() {
        let w: Vec<u64> = coeffs
            .iter()
            .zip(scales)
            .map(|(row, &m)| mul(row[i], m, p))
            .collect();
        let cs = vandermonde::solve(&bl.vals, &bl.master, &w, p);
        let mut pw = bl.vals.clone();
        for &wj in &w {
            let lhs = cs.iter().zip(&pw).fold(0, |acc, (&c, &v)| {
                zippel_interp::modp::add(acc, mul(c, v, p), p)
            });
            if lhs != wj {
                return Err(Failure::Abort);
            }
            pw.iter_mut()
                .zip(&bl.vals)
                .for_each(|(x, &v)| *x = mul(*x, v, p));
        }
        terms.extend(
            bl.monos
                .iter()
                .zip(cs)
                .filter(|t| t.1 != 0)
                .map(|(e, c)| (e.clone(), c)),
        );
    }
    let n = blocks[0].monos[0].len();
    Ok(ModPoly { n, terms })
}
