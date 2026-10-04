//! Polynomial relation discovery and block-triangular evaluation from numeric reductions.
//!
//! At each sample, reductions to masters impose linear constraints on a bounded
//! polynomial ansatz. Its null space supplies relations between a block and simpler
//! integrals. Search is performed once; other primes refit the learned ansatz and
//! pivot structure. Every candidate is checked on independent probes. This is a
//! bounded Monte Carlo learner, not the full adaptive/weighted search used by Blade.

use polycore::modp::{add, mul};
use polycore::modp_echelon::{Echelon, Row};
use polycore::sample::{hash, point};
use polycore::Lead;
use std::collections::HashMap;
use zippel_interp::{Exps, ModPoly};

type Oracle<'a> = &'a dyn Fn(&[u64], u64) -> Option<Vec<u64>>;

#[derive(Clone, Copy, Debug)]
pub struct Search {
    pub block_size: usize,
    pub max_degree: u32,
    pub max_unknowns: usize,
    /// Hard budget on distinct oracle calls, including rejected probes and validation.
    pub max_probes: usize,
}
impl Default for Search {
    fn default() -> Self {
        Self {
            block_size: 1,
            max_degree: 3,
            max_unknowns: 2048,
            max_probes: 4096,
        }
    }
}

#[derive(Clone, Debug)]
struct Shape {
    start: usize,
    end: usize,
    monomials: Vec<Exps>,
    pivots: Vec<usize>,
    selected: Vec<usize>,
}

/// Learned relation supports and pivot pattern, reusable across primes.
#[derive(Clone, Debug)]
pub struct BlockPlan {
    n: usize,
    targets: usize,
    masters: usize,
    blocks: Vec<Shape>,
    search: Search,
}

#[derive(Clone, Debug)]
struct Block {
    start: usize,
    end: usize,
    relations: Vec<Vec<ModPoly>>,
}

/// A fitted block-triangular form over one prime. Evaluations do not call the IBP oracle.
#[derive(Clone, Debug)]
pub struct BlockForm {
    n: usize,
    targets: usize,
    masters: usize,
    p: u64,
    blocks: Vec<Block>,
}

struct Samples<'a> {
    oracle: Oracle<'a>,
    p: u64,
    n: usize,
    width: usize,
    limit: usize,
    values: HashMap<Vec<u64>, Option<Vec<u64>>>,
}
impl Samples<'_> {
    fn get(&mut self, x: &[u64]) -> Option<Vec<u64>> {
        if let Some(v) = self.values.get(x) {
            return v.clone();
        }
        if self.values.len() >= self.limit {
            return None;
        }
        let v = (self.oracle)(x, self.p)
            .filter(|v| v.len() == self.width && v.iter().all(|&c| c < self.p));
        self.values.insert(x.to_vec(), v.clone());
        v
    }
    fn at(&mut self, seed: u64, i: u64) -> Option<(Vec<u64>, Vec<u64>)> {
        let x = (0..self.n as u64)
            .map(|j| point(&[seed, i, j], self.p))
            .collect::<Vec<_>>();
        let y = self.get(&x)?;
        Some((x, y))
    }
}

impl BlockPlan {
    /// Learn reductions of `targets` integrals in increasing complexity order.
    ///
    /// The oracle returns `targets * masters` coefficients, target-major. The caller
    /// must use one stable master basis at every sample and prime. `None` means the
    /// budgets/degree bound were insufficient, or samples/pivots were unlucky.
    pub fn learn(
        oracle: impl Fn(&[u64], u64) -> Option<Vec<u64>>,
        n: usize,
        targets: usize,
        masters: usize,
        p: u64,
        seed: u64,
        search: Search,
    ) -> Option<(Self, BlockForm)> {
        if search.block_size == 0 || targets == 0 || masters == 0 || n == 0 {
            return None;
        }
        let mut samples = Samples {
            oracle: &oracle,
            p,
            n,
            width: targets.checked_mul(masters)?,
            limit: search.max_probes,
            values: HashMap::new(),
        };
        let mut shapes = Vec::new();
        let mut blocks = Vec::new();
        for start in (0..targets).step_by(search.block_size) {
            let end = (start + search.block_size).min(targets);
            let mut found = None;
            for degree in 0..=search.max_degree {
                let monomials = monomials(n, degree, search.max_unknowns / (masters + end))?;
                let unknowns = monomials.len() * (masters + end);
                let Some((echelon, basis)) =
                    relations(&mut samples, &monomials, end, masters, seed, None)
                else {
                    continue;
                };
                if echelon.rank() == unknowns {
                    continue;
                }
                let shape = Shape {
                    start,
                    end,
                    monomials,
                    pivots: echelon.pivots(),
                    selected: vec![],
                };
                let all = polynomials(&basis, &shape, masters, n);
                let mut selected = None;
                for attempt in 0..8 {
                    let x = (0..n as u64)
                        .map(|j| point(&[seed, 101, attempt, j], p))
                        .collect::<Vec<_>>();
                    if let Some(indices) = choose(&all, &x, p, masters + start, masters + end) {
                        selected = Some(indices);
                        break;
                    }
                }
                let Some(selected) = selected else { continue };
                let shape = Shape { selected, ..shape };
                let block = Block {
                    start,
                    end,
                    relations: shape.selected.iter().map(|&i| all[i].clone()).collect(),
                };
                if !validate(&mut samples, &block, masters, hash(&[seed, 102])) {
                    continue;
                }
                found = Some((shape, block));
                break;
            }
            let (shape, block) = found?;
            shapes.push(shape);
            blocks.push(block);
        }
        let form = BlockForm {
            n,
            targets,
            masters,
            p,
            blocks,
        };
        let plan = Self {
            n,
            targets,
            masters,
            blocks: shapes,
            search,
        };
        Some((plan, form))
    }

    /// Fit the learned relations at another prime and verify against fresh IBP probes.
    pub fn fit(
        &self,
        oracle: impl Fn(&[u64], u64) -> Option<Vec<u64>>,
        p: u64,
        seed: u64,
    ) -> Option<BlockForm> {
        let mut samples = Samples {
            oracle: &oracle,
            p,
            n: self.n,
            width: self.targets * self.masters,
            limit: self.search.max_probes,
            values: HashMap::new(),
        };
        let mut blocks = Vec::new();
        for shape in &self.blocks {
            let (e, basis) = relations(
                &mut samples,
                &shape.monomials,
                shape.end,
                self.masters,
                seed,
                Some(&shape.pivots),
            )?;
            if e.pivots() != shape.pivots {
                return None;
            }
            let all = polynomials(&basis, shape, self.masters, self.n);
            let block = Block {
                start: shape.start,
                end: shape.end,
                relations: shape
                    .selected
                    .iter()
                    .map(|&i| all.get(i).cloned())
                    .collect::<Option<_>>()?,
            };
            if !validate(&mut samples, &block, self.masters, hash(&[seed, 102])) {
                return None;
            }
            blocks.push(block);
        }
        Some(BlockForm {
            n: self.n,
            targets: self.targets,
            masters: self.masters,
            p,
            blocks,
        })
    }
}

impl BlockForm {
    pub const fn blocks(&self) -> usize {
        self.blocks.len()
    }
    pub fn terms(&self) -> usize {
        self.blocks
            .iter()
            .flat_map(|b| &b.relations)
            .flatten()
            .map(|g| g.terms.len())
            .sum()
    }

    /// All target reductions. A singular block is a rejected probe (`None`).
    pub fn eval(&self, x: &[u64]) -> Option<Vec<u64>> {
        if x.len() != self.n || x.iter().any(|&v| v >= self.p) {
            return None;
        }
        let mut values = vec![vec![0; self.masters]; self.masters + self.targets];
        for (i, row) in values.iter_mut().take(self.masters).enumerate() {
            row[i] = 1;
        }
        for block in &self.blocks {
            let mut e = Echelon::new(Lead::High, self.p);
            for relation in &block.relations {
                let row: Row = relation
                    .iter()
                    .enumerate()
                    .map(|(i, g)| (i, g.eval(x, self.p)))
                    .collect();
                let pivot = e.insert(&row)?;
                if pivot < self.masters + block.start || pivot >= self.masters + block.end {
                    return None;
                }
            }
            for t in block.start..block.end {
                let mut value = vec![0; self.masters];
                for (j, c) in e.solve(self.masters + t).0 {
                    if j >= self.masters + block.start {
                        return None;
                    }
                    for (out, &v) in value.iter_mut().zip(&values[j]) {
                        *out = add(*out, mul(c, v, self.p), self.p);
                    }
                }
                values[self.masters + t] = value;
            }
        }
        Some(values.into_iter().skip(self.masters).flatten().collect())
    }
}

fn monomials(n: usize, degree: u32, limit: usize) -> Option<Vec<Exps>> {
    let mut out = vec![vec![]];
    for _ in 0..n {
        let mut next = Vec::new();
        for e in out {
            for d in 0..=degree - e.iter().sum::<u32>() {
                if next.len() >= limit {
                    return None;
                }
                let mut v = e.clone();
                v.push(d);
                next.push(v);
            }
        }
        out = next;
    }
    Some(out)
}

fn constraints(
    e: &mut Echelon,
    x: &[u64],
    y: &[u64],
    monomials: &[Exps],
    end: usize,
    masters: usize,
    p: u64,
) {
    let powers: Vec<_> = monomials
        .iter()
        .map(|exps| {
            exps.iter().zip(x).fold(1, |acc, (&d, &v)| {
                mul(acc, polycore::modp::pow(v, u64::from(d), p), p)
            })
        })
        .collect();
    for m in 0..masters {
        let mut row = Vec::new();
        for j in 0..masters + end {
            let c = if j < masters {
                u64::from(j == m)
            } else {
                y[(j - masters) * masters + m]
            };
            if c == 0 {
                continue;
            }
            row.extend(
                powers
                    .iter()
                    .enumerate()
                    .map(|(k, &v)| (j * powers.len() + k, mul(c, v, p))),
            );
        }
        e.insert(&row);
    }
}

fn relations(
    samples: &mut Samples<'_>,
    monomials: &[Exps],
    end: usize,
    masters: usize,
    seed: u64,
    expected: Option<&[usize]>,
) -> Option<(Echelon, Vec<Row>)> {
    let unknowns = monomials.len() * (masters + end);
    let mut e = Echelon::new(Lead::Low, samples.p);
    let mut stable = 0;
    for i in 0..samples.limit as u64 {
        let Some((x, y)) = samples.at(hash(&[seed, 100]), i) else {
            continue;
        };
        let before = e.rank();
        constraints(&mut e, &x, &y, monomials, end, masters, samples.p);
        if e.rank() == unknowns {
            return None;
        }
        if let Some(pivots) = expected {
            if e.rank() > pivots.len() {
                return None;
            }
            if e.rank() == pivots.len() {
                return Some((e.clone(), e.nullspace(unknowns)));
            }
        } else {
            stable = if e.rank() == before { stable + 1 } else { 0 };
            if stable == 3 {
                return Some((e.clone(), e.nullspace(unknowns)));
            }
        }
    }
    None
}

fn polynomials(basis: &[Row], shape: &Shape, masters: usize, n: usize) -> Vec<Vec<ModPoly>> {
    basis
        .iter()
        .map(|r| {
            let mut out = vec![ModPoly { n, terms: vec![] }; masters + shape.end];
            for &(j, c) in r {
                out[j / shape.monomials.len()]
                    .terms
                    .push((shape.monomials[j % shape.monomials.len()].clone(), c));
            }
            out
        })
        .collect()
}

fn choose(
    relations: &[Vec<ModPoly>],
    x: &[u64],
    p: u64,
    start: usize,
    end: usize,
) -> Option<Vec<usize>> {
    let mut e = Echelon::new(Lead::High, p);
    let mut out = Vec::new();
    for (i, r) in relations.iter().enumerate() {
        let row: Row = r
            .iter()
            .enumerate()
            .map(|(j, g)| (j, g.eval(x, p)))
            .collect();
        let reduced = e.reduce(&row, false).0;
        if reduced.first().is_some_and(|&(j, _)| j >= start && j < end) {
            e.insert(&row);
            out.push(i);
        }
        if out.len() == end - start {
            return Some(out);
        }
    }
    None
}

fn validate(samples: &mut Samples<'_>, block: &Block, masters: usize, seed: u64) -> bool {
    let mut good = 0;
    for i in 0..24 {
        let Some((x, y)) = samples.at(seed, i) else {
            continue;
        };
        if choose(
            &block.relations,
            &x,
            samples.p,
            masters + block.start,
            masters + block.end,
        )
        .is_none()
        {
            continue;
        }
        for r in &block.relations {
            for m in 0..masters {
                let sum = r.iter().enumerate().fold(0, |acc, (j, g)| {
                    let v = if j < masters {
                        u64::from(j == m)
                    } else {
                        y[(j - masters) * masters + m]
                    };
                    add(acc, mul(g.eval(&x, samples.p), v, samples.p), samples.p)
                });
                if sum != 0 {
                    return false;
                }
            }
        }
        good += 1;
        if good == 3 {
            return true;
        }
    }
    false
}

/// Learn once, refit once per prime, and lift coefficients using only block probes.
/// The oracle contract and target ordering are the same as [`BlockPlan::learn`].
pub fn lift(
    oracle: impl Fn(&[u64], u64) -> Option<Vec<u64>> + Sync,
    n: usize,
    targets: usize,
    masters: usize,
    seed: u64,
    search: Search,
) -> Option<Vec<zippel_lift::Fraction>> {
    use std::sync::{Arc, Mutex};
    let p = polycore::modp::Primes::new().next()?;
    let (plan, first) = BlockPlan::learn(&oracle, n, targets, masters, p, seed, search)?;
    let forms = Mutex::new(HashMap::from([(p, Some(Arc::new(first)))]));
    zippel_lift::lift(
        |x: &[u64], p| {
            let form = {
                let mut cache = forms.lock().unwrap();
                cache
                    .entry(p)
                    .or_insert_with(|| plan.fit(&oracle, p, seed).map(Arc::new))
                    .clone()?
            };
            form.eval(x)
        },
        n,
        seed,
    )
}
