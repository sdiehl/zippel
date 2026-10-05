//! Polynomial relation discovery and block-triangular evaluation from numeric reductions.
//!
//! At each sample, reductions to masters impose linear constraints on a bounded
//! polynomial ansatz. Its null space supplies relations between a block and simpler
//! integrals. Search is performed once; other primes refit the learned ansatz and
//! pivot structure. Every candidate is checked on independent probes. This is a
//! bounded Monte Carlo learner with adaptive blocks, degree weights, and intermediates.

use polycore::modp::{add, mul, sub, try_inv};
use polycore::modp_echelon::{Echelon, Row};
use polycore::sample::{hash, point};
use polycore::Lead;
use std::collections::{HashMap, HashSet};
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
    ansatz: Vec<(usize, Exps)>,
    pivots: Vec<usize>,
    selected: Vec<usize>,
}

/// Learned relation supports and pivot pattern, reusable across primes.
#[derive(Clone, Debug)]
pub struct BlockPlan {
    n: usize,
    target_order: Vec<usize>,
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
    degrees: Vec<usize>,
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
    fn at_slice(
        &mut self,
        seed: u64,
        i: u64,
        slice: Option<usize>,
    ) -> Option<(Vec<u64>, Vec<u64>)> {
        let x = (0..self.n)
            .map(|j| {
                if slice.is_some_and(|k| k != j) {
                    point(&[103, j as u64], self.p)
                } else {
                    point(&[seed, i, j as u64], self.p)
                }
            })
            .collect::<Vec<_>>();
        let y = self.get(&x)?;
        Some((x, y))
    }
}

impl BlockPlan {
    /// Integral IDs in oracle/output order, strictly increasing by complexity.
    pub fn targets(&self) -> &[usize] {
        &self.target_order
    }
    /// Learn reductions of integral IDs in strictly increasing complexity order.
    ///
    /// The oracle returns `targets.len() * masters` coefficients, target-major. The caller
    /// must use one stable master basis at every sample and prime. `None` means the
    /// budgets/degree bound were insufficient, or samples/pivots were unlucky.
    pub fn learn(
        oracle: impl Fn(&[u64], u64) -> Option<Vec<u64>>,
        n: usize,
        targets: &[usize],
        masters: usize,
        p: u64,
        seed: u64,
        search: Search,
    ) -> Option<(Self, BlockForm)> {
        if targets.windows(2).any(|w| w[0] >= w[1]) {
            return None;
        }
        let target_order = targets.to_vec();
        let targets = targets.len();

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
                let ansatz: Vec<_> = (0..masters + end)
                    .flat_map(|j| monomials.iter().cloned().map(move |e| (j, e)))
                    .collect();
                let unknowns = ansatz.len();
                let Some((echelon, basis)) =
                    relations(&mut samples, &ansatz, masters, seed, None, None)
                else {
                    continue;
                };
                if echelon.rank() == unknowns {
                    continue;
                }
                let shape = Shape {
                    start,
                    end,
                    ansatz,
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
            degrees: degree_bounds(&blocks, n),
            blocks,
        };
        let plan = Self {
            n,
            target_order,
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
                &shape.ansatz,
                self.masters,
                seed,
                Some(&shape.pivots),
                None,
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
            degrees: degree_bounds(&blocks, self.n),
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

    pub const fn prime(&self) -> u64 {
        self.p
    }

    /// All target reductions. Powers and elimination scratch are shared by all blocks.
    /// A singular block is a rejected probe (`None`).
    pub fn eval(&self, x: &[u64]) -> Option<Vec<u64>> {
        if x.len() != self.n || x.iter().any(|&v| v >= self.p) {
            return None;
        }
        let powers = polycore::evaluation::power_table(x, &self.degrees, self.p);
        let evaluate = |g: &ModPoly| {
            g.terms.iter().fold(0, |sum, (e, c)| {
                let v = e
                    .iter()
                    .zip(&powers)
                    .fold(*c, |v, (&d, powers)| mul(v, powers[d as usize], self.p));
                add(sum, v, self.p)
            })
        };
        let mut values = vec![vec![0; self.masters]; self.masters + self.targets];
        for (i, row) in values.iter_mut().take(self.masters).enumerate() {
            row[i] = 1;
        }
        let max = self
            .blocks
            .iter()
            .map(|b| b.end - b.start)
            .max()
            .unwrap_or(0);
        let mut scratch = vec![0; max * (max + self.masters)];
        for block in &self.blocks {
            let size = block.end - block.start;
            let width = size + self.masters;
            scratch.fill(0);
            for (i, relation) in block.relations.iter().enumerate() {
                for (j, g) in relation
                    .iter()
                    .enumerate()
                    .filter(|(_, g)| !g.terms.is_empty())
                {
                    let c = evaluate(g);
                    if j >= self.masters + block.start {
                        scratch[i * width + j - self.masters - block.start] = c;
                    } else {
                        for (k, &v) in values[j].iter().enumerate() {
                            let at = i * width + size + k;
                            scratch[at] = sub(scratch[at], mul(c, v, self.p), self.p);
                        }
                    }
                }
            }
            for k in 0..size {
                let pivot = (k..size).find(|&i| scratch[i * width + k] != 0)?;
                for j in 0..width {
                    scratch.swap(k * width + j, pivot * width + j);
                }
                let inverse = try_inv(scratch[k * width + k], self.p)?;
                for j in k..width {
                    scratch[k * width + j] = mul(scratch[k * width + j], inverse, self.p);
                }
                for i in 0..size {
                    if i != k {
                        let c = scratch[i * width + k];
                        if c != 0 {
                            for j in k..width {
                                scratch[i * width + j] = sub(
                                    scratch[i * width + j],
                                    mul(c, scratch[k * width + j], self.p),
                                    self.p,
                                );
                            }
                        }
                    }
                }
            }
            for i in 0..size {
                values[self.masters + block.start + i]
                    .copy_from_slice(&scratch[i * width + size..(i + 1) * width]);
            }
        }
        Some(values.into_iter().skip(self.masters).flatten().collect())
    }
}

fn degree_bounds(blocks: &[Block], n: usize) -> Vec<usize> {
    let mut degrees = vec![0; n];
    for (e, _) in blocks
        .iter()
        .flat_map(|b| &b.relations)
        .flatten()
        .flat_map(|g| &g.terms)
    {
        for (d, &v) in degrees.iter_mut().zip(e) {
            *d = (*d).max(v as usize);
        }
    }
    degrees
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
    ansatz: &[(usize, Exps)],
    masters: usize,
    p: u64,
) {
    let powers: Vec<_> = ansatz
        .iter()
        .map(|(_, exps)| {
            exps.iter().zip(x).fold(1, |acc, (&d, &v)| {
                mul(acc, polycore::modp::pow(v, u64::from(d), p), p)
            })
        })
        .collect();
    for m in 0..masters {
        let row: Row = ansatz
            .iter()
            .zip(&powers)
            .enumerate()
            .filter_map(|(k, (&(j, _), &power))| {
                let c = if j < masters {
                    u64::from(j == m)
                } else {
                    y[(j - masters) * masters + m]
                };
                (c != 0).then(|| (k, mul(c, power, p)))
            })
            .collect();
        e.insert(&row);
    }
}

fn relations(
    samples: &mut Samples<'_>,
    ansatz: &[(usize, Exps)],
    masters: usize,
    seed: u64,
    expected: Option<&[usize]>,
    slice: Option<usize>,
) -> Option<(Echelon, Vec<Row>)> {
    let unknowns = ansatz.len();
    let mut e = Echelon::new(Lead::Low, samples.p);
    let mut stable = 0;
    let mut distinct = HashSet::new();
    for i in 0..samples.limit as u64 {
        let Some((x, y)) = samples.at_slice(hash(&[seed, 100]), i, slice) else {
            continue;
        };
        if !distinct.insert(x.clone()) {
            continue;
        }
        let before = e.rank();
        constraints(&mut e, &x, &y, ansatz, masters, samples.p);
        if e.rank() == unknowns {
            return None;
        }
        stable = if e.rank() == before { stable + 1 } else { 0 };
        if let Some(pivots) = expected {
            if stable >= 3 || e.rank() > pivots.len() {
                return None;
            }
            if e.rank() == pivots.len() {
                return Some((e.clone(), e.nullspace(unknowns)));
            }
        } else if stable == 3 {
            return Some((e.clone(), e.nullspace(unknowns)));
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
                out[shape.ansatz[j].0]
                    .terms
                    .push((shape.ansatz[j].1.clone(), c));
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

fn validate_slice(
    samples: &mut Samples<'_>,
    block: &Block,
    masters: usize,
    seed: u64,
    slice: Option<usize>,
) -> bool {
    let mut good = 0;
    for i in 0..24 {
        let Some((x, y)) = samples.at_slice(seed, i, slice) else {
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
    targets: &[usize],
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

/// Adaptive polynomial ansätze. Empty weights are inferred on coordinate slices.
#[derive(Clone, Debug)]
pub struct AdaptiveSearch {
    pub limits: Search,
    pub variable_weights: Vec<u32>,
    /// Offsets in the weighted degree bound, masters first then oracle integrals.
    pub integral_weights: Vec<i32>,
    /// Additional independent degree bounds on groups of variable indices.
    pub groups: Vec<(Vec<usize>, u32)>,
    pub max_block_size: usize,
    /// Number of replay dependencies admitted as automatic intermediate integrals.
    pub max_intermediates: usize,
}
impl Default for AdaptiveSearch {
    fn default() -> Self {
        Self {
            limits: Search {
                max_degree: 8,
                ..Search::default()
            },
            variable_weights: vec![],
            integral_weights: vec![],
            groups: vec![],
            max_block_size: 4,
            max_intermediates: 32,
        }
    }
}

/// Search costs include degree inference, unsuccessful ansätze, and validation.
#[derive(Clone, Debug, Default)]
pub struct SearchReport {
    pub oracle_probes: usize,
    pub variable_weights: Vec<u32>,
    pub block_sizes: Vec<usize>,
    pub unknowns: Vec<usize>,
}

impl BlockPlan {
    /// Infer variable weights, raise relation degrees, and enlarge blocks as needed.
    ///
    /// Both inferred/user weights and uniform weights are tried, so an unlucky
    /// slice cannot exclude a relation that fits the configured uniform bounds.
    pub fn learn_adaptive(
        oracle: impl Fn(&[u64], u64) -> Option<Vec<u64>>,
        n: usize,
        targets: &[usize],
        masters: usize,
        p: u64,
        seed: u64,
        search: &AdaptiveSearch,
    ) -> Option<(Self, BlockForm, SearchReport)> {
        if targets.windows(2).any(|w| w[0] >= w[1]) {
            return None;
        }
        let target_order = targets.to_vec();
        let targets = targets.len();

        if n == 0
            || targets == 0
            || masters == 0
            || search.max_block_size == 0
            || (!search.variable_weights.is_empty()
                && (search.variable_weights.len() != n || search.variable_weights.contains(&0)))
            || (!search.integral_weights.is_empty()
                && search.integral_weights.len() != masters + targets)
            || search.groups.iter().any(|(g, _)| g.iter().any(|&k| k >= n))
        {
            return None;
        }
        let mut samples = Samples {
            oracle: &oracle,
            p,
            n,
            width: targets.checked_mul(masters)?,
            limit: search.limits.max_probes,
            values: HashMap::new(),
        };
        let weights = if search.variable_weights.is_empty() {
            infer_weights(&mut samples, masters, seed, search)
        } else {
            search.variable_weights.clone()
        };
        let mut shapes = Vec::new();
        let mut blocks = Vec::new();
        let mut start = 0;
        let mut report = SearchReport {
            variable_weights: weights.clone(),
            ..SearchReport::default()
        };
        while start < targets {
            let mut found = None;
            for size in 1..=search.max_block_size.min(targets - start) {
                let end = start + size;
                if let Some(pair) = search_block(
                    &mut samples,
                    start,
                    end,
                    masters,
                    seed,
                    search,
                    &weights,
                    None,
                ) {
                    found = Some(pair);
                    break;
                }
            }
            let (shape, block) = found?;
            let (shape, block) = compact(&mut samples, shape, block, masters, seed);
            report.block_sizes.push(block.end - block.start);
            report.unknowns.push(shape.ansatz.len());
            start = block.end;
            shapes.push(shape);
            blocks.push(block);
        }
        report.oracle_probes = samples.values.len();
        Some((
            Self {
                n,
                target_order,
                targets,
                masters,
                blocks: shapes,
                search: search.limits,
            },
            BlockForm {
                n,
                targets,
                masters,
                p,
                degrees: degree_bounds(&blocks, n),
                blocks,
            },
            report,
        ))
    }
}

fn infer_weights(
    samples: &mut Samples<'_>,
    masters: usize,
    seed: u64,
    search: &AdaptiveSearch,
) -> Vec<u32> {
    let mut degrees = vec![1; samples.n];
    for (k, degree) in degrees.iter_mut().enumerate() {
        // The first integral gives cheap lower bounds; uniform fallback protects
        // against later integrals having a different degree profile.
        let mut slice_search = search.clone();
        slice_search.groups.clear();
        slice_search.integral_weights.clear();
        if let Some((_, block)) = search_block(
            samples,
            0,
            1,
            masters,
            seed,
            &slice_search,
            &vec![1; samples.n],
            Some(k),
        ) {
            *degree = block
                .relations
                .iter()
                .flatten()
                .flat_map(|g| &g.terms)
                .map(|(e, _)| e[k])
                .max()
                .unwrap_or(0);
        }
    }
    let max = degrees.iter().copied().max().unwrap_or(1).max(1);
    degrees
        .iter()
        .map(|&d| {
            max.checked_div(d)
                .map_or_else(|| search.limits.max_degree.saturating_add(1), |w| w.max(1))
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn search_block(
    samples: &mut Samples<'_>,
    start: usize,
    end: usize,
    masters: usize,
    seed: u64,
    search: &AdaptiveSearch,
    weights: &[u32],
    slice: Option<usize>,
) -> Option<(Shape, Block)> {
    let uniform = vec![1; samples.n];
    for weights in
        std::iter::once(weights).chain((weights != uniform).then_some(uniform.as_slice()))
    {
        for degree in 0..=search.limits.max_degree {
            let Some(ansatz) =
                weighted_ansatz(samples.n, masters + end, degree, weights, search, slice)
            else {
                break;
            };
            if ansatz.is_empty() {
                continue;
            }
            let Some((e, basis)) = relations(samples, &ansatz, masters, seed, None, slice) else {
                continue;
            };
            let mut shape = Shape {
                start,
                end,
                ansatz,
                pivots: e.pivots(),
                selected: vec![],
            };
            let all = polynomials(&basis, &shape, masters, samples.n);
            for i in 0..8 {
                let Some((x, _)) = samples.at_slice(hash(&[seed, 104]), i, slice) else {
                    continue;
                };
                if let Some(selected) = choose(&all, &x, samples.p, masters + start, masters + end)
                {
                    shape.selected = selected;
                    break;
                }
            }
            if shape.selected.is_empty() {
                continue;
            }
            let block = Block {
                start,
                end,
                relations: shape.selected.iter().map(|&i| all[i].clone()).collect(),
            };
            if validate_slice(samples, &block, masters, hash(&[seed, 102]), slice) {
                return Some((shape, block));
            }
        }
    }
    None
}

// Refit only the monomials that actually occur in the selected relations.
// Discovering a dense degree envelope must not make every later prime pay for it.
fn compact(
    samples: &mut Samples<'_>,
    shape: Shape,
    block: Block,
    masters: usize,
    seed: u64,
) -> (Shape, Block) {
    let ansatz: Vec<_> = shape
        .ansatz
        .iter()
        .filter(|(j, e)| {
            block
                .relations
                .iter()
                .any(|r| r[*j].terms.iter().any(|(v, c)| v == e && *c != 0))
        })
        .cloned()
        .collect();
    if ansatz.len() == shape.ansatz.len() {
        return (shape, block);
    }
    let Some((e, basis)) = relations(samples, &ansatz, masters, seed, None, None) else {
        return (shape, block);
    };
    let mut sparse = Shape {
        ansatz,
        pivots: e.pivots(),
        selected: vec![],
        ..shape
    };
    let all = polynomials(&basis, &sparse, masters, samples.n);
    for i in 0..8 {
        let Some((x, _)) = samples.at_slice(hash(&[seed, 104]), i, None) else {
            continue;
        };
        if let Some(selected) = choose(
            &all,
            &x,
            samples.p,
            masters + shape.start,
            masters + shape.end,
        ) {
            sparse.selected = selected;
            break;
        }
    }
    if sparse.selected.is_empty() {
        return (shape, block);
    }
    let fitted = Block {
        start: shape.start,
        end: shape.end,
        relations: sparse.selected.iter().map(|&i| all[i].clone()).collect(),
    };
    if validate(samples, &fitted, masters, hash(&[seed, 105])) {
        (sparse, fitted)
    } else {
        (shape, block)
    }
}

fn weighted_ansatz(
    n: usize,
    integrals: usize,
    degree: u32,
    weights: &[u32],
    search: &AdaptiveSearch,
    slice: Option<usize>,
) -> Option<Vec<(usize, Exps)>> {
    let mut out = Vec::new();
    for j in 0..integrals {
        let budget =
            i64::from(degree) - i64::from(search.integral_weights.get(j).copied().unwrap_or(0));
        if budget < 0 {
            continue;
        }
        let budget = u32::try_from(budget).ok()?;
        let mut terms = vec![(vec![], 0u32)];
        for (k, &w) in weights.iter().enumerate().take(n) {
            let mut next = Vec::new();
            for (e, used) in terms {
                let max = if slice.is_some_and(|v| v != k) {
                    0
                } else {
                    (budget - used) / w
                };
                for d in 0..=max {
                    if out.len() + next.len() >= search.limits.max_unknowns {
                        return None;
                    }
                    let mut v = e.clone();
                    v.push(d);
                    next.push((v, used + d * w));
                }
            }
            terms = next;
        }
        out.extend(
            terms
                .into_iter()
                .filter(|(e, _)| {
                    search
                        .groups
                        .iter()
                        .all(|(g, b)| g.iter().map(|&k| e[k]).sum::<u32>() <= *b)
                })
                .map(|(e, _)| (j, e)),
        );
    }
    Some(out)
}

fn validate(samples: &mut Samples<'_>, block: &Block, masters: usize, seed: u64) -> bool {
    validate_slice(samples, block, masters, seed, None)
}

/// A block reduction including automatically chosen replay dependencies.
#[derive(Clone, Debug)]
pub struct Reduction {
    oracle_plan: crate::Plan,
    output: Vec<usize>,
    intermediates: usize,
    pub blocks: BlockPlan,
    pub first: BlockForm,
    pub report: SearchReport,
}
impl Reduction {
    /// Replay/export plan for just the requested outputs in the shared master basis.
    pub fn output_plan(&self) -> crate::Plan {
        let mut plan = self.oracle_plan.clone();
        plan.targets = self
            .output
            .iter()
            .map(|&i| self.oracle_plan.targets[i])
            .collect();
        plan
    }
    pub fn masters(&self) -> &[usize] {
        &self.oracle_plan.masters
    }
    pub const fn intermediates(&self) -> usize {
        self.intermediates
    }
    fn select(&self, all: &[u64]) -> Vec<u64> {
        let m = self.masters().len();
        self.output
            .iter()
            .flat_map(|&i| all[i * m..(i + 1) * m].iter().copied())
            .collect()
    }
    /// Evaluate the learned form only at its explicit prime; refit for other primes.
    pub fn eval(&self, x: &[u64], p: u64) -> Option<Vec<u64>> {
        if p != self.first.p {
            return None;
        }
        Some(self.select(&self.first.eval(x)?))
    }
    pub fn fit(&self, system: &crate::ibp::System, p: u64, seed: u64) -> Option<BlockForm> {
        self.blocks.fit(
            |x, p| self.oracle_plan.replay(|e| system.row(e, x, p), p),
            p,
            seed,
        )
    }
    /// Lift using fitted blocks and optional denominator candidates together.
    pub fn lift(
        &self,
        system: &crate::ibp::System,
        candidates: &[zippel_lift::Poly],
        seed: u64,
    ) -> Option<Vec<zippel_lift::Fraction>> {
        use std::sync::{Arc, Mutex};
        let forms = Mutex::new(HashMap::from([(
            self.first.p,
            Some(Arc::new(self.first.clone())),
        )]));
        let f = |x: &[u64], p| {
            let form = {
                let mut cache = forms.lock().unwrap();
                cache
                    .entry(p)
                    .or_insert_with(|| self.fit(system, p, seed).map(Arc::new))
                    .clone()?
            };
            Some(self.select(&form.eval(x)?))
        };
        zippel_lift::lift_with_discovered_factors(f, system.vars.len(), candidates, seed)
    }
}

impl crate::ibp::System {
    /// Select frequently used simpler pivot integrals as candidate intermediates.
    ///
    /// They are reduced in one shared master basis with the requested targets.
    /// The returned master list is explicit because admitting intermediates may
    /// expose additional masters hidden by cancellations in the requested outputs.
    pub fn learn_blocks(
        &self,
        plan: &crate::Plan,
        seed: u64,
        search: &AdaptiveSearch,
    ) -> Option<Reduction> {
        let p = polycore::modp::Primes::new().next()?;
        let mut frequency = HashMap::new();
        let max = plan.targets.iter().copied().max()?;
        let x = (0..self.vars.len() as u64)
            .map(|j| point(&[seed, 110, j], p))
            .collect::<Vec<_>>();
        for &e in &plan.eqs {
            for (j, c) in self.row(e, &x, p) {
                if c != 0 && j < max && plan.pivots.contains(&j) && !plan.targets.contains(&j) {
                    *frequency.entry(j).or_insert(0usize) += 1;
                }
            }
        }
        let mut intermediates: Vec<_> = frequency.into_iter().collect();
        intermediates.sort_by_key(|&(j, f)| (std::cmp::Reverse(f), std::cmp::Reverse(j)));
        intermediates.truncate(search.max_intermediates);
        let mut targets = plan.targets.clone();
        targets.extend(intermediates.iter().map(|&(j, _)| j));
        targets.sort_unstable();
        targets.dedup();
        let indices = targets
            .iter()
            .map(|&j| self.integrals[j].clone())
            .collect::<Vec<_>>();
        let oracle_plan = self.learn(&indices, seed)?;
        let output = plan
            .targets
            .iter()
            .map(|t| targets.binary_search(t).ok())
            .collect::<Option<Vec<_>>>()?;
        let (blocks, first, report) = BlockPlan::learn_adaptive(
            |x, p| oracle_plan.replay(|e| self.row(e, x, p), p),
            self.vars.len(),
            &targets,
            oracle_plan.masters.len(),
            p,
            seed,
            search,
        )?;
        Some(Reduction {
            oracle_plan,
            output,
            intermediates: intermediates.len(),
            blocks,
            first,
            report,
        })
    }
}
