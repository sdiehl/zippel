//! Exact parametric tests at fixed kinematics; no numerical symmetry guesses.
use super::{PolynomialFamily, Q};
use num_traits::{One, Zero};
use polycore::{Echelon, Lead, Monomial, Order, Poly};
use std::collections::{BTreeMap, BTreeSet};

type Terms = Vec<(Vec<u32>, Q)>;
type Mapping = Vec<(usize, usize)>;

/// Certified scaleless sectors and index maps from equal Symanzik polynomials.
/// The scaling test is sufficient, not a complete classification of zero sectors.
#[derive(Clone, Debug, Default)]
pub struct FamilyAnalysis {
    width: usize,
    pub zero_sectors: Vec<u32>,
    sector_maps: BTreeMap<u32, Vec<Mapping>>,
    permutations: Vec<Mapping>,
}
impl FamilyAnalysis {
    /// Choose an equivalent index under certified maps, or `None` for a zero sector.
    pub fn canonical_index(&self, a: &[i32]) -> Option<Vec<i32>> {
        if a.len() != self.width {
            return None;
        }
        let sector = a
            .iter()
            .enumerate()
            .filter(|(_, v)| **v > 0)
            .fold(0, |s, (i, _)| s | (1 << i));
        if sector == 0 || self.zero_sectors.iter().any(|&s| sector & !s == 0) {
            return None;
        }
        let mut best = a.to_vec();
        loop {
            let next = self
                .images(&best)
                .into_iter()
                .chain([best.clone()])
                .min_by_key(|a| super::key(a, a.len()))?;
            if next == best {
                return Some(best);
            }
            best = next;
        }
    }
    pub fn symmetry_count(&self) -> usize {
        self.permutations.len() + self.sector_maps.values().map(Vec::len).sum::<usize>()
    }
    pub(super) fn images(&self, a: &[i32]) -> Vec<Vec<i32>> {
        let sector = a
            .iter()
            .enumerate()
            .filter(|(_, v)| **v > 0)
            .fold(0, |s, (i, _)| s | (1 << i));
        let sector_maps = self
            .sector_maps
            .get(&sector)
            .into_iter()
            .flatten()
            .filter(|_| a.iter().all(|&v| v >= 0));
        self.permutations
            .iter()
            .chain(sector_maps)
            .map(|map| {
                let mut b = vec![0; a.len()];
                for &(i, j) in map {
                    b[j] = a[i];
                }
                b
            })
            .filter(|b| b != a)
            .collect()
    }
}

impl PolynomialFamily {
    /// Detect scaleless sectors by a weighted scaling of U+F and exact
    /// permutations of its parameters. Equality includes all kinematic coefficients.
    /// Symmetry enumeration is bounded; omitted maps only affect efficiency.
    pub fn analyze(&self) -> Option<FamilyAnalysis> {
        if !self.valid() {
            return None;
        }
        let n = self.propagators.len();
        let g = self.parametric()?;
        let mut sectors = BTreeSet::new();
        for &top in &self.top_sectors {
            let mut s = top;
            while s != 0 {
                sectors.insert(s);
                s = (s - 1) & top;
            }
        }
        let mut result = FamilyAnalysis {
            width: n,
            zero_sectors: self.zero_sectors.clone(),
            ..FamilyAnalysis::default()
        };
        let mut representatives: BTreeMap<Terms, Vec<usize>> = BTreeMap::new();
        for sector in sectors {
            let restricted = restrict(&g, n, sector);
            if scaleless(&restricted, n) {
                result.zero_sectors.push(sector);
                continue;
            }
            if self.zero_sectors.iter().any(|&s| sector & !s == 0)
                || self.cuts.iter().any(|&i| sector >> i & 1 == 0)
            {
                continue;
            }
            let active = (0..n).filter(|&i| sector >> i & 1 == 1).collect::<Vec<_>>();
            if let Some((key, orders)) = canonical(&restricted, n, &active, &self.cuts) {
                let representative = representatives
                    .entry(key)
                    .or_insert_with(|| orders[0].clone());
                let maps = orders
                    .into_iter()
                    .map(|order| {
                        order
                            .into_iter()
                            .zip(representative.iter().copied())
                            .collect()
                    })
                    .collect();
                result.sector_maps.insert(sector, maps);
            }
        }
        // Full-family maps also transform numerator indices, where sector-only
        // permutations would not justify an identity.
        if let Some((_, orders)) = canonical(&g, n, &(0..n).collect::<Vec<_>>(), &self.cuts) {
            let first = &orders[0];
            for order in &orders {
                let map: Mapping = order.iter().copied().zip(first.iter().copied()).collect();
                let mapped = |s: u32| {
                    map.iter()
                        .filter(|(i, _)| s >> i & 1 != 0)
                        .fold(0, |s, (_, j)| s | (1 << j))
                };
                if self
                    .top_sectors
                    .iter()
                    .all(|&s| self.top_sectors.contains(&mapped(s)))
                    && self
                        .zero_sectors
                        .iter()
                        .all(|&s| self.zero_sectors.contains(&mapped(s)))
                {
                    result.permutations.push(map);
                }
            }
        }
        Some(result)
    }

    fn parametric(&self) -> Option<Poly<Q>> {
        let l = self.loops;
        let e = self.products.len();
        let n = self.propagators.len();
        let nv = self.vars.len();
        // Determinants by subsets stay small for the intended loop counts.
        if l > 6 {
            return None;
        }
        let zero = Poly::constant(Q::zero(), n + nv, Order::Lex);
        let embed = |g: &Poly<Q>| {
            Poly::new(
                g.terms
                    .iter()
                    .map(|(m, c)| {
                        let mut exps = vec![0; n];
                        exps.extend(m.exps());
                        (Monomial::new(exps), c.clone())
                    })
                    .collect(),
                n + nv,
                Order::Lex,
            )
        };
        let mut a = vec![vec![zero.clone(); l]; l];
        let mut b = vec![vec![zero.clone(); e]; l];
        let mut c = zero.clone();
        let half = Q::new(1.into(), 2.into());
        for (k, g) in self.propagators.iter().enumerate() {
            for (m, coefficient) in &g.terms {
                let factors = m.exps()[..l + e]
                    .iter()
                    .enumerate()
                    .flat_map(|(i, &v)| std::iter::repeat_n(i, v as usize))
                    .collect::<Vec<_>>();
                let mut exps = vec![0; n];
                exps[k] = 1;
                exps.extend(&m.exps()[l + e..]);
                let value = Poly::new(
                    vec![(Monomial::new(exps), coefficient.clone())],
                    n + nv,
                    Order::Lex,
                );
                match factors.as_slice() {
                    [] => c = &c + &value,
                    &[i, j] if j < l => {
                        let value = if i == j { value } else { value.scale(&half) };
                        a[i][j] = &a[i][j] + &value;
                        if i != j {
                            a[j][i] = &a[j][i] + &value;
                        }
                    }
                    &[i, j] if i < l => b[i][j - l] = &b[i][j - l] + &value.scale(&half),
                    &[i, j] => c = &c + &(&value * &embed(&self.products[i - l][j - l])),
                    _ => return None,
                }
            }
        }
        let u = determinant(&a, &zero);
        let mut f = &c * &u;
        // det([[A,B_j],[B_k^T,0]]) = -B_k^T adj(A) B_j.
        for j in 0..e {
            for k in 0..e {
                if self.products[j][k].is_zero() {
                    continue;
                }
                let mut augmented = a.clone();
                for i in 0..l {
                    augmented[i].push(b[i][j].clone());
                }
                let mut last = (0..l).map(|i| b[i][k].clone()).collect::<Vec<_>>();
                last.push(zero.clone());
                augmented.push(last);
                f = &f + &(&determinant(&augmented, &zero) * &embed(&self.products[j][k]));
            }
        }
        Some(&u + &f)
    }
}

fn determinant(a: &[Vec<Poly<Q>>], zero: &Poly<Q>) -> Poly<Q> {
    let n = a.len();
    let mut det = vec![zero.clone(); 1 << n];
    det[0] = Poly::constant(Q::one(), zero.nvars, Order::Lex);
    for mask in 1usize..1 << n {
        let row = mask.count_ones() as usize - 1;
        let mut pos = 0;
        for j in 0..n {
            if mask >> j & 1 != 0 {
                let t = &det[mask ^ (1 << j)] * &a[row][j];
                det[mask] = if (row + pos).is_multiple_of(2) {
                    &det[mask] + &t
                } else {
                    &det[mask] - &t
                };
                pos += 1;
            }
        }
    }
    det.pop().unwrap()
}
fn restrict(g: &Poly<Q>, n: usize, sector: u32) -> Poly<Q> {
    Poly::new(
        g.terms
            .iter()
            .filter(|(m, _)| {
                m.exps()[..n]
                    .iter()
                    .enumerate()
                    .all(|(i, &v)| v == 0 || sector >> i & 1 != 0)
            })
            .cloned()
            .collect(),
        g.nvars,
        Order::Lex,
    )
}
fn scaleless(g: &Poly<Q>, n: usize) -> bool {
    let mut e = Echelon::new(Lead::Low);
    for (m, _) in &g.terms {
        let mut row = m.exps()[..n]
            .iter()
            .enumerate()
            .filter(|(_, v)| **v != 0)
            .map(|(i, &v)| (i, Q::from_integer(v.into())))
            .collect::<Vec<_>>();
        row.push((n, Q::one()));
        if e.insert(&row) == Some(n) {
            return false;
        }
    }
    true
}

fn canonical(
    g: &Poly<Q>,
    n: usize,
    active: &[usize],
    cuts: &[usize],
) -> Option<(Terms, Vec<Vec<usize>>)> {
    let mut groups = BTreeMap::new();
    for &i in active {
        let mut signature = g
            .terms
            .iter()
            .map(|(m, c)| {
                let mut powers = m.exps()[..n].to_vec();
                powers.sort_unstable();
                (m.exps()[i], powers, m.exps()[n..].to_vec(), c.clone())
            })
            .collect::<Vec<_>>();
        signature.sort();
        groups
            .entry((cuts.contains(&i), signature))
            .or_insert_with(Vec::new)
            .push(i);
    }
    let mut count = 1usize;
    for group in groups.values() {
        for i in 2..=group.len() {
            count = count.checked_mul(i)?;
            if count > 100_000 {
                return None;
            }
        }
    }
    let mut orders = vec![vec![]];
    for group in groups.into_values() {
        let mut permutations = vec![];
        permute(&group, &mut vec![], &mut permutations);
        orders = orders
            .into_iter()
            .flat_map(|prefix| {
                permutations.iter().map(move |suffix| {
                    let mut v = prefix.clone();
                    v.extend(suffix);
                    v
                })
            })
            .collect();
    }
    let mut best = None;
    let mut matches = vec![];
    for order in orders {
        let mut terms = g
            .terms
            .iter()
            .map(|(m, c)| {
                let mut exps = order.iter().map(|&i| m.exps()[i]).collect::<Vec<_>>();
                exps.extend(&m.exps()[n..]);
                (exps, c.clone())
            })
            .collect::<Terms>();
        terms.sort();
        if best.as_ref().is_none_or(|b| terms < *b) {
            best = Some(terms.clone());
            matches.clear();
        }
        if best.as_ref() == Some(&terms) {
            matches.push(order);
        }
    }
    Some((best?, matches))
}
fn permute(remaining: &[usize], prefix: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
    if remaining.is_empty() {
        out.push(prefix.clone());
        return;
    }
    for (k, &i) in remaining.iter().enumerate() {
        let mut rest = remaining.to_vec();
        rest.remove(k);
        prefix.push(i);
        permute(&rest, prefix, out);
        prefix.pop();
    }
}
