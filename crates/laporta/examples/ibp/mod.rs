//! IBP identities for the one-loop bubble with two equal masses,
//! `I(a, b) = int d^dk / (D1^a D2^b)` with `D1 = k^2 - m2`, `D2 = (k + q)^2 - m2`, `q^2 = s`.
//! Kinematics are the variables `(d, s, m2)`.

#![allow(unreachable_pub, dead_code)]
use std::fmt::Write;

use groebner::{MonomialOrder, PolynomialRing};
use num_rational::BigRational;
use zippel_interp::modp::{add, mul};
use zippel_laporta::{Plan, Row};
use zippel_lift::{lift, Fraction};

pub type Index = (i32, i32);

/// `c0 + c1*d + c2*s + c3*m2`, which is all IBP coefficients need here.
type Coef = [i64; 4];

/// Integrating `d/dk . k` and `d/dk . q` by parts at the seed `(a, b)`.
fn identities(a: i32, b: i32) -> [Vec<(Index, Coef)>; 2] {
    let (x, y) = (i64::from(a), i64::from(b));
    [
        vec![
            ((a, b), [-2 * x - y, 1, 0, 0]),
            ((a + 1, b), [0, 0, 0, -2 * x]),
            ((a - 1, b + 1), [-y, 0, 0, 0]),
            ((a, b + 1), [0, 0, y, -2 * y]),
        ],
        vec![
            ((a, b), [x - y, 0, 0, 0]),
            ((a + 1, b), [0, 0, x, 0]),
            ((a + 1, b - 1), [-x, 0, 0, 0]),
            ((a - 1, b + 1), [y, 0, 0, 0]),
            ((a, b + 1), [0, 0, -y, 0]),
        ],
    ]
}

/// Laporta's order: propagators, then dots, then numerators, then which line.
fn weight((a, b): Index) -> (i32, i32, i32, i32) {
    let (pos, neg) = [a, b].iter().fold(
        (0, 0),
        |(u, v), &i| {
            if i > 0 {
                (u + i, v)
            } else {
                (u, v - i)
            }
        },
    );
    (i32::from(a > 0) + i32::from(b > 0), pos, neg, b)
}

#[derive(Debug)]
pub struct Bubble {
    pub integrals: Vec<Index>,
    eqs: Vec<Vec<(usize, Coef)>>,
}

impl Bubble {
    /// All identities seeded with `1 <= a <= r`, `-1 <= b <= r`, plus the `a <-> b` symmetry.
    pub fn new(r: i32) -> Self {
        let mut raw: Vec<Vec<(Index, Coef)>> = Vec::new();
        for a in 1..=r {
            for b in -1..=r {
                raw.extend(identities(a, b));
            }
        }
        for eq in &mut raw {
            eq.retain(|((a, b), c)| (*a > 0 || *b > 0) && c.iter().any(|&v| v != 0));
        }
        let mut integrals: Vec<Index> = raw.iter().flatten().map(|t| t.0).collect();
        integrals.extend(integrals.clone().iter().map(|&(a, b)| (b, a)));
        integrals.sort_by_key(|i| (weight(*i), *i));
        integrals.dedup();
        for &(a, b) in &integrals {
            if a < b {
                raw.push(vec![((a, b), [1, 0, 0, 0]), ((b, a), [-1, 0, 0, 0])]);
            }
        }
        let column = |i: &Index| {
            integrals
                .binary_search_by_key(&(weight(*i), *i), |j| (weight(*j), *j))
                .unwrap()
        };
        let eqs = raw
            .iter()
            .map(|eq| eq.iter().map(|(i, c)| (column(i), *c)).collect())
            .collect();
        Self { integrals, eqs }
    }

    pub const fn len(&self) -> usize {
        self.eqs.len()
    }

    pub fn column(&self, i: Index) -> usize {
        self.integrals.iter().position(|&j| j == i).unwrap()
    }

    /// Equation `e` at the numeric point `x = (d, s, m2)`.
    pub fn row(&self, e: usize, x: &[u64], p: u64) -> Row {
        let residue = |v: i64| u64::try_from(i128::from(v).rem_euclid(i128::from(p))).unwrap();
        self.eqs[e]
            .iter()
            .map(|(j, c)| {
                let v = c[1..].iter().zip(x).fold(residue(c[0]), |acc, (&ci, &xi)| {
                    add(acc, mul(residue(ci), xi, p), p)
                });
                (*j, v)
            })
            .collect()
    }

    /// Learns at one point mod `p`, then lifts every coefficient to Q from replays.
    pub fn reduce(&self, targets: &[Index], p: u64) -> (Plan, Vec<Fraction>) {
        let targets: Vec<usize> = targets.iter().map(|&t| self.column(t)).collect();
        let x0 = [7, 11, 13].map(|v| v * 0x9e37_79b9 % p);
        let rows: Vec<Row> = (0..self.len()).map(|e| self.row(e, &x0, p)).collect();
        let plan = Plan::learn(&rows, &targets, p);
        let coefficients =
            lift(|x: &[u64], p| plan.replay(|e| self.row(e, x, p), p), 3, 1).unwrap();
        (plan, coefficients)
    }

    /// `I(a,b) = c * I(..) + ...` per target.
    pub fn render(&self, plan: &Plan, coefficients: &[Fraction]) -> String {
        let ring =
            PolynomialRing::<BigRational>::new(["d", "s", "m2"], MonomialOrder::Lex).unwrap();
        let show = |f: &Fraction| {
            let num = ring.format(&f.num).unwrap();
            match ring.format(&f.den).unwrap().as_str() {
                "1" => format!("({num})"),
                den => format!("({num}) / ({den})"),
            }
        };
        let name = |j: usize| format!("I({},{})", self.integrals[j].0, self.integrals[j].1);
        let mut out = String::new();
        for (t, row) in plan
            .targets
            .iter()
            .zip(coefficients.chunks(plan.masters.len()))
        {
            let terms: Vec<String> = row
                .iter()
                .zip(&plan.masters)
                .filter(|(f, _)| !f.num.is_zero())
                .map(|(f, &m)| format!("{} * {}", show(f), name(m)))
                .collect();
            writeln!(out, "{} = {}", name(*t), terms.join("\n    + ")).unwrap();
        }
        out
    }
}
