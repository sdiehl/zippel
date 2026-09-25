//! Integral families: propagator momenta over the loop then external momenta, and kinematics as
//! linear forms `[1, d, ..]` in the variables.

#![allow(dead_code, unreachable_pub)]

use zippel_laporta::ibp::Family;

/// One loop, two equal masses: `D1 = k^2 - m2`, `D2 = (k + q)^2 - m2`, `q^2 = s`.
pub fn bubble() -> Family {
    let m2 = vec![0, 0, 0, 1];
    Family {
        vars: vec!["d", "s", "m2"],
        loops: 1,
        props: vec![(vec![1, 0], m2.clone()), (vec![1, 1], m2)],
        lines: 2,
        legs: vec![vec![vec![0, 0, 2, 0]]],
        symmetries: vec![vec![1, 0]],
    }
}

/// Massless on-shell legs, `s = (p1 + p2)^2`, `t = (p2 + p3)^2`, `u = -s - t`.
fn massless_legs() -> Vec<Vec<Vec<i64>>> {
    let (z, s, t, u) = (
        vec![],
        vec![0, 0, 1, 0],
        vec![0, 0, 0, 1],
        vec![0, 0, -1, -1],
    );
    vec![
        vec![z.clone(), s.clone(), u.clone()],
        vec![s, z.clone(), t.clone()],
        vec![u, t, z],
    ]
}

/// The massless one-loop box: `k`, `k + p1`, `k + p1 + p2`, `k - p4`.
pub fn one_loop_box() -> Family {
    let props = [[1, 0, 0, 0], [1, 1, 0, 0], [1, 1, 1, 0], [1, 1, 1, 1]];
    Family {
        vars: vec!["d", "s", "t"],
        loops: 1,
        props: props.iter().map(|q| (q.to_vec(), vec![])).collect(),
        lines: 4,
        legs: massless_legs(),
        symmetries: vec![vec![2, 1, 0, 3], vec![0, 3, 2, 1]],
    }
}

/// The massless planar double box: lines `k1`, `k1 + p1`, `k1 + p1 + p2`, `k2 + p1 + p2`,
/// `k2 - p4`, `k2`, `k1 - k2`, and numerators `(k1 - p4)^2`, `(k2 + p1)^2`.
pub fn double_box() -> Family {
    let props = [
        [1, 0, 0, 0, 0],
        [1, 0, 1, 0, 0],
        [1, 0, 1, 1, 0],
        [0, 1, 1, 1, 0],
        [0, 1, 1, 1, 1],
        [0, 1, 0, 0, 0],
        [1, -1, 0, 0, 0],
        [1, 0, 1, 1, 1],
        [0, 1, 1, 0, 0],
    ];
    Family {
        vars: vec!["d", "s", "t"],
        loops: 2,
        props: props.iter().map(|q| (q.to_vec(), vec![])).collect(),
        lines: 7,
        legs: massless_legs(),
        symmetries: vec![
            vec![2, 1, 0, 5, 4, 3, 6, 7, 8],
            vec![5, 4, 3, 2, 1, 0, 6, 8, 7],
        ],
    }
}
