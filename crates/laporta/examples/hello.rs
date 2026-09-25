//! The massless one-loop bubble, `D1 = k^2`, `D2 = (k + q)^2`, `q^2 = s`, with dots reduced to
//! its one master integral.

use zippel_laporta::ibp::Family;

fn main() {
    let bubble = Family {
        vars: vec!["d", "s"],
        loops: 1,
        props: vec![(vec![1, 0], vec![]), (vec![1, 1], vec![])], // k^2, (k + q)^2
        lines: 2,
        legs: vec![vec![vec![0, 0, 2]]], // 2 q.q = 2s
        symmetries: vec![vec![1, 0]],
    };
    let system = bubble.system(2, 0);
    let (plan, coefficients) = system.reduce(&[vec![2, 1], vec![2, 2]], 1).unwrap();
    print!("{}", system.render(&plan, &coefficients));
}
