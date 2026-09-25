//! Reduce the massless planar double box, with a numerator and with a dot, to its masters.

mod families;

use std::time::Instant;

fn main() {
    let start = Instant::now();
    let targets = [
        [1, 1, 1, 1, 1, 1, 1, -1, 0],
        [1, 1, 1, 1, 1, 1, 1, 0, -1],
        [2, 1, 1, 1, 1, 1, 1, 0, 0],
    ]
    .map(|a| a.to_vec());
    let (system, plan, coefficients) = families::double_box().reduce(&targets, 1, 1).unwrap();
    print!("{}", system.render(&plan, &coefficients));
    println!(
        "\n{} equations in {} integrals, {} masters, {} kept after learning, {:.0?}",
        system.len(),
        system.integrals.len(),
        plan.masters.len(),
        plan.len(),
        start.elapsed()
    );
}
