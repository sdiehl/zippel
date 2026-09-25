//! Reduce the massless planar double box, with a numerator and with a dot, to its masters.

mod families;

use std::time::Instant;

fn main() {
    let start = Instant::now();
    let system = families::double_box().system(1, 1);
    let targets = [
        [1, 1, 1, 1, 1, 1, 1, -1, 0],
        [1, 1, 1, 1, 1, 1, 1, 0, -1],
        [2, 1, 1, 1, 1, 1, 1, 0, 0],
    ]
    .map(|a| a.to_vec());
    println!(
        "{} equations in {} integrals, {:.0?}",
        system.len(),
        system.integrals.len(),
        start.elapsed()
    );
    let (plan, coefficients) = system.reduce(&targets, 1).unwrap();
    print!("{}", system.render(&plan, &coefficients));
    println!(
        "\n{} masters, {} kept after learning, {:.0?}",
        plan.masters.len(),
        plan.len(),
        start.elapsed()
    );
}
