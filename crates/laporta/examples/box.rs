//! Reduce the massless one-loop box with dots and numerators to its masters.

mod families;

use std::time::Instant;

fn main() {
    let targets = [
        [2, 1, 1, 1],
        [1, 1, 1, 1],
        [2, 2, 1, 1],
        [1, 1, 2, 1],
        [2, 1, 1, 0],
    ]
    .map(|a| a.to_vec());
    let start = Instant::now();
    let (system, plan, coefficients) = families::one_loop_box().reduce(&targets, 2, 2).unwrap();
    print!("{}", system.render(&plan, &coefficients));
    println!(
        "\n{} equations in {} integrals, {} kept after learning, {:.0?}",
        system.len(),
        system.integrals.len(),
        plan.len(),
        start.elapsed()
    );
}
