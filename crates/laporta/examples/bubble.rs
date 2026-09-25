//! Reduce integrals of the one-loop equal-mass bubble to masters: Laporta elimination learned at
//! one point, replayed modulo many primes, and lifted to exact rational functions of `(d, s, m2)`.

mod families;

use std::time::Instant;

fn main() {
    let targets = [[2, 1], [2, 2], [3, 1], [3, 0], [1, -1]].map(|a| a.to_vec());
    let start = Instant::now();
    let (system, plan, coefficients) = families::bubble().reduce(&targets, 2, 1).unwrap();
    print!("{}", system.render(&plan, &coefficients));
    println!(
        "\n{} equations in {} integrals, {} kept after learning, {:.0?}",
        system.len(),
        system.integrals.len(),
        plan.len(),
        start.elapsed()
    );
}
