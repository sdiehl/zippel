//! Reduce integrals of the one-loop equal-mass bubble to masters: Laporta elimination learned at
//! one point, replayed modulo many primes, and lifted to exact rational functions of `(d, s, m2)`.

#[path = "ibp/mod.rs"]
mod ibp;

use std::time::Instant;
use zippel_interp::Primes;

fn main() {
    let bubble = ibp::Bubble::new(3);
    let start = Instant::now();
    let (plan, coefficients) = bubble.reduce(
        &[(2, 1), (2, 2), (3, 1), (3, 0), (1, -1)],
        Primes::new().next().unwrap(),
    );
    print!("{}", bubble.render(&plan, &coefficients));
    println!(
        "\n{} equations in {} integrals, {} kept after learning, {:.0?}",
        bubble.len(),
        bubble.integrals.len(),
        plan.len(),
        start.elapsed()
    );
}
