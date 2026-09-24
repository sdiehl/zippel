# tiny-zippel

Sparse polynomial interpolation over finite fields, building toward Feynman integral reduction. Everything works modulo word-sized primes, and polynomials are recovered from evaluations alone.

- `zippel-interp`: prime field arithmetic, sparse polynomials, and Zippel's black-box interpolation. Each variable is lifted in turn against the known support by transposed Vandermonde solves. Newton interpolation with early termination means no degree bounds are needed, and evaluation points the black box rejects are skipped.
- `zippel-interp::rational`: black-box rational functions, after FireFly (Klappert, Lange 2019). Along rays `x = t*z + s`, Thiele's continued fraction recovers the univariate function of `t`, normalized so the denominator is 1 at `t = 0`. The `t^r` coefficients are the homogeneous parts of the numerator and denominator, recovered top down by Zippel with the shift's spill from higher parts subtracted. Points are named by key rather than drawn from a stream, so all parts share their rays.
- `zippel-gcd`: multivariate GCD, cofactors and LCM of [groebner](https://crates.io/crates/groebner) polynomials over Q. It works by content removal and Zippel's recursion modulo primes with LINZIP sparse interpolation (de Kleine, Monagan, Wittkopf 2005), combining primes by Garner CRT and certifying the result by trial division.

```bash
cargo build
cargo test
```

```bash
cargo run -p zippel-interp --example determinant
cargo run -p zippel-interp --example linsolve
cargo run -p zippel-gcd --example demo
```

## License

MIT
