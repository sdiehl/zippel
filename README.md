# tiny-zippel

Sparse multivariate polynomial GCD by Zippel's modular interpolation. Works on [groebner](https://crates.io/crates/groebner) polynomials over the rationals and exposes `gcd`, `cofactors` and `lcm`.

- Integer driver: integer, monomial and main-variable content removal, images modulo 62-bit primes, Garner CRT, and a final trial division that certifies the result.
- Zippel's recursion over GF(p): dense Newton interpolation in the last variable, with the first image fixing the monomial skeleton for every later one.
- LINZIP sparse interpolation (de Kleine, Monagan, Wittkopf 2005): evaluation at powers turns each skeleton coefficient into a transposed Vandermonde solve, with the unknown image scalings recovered alongside.
- Unlucky primes and evaluation points are detected by comparing leading monomials and skeleton support.
- A coprime fast path from a single univariate image.

```bash
cargo build
cargo test
```

```bash
cargo run --example demo
```

## License

Released under the MIT License. See LICENSE for details.
