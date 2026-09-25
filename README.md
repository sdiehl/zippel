# tiny-zippel

Sparse polynomial interpolation over finite fields, building toward Feynman integral reduction. Everything works modulo word-sized primes, and polynomials are recovered from evaluations alone.

- [`zippel-interp`](crates/interp): Zippel's black-box sparse polynomial interpolation over prime fields.
- [`zippel-interp::rational`](crates/interp/src/rational.rs): FireFly-style rational function reconstruction via Thiele and Zippel.
- [`zippel-lift`](crates/lift): Lifts to Q by CRT and rational number reconstruction.
- [`zippel-laporta`](crates/laporta): Laporta elimination over GF(p), learned once and replayed.
- [`zippel-gcd`](crates/gcd): Multivariate GCD, cofactors and LCM of [groebner](https://crates.io/crates/groebner) polynomials via LINZIP.

```bash
cargo build
cargo test
```

```bash
cargo run -p zippel-interp --example determinant
cargo run -p zippel-interp --example linsolve
cargo run -p zippel-laporta --example bubble
cargo run -p zippel-gcd --example demo
```

## References

- R. Zippel, _Probabilistic algorithms for sparse polynomials_, EUROSAM 1979.
- R. Zippel, _Interpolating polynomials from their values_, J. Symbolic Comput. 9 (1990).
- J. de Kleine, M. Monagan, A. Wittkopf, _Algorithms for the non-monic case of the sparse modular GCD algorithm_, ISSAC 2005.
- P. S. Wang, M. J. T. Guy, J. H. Davenport, _P-adic reconstruction of rational numbers_, SIGSAM Bull. 16 (1982).
- S. Laporta, _High-precision calculation of multi-loop Feynman integrals by difference equations_, Int. J. Mod. Phys. A 15 (2000), [arXiv:hep-ph/0102033](https://arxiv.org/abs/hep-ph/0102033).
- J. Klappert, F. Lange, _Reconstructing rational functions with FireFly_, Comput. Phys. Commun. 247 (2020), [arXiv:1904.00009](https://arxiv.org/abs/1904.00009).

## License

MIT
