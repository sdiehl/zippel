# Changelog

## 0.2.0

Refactor onto `polycore` for shared polynomial and modular arithmetic.

- Polynomials are now `polycore::Poly` in place of `groebner::Polynomial`.
- `BlackBox`, `Primes` and `Rng` move to `polycore`.
- Drop the `modp`, `thiele`, `univariate` and `vandermonde` modules.
- Export `Dense` and `dense` from `zippel-interp`.
- Berlekamp-Massey, root finding, Vandermonde solves and CRT come from `polycore`.
- Bump `groebner` to 0.4.
- MSRV is now Rust 1.98, required by `polycore`.
- Optional `parallel` feature, on by default, gates rayon.

## 0.1.0

Initial release.

- Zippel sparse interpolation over word-sized primes.
- Ben-Or/Tiwari interpolation for very sparse polynomials.
- Rational function reconstruction via Thiele and Zippel.
- Lifting to Q by CRT and rational reconstruction.
- Multivariate GCD, cofactors and LCM via LINZIP.
- IBP generation from propagators, legs and symmetries.
- Laporta elimination learned once, replayed per sample.
- Per-sector seed trimming, widened when masters look wrong.
- Fix a kinematic scale to cut reconstruction variables.
- Parallel sampling with rayon.
- Bubble, box and double box examples.
