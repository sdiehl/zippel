# Changelog

## 0.2.0 (2026-09-27)

- Move polynomials onto `polycore::Poly`.
- Move `BlackBox`, `Primes` and `Rng` to `polycore`.
- Take Berlekamp-Massey, roots, Vandermonde and CRT from `polycore`.
- Drop the `modp`, `thiele`, `univariate` and `vandermonde` modules.
- Export `Dense` and `dense` from `zippel-interp`.
- Add `Frac`, rational functions in lowest terms.
- Add an optional `parallel` feature gating rayon.
- Bump `groebner` to 0.4.
- Set MSRV to Rust 1.88.

## 0.1.0 (2026-09-25)

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
