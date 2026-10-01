# Changelog

## Unreleased

- Expose each GCD algorithm through a named function.
- Return `Option` from named sparse GCD methods without automatic fallback.

## 0.3.0 (2026-10-01)

- Add Hu–Monagan GCD with univariate and bivariate images.
- Add Huang–Monagan GCD by prime substitution.
- Use Huang–Gao separated Hensel lifting by default.
- Add GCD algorithm selection with Zippel fallback.
- Add NTT convolution, Newton division and half-GCD.
- Bump `groebner` to 0.5.
- Bump `polycore` to 0.1.4.
- Cache monomial powers for modular evaluation.
- Stream variable evaluations and leading coefficients.
- Stop content GCDs at constants.
- Add exact heap division on packed exponents.
- Try CRT candidates early and reconstruct normalized coefficients as fractions.

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
