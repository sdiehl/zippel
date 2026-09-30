# Changelog

## Unreleased

- Bump `groebner` to 0.5, with faster F4 and rational reconstruction.
- Bump `polycore` to 0.1.3 and take integer contents by Lehmer's gcd, about 3x faster gcd on large coefficients.
- Evaluate monomials from precomputed power tables in `linzip` and `ModPoly` evaluation.
- Evaluate a variable of `ModPoly` in one pass without dense intermediates, read the leading coefficient directly, and stop the content gcd once it is constant; about 2x faster recursive `pgcd`.
- Divide exactly by heap division on packed exponents, after cheap rejections by degrees, end terms and the value at one.
- Try the trial division as soon as the CRT residues settle inside the modulus instead of waiting for two equal candidates, and also reconstruct `h / lc(h)` as fractions, which needs far fewer primes when `gcd(lc f, lc g)` is much larger than `lc(h)`.

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
