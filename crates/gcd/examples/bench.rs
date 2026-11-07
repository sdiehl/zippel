//! Deterministic GCD comparison: `cargo run --release -p zippel-gcd --example bench -- 5`.
//! Optional second argument filters case names. Output is CSV; timings exclude parsing
//! and validation. Direct timing includes declined recoveries; chain timing includes
//! actual fallback costs. Each case/algorithm is warmed up once before measurement.

use std::{collections::BTreeMap, hint::black_box, time::Instant};

use num_rational::BigRational;
use polycore::{Order, Ring};
use zippel_gcd::{
    gcd_with_report, hu_monagan, hu_monagan_bivariate, huang_gao, huang_monagan, zippel,
    GcdAlgorithm,
};

type Poly = polycore::Poly<BigRational>;

fn direct(f: &Poly, g: &Poly, algorithm: GcdAlgorithm) -> Option<Poly> {
    match algorithm {
        GcdAlgorithm::HuangGao => huang_gao(f, g),
        GcdAlgorithm::HuangMonagan => huang_monagan(f, g),
        GcdAlgorithm::HuMonagan => hu_monagan(f, g),
        GcdAlgorithm::HuMonaganBivariate => hu_monagan_bivariate(f, g),
        GcdAlgorithm::Zippel => Some(zippel(f, g)),
    }
}

fn measure(
    name: &str,
    f: &Poly,
    g: &Poly,
    expected: &Poly,
    algorithm: GcdAlgorithm,
    repeats: usize,
) {
    assert_eq!(gcd_with_report(f, g, algorithm).0, *expected);
    black_box(direct(f, g, algorithm));
    let mut direct_times = Vec::new();
    let mut chain_times = Vec::new();
    let mut declines = 0;
    let mut fallbacks = 0;
    let mut recovered = BTreeMap::<String, usize>::new();
    for _ in 0..repeats {
        let start = Instant::now();
        let value = black_box(direct(black_box(f), black_box(g), algorithm));
        direct_times.push(start.elapsed().as_micros());
        if let Some(value) = value {
            assert_eq!(value, *expected);
        } else {
            declines += 1;
        }
        let start = Instant::now();
        let (value, report) = black_box(gcd_with_report(black_box(f), black_box(g), algorithm));
        chain_times.push(start.elapsed().as_micros());
        assert_eq!(value, *expected);
        fallbacks += usize::from(report.attempted.len() > 1);
        let backend = format!("{:?}", report.attempted.last().unwrap());
        *recovered.entry(backend).or_default() += 1;
    }
    direct_times.sort_unstable();
    chain_times.sort_unstable();
    let recovered = recovered
        .into_iter()
        .map(|(a, n)| format!("{a}:{n}"))
        .collect::<Vec<_>>()
        .join("|");
    println!(
        "{name},{algorithm:?},{repeats},{},{},{},{},{declines},{fallbacks},{recovered}",
        f.terms.len(),
        g.terms.len(),
        direct_times[repeats / 2],
        chain_times[repeats / 2]
    );
}

fn main() {
    let mut args = std::env::args().skip(1);
    let repeats = args
        .next()
        .map_or(5, |s| s.parse::<usize>().expect("positive repeat count"));
    assert!(repeats > 0, "repeat count must be positive");
    let filter = args.next().unwrap_or_default();
    let ring = Ring::new(["x", "y", "z"], Order::Lex);
    let cases = [
        ("sparse", "x^3*y+2*y*z^2+3*z+1"),
        ("dense", "(x+y+z+1)^5"),
        (
            "large-coefficients",
            "1606938044258990275541962092341162602522202993782792835301377*x+3*y+7*z+1",
        ),
        ("high-degree", "x^80+y^3+z+1"),
    ];
    let algorithms = [
        GcdAlgorithm::HuangGao,
        GcdAlgorithm::HuangMonagan,
        GcdAlgorithm::HuMonagan,
        GcdAlgorithm::HuMonaganBivariate,
        GcdAlgorithm::Zippel,
    ];
    println!("case,algorithm,runs,f_terms,g_terms,direct_median_us,chain_median_us,direct_declines,fallback_runs,recovered_by");
    for (name, source) in cases {
        if !name.contains(&filter) {
            continue;
        }
        let h = ring.parse(source).unwrap().primitive();
        let f = &h * &ring.parse("x+2").unwrap();
        let g = &h * &ring.parse("y+3").unwrap();
        for algorithm in algorithms {
            measure(name, &f, &g, &h, algorithm, repeats);
        }
    }
    // A sparse univariate case that exceeds every bounded backend's image degree.
    if "fallback".contains(&filter) {
        let ring = Ring::new(["x"], Order::Lex);
        let f = ring.parse("x^65537+1").unwrap();
        for algorithm in algorithms {
            measure("fallback", &f, &f, &f, algorithm, repeats);
        }
    }
}
