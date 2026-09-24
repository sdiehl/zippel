use groebner::{MonomialOrder, Polynomial, PolynomialRing};
use num_rational::BigRational;
use std::fmt::Write;

type Ring = PolynomialRing<BigRational>;

/// A product of `(p)`, `(p)^k` or bare factors joined by ` * ` at parenthesis depth zero.
fn product(ring: &Ring, src: &str) -> Polynomial<BigRational> {
    let (mut depth, mut start, mut factors) = (0, 0, Vec::new());
    for (i, ch) in src.char_indices() {
        match ch {
            '(' => depth += 1,
            ')' => depth -= 1,
            '*' if depth == 0 && src[i..].starts_with("* ") => {
                factors.push(src[start..i].trim());
                start = i + 1;
            }
            _ => {}
        }
    }
    factors.push(src[start..].trim());
    factors.iter().fold(ring.parse("1").unwrap(), |acc, f| {
        let (body, k) = match f.rsplit_once(")^") {
            Some((b, k)) => (&b[1..], k.parse().unwrap()),
            None => (f.trim_start_matches('(').trim_end_matches(')'), 1),
        };
        let p = ring.parse(body).unwrap();
        (0..k).fold(acc, |acc, _| acc.multiply(&p))
    })
}

#[test]
fn cases() {
    insta::glob!("cases/*.txt", |path| {
        let src = std::fs::read_to_string(path).unwrap();
        let lines: Vec<&str> = src.lines().collect();
        let ring = Ring::new(lines[0].split(", "), MonomialOrder::GRevLex).unwrap();
        let (f, g) = (product(&ring, lines[1]), product(&ring, lines[2]));
        let (h, cf, cg) = zippel_gcd::cofactors(&f, &g);
        assert_eq!(h.multiply(&cf), f);
        assert_eq!(h.multiply(&cg), g);
        let show = |p| ring.format(p).unwrap();
        let mut out = String::new();
        writeln!(out, "gcd = {}", show(&h)).unwrap();
        writeln!(out, "f/gcd = {}", show(&cf)).unwrap();
        writeln!(out, "g/gcd = {}", show(&cg)).unwrap();
        writeln!(out, "lcm = {}", show(&zippel_gcd::lcm(&f, &g))).unwrap();
        insta::assert_snapshot!(out);
    });
}
