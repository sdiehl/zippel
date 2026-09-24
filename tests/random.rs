#![allow(clippy::many_single_char_names, clippy::cast_possible_truncation)]

use groebner::{Ideal, Monomial, MonomialOrder, Polynomial, Term};
use num_bigint::BigInt;
use num_rational::BigRational;

type Poly = Polynomial<BigRational>;

struct Lcg(u64);

impl Lcg {
    const fn below(&mut self, n: u64) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 33) % n
    }

    fn poly(&mut self, n: usize, max_terms: u64, deg: u64) -> Poly {
        let terms = (0..=self.below(max_terms))
            .map(|_| {
                let c = BigRational::from(BigInt::from(self.below(19).cast_signed() - 9));
                Term::new(
                    c,
                    Monomial::new((0..n).map(|_| self.below(deg + 1) as u32).collect()),
                )
            })
            .collect();
        Polynomial::new(terms, n, MonomialOrder::GRevLex)
    }
}

fn monic(p: &Poly) -> Poly {
    p.reorder(MonomialOrder::Lex).make_monic()
}

#[test]
fn common_factor_is_found() {
    let mut rng = Lcg(7);
    for _ in 0..200 {
        let n = 1 + rng.below(6) as usize;
        let c = rng.poly(n, 8, 3);
        let a = rng.poly(n, 6, 3);
        let b = rng.poly(n, 6, 3);
        let (f, g) = (a.multiply(&c), b.multiply(&c));
        let (h, cf, cg) = tiny_zippel::cofactors(&f, &g);
        assert_eq!(h.multiply(&cf), f);
        assert_eq!(h.multiply(&cg), g);
        assert_eq!(monic(&tiny_zippel::gcd(&h, &c)), monic(&c));
        if !f.is_zero() && !g.is_zero() {
            assert!(tiny_zippel::gcd(&cf, &cg).is_constant());
        }
    }
}

/// `<t f, (1 - t) g>` meets `Q[x]` in `<lcm(f, g)>`.
fn lcm_by_elimination(f: &Poly, g: &Poly) -> Poly {
    let n = f.nvars;
    let order = MonomialOrder::elimination(1, n);
    let lift = |p: &Poly, t: u32, c: i64| {
        let terms = p.terms.iter().map(|term| {
            let e = std::iter::once(t)
                .chain(term.monomial.exponents().iter().copied())
                .collect();
            Term::new(
                &term.coefficient * BigRational::from(BigInt::from(c)),
                Monomial::new(e),
            )
        });
        Polynomial::new(terms.collect(), n + 1, order.clone())
    };
    let tf = lift(f, 1, 1);
    let one_minus_t_g = lift(g, 0, 1).add(&lift(g, 1, -1));
    let ideal = Ideal::new(vec![tf, one_minus_t_g])
        .unwrap()
        .eliminate(1)
        .unwrap();
    let l = &ideal.basis()[0];
    let terms = l.terms.iter().map(|t| {
        Term::new(
            t.coefficient.clone(),
            Monomial::new(t.monomial.exponents()[1..].to_vec()),
        )
    });
    Polynomial::new(terms.collect(), n, MonomialOrder::Lex)
}

#[test]
fn lcm_matches_ideal_intersection() {
    let mut rng = Lcg(11);
    for _ in 0..30 {
        let n = 2 + rng.below(2) as usize;
        let c = rng.poly(n, 3, 2);
        let (f, g) = (
            rng.poly(n, 3, 2).multiply(&c),
            rng.poly(n, 3, 2).multiply(&c),
        );
        if f.is_zero() || g.is_zero() {
            continue;
        }
        assert_eq!(
            monic(&tiny_zippel::lcm(&f, &g)),
            monic(&lcm_by_elimination(&f, &g))
        );
    }
}
