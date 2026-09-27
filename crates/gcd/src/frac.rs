//! The field `Q(x_1, ..., x_n)` of rational functions, in lowest terms by [`crate::cofactors`].

use std::fmt;
use std::ops::{Add, Div, Mul, Neg, Sub};

use num_rational::BigRational;
use num_traits::{One, Zero};
use polycore::{Monomial, Order};

use crate::{cofactors, Poly};

/// `num / den` with no common factor and `den` monic in its order. The constants `Zero` and
/// `One` build have no variables and widen to whatever they meet, like `Fp`'s unbound modulus.
#[derive(Clone, Debug)]
pub struct Frac {
    num: Poly,
    den: Poly,
}

impl Frac {
    /// `num / den` in lowest terms, or `None` if `den` is zero.
    #[must_use]
    pub fn new(num: Poly, den: Poly) -> Option<Self> {
        if den.is_zero() {
            return None;
        }
        let (num, den) = align(num, den);
        if num.is_zero() {
            return Some(Self::from(Poly::zero(den.nvars, den.order)));
        }
        let (num, den) = if num.is_constant() || den.is_constant() {
            (num, den)
        } else {
            let (_, n, d) = cofactors(&num, &den);
            (n, d)
        };
        let l = den.lc().map(BigRational::recip)?;
        Some(Self {
            num: num.scale(&l),
            den: den.scale(&l),
        })
    }

    #[must_use]
    pub const fn num(&self) -> &Poly {
        &self.num
    }

    #[must_use]
    pub const fn den(&self) -> &Poly {
        &self.den
    }

    /// The value at `x`, or `None` at a pole.
    #[must_use]
    pub fn eval(&self, x: &[BigRational]) -> Option<BigRational> {
        let d = self.den.eval(x);
        (!d.is_zero()).then(|| self.num.eval(x) / d)
    }

    fn with(num: Poly, den: Poly) -> Self {
        Self::new(num, den).expect("nonzero denominator")
    }
}

/// `p` over `n` variables, the new ones last.
fn widen(p: Poly, n: usize, order: &Order) -> Poly {
    if p.nvars == n && p.order == *order {
        return p;
    }
    let terms = p.terms.into_iter().map(|(m, c)| {
        let mut e = m.exps().to_vec();
        e.resize(n, 0);
        (Monomial::new(e), c)
    });
    Poly::new(terms.collect(), n, order.clone())
}

fn align(a: Poly, b: Poly) -> (Poly, Poly) {
    let (n, order) = if a.nvars >= b.nvars {
        (a.nvars, a.order.clone())
    } else {
        (b.nvars, b.order.clone())
    };
    (widen(a, n, &order), widen(b, n, &order))
}

impl From<Poly> for Frac {
    fn from(num: Poly) -> Self {
        let den = Poly::constant(BigRational::one(), num.nvars, num.order.clone());
        Self { num, den }
    }
}

impl From<BigRational> for Frac {
    fn from(c: BigRational) -> Self {
        Poly::constant(c, 0, Order::GRevLex).into()
    }
}

/// Both over the same variables and order.
fn lift(a: Frac, b: Frac) -> (Frac, Frac) {
    let n = a.num.nvars.max(b.num.nvars);
    let order = if a.num.nvars >= b.num.nvars {
        a.num.order.clone()
    } else {
        b.num.order.clone()
    };
    let w = |f: Frac| Frac {
        num: widen(f.num, n, &order),
        den: widen(f.den, n, &order),
    };
    (w(a), w(b))
}

impl PartialEq for Frac {
    fn eq(&self, o: &Self) -> bool {
        let (a, b) = lift(self.clone(), o.clone());
        a.num == b.num && a.den == b.den
    }
}

impl Eq for Frac {}

impl Zero for Frac {
    fn zero() -> Self {
        BigRational::zero().into()
    }

    fn is_zero(&self) -> bool {
        self.num.is_zero()
    }
}

impl One for Frac {
    fn one() -> Self {
        BigRational::one().into()
    }
}

impl Neg for Frac {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            num: -&self.num,
            den: self.den,
        }
    }
}

impl Add for Frac {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        let (a, b) = lift(self, o);
        if a.den == b.den {
            return Self::with(&a.num + &b.num, a.den);
        }
        Self::with(&(&a.num * &b.den) + &(&b.num * &a.den), &a.den * &b.den)
    }
}

impl Sub for Frac {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        self + -o
    }
}

impl Mul for Frac {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        let (a, b) = lift(self, o);
        Self::with(&a.num * &b.num, &a.den * &b.den)
    }
}

impl Div for Frac {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        assert!(!o.is_zero(), "division by zero");
        let (a, b) = lift(self, o);
        Self::with(&a.num * &b.den, &a.den * &b.num)
    }
}

impl fmt::Display for Frac {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.den.is_constant() {
            write!(f, "{}", self.num)
        } else {
            write!(f, "({})/({})", self.num, self.den)
        }
    }
}
