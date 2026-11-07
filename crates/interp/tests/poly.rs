use polycore::{Fp, Order, Ring};
use zippel_interp::ModPoly;

#[test]
fn substitution_preserves_other_variables_order_and_cancellation() {
    let p = 101;
    let ring = Ring::new(["x", "y", "z"], Order::Lex);
    for source in ["x+y", "x*y^2+x^2*z+y^2-z", "x*y-y+x*z-z", "0"] {
        let f = ring
            .parse(source)
            .unwrap()
            .map(|c| Fp::from_rational(c, p).unwrap());
        let modular = ModPoly::from_poly(&f, p);
        for k in 0..3 {
            for a in [0, 1, 2, 100] {
                let expected = ModPoly::from_poly(&f.eval_var(k, &Fp::new(a, p)), p);
                assert_eq!(
                    modular.eval_var(k, a, p),
                    expected,
                    "{source}, k={k}, a={a}"
                );
            }
        }
    }
}
