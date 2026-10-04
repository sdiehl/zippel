#[path = "../examples/families/mod.rs"]
mod families;
use zippel_laporta::formats::{read_family_mathematica, read_family_yaml};
const FAMILY: &str = include_str!("fixtures/box.yaml");
const KIN: &str = include_str!("fixtures/kinematics.yaml");
const MATHEMATICA: &str = include_str!("fixtures/box.m");

#[test]
fn imports_match_native_family_and_export_signs() {
    let yaml = read_family_yaml(FAMILY, KIN, "box").unwrap();
    let mathematica = read_family_mathematica(
        &MATHEMATICA.replace("p1*p2", "p1 p2"),
        "box",
        &["s", "t"],
        4,
    )
    .unwrap();
    assert_eq!(yaml.family.products, mathematica.family.products);
    let target = vec![2, 1, 1, 1];
    let (system, plan, cs) = yaml.family.reduce(&[target], 1, 1).unwrap();
    let rules = yaml.rules(&system, &plan, &cs, None).unwrap();
    assert!(rules.contains("box[2,1,1,1] ->"));
    let (fs, fp, fc) = mathematica
        .family
        .reduce(&[vec![2, 1, 1, 1]], 1, 1)
        .unwrap();
    let tables = mathematica.reduction_tables(&fs, &fp, &fc, 1).unwrap();
    assert!(tables.contains("{1,{2,1,1,1}}"));
    insta::assert_snapshot!(rules);
    insta::assert_snapshot!(tables);
}

#[test]
fn rejects_unsupported_or_underdetermined_input() {
    assert!(read_family_yaml(&FAMILY.replace("[15]", "[7,15]"), KIN, "box").is_ok());
    assert!(read_family_yaml(FAMILY, &KIN.replace("- [[p1+p3,p1+p3],-s-t]", ""), "box").is_err());
    assert!(read_family_yaml(
        &FAMILY.replace("name: box", "name: box\n    cut_propagators: [1]"),
        KIN,
        "box"
    )
    .is_ok());
    assert!(read_family_mathematica(
        &format!("{MATHEMATICA}\nRun[\"anything\"];"),
        "box",
        &["s", "t"],
        4
    )
    .is_err());
    assert!(read_family_mathematica(
        &MATHEMATICA.replace("-(k+p1)^2", "-k^2"),
        "box",
        &["s", "t"],
        4
    )
    .is_err());
}

#[test]
fn reads_upstream_massive_box_with_redundant_kinematic_rules() {
    let imported = read_family_yaml(
        include_str!("fixtures/massive-integralfamilies.yaml"),
        include_str!("fixtures/massive-kinematics.yaml"),
        "box",
    )
    .unwrap();
    assert_eq!(imported.family.vars, vec!["d", "s", "t"]);
    let system = imported.family.system(1, 1).unwrap();
    let plan = system.learn(&[vec![2, 1, 1, 1]], 7);
    let p = polycore::modp::Primes::new().next().unwrap();
    assert!(plan
        .replay(|e| system.row(e, &[13, 19, 23], p), p)
        .is_some());
}

#[test]
fn polynomial_backend_agrees_with_native_rows_and_reductions() {
    use std::collections::BTreeMap;
    let imported = read_family_yaml(FAMILY, KIN, "box").unwrap();
    let mut native = families::one_loop_box();
    native.symmetries.clear();
    let a = imported.family.system(1, 1).unwrap();
    let b = native.system(1, 1);
    let targets = [vec![2, 1, 1, 1], vec![1, 2, 1, 1]];
    let ap = a.learn(&targets, 17);
    let bp = b.learn(&targets, 17);
    let mut rng = polycore::sample::Rng::new(27);
    for p in polycore::modp::Primes::new().take(3) {
        for _ in 0..6 {
            let x = (0..3).map(|_| rng.nonzero(p)).collect::<Vec<_>>();
            let ar = ap.replay(|e| a.row(e, &x, p), p).unwrap();
            let br = bp.replay(|e| b.row(e, &x, p), p).unwrap();
            for t in 0..targets.len() {
                let am: BTreeMap<_, _> = ap
                    .masters
                    .iter()
                    .enumerate()
                    .filter_map(|(i, &j)| {
                        let c = ar[t * ap.masters.len() + i];
                        (c != 0).then(|| (a.integrals[j].clone(), c))
                    })
                    .collect();
                let bm: BTreeMap<_, _> = bp
                    .masters
                    .iter()
                    .enumerate()
                    .filter_map(|(i, &j)| {
                        let c = br[t * bp.masters.len() + i];
                        (c != 0).then(|| (b.integrals[j].clone(), c))
                    })
                    .collect();
                assert_eq!(am, bm);
            }
        }
    }
}

#[test]
fn polynomial_masses_fractional_momenta_sector_unions_and_cuts() {
    let src =
        "Internal={k};External={q};Propagators={k^2-m^2,(k+q/2)^2-m^2};Replacements={q^2->s};";
    let imported = read_family_mathematica(src, "bubble", &["s", "m"], 2).unwrap();
    let system = imported.family.system(1, 1).unwrap();
    let plan = system.learn(&[vec![2, 1]], 11);
    let p = polycore::modp::Primes::new().next().unwrap();
    assert!(plan
        .replay(|e| system.row(e, &[13, 17, 19], p), p)
        .is_some());
    let mut cuts = imported.family.clone();
    cuts.cuts = vec![0];
    assert!(cuts
        .system(1, 1)
        .unwrap()
        .integrals
        .iter()
        .all(|a| a[0] > 0));
    let mut sectors = imported.family.clone();
    sectors.top_sectors = vec![1, 2];
    assert!(sectors
        .system(1, 1)
        .unwrap()
        .integrals
        .iter()
        .all(|a| a[0] <= 0 || a[1] <= 0));
    sectors.zero_sectors = vec![1];
    assert!(sectors
        .system(1, 1)
        .unwrap()
        .integrals
        .iter()
        .all(|a| a[1] > 0));
    let ring = polycore::Ring::new(["d", "s", "m"], polycore::Order::Lex);
    let candidates = imported.family.denominator_candidates();
    assert!(candidates.contains(&ring.parse("m^2").unwrap()));
    assert!(candidates.contains(&ring.parse("s").unwrap()));
    // An explicitly expanded quadratic has the same meaning as the square.
    let expanded = read_family_mathematica(
        &src.replace("(k+q/2)^2", "k^2+k*q+q^2/4"),
        "bubble",
        &["s", "m"],
        2,
    )
    .unwrap();
    assert_eq!(expanded.family.propagators, imported.family.propagators);
}

#[test]
fn numeric_reference_and_table_roundtrip() {
    use zippel_laporta::formats::read_reduction_tables;
    // Reference version 7.1 d132e5365dd2a13db9cd9dbaf5c200b53d489cfd, tests/box/box.start;
    // reference executable, d=13, s=t=1, prime index 1 (18446744073709551557).
    let table = read_reduction_tables(include_str!("fixtures/reference-box-d13.tables")).unwrap();
    let imported = read_family_yaml(FAMILY, KIN, "box").unwrap();
    let system = imported.family.system(1, 1).unwrap();
    let target = vec![2, 1, 1, 1];
    let plan = system.learn(std::slice::from_ref(&target), 7);
    let p = 18_446_744_073_709_551_557;
    let x = [13, 1, 1];
    let values = plan.replay(|e| system.row(e, &x, p), p).unwrap();
    let got = plan
        .masters
        .iter()
        .zip(values)
        .filter(|(_, c)| *c != 0)
        .map(|(&j, c)| ((1, system.integrals[j].clone()), c))
        .collect();
    assert_eq!(table.evaluate(1, &target, &system.vars, &x, p), Some(got));
    let cs = system.lift(&plan, 7).unwrap();
    let exported = imported.reduction_tables(&system, &plan, &cs, 1).unwrap();
    let roundtrip = read_reduction_tables(&exported).unwrap();
    assert_eq!(
        roundtrip.evaluate(1, &target, &system.vars, &x, p),
        table.evaluate(1, &target, &system.vars, &x, p)
    );
    assert!(read_reduction_tables("{{{1,{{2,\"1\"}}}},{{1,{1,{1}}}}}").is_err());
}

#[test]
fn imports_nonplanar_five_point_family_with_indexed_momenta() {
    let imported = read_family_mathematica(
        include_str!("fixtures/nonplanar-double-pentagon.m"),
        "nonplanar",
        &["s12", "s23", "s34", "s45", "s51"],
        8,
    )
    .unwrap();
    assert_eq!(imported.family.propagators.len(), 11);
    assert_eq!(imported.family.vars.len(), 6);
    let candidates = imported.family.denominator_candidates();
    assert!(candidates.iter().any(|g| g.total_degree() == 4));
    // The maximal cut keeps this regression small while exercising all 12 IBPs.
    let mut family = imported.family;
    family.cuts = (0..8).collect();
    let system = family.system(1, 1).unwrap();
    let plan = system.learn(&[vec![2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0]], 21);
    for p in polycore::modp::Primes::new().take(2) {
        assert!(plan
            .replay(|e| system.row(e, &[13, 17, 19, 23, 29, 31], p), p)
            .is_some());
    }
}

#[test]
fn doublebox_agrees_with_external_table_after_basis_conversion() {
    use polycore::modp::{add, mul};
    use std::collections::BTreeMap;
    let imported = read_family_mathematica(
        include_str!("fixtures/doublebox.m"),
        "doublebox",
        &["s", "t"],
        7,
    )
    .unwrap();
    let table = zippel_laporta::formats::read_reduction_tables(include_str!(
        "fixtures/reference-doublebox-d13-s17-t19.tables"
    ))
    .unwrap();
    let system = imported.family.system(1, 1).unwrap();
    let targets = [
        vec![2, 1, 1, 1, 1, 1, 1, 0, 0],
        vec![1, 2, 1, 1, 1, 1, 1, 0, 0],
    ];
    let plan = system.learn(&targets, 7);
    let p = 18_446_744_073_709_551_557;
    let x = [13, 17, 19];
    let cs = plan.replay(|e| system.row(e, &x, p), p).unwrap();
    for (i, target) in targets.iter().enumerate() {
        let mut composed = BTreeMap::new();
        for (j, &m) in plan.masters.iter().enumerate() {
            let c = cs[i * plan.masters.len() + j];
            if c == 0 {
                continue;
            }
            for (master, v) in table
                .evaluate(1, &system.integrals[m], &system.vars, &x, p)
                .unwrap()
            {
                let entry = composed.entry(master).or_insert(0);
                *entry = add(*entry, mul(c, v, p), p);
            }
        }
        composed.retain(|_, c| *c != 0);
        assert_eq!(
            Some(composed),
            table.evaluate(1, target, &system.vars, &x, p)
        );
    }
}

#[test]
fn mapped_negative_momenta_preserve_square_sign() {
    let a = read_family_mathematica(
        "Internal={k};External={q};Propagators=(-#)^2 & /@ {k,k+q};Replacements={q^2->s};",
        "bubble",
        &["s"],
        2,
    )
    .unwrap();
    let b = read_family_mathematica(
        "Internal={k};External={q};Propagators=#^2 & /@ {k,k+q};Replacements={q^2->s};",
        "bubble",
        &["s"],
        2,
    )
    .unwrap();
    assert_eq!(a.family.propagators, b.family.propagators);
    assert!(read_family_yaml(&FAMILY.replace("[15]", "[0]"), KIN, "box").is_err());
}
