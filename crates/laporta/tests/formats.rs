#[path = "../examples/families/mod.rs"]
mod families;
use zippel_laporta::formats::{read_fire, read_kira};
const FAMILY: &str = include_str!("fixtures/box.yaml");
const KIN: &str = include_str!("fixtures/kinematics.yaml");
const FIRE: &str = include_str!("fixtures/box.m");

#[test]
fn imports_match_native_family_and_export_signs() {
    let kira = read_kira(FAMILY, KIN, "box").unwrap();
    let fire = read_fire(&FIRE.replace("p1*p2", "p1 p2"), "box", &["s", "t"], 4).unwrap();
    let mut native = families::one_loop_box();
    for (_, mass) in &mut native.props {
        mass.resize(4, 0);
    }
    assert_eq!(kira.family.props, native.props);
    assert_eq!(kira.family.legs, fire.family.legs);
    assert_eq!(kira.family.props, fire.family.props);
    let target = vec![2, 1, 1, 1];
    let (system, plan, cs) = kira.family.reduce(&[target], 1, 1).unwrap();
    let rules = kira.rules(&system, &plan, &cs, None).unwrap();
    assert!(rules.contains("box[2,1,1,1] ->"));
    let tables = fire.fire_tables(&system, &plan, &cs, 1).unwrap();
    assert!(tables.contains("{1,{2,1,1,1}}"));
    insta::assert_snapshot!(rules);
    insta::assert_snapshot!(tables);
}

#[test]
fn rejects_unsupported_or_underdetermined_input() {
    assert!(read_kira(&FAMILY.replace("[15]", "[7,15]"), KIN, "box").is_err());
    assert!(read_kira(FAMILY, &KIN.replace("- [[p1+p3,p1+p3],-s-t]", ""), "box").is_err());
    assert!(read_kira(
        &FAMILY.replace("name: box", "name: box\n    cut_propagators: [1]"),
        KIN,
        "box"
    )
    .is_err());
    assert!(read_fire(
        &format!("{FIRE}\nRun[\"anything\"];"),
        "box",
        &["s", "t"],
        4
    )
    .is_err());
    assert!(read_fire(&FIRE.replace("-(k+p1)^2", "-k^2"), "box", &["s", "t"], 4).is_err());
}

#[test]
fn reads_upstream_massive_box_with_redundant_kinematic_rules() {
    let imported = read_kira(
        include_str!("fixtures/massive-integralfamilies.yaml"),
        include_str!("fixtures/massive-kinematics.yaml"),
        "box",
    )
    .unwrap();
    assert_eq!(imported.family.vars, vec!["d", "s", "t"]);
    assert_eq!(imported.family.props[0].1[0], 1);
    let system = imported.family.system(1, 1);
    let plan = system.learn(&[vec![2, 1, 1, 1]], 7);
    let p = polycore::modp::Primes::new().next().unwrap();
    assert!(plan
        .replay(|e| system.row(e, &[13, 19, 23], p), p)
        .is_some());
}
