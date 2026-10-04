//! Import a declarative family, reduce a target list, and export replacement tables.
use std::{env, error::Error, fs};
use zippel_laporta::formats::{read_fire, read_kira};

fn main() -> Result<(), Box<dyn Error>> {
    let a: Vec<_> = env::args().skip(1).collect();
    let (imported,offset,fire)=match a.first().map(String::as_str) {
        Some("kira") if a.len()==9=>(read_kira(&fs::read_to_string(&a[1])?,&fs::read_to_string(&a[2])?,&a[3])?,4,false),
        Some("fire") if a.len()==10=>(read_fire(&fs::read_to_string(&a[1])?,&a[2],&a[3].split(',').filter(|s|!s.is_empty()).collect::<Vec<_>>(),a[4].parse()?)?,5,true),
        _=>return Err("usage: reduce_file kira families.yaml kinematics.yaml name targets.txt output dots numerators seed\n   or: reduce_file fire family.m name s,t lines targets.txt output dots numerators seed".into()),
    };
    let targets = fs::read_to_string(&a[offset])?
        .lines()
        .filter(|s| !s.trim().is_empty() && !s.trim().starts_with('#'))
        .map(|s| {
            s.trim()
                .split(',')
                .map(|v| v.trim().parse())
                .collect::<Result<Vec<i32>, _>>()
        })
        .collect::<Result<Vec<_>, _>>()?;
    if targets
        .iter()
        .any(|t| t.len() != imported.family.props.len())
    {
        return Err("target index dimension".into());
    }
    let dots = a[offset + 2].parse()?;
    let numerators = a[offset + 3].parse()?;
    let seed = a[offset + 4].parse()?;
    let system = imported.family.system(dots, numerators);
    if targets.iter().any(|t| !system.integrals.contains(t)) {
        return Err("target outside generated system; increase dots/numerators".into());
    }
    let plan = system.learn(&targets, seed);
    let cs = system.lift(&plan, seed).ok_or("reconstruction failed")?;
    let output = if fire {
        imported.fire_tables(&system, &plan, &cs, 1)?
    } else {
        imported.rules(&system, &plan, &cs, None)?
    };
    fs::write(&a[offset + 1], output)?;
    eprintln!(
        "{} equations, {} retained, {} masters",
        system.len(),
        plan.len(),
        plan.masters.len()
    );
    Ok(())
}
