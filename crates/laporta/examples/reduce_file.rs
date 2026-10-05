//! Import a declarative family, reduce a target list, and export replacement tables.
use std::{env, error::Error, fs};
use zippel_laporta::formats::{read_family_mathematica, read_family_yaml};

fn main() -> Result<(), Box<dyn Error>> {
    let a: Vec<_> = env::args().skip(1).collect();
    let (imported, offset, mathematica) = match a.first().map(String::as_str) {
        Some("yaml") if a.len() == 9 => {
            let family = read_family_yaml(
                &fs::read_to_string(&a[1])?,
                &fs::read_to_string(&a[2])?,
                &a[3],
            )?;
            (family, 4, false)
        }
        Some("mathematica") if a.len() == 10 => {
            let variables = a[3]
                .split(',')
                .filter(|s| !s.is_empty())
                .collect::<Vec<_>>();
            let family = read_family_mathematica(
                &fs::read_to_string(&a[1])?,
                &a[2],
                &variables,
                a[4].parse()?,
            )?;
            (family, 5, true)
        }
        _ => {
            return Err(concat!(
                "usage: reduce_file yaml families.yaml kinematics.yaml name ",
                "targets.txt output dots numerators seed\n",
                "   or: reduce_file mathematica family.m name s,t lines ",
                "targets.txt output dots numerators seed",
            )
            .into())
        }
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
        .any(|t| t.len() != imported.family.propagators.len())
    {
        return Err("target index dimension".into());
    }
    let dots = a[offset + 2].parse()?;
    let numerators = a[offset + 3].parse()?;
    let seed = a[offset + 4].parse()?;
    let system = imported
        .family
        .system(dots, numerators)
        .ok_or("invalid family")?;
    if targets.iter().any(|t| !system.integrals.contains(t)) {
        return Err("target outside generated system; increase dots/numerators".into());
    }
    let plan = system
        .learn(&targets, seed)
        .ok_or("invalid targets or no nonsingular sample")?;
    let candidates = imported.family.denominator_candidates();
    let search = zippel_laporta::block::AdaptiveSearch::default();
    let (plan, cs) = if let Some(reduction) = system.learn_blocks(&plan, seed, &search) {
        eprintln!(
            "{} blocks, {} intermediates, {} discovery probes",
            reduction.first.blocks(),
            reduction.intermediates(),
            reduction.report.oracle_probes
        );
        let cs = reduction
            .lift(&system, &candidates, seed)
            .ok_or("block reconstruction failed")?;
        (reduction.output_plan(), cs)
    } else {
        eprintln!("block search exhausted its bounds; reconstructing the replay oracle");
        let cs = zippel_lift::lift_with_discovered_factors(
            |x, p| plan.replay(|e| system.row(e, x, p), p),
            system.vars.len(),
            &candidates,
            seed,
        )
        .ok_or("reconstruction failed")?;
        (plan, cs)
    };
    let output = if mathematica {
        imported.reduction_tables(&system, &plan, &cs, 1)?
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
