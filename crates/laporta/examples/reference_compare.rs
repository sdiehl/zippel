//! Compare against an installed reference reducer, including changes of master basis.
//! `cargo run --release -p zippel-laporta --example reference_compare -- /path/to/root /path/to/executable [double-box]`
use polycore::modp::{add, mul};
use std::{
    collections::BTreeMap, env, error::Error, fs, path::PathBuf, process::Command, time::Instant,
};
use zippel_laporta::formats::{read_family_mathematica, read_reduction_tables};

#[allow(clippy::too_many_lines)]
fn main() -> Result<(), Box<dyn Error>> {
    let mut args = env::args().skip(1);
    let root = fs::canonicalize(
        args.next()
            .ok_or("provide the reference reducer directory")?,
    )?;
    // Preserve the executable name: some programs locate helpers through argv[0].
    let executable =
        std::path::absolute(args.next().ok_or("provide the reference executable path")?)?;
    let double = args.next().is_some_and(|v| v == "double-box");
    let (src, start, lines, targets) = if double {
        (
            include_str!("../tests/fixtures/doublebox.m"),
            "examples/doublebox.start",
            7,
            vec![
                vec![2, 1, 1, 1, 1, 1, 1, 0, 0],
                vec![1, 2, 1, 1, 1, 1, 1, 0, 0],
            ],
        )
    } else {
        // The upstream box.start uses positive propagators and s=t=1.
        (
            include_str!("../tests/fixtures/box.m"),
            "tests/box/box.start",
            4,
            vec![vec![2, 1, 1, 1], vec![1, 2, 1, 1]],
        )
    };
    let src = if double {
        src.to_owned()
    } else {
        src.replace("-k^2", "k^2").replace("-(k+", "(k+")
    };
    let imported = read_family_mathematica(&src, "reference", &["s", "t"], lines)?;
    let analysis = imported.family.analyze().ok_or("invalid family")?;
    let system = imported.family.system(1, 1).ok_or("invalid family")?;
    let plan = system.learn(&targets, 7).unwrap();
    eprintln!(
        "{} equations, {} retained, {} masters",
        system.len(),
        plan.len(),
        plan.masters.len()
    );
    let folder = env::temp_dir().join(format!("zippel-reference-{}", std::process::id()));
    fs::create_dir(&folder)?;
    let mut all = targets.clone();
    all.extend(plan.masters.iter().map(|&i| system.integrals[i].clone()));
    all.sort();
    all.dedup();
    let integrals = all
        .iter()
        .map(|a| {
            format!(
                "{{1,{{{}}}}}",
                a.iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>()
                    .join(",")
            )
        })
        .collect::<Vec<_>>();
    fs::write(
        folder.join("targets.m"),
        format!("{{{}}}", integrals.join(",")),
    )?;
    let config = folder.join("reference");
    let variables = if double { "d,s,t" } else { "d" };
    let database = folder.join("db");
    let input = folder.join("targets.m");
    let output = folder.join("result.tables");
    let configuration = format!(
        concat!(
            "#threads 1\n#fthreads 1\n#variables {}\n#database {}\n",
            "#start\n#folder {}/\n#problem 1 {}\n#integrals {}\n#output {}\n",
        ),
        variables,
        database.display(),
        root.display(),
        start,
        input.display(),
        output.display(),
    );
    fs::write(config.with_extension("config"), configuration)?;
    for (index, p) in [
        (1, 18_446_744_073_709_551_557),
        (2, 18_446_744_073_709_551_533),
    ] {
        let x = if double {
            vec![13, 17, 19]
        } else {
            vec![13, 1, 1]
        };
        let values = if double {
            format!("13_17_19_{index}")
        } else {
            format!("13_{index}")
        };
        let now = Instant::now();
        let result = Command::new(&executable)
            .current_dir(&root)
            .args(["-c", config.to_str().ok_or("config path")?, "-v", &values])
            .output()?;
        fs::write(
            folder.join(format!("reference-{index}.log")),
            &result.stdout,
        )?;
        if !result.status.success() {
            return Err(format!(
                "Reference reducer failed: {}; logs in {}",
                String::from_utf8_lossy(&result.stderr),
                folder.display()
            )
            .into());
        }
        let reference_time = now.elapsed();
        let table = read_reduction_tables(&fs::read_to_string(
            folder.join(format!("result_{values}.tables")),
        )?)?;
        let now = Instant::now();
        let cs = plan
            .replay(|e| system.row(e, &x, p), p)
            .ok_or("replay failed")?;
        let replay_time = now.elapsed();
        for (i, target) in targets.iter().enumerate() {
            let mut composed = BTreeMap::new();
            for (j, &m) in plan.masters.iter().enumerate() {
                let c = cs[i * plan.masters.len() + j];
                if c == 0 {
                    continue;
                }
                for (master, v) in table
                    .evaluate(1, &system.integrals[m], &system.vars, &x, p)
                    .ok_or("missing master reduction")?
                {
                    let entry = composed.entry(master).or_insert(0);
                    *entry = add(*entry, mul(c, v, p), p);
                }
            }
            composed.retain(|_, v| *v != 0);
            let reference = table
                .evaluate(1, target, &system.vars, &x, p)
                .ok_or("missing target reduction")?;
            let canonicalize = |terms: BTreeMap<(u32, Vec<i32>), u64>| {
                let mut result = BTreeMap::new();
                for ((family, a), c) in terms {
                    if let Some(a) = analysis.canonical_index(&a) {
                        let value = result.entry((family, a)).or_insert(0);
                        *value = add(*value, c, p);
                    }
                }
                result.retain(|_, c| *c != 0);
                result
            };
            let composed = canonicalize(composed);
            let reference = canonicalize(reference);
            if composed != reference {
                return Err(format!(
                    "mismatch for {target:?} at prime {p}: {composed:?} versus {reference:?}"
                )
                .into());
            }
        }
        println!(
            "prime {p}: {} targets agree; reference process {:?}, learned replay {:?}",
            targets.len(),
            reference_time,
            replay_time
        );
    }
    println!(
        "Inputs, tables, and logs: {}",
        PathBuf::from(&folder).display()
    );
    Ok(())
}
