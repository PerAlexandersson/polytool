//! A fixed set of real workloads of the `ore-algebra` crate and of the
//! image recurrences built on it, with recorded output fingerprints.
//!
//! The regression test runs every workload and compares a fingerprint of
//! its complete output (the `Debug` form of every returned rule and
//! certificate) with `tests/fixtures/ore-algebra-workload-fingerprints.tsv`,
//! so an optimisation of the algebra cannot change a result unnoticed.
//! The Hadamard-factor workloads, whose time is mostly spent outside the
//! crate, are compared only by the ignored full test.  Re-bless only for an
//! intentional change of results: `POLYTOOL_BLESS_ORE_WORKLOADS=1 cargo test
//! --release --test ore_algebra_workloads -- --ignored every_`.
//!
//! The ignored timing test prints the median wall time of each workload:
//! `cargo test --release --test ore_algebra_workloads -- --ignored --nocapture time_`
//! (`ORE_WORKLOAD_REPEAT` sets the number of runs, `ORE_WORKLOAD_ONLY` a
//! name filter).

use num_bigint::BigInt;
use ore_algebra::{
    find_annihilator_system, find_difference_image_annihilator, push_forward_substitution,
    right_pseudo_division, AnnihilatorSystemBudget, BivariatePolynomial,
    DifferenceImageAnnihilatorBudget, DifferenceRowMapStep, DifferentialDifferenceOperator,
    OreOperator, Substitution,
};
use polytool::recurrence::{
    apply_basis_map_to_family, find_hadamard_factor_lifts, image_recurrence_of_operator, BasisMap,
    BasisMapImageOptions, BivarPoly, BivarPolyJson, DegreeSchedule, HadamardFactorLiftOptions,
    ImageOperatorSpecification, ImageRecurrenceOptions, ImageRowMap, ImplicitRelation,
    OperatorTermJson, PreserverParentLibrary, Recurrence, RecurrenceJson, RecurrenceTerm,
};
use polytool::{vec_poly, BigRational};
use std::collections::BTreeMap;
use std::time::Instant;

type Q = BigRational;

const FINGERPRINTS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/fixtures/ore-algebra-workload-fingerprints.tsv"
);

mod common;
use common::q;

/// 64-bit FNV-1a of the text (stable across platforms and toolchains).
fn fingerprint(text: &str) -> String {
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    for byte in text.bytes() {
        hash ^= u64::from(byte);
        hash = hash.wrapping_mul(0x0100_0000_01b3);
    }
    format!("{hash:016x}:{}", text.len())
}

// ---------------------------------------------------------------------------
// Inputs
// ---------------------------------------------------------------------------

fn parent(name: &str) -> RecurrenceJson {
    PreserverParentLibrary::named_family_records()
        .into_iter()
        .find(|(n, _)| n == name)
        .unwrap_or_else(|| panic!("{name} is a named record"))
        .1
}

/// `coeffs[a][b]` multiplies `n^a x^b`.
fn term(deriv_order: usize, coeffs: &[&[i64]]) -> OperatorTermJson {
    OperatorTermJson {
        deriv_order,
        coeff: BivarPolyJson {
            coeffs: coeffs
                .iter()
                .map(|row| row.iter().map(|v| v.to_string()).collect())
                .collect(),
        },
    }
}

fn weyl(terms: Vec<OperatorTermJson>) -> ImageOperatorSpecification {
    ImageOperatorSpecification::new(ImageRowMap::WeylOperator {
        terms,
        denominator: vec!["1".into()],
    })
}

fn multiplier(numerator: &[&str], denominator: &[&str]) -> ImageOperatorSpecification {
    ImageOperatorSpecification::new(ImageRowMap::CoefficientMultiplier {
        numerator: numerator.iter().map(|s| s.to_string()).collect(),
        denominator: denominator.iter().map(|s| s.to_string()).collect(),
    })
}

/// The image recurrences: (name, parent, map) for the preserver menu shapes.
fn image_cases() -> Vec<(&'static str, &'static str, ImageOperatorSpecification)> {
    let theta_theta_plus_one = || weyl(vec![term(1, &[&[0, 2]]), term(2, &[&[0, 0, 1]])]);
    let one_plus_d = || weyl(vec![term(0, &[&[1]]), term(1, &[&[1]])]);
    vec![
        ("image/eulerian/1+D", "eulerian", one_plus_d()),
        (
            "image/eulerian/D o x",
            "eulerian",
            weyl(vec![term(0, &[&[1]]), term(1, &[&[0, 1]])]),
        ),
        (
            "image/eulerian/x D",
            "eulerian",
            weyl(vec![term(1, &[&[0, 1]])]),
        ),
        (
            "image/eulerian/D^2",
            "eulerian",
            weyl(vec![term(2, &[&[1]])]),
        ),
        (
            "image/eulerian/theta(theta+1)",
            "eulerian",
            theta_theta_plus_one(),
        ),
        (
            "image/eulerian/affine 2x+1",
            "eulerian",
            ImageOperatorSpecification::new(ImageRowMap::Affine {
                slope: "2".into(),
                intercept: "1".into(),
            }),
        ),
        (
            "image/narayana/theta(theta+1)",
            "narayana",
            theta_theta_plus_one(),
        ),
        ("image/narayana/1+D", "narayana", one_plus_d()),
        (
            "image/narayana/reversal",
            "narayana",
            ImageOperatorSpecification::new(ImageRowMap::Reversal {
                schedule: DegreeSchedule::new(1, 0, 1),
            }),
        ),
        (
            "image/narayana/f-polynomial",
            "narayana",
            ImageOperatorSpecification::new(ImageRowMap::Mobius {
                a: "1".into(),
                b: "0".into(),
                c: "1".into(),
                e: "1".into(),
                schedule: DegreeSchedule::new(1, 0, 1),
            }),
        ),
        (
            "image/binomial/1/k!",
            "binomial",
            multiplier(&["1"], &["1", "1"]),
        ),
        (
            "image/binomial/1/(1/2)_k",
            "binomial",
            multiplier(&["1"], &["1/2", "1"]),
        ),
        ("image/binomial/1+D", "binomial", one_plus_d()),
        (
            "image/delannoy/theta(theta+1)",
            "delannoy",
            theta_theta_plus_one(),
        ),
    ]
}

fn image_workload(parent_name: &str, specification: &ImageOperatorSpecification) -> String {
    let record = parent(parent_name);
    let result = image_recurrence_of_operator(
        &record,
        record.first_index,
        specification,
        &ImageRecurrenceOptions::default(),
    );
    format!("{result:?}")
}

/// `bivariate(&[(i, j, c)])`: `sum c n^i x^j`.
fn bivariate(monomials: &[(usize, usize, i64)]) -> BivarPoly {
    let rows = monomials.iter().map(|m| m.0).max().unwrap_or(0) + 1;
    let width = monomials.iter().map(|m| m.1).max().unwrap_or(0) + 1;
    let mut coeffs = vec![vec![q(0); width]; rows];
    for (i, j, c) in monomials {
        coeffs[*i][*j] += q(*c);
    }
    BivarPoly::from_coefficients_normalized(coeffs)
}

/// `R_n = (u + c)^n`, `n >= 0`.
fn power_family(c: i64) -> RecurrenceJson {
    let rule = Recurrence {
        terms: vec![RecurrenceTerm::new(
            1,
            0,
            bivariate(&[(0, 0, c), (0, 1, 1)]),
        )],
        denominator: None,
        inhomogeneous: None,
    };
    RecurrenceJson::from_recurrence_rational(&rule, 0, &[vec![q(1)]], None)
}

fn basis_map_workload(map: &str, c: i64) -> String {
    let map = BasisMap::from_oeis(map).unwrap();
    let image = apply_basis_map_to_family(&map, &power_family(c), &BasisMapImageOptions::default());
    format!("{image:?}")
}

/// The descent polynomials of 321-avoiding permutations, `n >= 1`.
const AV321_DESCENTS: &str = r#"{"schema":"polytool.recurrence.v1","first_index":1,
  "initial_polynomials":[["1"],["1","1"],["1","4"]],
  "recurrence":{"terms":[
    {"offset":1,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["-1","-1"],["5","-1"]]}},
    {"offset":2,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["10","-10"],["-8","8"]]}},
    {"offset":3,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["-8","16","-8"],["4","-8","4"]]}}],
  "denominator":{"coeffs":[["1"],["1"]]}}}"#;

/// The peak polynomials of 123-avoiding permutations, `n >= 1`.
const AV123_PEAKS: &str = r#"{"schema":"polytool.recurrence.v1","first_index":1,
  "initial_polynomials":[["1"],["2"],["3","2"]],
  "recurrence":{"terms":[
    {"offset":1,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["0","-2"],["5","-1"]]}},
    {"offset":2,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["6","-6"],["-8","8"]]}},
    {"offset":3,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["-4","8","-4"],["4","-8","4"]]}}],
  "denominator":{"coeffs":[["1"],["1"]]}}}"#;

/// The relation `L` (base `n`) of a record and its rows from the first index.
fn relation_and_rows(record: &RecurrenceJson, count: usize) -> (OreOperator, usize, Vec<Vec<Q>>) {
    let (rule, first, seeds) = record.to_recurrence_parts().unwrap();
    let (relation, _) = ImplicitRelation::from_recurrence(&rule).unwrap();
    let mut op = OreOperator::default();
    for term in &relation.terms {
        op.add_term(
            (term.shift, term.deriv_order),
            BivariatePolynomial::new(term.coeff.coeffs.clone()),
        );
    }
    let rows = rule
        .generate_rows_rational(&seeds, first, count.max(seeds.len()))
        .unwrap()
        .into_iter()
        .take(count)
        .map(|row| vec_poly::trim(&row))
        .collect();
    (op, first, rows)
}

/// The Hadamard-factor lifts of a record (row-annihilator systems of the
/// parent through the crate's kernel and operator algebra).
fn hadamard_workload(record: &RecurrenceJson) -> String {
    let (rule, first, seeds) = record.to_recurrence_parts().unwrap();
    let report =
        find_hadamard_factor_lifts(&rule, &seeds, first, &HadamardFactorLiftOptions::default());
    format!("{report:?}")
}

fn annihilator_system_workload(record: &RecurrenceJson) -> String {
    let (relation, first, rows) = relation_and_rows(record, 32);
    let system =
        find_annihilator_system(&relation, &rows, first, &AnnihilatorSystemBudget::default());
    format!("{system:?}")
}

/// `rows[i][j]` is the coefficient of `n^i x^j`.
fn biv(rows: &[&[i64]]) -> BivariatePolynomial {
    BivariatePolynomial::new(
        rows.iter()
            .map(|row| row.iter().map(|v| q(*v)).collect())
            .collect(),
    )
}

/// `S - step` for a first-order rule `F_(n+1) = step F_n`.
fn first_order_parent(step: &OreOperator) -> OreOperator {
    OreOperator::single(1, 0, BivariatePolynomial::one()).subtract(step)
}

/// Eulerian: `A_(n+1) = (1 + n x) A_n + x (1 - x) A_n'`.
fn eulerian_step() -> OreOperator {
    let mut step = OreOperator::single(0, 0, biv(&[&[1], &[0, 1]]));
    step.add_term((0, 1), biv(&[&[0, 1, -1]]));
    step
}

/// Touchard: `T_(n+1) = x (T_n + T_n')`.
fn touchard_step() -> DifferentialDifferenceOperator {
    let mut step = DifferentialDifferenceOperator::multiplication(biv(&[&[0, 1]]));
    step.add_term((0, 1, 0), biv(&[&[0, 1]]));
    step
}

fn difference_workload(
    parent: &DifferentialDifferenceOperator,
    steps: &[DifferenceRowMapStep],
) -> String {
    let search = find_difference_image_annihilator(
        parent,
        steps,
        &DifferenceImageAnnihilatorBudget::default(),
    );
    format!("{search:?}")
}

fn touchard_parent() -> DifferentialDifferenceOperator {
    DifferentialDifferenceOperator::single(1, 0, 0, BivariatePolynomial::one())
        .subtract(&touchard_step())
}

fn substitutions() -> Vec<Substitution> {
    vec![
        Substitution::FPolynomial {
            degree: vec![q(-1), q(1)],
        },
        Substitution::Reversal {
            degree: vec![q(-1), q(1)],
        },
        Substitution::Mobius(Box::new(ore_algebra::MobiusMap {
            a: q(2),
            b: q(-1),
            c: q(3),
            e: q(5),
            scale: Q::new(BigInt::from(-1), BigInt::from(2)),
            degree: vec![q(1), q(1)],
        })),
        Substitution::Affine {
            slope: q(2),
            intercept: q(1),
        },
        Substitution::PowerGauge {
            power: vec![q(1), q(1)],
        },
    ]
}

fn substitution_workload() -> String {
    let parent = first_order_parent(&eulerian_step());
    let embedded = DifferentialDifferenceOperator::from_ore_operator(&parent);
    let mut out = String::new();
    for substitution in substitutions() {
        let forward = push_forward_substitution(&parent, &substitution);
        out.push_str(&format!("{forward:?}\n"));
        let steps = [DifferenceRowMapStep::Substitution(substitution)];
        out.push_str(&difference_workload(&embedded, &steps));
        out.push('\n');
    }
    out
}

/// A deterministic pseudo-random operator with small integer coefficients.
fn random_operator(seed: u64, shift: usize, derivative: usize, degree: usize) -> OreOperator {
    let mut state = seed;
    let mut next = move || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((state >> 33) % 7) as i64 - 3
    };
    let mut op = OreOperator::default();
    for i in 0..=shift {
        for d in 0..=derivative {
            let rows: Vec<Vec<Q>> = (0..=degree)
                .map(|_| (0..=degree).map(|_| q(next())).collect())
                .collect();
            op.add_term((i, d), BivariatePolynomial::new(rows));
        }
    }
    op
}

fn division_workload(divisor: &OreOperator, seed: u64, shift: usize) -> String {
    let dividend = random_operator(seed, shift, 3, 4);
    let division = right_pseudo_division(&dividend, divisor, shift);
    format!("{division:?}")
}

type Workload = (String, Box<dyn Fn() -> String>);

fn workloads() -> Vec<Workload> {
    let mut out: Vec<Workload> = Vec::new();
    for (name, parent_name, specification) in image_cases() {
        out.push((
            name.into(),
            Box::new(move || image_workload(parent_name, &specification)),
        ));
    }
    for c in [1, 2, 4] {
        out.push((
            format!("basis-map/A035607/(u+{c})^n"),
            Box::new(move || basis_map_workload("A035607", c)),
        ));
    }
    for (name, json) in [
        ("annihilator-system/av321", AV321_DESCENTS),
        ("annihilator-system/av123", AV123_PEAKS),
    ] {
        out.push((
            name.into(),
            Box::new(move || annihilator_system_workload(&serde_json::from_str(json).unwrap())),
        ));
    }
    for (name, json) in [
        ("hadamard-lift/av321", AV321_DESCENTS),
        ("hadamard-lift/av123", AV123_PEAKS),
    ] {
        out.push((
            name.into(),
            Box::new(move || hadamard_workload(&serde_json::from_str(json).unwrap())),
        ));
    }
    for name in ["narayana", "binomial", "delannoy"] {
        out.push((
            format!("annihilator-system/{name}"),
            Box::new(move || annihilator_system_workload(&parent(name))),
        ));
    }
    out.push((
        "difference/touchard/Delta".into(),
        Box::new(|| {
            difference_workload(
                &touchard_parent(),
                &[DifferenceRowMapStep::forward_difference(1)],
            )
        }),
    ));
    out.push((
        "difference/touchard/Delta^2".into(),
        Box::new(|| {
            difference_workload(
                &touchard_parent(),
                &[DifferenceRowMapStep::forward_difference(2)],
            )
        }),
    ));
    out.push((
        "difference/eulerian/Delta".into(),
        Box::new(|| {
            let parent = first_order_parent(&eulerian_step());
            difference_workload(
                &DifferentialDifferenceOperator::from_ore_operator(&parent),
                &[DifferenceRowMapStep::forward_difference(1)],
            )
        }),
    ));
    out.push((
        "substitution/eulerian/five maps".into(),
        Box::new(substitution_workload),
    ));
    out.push((
        "division/eulerian/shift 12".into(),
        Box::new(|| division_workload(&first_order_parent(&eulerian_step()), 7, 12)),
    ));
    out.push((
        "division/av123/shift 12".into(),
        Box::new(|| {
            let (relation, _, _) =
                relation_and_rows(&serde_json::from_str(AV123_PEAKS).unwrap(), 4);
            division_workload(&relation, 11, 12)
        }),
    ));
    out
}

fn recorded() -> BTreeMap<String, String> {
    std::fs::read_to_string(FINGERPRINTS)
        .unwrap_or_default()
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .filter_map(|line| {
            let (name, value) = line.split_once('\t')?;
            Some((name.to_string(), value.to_string()))
        })
        .collect()
}

/// Workloads dominated by code outside the crate (the recurrence fits of
/// the Hadamard-factor search); checked by the ignored full comparison.
fn is_slow(name: &str) -> bool {
    name.starts_with("hadamard-lift/")
}

fn compare_fingerprints(include_slow: bool) {
    let computed: BTreeMap<String, String> = workloads()
        .into_iter()
        .filter(|(name, _)| include_slow || !is_slow(name))
        .map(|(name, run)| (name, fingerprint(&run())))
        .collect();
    if include_slow && std::env::var_os("POLYTOOL_BLESS_ORE_WORKLOADS").is_some() {
        let mut text = String::from(
            "# Fingerprints (FNV-1a 64 of the Debug output, and its length) of\n\
             # tests/ore_algebra_workloads.rs; re-bless only for intentional changes.\n",
        );
        for (name, value) in &computed {
            text.push_str(&format!("{name}\t{value}\n"));
        }
        std::fs::write(FINGERPRINTS, text).unwrap();
        return;
    }
    let expected: BTreeMap<String, String> = recorded()
        .into_iter()
        .filter(|(name, _)| include_slow || !is_slow(name))
        .collect();
    assert_eq!(computed, expected);
}

#[test]
fn ore_algebra_workload_outputs_are_unchanged() {
    compare_fingerprints(false);
}

/// Also the slow workloads; blesses with `POLYTOOL_BLESS_ORE_WORKLOADS=1`.
#[test]
#[ignore = "slow in debug builds: cargo test --release --test ore_algebra_workloads -- --ignored"]
fn every_ore_algebra_workload_output_is_unchanged() {
    compare_fingerprints(true);
}

#[test]
#[ignore = "timing run: cargo test --release --test ore_algebra_workloads -- --ignored --nocapture"]
fn time_ore_algebra_workloads() {
    let repeat: usize = std::env::var("ORE_WORKLOAD_REPEAT")
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(3)
        .max(1);
    let only = std::env::var("ORE_WORKLOAD_ONLY").unwrap_or_default();
    let mut total = 0.0;
    println!("workload\tmedian_ms\tmin_ms\tfingerprint");
    for (name, run) in workloads() {
        if !name.contains(&only) {
            continue;
        }
        let mut times = Vec::with_capacity(repeat);
        let mut print = String::new();
        for _ in 0..repeat {
            let start = Instant::now();
            print = fingerprint(&run());
            times.push(start.elapsed().as_secs_f64() * 1000.0);
        }
        times.sort_by(f64::total_cmp);
        let median = times[times.len() / 2];
        total += median;
        println!("{name}\t{median:.1}\t{:.1}\t{print}", times[0]);
    }
    println!("total\t{total:.1}");
}
