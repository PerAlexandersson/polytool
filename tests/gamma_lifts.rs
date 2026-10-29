//! Gamma lifts: peak polynomials proved real-rooted from the Eulerian
//! polynomials, gamma rows of other palindromic parents, and near misses.

use num_bigint::BigInt;
use polytool::oeis;
use polytool::recurrence::{
    certify_gamma_lift, find_gamma_lifts, transport_gamma, BigRational, DegreeSchedule,
    GammaLiftCandidate, GammaLiftCertification, GammaLiftIdentityStatus, GammaLiftOptions,
    GammaLiftParent, GammaLiftReport, GammaLiftSearchOptions, GammaParentSymmetryStatus,
    LiftIdentityStatus, LiftOperation, LiftVerdict, Recurrence, RecurrenceJson,
    RecurrenceTransformation, VeroneseSectionGauge,
};
use polytool::sequences::{self, literature};

type Q = BigRational;

fn q(value: i64) -> Q {
    Q::from_integer(BigInt::from(value))
}

fn rational_rows(rows: &[Vec<BigInt>]) -> Vec<Vec<Q>> {
    rows.iter()
        .map(|row| row.iter().cloned().map(Q::from_integer).collect())
        .collect()
}

fn doubling_gauge() -> VeroneseSectionGauge {
    VeroneseSectionGauge::Ratio {
        numerator: vec!["2".into()],
        denominator: vec!["1".into()],
        initial_value: "1".into(),
    }
}

fn assert_proved_nonpositive(certification: &GammaLiftCertification) {
    assert!(certification.is_proved(), "{certification:#?}");
    certification.verify().expect("the certification replays");
    let summary = &certification.lift_summary;
    assert_eq!(summary.verdict, LiftVerdict::Proved);
    assert_eq!(
        summary.root_transfer.input_root_region.to_string(),
        "(-inf, 0]"
    );
    assert!(matches!(
        summary.identity,
        LiftIdentityStatus::Certified { .. }
    ));
    assert!(certification.parent_symmetry.is_certified());
}

/// The catalog record A008303 (row `m` is `W_(m+1)`) is found as
/// `W_(n+1) = 2^n Gamma(A_(n+1))(z/4)`, certified for every `n` and proved.
#[test]
fn peak_polynomials_from_the_catalog_are_proved_from_eulerian() {
    let entry = oeis::by_id("A008303").expect("catalog entry");
    let (rule, seeds) = entry.recurrence_parts().expect("catalog recurrence");
    let first_index = entry.recurrence_first_index;
    let report = find_gamma_lifts(
        &rule,
        &seeds,
        first_index,
        &GammaLiftSearchOptions::default(),
    )
    .expect("the search runs");
    report.verify().expect("the report replays");
    let proved: Vec<_> = report.proved_lifts().collect();
    assert_eq!(proved.len(), 1, "{report:#?}");
    let lift = proved[0];
    assert_eq!(lift.candidate.parent.name, "eulerian");
    assert_eq!(lift.candidate.index_shift, 1);
    assert_eq!(lift.candidate.scale, "1/4");
    assert_eq!(lift.candidate.gauge, doubling_gauge());
    assert_proved_nonpositive(lift);
    assert!(matches!(
        lift.lift_summary.operation,
        LiftOperation::GammaPolynomial { .. }
    ));
    // The catalog rows are the generator's rows W_1, W_2, ... .
    let generated = rule
        .generate_rows_rational(&seeds, first_index, 20)
        .unwrap();
    let peaks = rational_rows(&literature::peak_polynomials_bigint(20));
    assert_eq!(generated, peaks[1..]);
    // The JSON form round-trips and replays.
    let json = serde_json::to_string(&report).unwrap();
    let back: GammaLiftReport = serde_json::from_str(&json).unwrap();
    assert_eq!(back, report);
    back.verify().unwrap();
}

/// The library recurrence of the peak (= valley) polynomials, with the
/// candidate supplied by the caller.
#[test]
fn peak_polynomials_from_the_library_recurrence_are_proved() {
    let (rule, seeds, first_index) = literature::peak_polynomials_recurrence();
    let candidate = GammaLiftCandidate {
        parent: GammaLiftParent::eulerian(),
        index_shift: 0,
        scale: "1/4".into(),
        gauge: doubling_gauge(),
    };
    let certification = certify_gamma_lift(
        &rule,
        &seeds,
        first_index,
        &candidate,
        &GammaLiftOptions::default(),
    )
    .unwrap();
    assert_proved_nonpositive(&certification);
    let GammaLiftIdentityStatus::Certified { certificate } = &certification.identity else {
        unreachable!()
    };
    // A tampered certificate is rejected.
    let mut tampered = (**certificate).clone();
    tampered.multiplier = vec!["2".into()];
    assert!(tampered.verify().is_err());
    let mut tampered = (**certificate).clone();
    tampered.index_shift = 1;
    assert!(tampered.verify().is_err());
    let mut tampered = (**certificate).clone();
    tampered.scale = "1/2".into();
    assert!(tampered.verify().is_err());
}

/// The input is the gamma rows of `parent` themselves: its rule is the
/// transported gamma rule and its seeds are gamma rows.
fn gamma_input(parent: &GammaLiftParent, seed_count: usize) -> (Recurrence, Vec<Vec<Q>>, usize) {
    let (rule, first_index, seeds) = parent.recurrence.to_recurrence_parts().unwrap();
    let schedule = &parent.symmetry_schedule;
    let gamma = transport_gamma(&rule, schedule).unwrap().scalar.unwrap();
    let rows = rule
        .generate_rows_rational(&seeds, first_index, seed_count)
        .unwrap();
    let map = RecurrenceTransformation::Gamma {
        schedule: schedule.clone(),
    };
    let gamma_rows = rows
        .iter()
        .enumerate()
        .map(|(i, row)| map.map_row(first_index + i, row).unwrap().swap_remove(0))
        .collect();
    (gamma, gamma_rows, first_index)
}

#[test]
fn gamma_rows_of_type_b_eulerian_and_narayana_are_proved() {
    for (parent, seed_count) in [
        (GammaLiftParent::type_b_eulerian(), 1),
        (GammaLiftParent::narayana(), 2),
    ] {
        let (rule, seeds, first_index) = gamma_input(&parent, seed_count);
        let report = find_gamma_lifts(
            &rule,
            &seeds,
            first_index,
            &GammaLiftSearchOptions::default(),
        )
        .unwrap();
        report.verify().unwrap();
        let lift = report
            .proved_lifts()
            .find(|lift| lift.candidate.parent.name == parent.name)
            .unwrap_or_else(|| panic!("{report:#?}"));
        assert_eq!(lift.candidate.index_shift, 0);
        assert_eq!(lift.candidate.scale, "1");
        assert_eq!(lift.candidate.gauge, VeroneseSectionGauge::None);
        assert_proved_nonpositive(lift);
    }
}

/// The menu parent rules reproduce the library generators.
#[test]
fn menu_parents_match_the_generators() {
    let check = |parent: GammaLiftParent, generator: Vec<Vec<Q>>| {
        let (rule, first_index, seeds) = parent.recurrence.to_recurrence_parts().unwrap();
        let rows = rule
            .generate_rows_rational(&seeds, first_index, generator.len() - first_index)
            .unwrap();
        assert_eq!(rows, generator[first_index..], "{}", parent.name);
    };
    check(
        GammaLiftParent::eulerian(),
        rational_rows(&sequences::eulerian_polynomials_bigint(20)),
    );
    check(
        GammaLiftParent::type_b_eulerian(),
        rational_rows(&sequences::type_b_eulerian_polynomials_bigint(20)),
    );
    check(
        GammaLiftParent::narayana(),
        rational_rows(&sequences::narayana_polynomials_bigint(20)),
    );
    let k = q(2);
    let parts: Vec<Vec<Q>> = literature::one_over_k_eulerian_symmetric_parts_rational(20, &k)
        .into_iter()
        .map(|(a, _)| a)
        .collect();
    check(
        GammaLiftParent::one_over_k_eulerian_symmetric_part(&k),
        parts,
    );
}

/// The 1/2-Eulerian symmetric parts as a parent: the identity is certified
/// through the second-order rule; the verdict follows the parent provers.
#[test]
fn symmetric_part_parent_certifies_the_identity() {
    let parent = GammaLiftParent::one_over_k_eulerian_symmetric_part(&q(2));
    let (rule, seeds, first_index) = gamma_input(&parent, 2);
    let candidate = GammaLiftCandidate {
        parent,
        index_shift: 0,
        scale: "1".into(),
        gauge: VeroneseSectionGauge::None,
    };
    let certification = certify_gamma_lift(
        &rule,
        &seeds,
        first_index,
        &candidate,
        &GammaLiftOptions::default(),
    )
    .unwrap();
    certification.verify().unwrap();
    assert!(
        matches!(
            certification.identity,
            GammaLiftIdentityStatus::Certified { .. }
        ),
        "{certification:#?}"
    );
    assert!(certification.parent_symmetry.is_certified());
    if !certification.is_proved() {
        assert_eq!(certification.lift_summary.verdict, LiftVerdict::Open);
    }
}

#[test]
fn near_misses_are_reported_exactly() {
    let (rule, seeds, first_index) = literature::peak_polynomials_recurrence();
    let base = GammaLiftCandidate {
        parent: GammaLiftParent::eulerian(),
        index_shift: 0,
        scale: "1/4".into(),
        gauge: doubling_gauge(),
    };
    let options = GammaLiftOptions::default();
    // A wrong scale: refuted at the first row with two coefficients (n = 3).
    let wrong_scale = GammaLiftCandidate {
        scale: "1/2".into(),
        ..base.clone()
    };
    let certification =
        certify_gamma_lift(&rule, &seeds, first_index, &wrong_scale, &options).unwrap();
    certification.verify().unwrap();
    assert!(matches!(
        certification.identity,
        GammaLiftIdentityStatus::Refuted { index: 3, .. }
    ));
    assert_eq!(certification.lift_summary.verdict, LiftVerdict::RouteFailed);
    // A wrong symmetry schedule (d_n = n): A_1 = 1 is not symmetric at 1.
    let mut parent = GammaLiftParent::eulerian();
    parent.symmetry_schedule = DegreeSchedule::new(1, 1, 1);
    let wrong_schedule = GammaLiftCandidate {
        parent,
        ..base.clone()
    };
    let certification =
        certify_gamma_lift(&rule, &seeds, first_index, &wrong_schedule, &options).unwrap();
    certification.verify().unwrap();
    assert_eq!(
        certification.parent_symmetry,
        GammaParentSymmetryStatus::Refuted { index: 1 }
    );
    assert_eq!(certification.lift_summary.verdict, LiftVerdict::RouteFailed);
    // A wrong gauge: refuted at n = 2.
    let wrong_gauge = GammaLiftCandidate {
        gauge: VeroneseSectionGauge::None,
        ..base.clone()
    };
    let certification =
        certify_gamma_lift(&rule, &seeds, first_index, &wrong_gauge, &options).unwrap();
    assert!(matches!(
        certification.identity,
        GammaLiftIdentityStatus::Refuted { index: 2, .. }
    ));
    // A parent before its first index and a zero scale are invalid input.
    let early = GammaLiftCandidate {
        index_shift: -1,
        ..base.clone()
    };
    assert!(certify_gamma_lift(&rule, &seeds, first_index, &early, &options).is_err());
    let zero = GammaLiftCandidate {
        scale: "0".into(),
        ..base
    };
    assert!(certify_gamma_lift(&rule, &seeds, first_index, &zero, &options).is_err());
}

/// A parent whose rule is not invariant under reversal: the rows
/// `(1 + x)^n` are symmetric, but the induction does not apply, so the
/// identity is supported by rows only and the verdict is open.
#[test]
fn symmetry_without_an_invariant_rule_is_only_supported_by_rows() {
    let parent_json: RecurrenceJson = serde_json::from_str(
        r#"{"schema": "polytool.recurrence.v1", "first_index": 0,
            "initial_polynomials": [["1"]],
            "recurrence": {"terms": [
              {"offset": 1, "deriv_order": 0, "sign": "none",
               "coeff": {"coeffs": [["1", "2"], ["0", "-1"]]}},
              {"offset": 1, "deriv_order": 1, "sign": "none",
               "coeff": {"coeffs": [["0", "1", "1"]]}}]}}"#,
    )
    .unwrap();
    let parent = GammaLiftParent {
        name: "binomial".into(),
        recurrence: parent_json,
        symmetry_schedule: DegreeSchedule::new(0, 0, 1),
    };
    // Q_n = 1: Q_n = Q_(n-1).
    let input_json: RecurrenceJson = serde_json::from_str(
        r#"{"schema": "polytool.recurrence.v1", "first_index": 0,
            "initial_polynomials": [["1"]],
            "recurrence": {"terms": [
              {"offset": 1, "deriv_order": 0, "sign": "none",
               "coeff": {"coeffs": [["1"]]}}]}}"#,
    )
    .unwrap();
    let (rule, first_index, seeds) = input_json.to_recurrence_parts().unwrap();
    let candidate = GammaLiftCandidate {
        parent,
        index_shift: 0,
        scale: "1".into(),
        gauge: VeroneseSectionGauge::None,
    };
    let certification = certify_gamma_lift(
        &rule,
        &seeds,
        first_index,
        &candidate,
        &GammaLiftOptions::default(),
    )
    .unwrap();
    certification.verify().unwrap();
    assert!(matches!(
        certification.parent_symmetry,
        GammaParentSymmetryStatus::SupportedByRows { .. }
    ));
    assert!(matches!(
        certification.identity,
        GammaLiftIdentityStatus::SupportedByRows { .. }
    ));
    assert_eq!(certification.lift_summary.verdict, LiftVerdict::Open);
}
