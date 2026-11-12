//! Liu--Wang three-term theorem (`check_liu_wang_three_term_theorem`) on the
//! bundled OEIS catalog: the lag recurrences `P_n = A_n(x) P_(n-1) +
//! B_n(x) P_(n-2)` with `x`-dependent coefficients that satisfy Liu and
//! Wang's Corollary 2.4 (with `b_n = 0`), root windows, derivative-lag rules
//! with half degree growth, their parent proofs, negative cases and tampered
//! reports.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use polytool::recurrence::{
    check_liu_wang_three_term_theorem, prove_parent_real_rootedness, BigRational,
    ParentRecurrenceRecord, RealRootednessTheoremOptions, RealRootednessTheoremReport,
    RealRootednessTheoremStatus, Recurrence, RowInterlacing, VeroneseParentProofRoute,
    VeroneseParentProofStatus,
};

/// The catalog entries whose recurrence satisfies Liu--Wang Corollary 2.4
/// with `b_n = 0` (in the recurrence's own coordinates, after dividing by a
/// power of `x`, after the reversal or after a gauge `x^n`), all of them lag
/// rules without a derivative term that no other polytool route proved.
const COROLLARY_2_4_CLASS: [&str; 80] = [
    "A008306", "A008556", "A034807", "A049403", "A056241", "A061896", "A079510", "A084534",
    "A084950", "A098158", "A102413", "A103450", "A104029", "A105868", "A106828", "A107131",
    "A109983", "A110162", "A110552", "A111049", "A115951", "A119462", "A119808", "A119900",
    "A121314", "A127952", "A128908", "A128966", "A129710", "A129718", "A131222", "A131436",
    "A134239", "A134346", "A136394", "A146986", "A146987", "A146988", "A153516", "A153518",
    "A153520", "A153521", "A153648", "A153649", "A153650", "A153651", "A153652", "A153653",
    "A153654", "A153655", "A153656", "A153657", "A177267", "A180047", "A185263", "A199478",
    "A199881", "A201701", "A201730", "A202023", "A206306", "A208509", "A209414", "A228955",
    "A235706", "A249247", "A249253", "A253284", "A269951", "A281260", "A306364", "A308737",
    "A330381", "A330510", "A330609", "A350512", "A357367", "A360282", "A372873", "A386789",
];

/// Members of the class the check does not prove.  A235706 (`n! L_n^(3-n)(-x)`,
/// `P_n = (x - n + 5) P_(n-1) + (n - 2) x P_(n-2)`): `B_n <= 0` holds on
/// `x <= 0` only, and neither root location applies, because `P_n(0) = 0`
/// from `n = 5` on (the power of `x` grows with `n`, so no fixed power can be
/// divided out) and `[x^0] A_n = 5 - n` is negative from `n = 6` on.
const NOT_PROVED: [&str; 1] = ["A235706"];

/// Members whose rows have roots of both signs (no root region; the parent
/// statement "coefficients of one sign" is false for them).
const BOTH_SIGNS: [&str; 1] = ["A110162"];

/// Catalog entries proved with a root window `[L, U]`, `(-inf, U]` or
/// `[L, +inf)` (`B_n <= 0` only there), with the certified region of the
/// rows (`None`: real-rooted, the region of the reversed rows does not give
/// an interval).
const WINDOW_CLASS: [(&str, Option<&str>); 7] = [
    ("A155112", Some("(-inf, 0]")),
    ("A181289", Some("[-1, +inf)")),
    ("A201972", None),
    ("A209599", Some("(-inf, -1]")),
    ("A236076", Some("[-1, 0]")),
    ("A236471", Some("[-4, +inf)")),
    ("A245962", Some("(-inf, -1]")),
];

/// The catalog family: rule, first index and every initial row (the rule
/// applies after the last one).
fn catalog_family(id: &str) -> (Recurrence, usize, Vec<Vec<BigRational>>) {
    let entry = polytool::oeis::by_id(id).expect("catalog entry");
    let (rule, rows) = entry.recurrence_parts().expect("catalog recurrence");
    (rule, entry.recurrence_first_index, rows)
}

fn check(id: &str) -> RealRootednessTheoremReport {
    let (rule, first_index, seeds) = catalog_family(id);
    check_liu_wang_three_term_theorem(
        &rule,
        first_index,
        &seeds,
        &RealRootednessTheoremOptions::default(),
    )
    .unwrap_or_else(|error| panic!("{id}: {error}"))
}

fn integer_row(row: &[BigRational]) -> Vec<BigInt> {
    let common = row
        .iter()
        .fold(BigInt::one(), |accumulator, c| accumulator.lcm(c.denom()));
    let mut out: Vec<BigInt> = row
        .iter()
        .map(|c| c.numer() * (&common / c.denom()))
        .collect();
    while out.len() > 1 && out.last().is_some_and(Zero::is_zero) {
        out.pop();
    }
    out
}

/// Exact checks of a proved conclusion on the first rows (evidence that the
/// conclusion is consistent with the data; the proof is the report).
fn assert_conclusion_on_rows(id: &str, report: &RealRootednessTheoremReport) {
    let (rule, first_index, seeds) = catalog_family(id);
    let conclusion = report.conclusion.as_ref().expect("a proved report");
    let rows = rule
        .generate_rows_rational(&seeds, first_index, seeds.len() + 16)
        .expect("rows");
    let region = conclusion
        .root_region
        .as_ref()
        .map(|region| region.at(0).expect("an index-independent region"));
    for (offset, row) in rows.iter().enumerate() {
        let n = first_index + offset;
        if n < conclusion.first_index {
            continue;
        }
        let row = integer_row(row);
        assert!(row.iter().any(|c| !c.is_zero()), "{id}: P_{n} is zero");
        assert!(
            polytool::is_real_rooted_bigint_coeffs(&row),
            "{id}: P_{n} is not real-rooted"
        );
        if let Some(region) = &region {
            let in_region = |p: &[BigInt]| {
                let p: Vec<BigRational> =
                    p.iter().cloned().map(BigRational::from_integer).collect();
                let one_signed = |q: &[BigRational]| {
                    q.iter().all(|c| !c.is_negative()) || q.iter().all(|c| !c.is_positive())
                };
                let shifted = |a: &BigRational, b: i64| {
                    polytool::vec_poly::compose(
                        &p,
                        &[a.clone(), BigRational::from_integer(b.into())],
                    )
                };
                let upper = region.upper.as_ref().is_none_or(|u| {
                    let q = shifted(u, 1);
                    one_signed(&q) && (region.upper_closed || !q[0].is_zero())
                });
                let lower = region.lower.as_ref().is_none_or(|l| {
                    let q = shifted(l, -1);
                    one_signed(&q) && (region.lower_closed || !q[0].is_zero())
                });
                upper && lower
            };
            assert!(in_region(&row), "{id}: P_{n} has a root outside the region");
        }
        if conclusion.interlacing == Some(RowInterlacing::Weak)
            && conclusion.interlacing_from.is_some_and(|from| n >= from)
            && offset > 0
        {
            let previous = integer_row(&rows[offset - 1]);
            assert_eq!(
                polytool::check_weak_interlacing_bigint_coeffs(&previous, &row),
                Some(true),
                "{id}: P_{} does not interlace P_{n}",
                n - 1
            );
        }
    }
}

/// Every member of the class except the documented exception is proved from
/// the catalog recurrence alone, the report replays, and the conclusion is
/// consistent with the first rows.
#[test]
fn the_corollary_2_4_class_is_proved() {
    let mut unproved = Vec::new();
    for id in COROLLARY_2_4_CLASS {
        let report = check(id);
        if !report.is_proved() {
            unproved.push(id);
            continue;
        }
        report
            .verify()
            .unwrap_or_else(|error| panic!("{id}: {error}"));
        assert!(report.reference.contains(
            "L. L. Liu and Y. Wang, Adv. Appl. Math. 38 (2007) 542-560, Theorem 2.1 and \
             Corollary 2.4"
        ));
        assert_conclusion_on_rows(id, &report);
        let region = report
            .conclusion
            .as_ref()
            .and_then(|c| c.root_region.as_ref());
        assert_eq!(
            region.is_none(),
            BOTH_SIGNS.contains(&id),
            "{id}: {region:?}"
        );
    }
    assert_eq!(unproved, NOT_PROVED);
}

/// The parent provers prove every member with a root region in `(-inf, 0]`
/// (most of them by the new route; the others by earlier routes), and refute
/// the members with roots of both signs.
#[test]
fn the_parent_provers_prove_the_one_signed_members() {
    let mut by_new_route = 0;
    for id in COROLLARY_2_4_CLASS {
        if NOT_PROVED.contains(&id) {
            continue;
        }
        let (rule, first_index, seeds) = catalog_family(id);
        let record = ParentRecurrenceRecord {
            rule: &rule,
            seeds: &seeds,
            first_index,
        };
        let proof = prove_parent_real_rootedness(&record, first_index, &[], false).unwrap();
        if BOTH_SIGNS.contains(&id) {
            assert_eq!(proof.status, VeroneseParentProofStatus::Refuted, "{id}");
            continue;
        }
        assert_eq!(
            proof.status,
            VeroneseParentProofStatus::Proved,
            "{id}: {proof:#?}"
        );
        if proof.route == Some(VeroneseParentProofRoute::LiuWangThreeTerm) {
            by_new_route += 1;
            assert!(proof.conditional_on.is_empty());
        }
    }
    // 59: A121314, A209414 and A360282 are proved earlier in the cascade by
    // the derivative-lag check, whose degree law comes from the shared
    // leading-coefficient module where a lag multiplier is negative.
    assert!(by_new_route >= 59, "{by_new_route}");
}

/// Derivative-lag rules `P_n = A_n P_(n-1) + V_n P'_(n-1) + B_n P_(n-2)` with
/// half degree growth (`deg P_(2m) = deg P_(2m+1) = m`), Liu--Wang's Theorem
/// 2.3 and Corollary 2.4; `check_liu_wang_derivative_lag_theorem` (strict
/// interlacing, growth one at every step) does not apply to them.  A358623
/// needs the division by `x` (with `A_n + V_n / x`) and a start two rows
/// later (`P_2 = 0`).
#[test]
fn half_growth_derivative_lag_rules_are_proved() {
    for (id, region) in [
        ("A008299", "(-inf, 0)"),
        ("A124324", "(-inf, 0)"),
        ("A136011", "(-inf, 0)"),
        ("A358623", "(-inf, 0]"),
    ] {
        let (rule, first_index, seeds) = catalog_family(id);
        assert!(rule.terms.iter().any(|t| t.deriv_order == 1), "{id}");
        let strict = polytool::recurrence::check_liu_wang_derivative_lag_theorem(
            &rule,
            first_index,
            &seeds,
            &RealRootednessTheoremOptions::default(),
        )
        .unwrap();
        assert!(!strict.is_proved(), "{id}");
        let report = check(id);
        assert!(report.is_proved(), "{id}: {:#?}", report.open_hypotheses());
        report.verify().unwrap();
        assert!(report.reference.contains("Theorem 2.3"), "{id}");
        assert!(report
            .hypotheses
            .iter()
            .any(|h| h.name == "derivative_coefficient_nonpositive"));
        let conclusion = report.conclusion.as_ref().unwrap();
        assert_eq!(conclusion.interlacing, Some(RowInterlacing::Weak));
        assert_eq!(
            conclusion.root_region.as_ref().unwrap().to_string(),
            region,
            "{id}"
        );
        assert_conclusion_on_rows(id, &report);
        let record = ParentRecurrenceRecord {
            rule: &rule,
            seeds: &seeds,
            first_index,
        };
        let proof = prove_parent_real_rootedness(&record, first_index, &[], false).unwrap();
        assert_eq!(proof.status, VeroneseParentProofStatus::Proved, "{id}");
        assert_eq!(
            proof.route,
            Some(VeroneseParentProofRoute::LiuWangThreeTerm)
        );
    }
}

/// The root windows: real-rootedness, the window as root region (after the
/// root `0` of a divided power of `x`), and the RealRooted window routes.
#[test]
fn the_window_class_is_proved() {
    for (id, region) in WINDOW_CLASS {
        let report = check(id);
        assert!(report.is_proved(), "{id}: {:#?}", report.open_hypotheses());
        report.verify().unwrap();
        assert!(report
            .hypotheses
            .iter()
            .any(|h| h.name == "seed_roots_in_window"));
        let certified = report
            .conclusion
            .as_ref()
            .and_then(|c| c.root_region.as_ref())
            .map(ToString::to_string);
        assert_eq!(certified.as_deref(), region, "{id}");
        assert_conclusion_on_rows(id, &report);
    }
    let window = check("A236076");
    assert_eq!(window.lean_route.as_deref(), Some("window"));
    assert_eq!(
        window.lean_declarations,
        ["RealRooted.threeTerm_interlaces_of_roots_mem_Icc"]
    );
}

/// The coordinates: the reversal (A098158, whose rows gain a factor `x`
/// every other step), the gauge `x^n` (A109983) and the gauge `(x + 1)^n`
/// (A090452) give real-rootedness with a root region but no interlacing of
/// the original rows.
#[test]
fn the_reversal_and_gauge_coordinates() {
    for (id, marker) in [
        ("A098158", "P_n(1/x)"),
        ("A109983", "Q_n = P_n / x^(n"),
        ("A090452", "Q_n = P_n / (x + 1)^"),
    ] {
        let report = check(id);
        assert!(report.is_proved(), "{id}: {report:#?}");
        assert!(
            report.normal_form.contains(marker),
            "{id}: {}",
            report.normal_form
        );
        let conclusion = report.conclusion.as_ref().unwrap();
        assert_eq!(conclusion.interlacing, None);
        assert!(conclusion.root_region.is_some());
        report.verify().unwrap();
        assert_conclusion_on_rows(id, &report);
    }
}

/// Lucas polynomials: the half-growth degree law (degrees 0, 0, 1, 1, 2, ...),
/// `B_n = x <= 0` on `x <= 0` and nonnegative coefficients; the Narayana-type
/// square `B_n = -x^2` (A084534) is Liu--Wang's Proposition 3.3.
#[test]
fn half_growth_and_proposition_3_3() {
    let lucas = check("A034807");
    assert!(lucas.is_proved());
    assert!(lucas
        .hypotheses
        .iter()
        .any(|h| h.name.starts_with("degree_law_from_step_")));
    assert_eq!(lucas.lean_route, None);
    let square = check("A084534");
    assert!(square.is_proved());
    assert!(
        square.reference.contains("Proposition 3.3"),
        "{}",
        square.reference
    );
    assert_eq!(square.lean_route.as_deref(), Some("sign"));
}

/// A039692 (`P_n = (x + n) P_(n-1) + n (n - 1 + 2x) P_(n-2)`): `B_n` is
/// positive at `x = -1`, a root of every row, so the sign condition fails;
/// the check must not claim the family.
#[test]
fn a039692_is_not_claimed() {
    let report = check("A039692");
    assert_eq!(report.status, RealRootednessTheoremStatus::Unresolved);
    assert!(report.conclusion.is_none());
    let sign = report
        .hypotheses
        .iter()
        .find(|h| h.name == "lag_coefficient_nonpositive")
        .expect("the sign hypothesis");
    assert_ne!(
        sign.status,
        polytool::recurrence::LiftHypothesisStatus::Verified
    );
    report.verify().unwrap();
    let (rule, first_index, seeds) = catalog_family("A039692");
    let record = ParentRecurrenceRecord {
        rule: &rule,
        seeds: &seeds,
        first_index,
    };
    let proof = prove_parent_real_rootedness(&record, first_index, &[], false).unwrap();
    assert_ne!(proof.status, VeroneseParentProofStatus::Proved);
    assert!(proof
        .attempts
        .iter()
        .any(|a| a.route == VeroneseParentProofRoute::LiuWangThreeTerm && !a.success));
}

/// A tampered report fails `verify`: a changed conclusion, a hypothesis
/// marked verified, changed seeds, and an exchanged sign certificate.
#[test]
fn tampered_reports_fail_verification() {
    let report = check("A056241");
    assert!(report.is_proved());
    report.verify().unwrap();

    let mut tampered = report.clone();
    tampered.conclusion.as_mut().unwrap().root_region = None;
    assert!(tampered.verify().is_err());

    let mut tampered = report.clone();
    tampered.conclusion.as_mut().unwrap().interlacing = Some(RowInterlacing::Strict);
    assert!(tampered.verify().is_err());

    let mut tampered = report.clone();
    let seed = tampered.recurrence.initial_polynomials.last_mut().unwrap();
    seed[0] = "5".into();
    assert!(tampered.verify().is_err());

    let with_certificates: Vec<usize> = report
        .hypotheses
        .iter()
        .enumerate()
        .filter(|(_, h)| h.sign_decision.is_some())
        .map(|(i, _)| i)
        .collect();
    assert!(!with_certificates.is_empty(), "{report:#?}");
    let mut tampered = report.clone();
    tampered.hypotheses[with_certificates[0]].sign_decision = None;
    assert!(tampered.verify().is_err());
    // A certificate from another family's report.
    let other = check("A084534");
    let foreign = other
        .hypotheses
        .iter()
        .filter_map(|h| h.sign_decision.clone())
        .find(|decision| {
            Some(decision)
                != report.hypotheses[with_certificates[0]]
                    .sign_decision
                    .as_ref()
        })
        .expect("a certificate");
    let mut tampered = report.clone();
    tampered.hypotheses[with_certificates[0]].sign_decision = Some(foreign);
    assert!(tampered.verify().is_err());

    let failing = check("A039692");
    let mut tampered = failing.clone();
    for hypothesis in &mut tampered.hypotheses {
        hypothesis.status = polytool::recurrence::LiftHypothesisStatus::Verified;
    }
    tampered.status = RealRootednessTheoremStatus::Proved;
    assert!(tampered.verify().is_err());
}
