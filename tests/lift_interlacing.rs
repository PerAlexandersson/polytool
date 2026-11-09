//! Interlacing of consecutive input rows as a proved property of the lifts
//! (`LiftSummary::input_interlacing`): the Veronese-section rule
//! (Athanasiadis--Wagner), the fixed-preserver rule (Hermite--Kakeya--
//! Obreschkoff; `n`-dependent Hadamard factors get no claim), the common-factor
//! rule of row factors with a constant exponent, near misses (`n`-dependent
//! operators, parents without a proved interlacing) and replay.

use num_bigint::BigInt;
use polytool::check_weak_interlacing_bigint_coeffs;
use polytool::oeis;
use polytool::recurrence::*;

mod common;
use common::q;

fn integer(row: &[BigRational]) -> Vec<BigInt> {
    row.iter()
        .map(|c| {
            assert!(c.is_integer());
            c.to_integer()
        })
        .collect()
}

/// `coeffs[i][j]` = coefficient of `n^i t^j`.
fn bivar(coeffs: &[&[i64]]) -> BivarPoly {
    BivarPoly::from_coefficients_normalized(
        coeffs
            .iter()
            .map(|r| r.iter().map(|&v| q(v)).collect())
            .collect(),
    )
}

/// `P_n = sum c_k(n, t) D^d P_(n-k)` from `(k, d, c_k)`.
fn rule(terms: &[(usize, usize, &[&[i64]])]) -> Recurrence {
    Recurrence {
        terms: terms
            .iter()
            .map(|(offset, order, coeffs)| RecurrenceTerm::new(*offset, *order, bivar(coeffs)))
            .collect(),
        denominator: None,
        inhomogeneous: None,
    }
}

/// Check `Q_(n-1) << Q_n` exactly on the rows `n = from, ..., from + count -
/// 1` (evidence for the proved claim, not part of it).
fn assert_rows_interlace(rows: &[Vec<BigRational>], first_index: usize, from: usize) {
    for n in from..first_index + rows.len() {
        let previous = integer(&rows[n - 1 - first_index]);
        let current = integer(&rows[n - first_index]);
        assert_eq!(
            check_weak_interlacing_bigint_coeffs(&previous, &current),
            Some(true),
            "Q_{} << Q_{n}",
            n - 1
        );
    }
}

fn proved(claim: &LiftInterlacingClaim) -> &LiftInterlacing {
    claim
        .proved()
        .unwrap_or_else(|| panic!("interlacing not proved: {claim:?}"))
}

fn not_proved_reason(claim: &LiftInterlacingClaim) -> &str {
    match claim {
        LiftInterlacingClaim::NotProved { reason } => reason,
        other => panic!("unexpected claim {other:?}"),
    }
}

/// The Eulerian rule (A008292) with its seeds; the first row is `A_1 = 1`.
fn eulerian_parent() -> (Recurrence, Vec<Vec<BigRational>>) {
    let (rule, seeds, first_row) = oeis::by_id("A008292")
        .unwrap()
        .recurrence_export_parts()
        .unwrap();
    assert_eq!(first_row, 1);
    (rule, seeds)
}

// ---------------------------------------------------------------------------
// Veronese sections
// ---------------------------------------------------------------------------

/// The monic order-3 rule of the even sections of the Eulerian polynomials
/// (see `veronese_section_lifts_v2.rs`).
fn eulerian_even_section_rule() -> Recurrence {
    rule(&[
        (1, 0, &[&[-1], &[1]]),
        (1, 1, &[&[0, 4]]),
        (2, 0, &[&[2, 4], &[-1, -4], &[0, 1]]),
        (2, 1, &[&[0, 0, 12], &[0, -4, -4]]),
        (2, 2, &[&[0, 0, -4, 4]]),
        (3, 0, &[&[0, 12], &[0, -19], &[0, 8], &[0, -1]]),
        (3, 1, &[&[0, -6, 18], &[0, 6, -22], &[0, 0, 4]]),
        (3, 2, &[&[0, 0, -4, 4], &[0, 0, 4, -4]]),
    ])
}

fn even_sections(
    parent_rule: &Recurrence,
    parent_seeds: &[Vec<BigRational>],
    count: usize,
) -> Vec<Vec<BigRational>> {
    parent_rule
        .generate_rows_rational(parent_seeds, 1, count)
        .unwrap()
        .iter()
        .map(|row| row.iter().step_by(2).cloned().collect())
        .collect()
}

fn certify_eulerian_even_sections(try_root_domain: bool) -> VeroneseSectionCertification {
    let (parent_rule, parent_seeds) = eulerian_parent();
    let parent = RecurrenceJson::from_recurrence_rational(&parent_rule, 1, &parent_seeds, None);
    let certification = certify_veronese_section_lift(
        &eulerian_even_section_rule(),
        &even_sections(&parent_rule, &parent_seeds, 3),
        1,
        &VeroneseSectionCandidate {
            modulus: 2,
            residue: 0,
            parent,
            index_shift: 0,
            gauge: VeroneseSectionGauge::None,
        },
        &VeroneseSectionCertifyOptions {
            try_root_domain,
            ..VeroneseSectionCertifyOptions::default()
        },
    )
    .unwrap();
    certification.verify().unwrap();
    let json = serde_json::to_string(&certification).unwrap();
    let back: VeroneseSectionCertification = serde_json::from_str(&json).unwrap();
    assert_eq!(back, certification);
    back.verify().unwrap();
    certification
}

/// The even sections of the Eulerian polynomials: with the Liu--Wang parent
/// (root-domain route off) the sections interlace, directed, by
/// Athanasiadis--Wagner Corollary 5.6(b); with the root-domain parent, which
/// concludes no interlacing, nothing is claimed.  The verdict is `proved`
/// either way.
#[test]
fn veronese_sections_inherit_the_parent_interlacing() {
    let certification = certify_eulerian_even_sections(false);
    assert_eq!(
        certification.parent_proof.route,
        Some(VeroneseParentProofRoute::LiuWang)
    );
    let parent = certification
        .parent_proof
        .interlacing
        .as_ref()
        .expect("the Liu--Wang interlacing");
    let summary = certification.lift_summary.as_ref().unwrap();
    assert_eq!(summary.verdict, LiftVerdict::Proved);
    let claim = proved(&summary.input_interlacing);
    assert_eq!(claim.kind, RowInterlacing::Weak);
    assert_eq!(
        claim.orientation,
        InterlacingOrientation::PreviousInterlacesNext
    );
    assert_eq!(claim.from_index, parent.from_index.max(2));
    assert!(claim.justification.theorem.contains("arXiv:2404.12989"));
    assert!(claim
        .justification
        .premises
        .iter()
        .any(|p| p.contains("liu_wang")));
    let (parent_rule, parent_seeds) = eulerian_parent();
    assert_rows_interlace(
        &even_sections(&parent_rule, &parent_seeds, 14),
        1,
        claim.from_index,
    );
    // The JSON form.
    let json = serde_json::to_value(summary).unwrap();
    assert_eq!(json["input_interlacing"]["status"], "proved");
    assert_eq!(
        json["input_interlacing"]["interlacing"]["orientation"],
        "previous_interlaces_next"
    );
    assert_eq!(json["input_interlacing"]["interlacing"]["kind"], "weak");

    // Tampering with the claim is caught by the replay.
    let mut tampered = certification.clone();
    if let Some(summary) = tampered.lift_summary.as_mut() {
        summary.input_interlacing = LiftInterlacingClaim::not_proved("x");
    }
    assert!(tampered.verify().is_err());

    let certification = certify_eulerian_even_sections(true);
    assert_eq!(
        certification.parent_proof.route,
        Some(VeroneseParentProofRoute::RootDomain)
    );
    assert!(certification.parent_proof.interlacing.is_none());
    let summary = certification.lift_summary.as_ref().unwrap();
    assert_eq!(summary.verdict, LiftVerdict::Proved);
    assert!(not_proved_reason(&summary.input_interlacing).contains("root_domain"));
}

/// A periodic-product parent concludes no interlacing: the v1 lifts and the
/// certifications on such parents claim nothing.
#[test]
fn periodic_product_parents_give_no_claim() {
    // Q_n = V_(2,0)((1 + x)^(2n)), Q_(n+2) = 2 Q_(n+1) + (y - 1) Q_n.
    let input = rule(&[(1, 0, &[&[2]]), (2, 0, &[&[-1, 1]])]);
    let report = find_veronese_section_lifts(
        &input,
        &[vec![q(1)], vec![q(1), q(1)]],
        0,
        &VeroneseSectionLiftOptions::default(),
    )
    .unwrap();
    assert!(!report.lifts.is_empty());
    for lift in &report.lifts {
        let summary = lift.lift_summary();
        assert!(!summary.input_interlacing.is_proved());
        assert!(not_proved_reason(&summary.input_interlacing).contains("periodic-product"));
    }
}

// ---------------------------------------------------------------------------
// Hadamard factors
// ---------------------------------------------------------------------------

fn fitted(rows: &[Vec<BigRational>]) -> (Recurrence, Vec<Vec<BigRational>>) {
    let fit = find_recurrence_adaptive_rational(rows, &AdaptiveSearchOptions::default())
        .expect("a recurrence");
    let seeds = rows[..fit.recurrence.max_offset()].to_vec();
    (fit.recurrence, seeds)
}

// The positive Hadamard rule (`1/k!`, `1/(alpha)_k`) is tested on the rule
// itself in the `hadamard_factor` report module; an end-to-end search for an
// interlacing parent is too slow for a debug test.

/// `binom(n, k) ⊙ (1 + x)^n` (type-B Narayana): the factor depends on `n`,
/// so no claim, although the lift is proved.
#[test]
fn an_index_dependent_hadamard_factor_gives_no_claim() {
    let family = oeis::by_id("A008459").unwrap();
    let (rule, rows) = family.recurrence_parts().unwrap();
    let report = find_hadamard_factor_lifts(
        &rule,
        &rows,
        family.recurrence_first_index,
        &HadamardFactorLiftOptions::default(),
    )
    .unwrap();
    let lift = report.proved_lifts().next().expect("a proved lift");
    assert_eq!(lift.factor_text, "binom(n, k)");
    assert!(not_proved_reason(&lift.lift_summary.input_interlacing).contains("depends on n"));
}

// ---------------------------------------------------------------------------
// Row factors
// ---------------------------------------------------------------------------

fn multiply(a: &[BigRational], b: &[BigRational]) -> Vec<BigRational> {
    let mut out = vec![q(0); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            out[i + j] += x * y;
        }
    }
    out
}

fn row_factor_report(
    input: &Recurrence,
    seeds: &[Vec<BigRational>],
    first_index: usize,
) -> RowFactorLiftReport {
    let report = find_row_factor_lifts(
        input,
        seeds,
        first_index,
        &RowFactorLiftOptions {
            try_root_domain: false,
            ..RowFactorLiftOptions::default()
        },
    )
    .unwrap();
    report.verify().unwrap();
    report
}

/// `Q_n = (2 + x) S_n` with `S_n = n! L_n(-x)` (Laguerre rows,
/// `S_n = (2n - 1 + x) S_(n-1) - (n - 1)^2 S_(n-2)`, which the three-term
/// theorem proves with strict interlacing): the constant exponent keeps the
/// directed interlacing, weakly.  `Q_n = (1 + x)^(n-1) A_n` (Eulerian `A_n`)
/// has an `n`-dependent exponent and gets no claim.
#[test]
fn row_factors_with_a_constant_exponent_keep_the_interlacing() {
    let laguerre = rule(&[(1, 0, &[&[-1, 1], &[2]]), (2, 0, &[&[-1], &[2], &[-1]])]);
    let factor = [q(2), q(1)];
    let seeds = [multiply(&[q(1)], &factor), multiply(&[q(1), q(1)], &factor)];
    let report = row_factor_report(&laguerre, &seeds, 0);
    let lift = report
        .lifts
        .iter()
        .find(|lift| lift.factor == ["2", "1"])
        .unwrap_or_else(|| panic!("the factor 2 + x: {report:#?}"));
    assert_eq!((lift.exponent_slope, lift.exponent_offset), (0, 1));
    assert!(lift.lift_summary.is_proved(), "{lift:#?}");
    assert_eq!(
        lift.parent_proof.route,
        Some(VeroneseParentProofRoute::ThreeTermRecurrence)
    );
    let claim = proved(&lift.lift_summary.input_interlacing);
    assert_eq!(claim.kind, RowInterlacing::Weak);
    assert_eq!(
        claim.orientation,
        InterlacingOrientation::PreviousInterlacesNext
    );
    assert!(claim.justification.theorem.contains("mul_common_factor"));
    let rows = laguerre.generate_rows_rational(&seeds, 0, 12).unwrap();
    assert_rows_interlace(&rows, 0, claim.from_index);
    let back: RowFactorLiftReport =
        serde_json::from_value(serde_json::to_value(&report).unwrap()).unwrap();
    assert_eq!(back, report);
    back.verify().unwrap();

    let (parent_rule, parent_seeds) = eulerian_parent();
    let mut power = vec![q(1)];
    let rows: Vec<_> = parent_rule
        .generate_rows_rational(&parent_seeds, 1, 14)
        .unwrap()
        .iter()
        .map(|a| {
            let out = multiply(a, &power);
            power = multiply(&power, &[q(1), q(1)]);
            out
        })
        .collect();
    let (input, seeds) = fitted(&rows);
    let report = row_factor_report(&input, &seeds, 1);
    let lift = report
        .lifts
        .iter()
        .find(|lift| lift.factor == ["1", "1"])
        .expect("the (1 + x) factor");
    assert_eq!(lift.exponent_slope, 1);
    assert!(not_proved_reason(&lift.lift_summary.input_interlacing).contains("depends on n"));
}

// ---------------------------------------------------------------------------
// Basis maps
// ---------------------------------------------------------------------------

/// The Hermite map `T(u^k) = H_k` (`H_(k+1) = 2t H_k - H_k'`) preserves
/// real-rootedness on the whole line and does not depend on `n`: the
/// interlacing of the Laguerre rows `R_n = n! L_n(-u)` (three-term theorem,
/// root-domain route off) transfers to `Q_n = T(R_n)`, weakly and unordered.
/// The power family `(u + 1)^n` (a periodic product) gives no claim.
#[test]
fn a_whole_line_basis_map_transfers_interlacing() {
    let hermite = BasisMap::new(
        "hermite",
        rule(&[(1, 0, &[&[0, 2]]), (1, 1, &[&[-1]])]),
        0,
        vec![vec![q(1)]],
    )
    .unwrap();
    let laguerre = rule(&[(1, 0, &[&[-1, 1], &[2]]), (2, 0, &[&[-1], &[2], &[-1]])]);
    let options = BasisMapImageOptions {
        certification: BasisMapLiftOptions {
            try_root_domain: false,
            ..BasisMapLiftOptions::default()
        },
        ..BasisMapImageOptions::default()
    };
    let family = RecurrenceJson::from_recurrence_rational(
        &laguerre,
        0,
        &[vec![q(1)], vec![q(1), q(1)]],
        None,
    );
    let image = apply_basis_map_to_family(&hermite, &family, &options).unwrap();
    let summary = &image.certification.lift_summary;
    assert!(image.is_proved(), "{summary:#?}");
    let claim = proved(&summary.input_interlacing);
    assert_eq!(claim.kind, RowInterlacing::Weak);
    assert_eq!(claim.orientation, InterlacingOrientation::Unordered);
    assert!(claim
        .justification
        .theorem
        .starts_with("Hermite--Kakeya--Obreschkoff"));
    image.certification.verify().unwrap();

    let power = rule(&[(1, 0, &[&[1, 1]])]);
    let family = RecurrenceJson::from_recurrence_rational(&power, 0, &[vec![q(1)]], None);
    let image = apply_basis_map_to_family(&hermite, &family, &options).unwrap();
    assert!(image.is_proved());
    assert!(
        not_proved_reason(&image.certification.lift_summary.input_interlacing)
            .contains("periodic_product")
    );
}

// ---------------------------------------------------------------------------
// The verdict does not depend on the claim
// ---------------------------------------------------------------------------

#[test]
fn the_claim_is_not_part_of_the_verdict() {
    let certification = certify_eulerian_even_sections(false);
    let summary = certification.lift_summary.clone().unwrap();
    let cleared = summary
        .clone()
        .with_input_interlacing(LiftInterlacingClaim::not_proved("cleared"));
    assert_eq!(cleared.verdict, summary.verdict);
    assert!(cleared.is_consistent());
    assert_eq!(
        lift_verdict(
            &cleared.identity,
            &cleared.parent_proof,
            &cleared.hypotheses
        ),
        LiftVerdict::Proved
    );
}
