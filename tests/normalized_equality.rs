//! Normalized all-index equality (`certify_normalized_recurrence_equality`):
//! each normalization of the menu by hand (index shift, dropped zero and
//! exceptional rows, a power of `x` with an affine exponent, a
//! hypergeometric scalar, all at once, a derivative rule), the catalog
//! cases found by the atlas (A257142, a forward Narayana image; A263789
//! from the peak polynomials; A348576 from a Möbius image), negative cases
//! (a scalar that vanishes, a wrong shift), replay and tampering, and the
//! transfer of a proved root region.

use num_bigint::BigInt;
use polytool::oeis;
use polytool::real_rootedness::RealInterval;
use polytool::recurrence::{
    apply_basis_map_to_family, certify_normalized_recurrence_equality,
    certify_recurrence_normalization, find_basis_map_lifts,
    transfer_real_rootedness_along_normalization, BasisMap, BasisMapImageOptions, BasisMapLibrary,
    BasisMapLibraryOptions, BasisMapLiftCertification, BasisMapLiftSearchOptions, BigRational,
    BivarPoly, NormalizationTransferDirection, NormalizedEqualityOptions,
    NormalizedRecurrenceEqualityCertificate, Recurrence, RecurrenceJson, RecurrenceNormalization,
    RecurrenceTerm, RootEndpoint, RootInterval, NORMALIZED_RECURRENCE_EQUALITY_SCHEMA,
};
use std::sync::Arc;

type Q = BigRational;

fn q(value: i64) -> Q {
    Q::from_integer(BigInt::from(value))
}

/// `coeffs[i][j]` is the coefficient of `n^i x^j`.
fn bivar(coeffs: &[&[i64]]) -> BivarPoly {
    BivarPoly::from_coefficients_normalized(
        coeffs
            .iter()
            .map(|row| row.iter().map(|&c| q(c)).collect())
            .collect(),
    )
}

fn rows(values: &[&[i64]]) -> Vec<Vec<Q>> {
    values
        .iter()
        .map(|row| row.iter().map(|&c| q(c)).collect())
        .collect()
}

/// `f(n, x) P_n = c(n, x) P_(n-1)` from the seeds.
fn first_order(
    coefficient: &[&[i64]],
    denominator: Option<&[&[i64]]>,
    first_index: usize,
    seeds: &[&[i64]],
) -> RecurrenceJson {
    let rule = Recurrence {
        terms: vec![RecurrenceTerm::new(1, 0, bivar(coefficient))],
        denominator: denominator.map(bivar),
        inhomogeneous: None,
    };
    RecurrenceJson::from_recurrence_rational(&rule, first_index, &rows(seeds), None)
}

/// `(1 + x)^n`, `n >= 0`.
fn binomial_rows() -> RecurrenceJson {
    first_order(&[&[1, 1]], None, 0, &[&[1]])
}

fn certify(
    source: &RecurrenceJson,
    target: &RecurrenceJson,
) -> NormalizedRecurrenceEqualityCertificate {
    let certificate = certify_normalized_recurrence_equality(
        source,
        target,
        &NormalizedEqualityOptions::default(),
    )
    .unwrap();
    certificate.verify().unwrap();
    certificate
}

fn certified(
    source: &RecurrenceJson,
    target: &RecurrenceJson,
) -> (
    NormalizedRecurrenceEqualityCertificate,
    RecurrenceNormalization,
) {
    let certificate = certify(source, target);
    assert!(certificate.is_certified(), "{certificate:#?}");
    let normalization = certificate.normalization.clone().unwrap();
    (certificate, normalization)
}

fn ratio(normalization: &RecurrenceNormalization) -> (Vec<&str>, Vec<&str>) {
    (
        normalization
            .scalar_ratio_numerator
            .iter()
            .map(String::as_str)
            .collect(),
        normalization
            .scalar_ratio_denominator
            .iter()
            .map(String::as_str)
            .collect(),
    )
}

#[test]
fn plain_equality_is_the_identity_normalization() {
    let (certificate, normalization) = certified(&binomial_rows(), &binomial_rows());
    assert_eq!(normalization, RecurrenceNormalization::identity());
    assert_eq!(certificate.identity_first_index, Some(0));
    assert!(certificate.dropped_target_rows.is_empty());
}

#[test]
fn an_index_shift_is_found() {
    // T_n = (1 + x)^(n + 2), n >= 0.
    let target = first_order(&[&[1, 1]], None, 0, &[&[1, 2, 1]]);
    let (_, normalization) = certified(&binomial_rows(), &target);
    assert_eq!(normalization.index_shift, 2);
    assert!(normalization.has_trivial_x_power());
    // And backwards: S_n = T_(n - 2) from n = 2 on, the first two rows of
    // the binomial family are exceptional.
    let (certificate, normalization) = certified(&target, &binomial_rows());
    assert_eq!(normalization.index_shift, -2);
    assert_eq!(normalization.dropped_leading_rows, 2);
    assert_eq!(
        certificate.dropped_target_rows,
        vec![
            vec!["1".to_string()],
            vec!["1".to_string(), "1".to_string()]
        ]
    );
}

#[test]
fn an_affine_power_of_x_is_found() {
    // T_n = x^(2n + 1) (1 + x)^n: T_n = x^2 (1 + x) T_(n-1), T_0 = x.
    let target = first_order(&[&[0, 0, 1, 1]], None, 0, &[&[0, 1]]);
    let (_, normalization) = certified(&binomial_rows(), &target);
    assert_eq!(
        (normalization.x_power_slope, normalization.x_power_offset),
        (2, 1)
    );
    // The reverse direction divides by x^(2n + 1) exactly.
    let (_, normalization) = certified(&target, &binomial_rows());
    assert_eq!(
        (normalization.x_power_slope, normalization.x_power_offset),
        (-2, -1)
    );
}

#[test]
fn a_hypergeometric_scalar_is_found() {
    // T_n = (-2)^n n! (1 + x)^n: T_n = -2 n (1 + x) T_(n-1), T_0 = 1.
    let target = first_order(&[&[0, 0], &[-2, -2]], None, 0, &[&[1]]);
    let (certificate, normalization) = certified(&binomial_rows(), &target);
    assert_eq!(normalization.scalar_first_value, "1");
    assert_eq!(ratio(&normalization), (vec!["0", "-2"], vec!["1"]));
    assert!(certificate.scalar_nonzero.as_ref().unwrap().holds());
}

#[test]
fn zero_and_exceptional_leading_rows_are_dropped() {
    // A zero first row, then (1 + x)^n from n = 1 on.
    let zero_first = first_order(&[&[1, 1]], None, 0, &[&[], &[1, 1]]);
    let (certificate, normalization) = certified(&binomial_rows(), &zero_first);
    assert_eq!(normalization.dropped_leading_rows, 1);
    assert_eq!(certificate.dropped_target_rows, vec![Vec::<String>::new()]);
    // An exceptional first row 5.
    let exceptional = first_order(&[&[1, 1]], None, 0, &[&[5], &[1, 1]]);
    let (certificate, normalization) = certified(&binomial_rows(), &exceptional);
    assert_eq!(normalization.dropped_leading_rows, 1);
    assert_eq!(certificate.dropped_target_rows, vec![vec!["5".to_string()]]);
}

/// `T_n = n x (1 + x)^(n - 1)`, `n >= 1`: `(n - 1) T_n = n (1 + x) T_(n-1)`
/// from `n = 2` on.
fn derivative_like_rows() -> RecurrenceJson {
    first_order(&[&[0, 0], &[1, 1]], Some(&[&[-1], &[1]]), 1, &[&[0, 1]])
}

#[test]
fn every_normalization_composes() {
    let (certificate, normalization) = certified(&binomial_rows(), &derivative_like_rows());
    assert_eq!(normalization.index_shift, -1);
    assert_eq!(normalization.dropped_leading_rows, 0);
    assert_eq!(
        (normalization.x_power_slope, normalization.x_power_offset),
        (0, 1)
    );
    assert_eq!(normalization.scalar_first_value, "1");
    assert_eq!(ratio(&normalization), (vec!["0", "1"], vec!["-1", "1"]));
    assert_eq!(certificate.identity_first_index, Some(1));
}

#[test]
fn derivative_rules_are_normalized() {
    // The Eulerian polynomials A_n = (1 + (n - 1) x) A_(n-1) + x (1 - x) A'_(n-1),
    // A_1 = 1, and T_n = (n + 1) x^n A_(n+1), n >= 0.
    let eulerian_rule = Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, bivar(&[&[1, -1], &[0, 1]])),
            RecurrenceTerm::new(1, 1, bivar(&[&[0, 1, -1]])),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let eulerian =
        RecurrenceJson::from_recurrence_rational(&eulerian_rule, 1, &rows(&[&[1]]), None);
    let source_rows = eulerian_rule
        .generate_rows_rational(&rows(&[&[1]]), 1, 4)
        .unwrap();
    // The target's own rule is the normalized one, written by hand would be
    // long: generate it from the certificate of a first candidate and use
    // the plain rows as its seeds instead.
    let normalization = RecurrenceNormalization {
        index_shift: 1,
        dropped_leading_rows: 0,
        x_power_slope: 1,
        x_power_offset: 0,
        scalar_first_value: "1".into(),
        scalar_ratio_numerator: vec!["1".into(), "1".into()],
        scalar_ratio_denominator: vec!["0".into(), "1".into()],
    };
    // Any target starting at index 0 fixes n1 = 0.
    let probe =
        certify_recurrence_normalization(&eulerian, &binomial_rows(), &normalization).unwrap();
    let normalized = probe.normalized_source.expect("the normalized rule");
    let (rule, first, seeds) = normalized.to_recurrence_parts().unwrap();
    assert_eq!(first, 0);
    // T_0 = 1 A_1, T_1 = 2 x A_2.
    assert_eq!(seeds[0], vec![q(1)]);
    let generated = rule.generate_rows_rational(&seeds, 0, 3).unwrap();
    let expected: Vec<Q> = std::iter::once(q(0))
        .chain(source_rows[1].iter().map(|c| c * q(2)))
        .collect();
    assert_eq!(generated[1], expected);
    let (certificate, found) = certified(&eulerian, &normalized);
    assert_eq!(found, normalization);
    assert!(certificate.equality.as_ref().unwrap().is_certified());
}

#[test]
fn a_vanishing_scalar_is_refused() {
    // T_n = (n - 3) (1 + x)^n has the zero row T_3: no candidate matches.
    let target = first_order(
        &[&[-3, -3], &[1, 1]],
        Some(&[&[-4], &[1]]),
        0,
        &[&[-3], &[-2, -2], &[-1, -2, -1], &[], &[1, 4, 6, 4, 1]],
    );
    let certificate = certify(&binomial_rows(), &target);
    assert!(!certificate.is_certified());
    assert!(certificate.normalization.is_none());
    // The recorded normalization c(n) = n - 3 is refused by the sign
    // decision: the ratio (n - 3) / (n - 4) vanishes at n = 3.
    let normalization = RecurrenceNormalization {
        index_shift: 0,
        dropped_leading_rows: 0,
        x_power_slope: 0,
        x_power_offset: 0,
        scalar_first_value: "-3".into(),
        scalar_ratio_numerator: vec!["-3".into(), "1".into()],
        scalar_ratio_denominator: vec!["-4".into(), "1".into()],
    };
    let refused =
        certify_recurrence_normalization(&binomial_rows(), &target, &normalization).unwrap();
    assert!(!refused.is_certified());
    assert!(
        refused.reasons[0].contains("not proved nonzero"),
        "{:?}",
        refused.reasons
    );
    assert_eq!(
        refused.scalar_nonzero.as_ref().unwrap().failing_index(),
        Some(&BigInt::from(3))
    );
    refused.verify().unwrap();
}

#[test]
fn a_wrong_shift_is_not_certified() {
    let target = first_order(&[&[1, 1]], None, 0, &[&[1, 2, 1]]);
    let wrong = RecurrenceNormalization {
        index_shift: 1,
        ..RecurrenceNormalization::identity()
    };
    let certificate = certify_recurrence_normalization(&binomial_rows(), &target, &wrong).unwrap();
    assert!(!certificate.is_certified());
    assert!(certificate.verify().is_ok());
    // Unrelated rows: nothing matches.
    let unrelated = first_order(&[&[2, 1]], None, 0, &[&[1]]);
    assert!(!certify(&binomial_rows(), &unrelated).is_certified());
}

#[test]
fn tampered_certificates_are_rejected() {
    let (certificate, _) = certified(&binomial_rows(), &derivative_like_rows());
    let mut tampered = certificate.clone();
    tampered.normalization.as_mut().unwrap().x_power_offset = 2;
    assert!(tampered.verify().is_err());
    let mut tampered = certificate.clone();
    tampered.normalization.as_mut().unwrap().scalar_first_value = "2".into();
    assert!(tampered.verify().is_err());
    let mut tampered = certificate.clone();
    tampered.statement.push('!');
    assert!(tampered.verify().is_err());
    let mut tampered = certificate.clone();
    tampered.target = binomial_rows();
    assert!(tampered.verify().is_err());
    let mut tampered = certificate.clone();
    tampered.schema = "polytool.recurrence-equality.v1".into();
    assert!(tampered.verify().is_err());
    // Claiming without a normalization.
    let mut tampered = certificate.clone();
    tampered.normalization = None;
    assert!(tampered.verify().is_err());
    // The JSON form round-trips.
    let text = serde_json::to_string(&certificate).unwrap();
    assert!(text.contains(NORMALIZED_RECURRENCE_EQUALITY_SCHEMA));
    let back: NormalizedRecurrenceEqualityCertificate = serde_json::from_str(&text).unwrap();
    assert_eq!(back, certificate);
}

fn closed(lower: i64, upper: i64) -> RootInterval {
    RootInterval::from_real(&RealInterval {
        lower: Some(q(lower)),
        lower_closed: true,
        upper: Some(q(upper)),
        upper_closed: true,
    })
}

#[test]
fn root_regions_transfer_along_the_normalization() {
    // (1 + x)^n has its roots in [-1, -1]; T_n = n x (1 + x)^(n - 1) gains
    // the root 0, so the region becomes [-1, 0] from T_1 on.
    let (certificate, _) = certified(&binomial_rows(), &derivative_like_rows());
    let transfer = transfer_real_rootedness_along_normalization(
        &certificate,
        NormalizationTransferDirection::SourceToTarget,
        0,
        &closed(-1, -1),
    )
    .unwrap();
    assert_eq!(transfer.conclusion_from_index, 1);
    assert_eq!(transfer.conclusion_region, closed(-1, 0));
    transfer.verify().unwrap();
    // Backwards the region also contains 0 (a power of x is removed).
    let transfer = transfer_real_rootedness_along_normalization(
        &certificate,
        NormalizationTransferDirection::TargetToSource,
        1,
        &closed(-1, 0),
    )
    .unwrap();
    assert_eq!(transfer.conclusion_from_index, 0);
    assert_eq!(transfer.conclusion_region, closed(-1, 0));
    transfer.verify().unwrap();
    // A pure shift keeps the region; the premise index moves.
    let shifted = first_order(&[&[1, 1]], None, 0, &[&[1, 2, 1]]);
    let (certificate, _) = certified(&binomial_rows(), &shifted);
    let transfer = transfer_real_rootedness_along_normalization(
        &certificate,
        NormalizationTransferDirection::SourceToTarget,
        5,
        &closed(-1, -1),
    )
    .unwrap();
    // T_0, T_1, T_2 (premise rows S_2, S_3, S_4 before 5) are checked.
    assert_eq!(transfer.conclusion_from_index, 0);
    assert_eq!(transfer.conclusion_region, closed(-1, -1));
    assert_eq!(transfer.checked_rows, vec![0, 1, 2]);
    // Tampering with the conclusion.
    let mut tampered = transfer.clone();
    tampered.conclusion_region = closed(-1, -1);
    tampered.conclusion_from_index = 0;
    tampered.premise_from_index = 9;
    assert!(tampered.verify().is_err());
    let mut tampered = transfer;
    tampered.conclusion_region = closed(-2, -1);
    assert!(tampered.verify().is_err());
}

#[test]
fn exceptional_rows_are_checked_exactly_in_a_transfer() {
    // An exceptional first row 1 + 3 x (root -1/3) outside [-1, -1]: the
    // conclusion starts after it; the row 5 has no root and is covered.
    let bad = first_order(&[&[1, 1]], None, 0, &[&[1, 3], &[1, 1]]);
    let (certificate, _) = certified(&binomial_rows(), &bad);
    let transfer = transfer_real_rootedness_along_normalization(
        &certificate,
        NormalizationTransferDirection::SourceToTarget,
        0,
        &closed(-1, -1),
    )
    .unwrap();
    assert_eq!(transfer.conclusion_from_index, 1);
    let good = first_order(&[&[1, 1]], None, 0, &[&[5], &[1, 1]]);
    let (certificate, _) = certified(&binomial_rows(), &good);
    let transfer = transfer_real_rootedness_along_normalization(
        &certificate,
        NormalizationTransferDirection::SourceToTarget,
        0,
        &closed(-1, -1),
    )
    .unwrap();
    assert_eq!(transfer.conclusion_from_index, 0);
    assert_eq!(transfer.checked_rows, vec![0]);
    transfer.verify().unwrap();
    // A region with an index-dependent end is refused.
    let mut moving = closed(-1, -1);
    moving.lower = Some(RootEndpoint {
        per_index: Some("1".into()),
        ..RootEndpoint::rational(&q(-1), true)
    });
    assert!(transfer_real_rootedness_along_normalization(
        &certificate,
        NormalizationTransferDirection::SourceToTarget,
        0,
        &moving,
    )
    .is_err());
}

// ---------------------------------------------------------------------------
// Catalog cases (the atlas definitions: a catalog entry at its recurrence's
// first index; the atlas nodes as recorded).
// ---------------------------------------------------------------------------

fn catalog_record(id: &str) -> RecurrenceJson {
    let entry = oeis::by_id(id).expect("catalog entry");
    let (rule, seeds) = entry.recurrence_parts().expect("catalog recurrence");
    RecurrenceJson::from_recurrence_rational(&rule, entry.recurrence_first_index, &seeds, None)
}

fn record(text: &str) -> RecurrenceJson {
    serde_json::from_str(text).expect("a recurrence record")
}

/// The atlas node `peak` (named family; first index 1).
const PEAK: &str = r#"{"schema": "polytool.recurrence.v1", "first_index": 1,
 "initial_polynomials": [["1"]], "recurrence": {"terms": [
 {"offset": 1, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["2", "-2"], ["0", "1"]]}},
 {"offset": 1, "deriv_order": 1, "sign": "none", "coeff": {"coeffs": [["0", "2", "-2"]]}}]}}"#;

/// The atlas node `mobius_x_over_x_plus_one of T(n,k) = k E(n,k)`: its
/// left-hand factor `(1 + x)^3` divides every coefficient.
const MOBIUS_OF_K_TIMES_EULERIAN: &str = r#"{"schema": "polytool.recurrence.v1",
 "first_index": 1, "initial_polynomials": [["1", "10", "12"], ["1", "25", "80", "60"]],
 "recurrence": {"terms": [
 {"offset": 2, "deriv_order": 0, "sign": "none", "coeff": {"coeffs":
   [["0", "0", "-2", "-6", "-6", "-2"], ["0", "0", "-1", "-3", "-3", "-1"]]}},
 {"offset": 2, "deriv_order": 1, "sign": "none", "coeff": {"coeffs":
   [["0", "0", "-2", "-8", "-12", "-8", "-2"], ["0", "0", "-1", "-4", "-6", "-4", "-1"]]}},
 {"offset": 1, "deriv_order": 0, "sign": "none", "coeff": {"coeffs":
   [["1", "6", "12", "10", "3"], ["0", "1", "3", "3", "1"]]}},
 {"offset": 1, "deriv_order": 1, "sign": "none", "coeff": {"coeffs":
   [["0", "1", "4", "6", "4", "1"]]}}],
 "denominator": {"coeffs": [["1", "3", "3", "1"]]}}}"#;

#[test]
fn a263789_is_a_normalized_peak_family() {
    let (certificate, normalization) = certified(&record(PEAK), &catalog_record("A263789"));
    // T_n = n x peak_(n - 1) for n >= 2; T_0 = T_1 = 1 are exceptional.
    assert_eq!(normalization.index_shift, -1, "{normalization:?}");
    assert_eq!(normalization.dropped_leading_rows, 2);
    assert_eq!(
        (normalization.x_power_slope, normalization.x_power_offset),
        (0, 1)
    );
    assert_eq!(normalization.scalar_first_value, "2");
    assert_eq!(ratio(&normalization), (vec!["0", "1"], vec!["-1", "1"]));
    // The atlas region (-inf, 0] of the peak family transfers unchanged.
    let ray = RootInterval::from_real(&RealInterval::at_most(q(0)));
    let transfer = transfer_real_rootedness_along_normalization(
        &certificate,
        NormalizationTransferDirection::SourceToTarget,
        1,
        &ray,
    )
    .unwrap();
    assert_eq!(transfer.conclusion_region, ray);
    assert_eq!(transfer.conclusion_from_index, 0);
    transfer.verify().unwrap();
}

#[test]
fn a348576_is_a_shifted_mobius_image() {
    let (certificate, normalization) = certified(
        &record(MOBIUS_OF_K_TIMES_EULERIAN),
        &catalog_record("A348576"),
    );
    assert_eq!(normalization.index_shift, -1, "{normalization:?}");
    assert!(normalization.has_trivial_x_power());
    assert_eq!(
        certificate.dropped_target_rows,
        vec![vec!["1".to_string(), "3".to_string()]]
    );
    let region = closed(-1, 0);
    let transfer = transfer_real_rootedness_along_normalization(
        &certificate,
        NormalizationTransferDirection::SourceToTarget,
        1,
        &region,
    )
    .unwrap();
    // The exceptional row 1 + 3 x (root -1/3) lies in [-1, 0].
    assert_eq!(transfer.conclusion_from_index, 1);
    assert_eq!(transfer.conclusion_region, region);
    transfer.verify().unwrap();
}

fn narayana_library() -> Arc<BasisMapLibrary> {
    Arc::new(
        BasisMapLibrary::build(&BasisMapLibraryOptions {
            maps: vec![BasisMap::narayana_transformation(1)],
            ..BasisMapLibraryOptions::default()
        })
        .unwrap(),
    )
}

#[test]
fn a257142_is_a_shifted_forward_narayana_image() {
    // The coefficient family of the backward lift, mapped forward by N_1:
    // the image starts at index 0, the catalog rule at index 1.
    let target = catalog_record("A257142");
    // The backward lift from the displayed indices (as in basis_map_published).
    let entry = oeis::by_id("A257142").unwrap();
    let (rule, seeds, first) = entry.recurrence_export_parts().unwrap();
    let first = usize::try_from(first).unwrap();
    let target_rows = rule.generate_rows_rational(&seeds, first, 16).unwrap();
    let lift = find_basis_map_lifts(
        &rule,
        &target_rows,
        first,
        &BasisMapLiftSearchOptions {
            library: Some(narayana_library()),
            ..BasisMapLiftSearchOptions::default()
        },
    )
    .unwrap()
    .lifts
    .into_iter()
    .find(BasisMapLiftCertification::is_proved)
    .expect("a proved backward lift");
    // The coefficient family as found (first index 1) and relabelled from 0:
    // T_n = N_1(R)_(n + 2) = N_1(R')_(n + 1).
    let family = lift.candidate.coefficient_family.clone();
    for (family, shift) in [(family.clone(), 2), (family.reindexed(0).unwrap(), 1)] {
        let image = apply_basis_map_to_family(
            &BasisMap::narayana_transformation(1),
            &family,
            &BasisMapImageOptions::default(),
        )
        .unwrap();
        assert!(image.is_proved());
        assert_eq!(image.image.first_index, family.first_index);
        let (_, normalization) = certified(&image.image, &target);
        assert_eq!(normalization.index_shift, shift, "{normalization:?}");
        assert!(normalization.has_trivial_x_power());
        assert_eq!(normalization.scalar_first_value, "1");
    }
}
