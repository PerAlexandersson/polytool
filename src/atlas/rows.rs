//! Canonical rows, stable hashes, fingerprints, observations and the exact
//! check of a claim on the stored rows.
//!
//! - Rows are ascending exact integer coefficient vectors with trailing zeros
//!   removed (the zero row is empty).
//! - [`atlas_hash`] is the 64-bit FNV-1a hash of a string, printed as 16
//!   hexadecimal digits; it is stable across platforms and versions (ids are
//!   derived from it).
//! - [`fingerprint_rows`] normalizes every row by its content, its sign, a
//!   power of `x` and reversal, drops leading rows of degree below two and
//!   hashes a window of at most [`FINGERPRINT_WINDOW`] rows: two sequences
//!   that differ by a row-wise scalar (for example a gauge `lambda(n)`), by
//!   `x -> 1/x` (with the row degree), by a power of `x` per row or by a
//!   small index shift get the same fingerprint.  A fingerprint is a search
//!   key, never an identification.

use super::types::{AtlasClaim, AtlasObservations, DegreePattern};
use crate::real_rootedness::{isolate_real_roots_rational_coeffs, RealInterval};
use crate::recurrence::{InterlacingOrientation, RootInterval};
use crate::{
    check_weak_interlacing_bigint_coeffs, interlaces_either_order_bigint_coeffs,
    is_gamma_positive_bigint_coeffs, is_palindromic_bigint_coeffs, BigRational,
};
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

/// Number of normalized rows hashed by [`fingerprint_rows`].
pub const FINGERPRINT_WINDOW: usize = 6;

/// The 64-bit FNV-1a hash of `text` as 16 lowercase hexadecimal digits.
pub fn atlas_hash(text: &str) -> String {
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    for byte in text.bytes() {
        hash ^= u64::from(byte);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    format!("{hash:016x}")
}

/// Exact integer rows (trailing zeros removed) from rational rows, or
/// `None` when some coefficient is not an integer.
pub fn integer_rows(rows: &[Vec<BigRational>]) -> Option<Vec<Vec<BigInt>>> {
    rows.iter()
        .map(|row| {
            let mut out = Vec::with_capacity(row.len());
            for c in row {
                if !c.is_integer() {
                    return None;
                }
                out.push(c.to_integer());
            }
            while out.last().is_some_and(Zero::is_zero) {
                out.pop();
            }
            Some(out)
        })
        .collect()
}

/// Rows as exact decimal strings.
pub(crate) fn rows_to_strings(rows: &[Vec<BigInt>]) -> Vec<Vec<String>> {
    rows.iter()
        .map(|row| row.iter().map(ToString::to_string).collect())
        .collect()
}

/// Rows from exact decimal strings.
pub(crate) fn rows_from_strings(rows: &[Vec<String>]) -> Result<Vec<Vec<BigInt>>, String> {
    rows.iter()
        .map(|row| {
            row.iter()
                .map(|c| {
                    c.parse::<BigInt>()
                        .map_err(|_| format!("`{c}` is not an integer"))
                })
                .collect()
        })
        .collect()
}

/// The row with its lowest zero coefficients removed.
pub(crate) fn strip_low_zeros(row: &[BigInt]) -> &[BigInt] {
    match row.iter().position(|c| !c.is_zero()) {
        Some(start) => &row[start..],
        None => &[],
    }
}

/// `row` divided by its content, with a positive leading coefficient.
fn primitive_positive(row: &[BigInt]) -> Vec<BigInt> {
    let mut content = BigInt::zero();
    for c in row {
        content = content.gcd(c);
    }
    if content.is_zero() {
        return Vec::new();
    }
    if row.last().is_some_and(Signed::is_negative) {
        content = -content;
    }
    row.iter().map(|c| c / &content).collect()
}

/// The normal form of one row under content, sign, powers of `x` and
/// reversal.
fn normalized_row(row: &[BigInt]) -> Vec<BigInt> {
    let stripped = strip_low_zeros(row);
    let forward = primitive_positive(stripped);
    let reversed: Vec<BigInt> = stripped.iter().rev().cloned().collect();
    let backward = primitive_positive(&reversed);
    forward.min(backward)
}

/// The fingerprint of a sequence of rows (module documentation).
pub fn fingerprint_rows(rows: &[Vec<BigInt>]) -> String {
    let normalized: Vec<Vec<BigInt>> = rows.iter().map(|row| normalized_row(row)).collect();
    let start = normalized
        .iter()
        .position(|row| row.len() >= 3)
        .or_else(|| normalized.iter().position(|row| !row.is_empty()))
        .unwrap_or(0);
    let window = &normalized[start..normalized.len().min(start + FINGERPRINT_WINDOW)];
    let text = window
        .iter()
        .map(|row| {
            row.iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(",")
        })
        .collect::<Vec<_>>()
        .join(";");
    format!("f{}-{}", window.len(), atlas_hash(&text))
}

/// The degree pattern of the nonzero rows `degrees[i] = deg P_(first_index + i)`
/// (`-1` for a zero row).
pub fn degree_pattern(first_index: usize, degrees: &[i64]) -> DegreePattern {
    let points: Vec<(i64, i64)> = degrees
        .iter()
        .enumerate()
        .filter(|(_, d)| **d >= 0)
        .map(|(i, d)| ((first_index + i) as i64, *d))
        .collect();
    if points.len() <= 1 {
        let intercept = points.first().map_or(0, |p| p.1);
        return DegreePattern::Affine {
            slope: 0,
            intercept,
        };
    }
    let (n0, d0) = points[0];
    let (n1, d1) = points[1];
    if (d1 - d0) % (n1 - n0) == 0 {
        let slope = (d1 - d0) / (n1 - n0);
        let intercept = d0 - slope * n0;
        if points.iter().all(|&(n, d)| d == slope * n + intercept) {
            return DegreePattern::Affine { slope, intercept };
        }
    }
    for slope in 1..=4i64 {
        for intercept in -2 * n0 - 8..=2 * n0 + 8 {
            if points
                .iter()
                .all(|&(n, d)| (slope * n + intercept).div_euclid(2) == d)
            {
                return DegreePattern::HalfAffine { slope, intercept };
            }
        }
    }
    DegreePattern::Irregular
}

/// The observations of [`AtlasObservations`] on the stored rows.
pub fn observe_rows(first_index: usize, rows: &[Vec<BigInt>]) -> AtlasObservations {
    let degrees: Vec<i64> = rows
        .iter()
        .map(|row| {
            row.iter()
                .rposition(|c| !c.is_zero())
                .map_or(-1, |d| d as i64)
        })
        .collect();
    let nonzero: Vec<&[BigInt]> = rows
        .iter()
        .map(|row| strip_low_zeros(row))
        .filter(|row| !row.is_empty())
        .collect();
    let palindromic = nonzero.iter().all(|row| is_palindromic_bigint_coeffs(row));
    let gamma_positive = palindromic
        && nonzero
            .iter()
            .all(|row| row.len() == 1 || is_gamma_positive_bigint_coeffs(row));
    AtlasObservations {
        rows_observed: rows.len(),
        degree_pattern: degree_pattern(first_index, &degrees),
        degrees,
        palindromic,
        gamma_positive,
        nonnegative_coefficients: rows.iter().flatten().all(|c| !c.is_negative()),
        row_sums: rows
            .iter()
            .map(|row| row.iter().fold(BigInt::zero(), |s, c| s + c).to_string())
            .collect(),
    }
}

/// The interval of an index-independent root region.
pub(crate) fn region_interval(region: &RootInterval) -> Result<RealInterval, String> {
    if !region.is_index_independent() {
        return Err(format!("the region {region} depends on the index"));
    }
    region.at(0).map_err(|e| e.to_string())
}

/// The ray `(-inf, 0]`.
pub(crate) fn nonpositive_ray() -> RootInterval {
    RootInterval::from_real(&RealInterval::at_most(BigRational::zero()))
}

/// The whole line.
pub(crate) fn whole_line() -> RootInterval {
    RootInterval {
        lower: None,
        upper: None,
    }
}

fn rational_row(row: &[BigInt]) -> Vec<BigRational> {
    row.iter().cloned().map(BigRational::from_integer).collect()
}

/// Check `claim` exactly on the stored rows `rows[i] = P_(first_index + i)`:
/// a real-rootedness claim on every covered row (zero, or real-rooted with
/// every root in the region, by exact root isolation), one-signed
/// coefficients, and the interlacing of covered consecutive nonzero rows
/// (weakly, in the claimed orientation).  Map claims and Hadamard kernels
/// are checked through the real-rootedness of the rows they rest on, so
/// they pass here.
///
/// # Errors
///
/// A description of the first stored row that contradicts the claim.
pub fn rows_satisfy_claim(
    first_index: usize,
    rows: &[Vec<BigInt>],
    claim: &AtlasClaim,
) -> Result<(), String> {
    let covered = |from: usize| {
        rows.iter()
            .enumerate()
            .filter(move |(i, _)| first_index + i >= from)
            .map(move |(i, row)| (first_index + i, row))
    };
    match claim {
        AtlasClaim::RealRooted { from_index, region } => {
            let interval = region_interval(region)?;
            for (n, row) in covered(*from_index) {
                if row.iter().all(Zero::is_zero) {
                    continue;
                }
                let isolation = isolate_real_roots_rational_coeffs(&rational_row(row));
                if !isolation.real_rooted_in(&interval) {
                    return Err(format!(
                        "row P_{n} is not real-rooted with every root in {region}"
                    ));
                }
            }
            Ok(())
        }
        AtlasClaim::OneSignedCoefficients {
            from_index,
            alternating,
        } => {
            for (n, row) in covered(*from_index) {
                let signs: Vec<bool> = row
                    .iter()
                    .enumerate()
                    .filter(|(_, c)| !c.is_zero())
                    .map(|(k, c)| c.is_positive() != (*alternating && k % 2 == 1))
                    .collect();
                if signs.windows(2).any(|w| w[0] != w[1]) {
                    return Err(format!("row P_{n} has coefficients of both signs"));
                }
            }
            Ok(())
        }
        AtlasClaim::Interlacing {
            from_index,
            orientation,
            ..
        } => {
            for (n, row) in covered(*from_index) {
                if n == first_index {
                    continue;
                }
                let previous = &rows[n - first_index - 1];
                if row.is_empty() || previous.is_empty() {
                    continue;
                }
                let holds = match orientation {
                    InterlacingOrientation::PreviousInterlacesNext => {
                        check_weak_interlacing_bigint_coeffs(previous, row) != Some(false)
                    }
                    InterlacingOrientation::NextInterlacesPrevious => {
                        check_weak_interlacing_bigint_coeffs(row, previous) != Some(false)
                    }
                    InterlacingOrientation::Unordered => {
                        previous.len() <= 1
                            || row.len() <= 1
                            || interlaces_either_order_bigint_coeffs(previous, row)
                    }
                };
                if !holds {
                    return Err(format!("rows P_{} and P_{n} do not interlace", n - 1));
                }
            }
            Ok(())
        }
        AtlasClaim::HadamardKernel { .. } | AtlasClaim::BasisMapPreserver { .. } => Ok(()),
    }
}

/// The intersection of two intervals.
pub(crate) fn intersect_intervals(a: &RealInterval, b: &RealInterval) -> RealInterval {
    let (lower, lower_closed) = match (&a.lower, &b.lower) {
        (None, None) => (None, false),
        (Some(x), None) => (Some(x.clone()), a.lower_closed),
        (None, Some(y)) => (Some(y.clone()), b.lower_closed),
        (Some(x), Some(y)) => match x.cmp(y) {
            std::cmp::Ordering::Greater => (Some(x.clone()), a.lower_closed),
            std::cmp::Ordering::Less => (Some(y.clone()), b.lower_closed),
            std::cmp::Ordering::Equal => (Some(x.clone()), a.lower_closed && b.lower_closed),
        },
    };
    let (upper, upper_closed) = match (&a.upper, &b.upper) {
        (None, None) => (None, false),
        (Some(x), None) => (Some(x.clone()), a.upper_closed),
        (None, Some(y)) => (Some(y.clone()), b.upper_closed),
        (Some(x), Some(y)) => match x.cmp(y) {
            std::cmp::Ordering::Less => (Some(x.clone()), a.upper_closed),
            std::cmp::Ordering::Greater => (Some(y.clone()), b.upper_closed),
            std::cmp::Ordering::Equal => (Some(x.clone()), a.upper_closed && b.upper_closed),
        },
    };
    RealInterval {
        lower,
        lower_closed,
        upper,
        upper_closed,
    }
}

/// `true` when `inner` is contained in `outer` (both index-independent).
pub(crate) fn region_within(inner: &RootInterval, outer: &RootInterval) -> bool {
    match (region_interval(inner), region_interval(outer)) {
        (Ok(a), Ok(b)) => a.is_empty() || a.is_subset_of(&b),
        _ => false,
    }
}

/// `1` as a rational (used by the row helpers of other modules).
pub(crate) fn one() -> BigRational {
    BigRational::one()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rows(values: &[&[i64]]) -> Vec<Vec<BigInt>> {
        values
            .iter()
            .map(|row| row.iter().map(|&c| BigInt::from(c)).collect())
            .collect()
    }

    #[test]
    fn hash_is_stable() {
        assert_eq!(atlas_hash(""), "cbf29ce484222325");
        assert_eq!(atlas_hash("a"), "af63dc4c8601ec8c");
    }

    #[test]
    fn fingerprint_ignores_content_sign_reversal_and_shift() {
        let eulerian = rows(&[
            &[1],
            &[1, 1],
            &[1, 4, 1],
            &[1, 11, 11, 1],
            &[1, 26, 66, 26, 1],
        ]);
        let scaled = rows(&[
            &[1, 1],
            &[-2, -8, -2],
            &[0, 3, 33, 33, 3],
            &[1, 26, 66, 26, 1],
        ]);
        assert_eq!(fingerprint_rows(&eulerian), fingerprint_rows(&scaled));
        let narayana = rows(&[&[1], &[1, 1], &[1, 3, 1], &[1, 6, 6, 1]]);
        assert_ne!(fingerprint_rows(&eulerian), fingerprint_rows(&narayana));
        let a = rows(&[&[1, 2, 3], &[1, 3, 3, 1]]);
        let b = rows(&[&[3, 2, 1], &[1, 3, 3, 1]]);
        assert_eq!(fingerprint_rows(&a), fingerprint_rows(&b));
    }

    #[test]
    fn intersections() {
        let q = |v: i64| BigRational::from_integer(v.into());
        let a = RealInterval::at_most(q(1));
        let b = RealInterval::at_least(q(-2));
        let i = intersect_intervals(&a, &b);
        assert!(i.is_subset_of(&RealInterval::closed(q(-2), q(1))));
        assert!(!i.is_subset_of(&RealInterval::closed(q(-1), q(1))));
    }

    #[test]
    fn degree_patterns() {
        assert_eq!(
            degree_pattern(1, &[0, 1, 2, 3]),
            DegreePattern::Affine {
                slope: 1,
                intercept: -1
            }
        );
        assert_eq!(
            degree_pattern(0, &[0, 0, 1, 1, 2, 2]),
            DegreePattern::HalfAffine {
                slope: 1,
                intercept: 0
            }
        );
    }

    #[test]
    fn claims_on_rows() {
        let eulerian = rows(&[&[1], &[1, 1], &[1, 4, 1], &[1, 11, 11, 1]]);
        let claim = AtlasClaim::RealRooted {
            from_index: 1,
            region: nonpositive_ray(),
        };
        assert!(rows_satisfy_claim(1, &eulerian, &claim).is_ok());
        let bad = rows(&[&[1, 1], &[1, 1, 1]]);
        assert!(rows_satisfy_claim(1, &bad, &claim).is_err());
        let interlacing = AtlasClaim::Interlacing {
            from_index: 2,
            kind: crate::recurrence::RowInterlacing::Weak,
            orientation: InterlacingOrientation::PreviousInterlacesNext,
        };
        assert!(rows_satisfy_claim(1, &eulerian, &interlacing).is_ok());
    }
}
