//! In-process finite checks and exact, replayable failure records.
use super::*;
use crate::{check_weak_interlacing_bigint_coeffs, is_real_rooted_bigint_coeffs, root_count};
use exact_algebra::denominators::clear_row_denominators;
use num_bigint::BigInt;

/// Finite property to check; first-failure stops at the first non-real-rooted row.
#[derive(Clone, Copy, Debug, Default, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum ExponentialRowCheck {
    /// Generate only.
    #[default]
    None,
    /// Check every selected row.
    RealRooted,
    /// Check rows and both directed weak interlacings of consecutive selected ranks.
    Interlacing,
    /// Stop at the first non-real-rooted selected row.
    FirstFailure,
}

/// Request shared by the library, CLI and MCP adapters.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ExponentialRowsRequest {
    /// Canonical labelled-block family.
    pub family: ExponentialFamily,
    /// Strictly increasing nonempty selected ranks. Gaps are allowed.
    pub rows: Vec<usize>,
    /// Finite check, default generation only.
    #[serde(default)]
    pub check: ExponentialRowCheck,
    /// Suppress passing coefficient rows when true; failures retain certificates.
    #[serde(default)]
    pub summary: bool,
    /// Optional per-rank reference certificates.
    #[serde(default)]
    pub finite_rank: Vec<FiniteRankRequest>,
}

/// One generated row and optional finite verdicts.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct ExponentialRowResult {
    /// Absolute rank n.
    pub index: usize,
    /// Ascending exact rational coefficients; absent in summary mode.
    pub coefficients: Option<Vec<String>>,
    /// Actual degree; null for the zero polynomial.
    pub degree: Option<usize>,
    /// Exact verdict; null when not requested.
    pub real_rooted: Option<bool>,
}
/// Both orientations of one adjacent-rank pair. Gaps are never silently checked.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct ExponentialInterlacingResult {
    /// Lower rank.
    pub left_index: usize,
    /// Upper rank (left + 1).
    pub right_index: usize,
    /// `F_left << F_right`; null for an inadmissible degree relation.
    pub forward: Option<bool>,
    /// `F_right << F_left`; null for an inadmissible degree relation.
    pub reverse: Option<bool>,
}
/// Exact failure certificate, replayed by the existing PRS/interlacing engine.
/// Counts are distinct real roots of the square-free part, not multiplicities.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct ExponentialFailureCertificate {
    /// Failed absolute rank (right rank for a pair).
    pub index: usize,
    /// Positive-denominator clearing of the failed row, ascending integers.
    pub integer_coefficients: Vec<String>,
    /// Square-free part, ascending primitive integer coefficients.
    pub squarefree_coefficients: Vec<String>,
    /// Original degree, null for zero.
    pub degree: Option<usize>,
    /// Degree of the square-free part.
    pub squarefree_degree: usize,
    /// Exact distinct real-root count, from signed PRS.
    pub distinct_real_roots: usize,
    /// Preceding row for an interlacing failure; otherwise null.
    pub previous_integer_coefficients: Option<Vec<String>>,
    /// Failure category: real_rooted or interlacing_either_orientation.
    pub property: String,
}
impl ExponentialFailureCertificate {
    /// Replay normalization, square-free part, root counts and any failed pair.
    pub fn verify(&self) -> Result<(), ExponentialRowError> {
        let p = parse_integers(&self.integer_coefficients)?;
        let q = self
            .previous_integer_coefficients
            .as_ref()
            .map(|r| parse_integers(r))
            .transpose()?;
        let regenerated = failure(self.index, &p, q.as_deref());
        let failed = if let Some(q) = q {
            check_weak_interlacing_bigint_coeffs(&q, &p) != Some(true)
                && check_weak_interlacing_bigint_coeffs(&p, &q) != Some(true)
        } else {
            !is_real_rooted_bigint_coeffs(&p)
        };
        if !failed || regenerated != *self {
            return Err(ExponentialRowError(
                "failure certificate does not replay".into(),
            ));
        }
        Ok(())
    }
}
fn parse_integers(row: &[String]) -> Result<Vec<BigInt>, ExponentialRowError> {
    row.iter()
        .map(|s| {
            s.parse()
                .map_err(|_| ExponentialRowError("invalid integer certificate".into()))
        })
        .collect()
}
fn failure(
    index: usize,
    p: &[BigInt],
    previous: Option<&[BigInt]>,
) -> ExponentialFailureCertificate {
    let sf = root_count::squarefree_part_primitive_bigint_coeffs(p);
    ExponentialFailureCertificate {
        index,
        integer_coefficients: p.iter().map(ToString::to_string).collect(),
        squarefree_degree: vec_poly::degree(&sf).unwrap_or(0),
        distinct_real_roots: root_count::count_real_roots_prs_bigint_coeffs(p),
        squarefree_coefficients: sf.iter().map(ToString::to_string).collect(),
        degree: vec_poly::degree(p),
        previous_integer_coefficients: previous
            .map(|r| r.iter().map(ToString::to_string).collect()),
        property: if previous.is_some() {
            "interlacing_either_orientation"
        } else {
            "real_rooted"
        }
        .into(),
    }
}

/// Finite report; only theorem-backed route entries make all-rank statements.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct ExponentialRowsReport {
    /// Stable schema identifier.
    pub schema: String,
    /// Embedded exact request, for reproduction.
    pub request: ExponentialRowsRequest,
    /// Selected rows actually visited, in increasing rank order.
    pub rows: Vec<ExponentialRowResult>,
    /// Checked adjacent-rank pairs, in both orientations.
    pub interlacing: Vec<ExponentialInterlacingResult>,
    /// First finite failure, with coefficients and exact counts for replay.
    pub first_failure: Option<ExponentialFailureCertificate>,
    /// Exact symbolic routes; never inferred from finite success.
    pub routes: Vec<ExponentialRouteReport>,
    /// Reference-symbol or pencil certificates for the requested ranks.
    pub finite_rank: Vec<FiniteRankCertificate>,
}

/// Generate/check selected rows in process. Dense ranges stream the insertion
/// iterator; sparse selections share scalar extraction tables. `first_failure`
/// refers only to the requested ranks, and stops without generating later rows
/// on a dense selection. Interlacing accepts either orientation and records both.
pub fn check_exponential_rows(
    request: &ExponentialRowsRequest,
) -> Result<ExponentialRowsReport, ExponentialRowError> {
    if request.rows.is_empty()
        || request.rows.windows(2).any(|w| w[0] >= w[1])
        || request.rows.last() == Some(&usize::MAX)
    {
        return Err(ExponentialRowError(
            "rows must be nonempty, strictly increasing, and below usize::MAX".into(),
        ));
    }
    if request
        .finite_rank
        .iter()
        .any(|r| !request.rows.contains(&r.index))
    {
        return Err(ExponentialRowError(
            "finite-rank reference index is not selected".into(),
        ));
    }
    request.family.unmarked.validate(0)?;
    request.family.marked.validate(0)?;
    let max_n = *request.rows.last().expect("validated rows");
    let dense =
        request.rows[0] <= request.rows.len() && request.rows.windows(2).all(|w| w[1] == w[0] + 1);
    let input: Box<dyn Iterator<Item = (usize, Vec<Q>)>> = if dense {
        let start = request.rows[0];
        Box::new(request.family.rows_iter(max_n)?.enumerate().skip(start))
    } else {
        Box::new(
            request
                .rows
                .clone()
                .into_iter()
                .zip(request.family.selected_rows(&request.rows)?),
        )
    };
    let mut report = ExponentialRowsReport {
        schema: "polytool.exp-rows.v1".into(),
        request: request.clone(),
        rows: vec![],
        interlacing: vec![],
        first_failure: None,
        routes: exponential_row_routes(&request.family)?,
        finite_rank: vec![],
    };
    let mut previous: Option<(usize, Vec<BigInt>)> = None;
    for (index, row) in input {
        let p = clear_row_denominators(&row).0;
        let rooted =
            (request.check != ExponentialRowCheck::None).then(|| is_real_rooted_bigint_coeffs(&p));
        if rooted == Some(false) && report.first_failure.is_none() {
            report.first_failure = Some(failure(index, &p, None));
        }
        if request.check == ExponentialRowCheck::Interlacing {
            if let Some((left, q)) = &previous {
                if *left + 1 == index {
                    let forward = check_weak_interlacing_bigint_coeffs(q, &p);
                    let reverse = check_weak_interlacing_bigint_coeffs(&p, q);
                    if forward != Some(true)
                        && reverse != Some(true)
                        && report.first_failure.is_none()
                    {
                        report.first_failure = Some(failure(index, &p, Some(q)));
                    }
                    report.interlacing.push(ExponentialInterlacingResult {
                        left_index: *left,
                        right_index: index,
                        forward,
                        reverse,
                    });
                }
            }
        }
        for reference in request.finite_rank.iter().filter(|r| r.index == index) {
            report
                .finite_rank
                .push(certify_finite_rank(&row, reference)?);
        }
        report.rows.push(ExponentialRowResult {
            index,
            degree: vec_poly::degree(&row),
            coefficients: (!request.summary).then(|| strings(&row)),
            real_rooted: rooted,
        });
        previous = Some((index, p));
        if rooted == Some(false) && request.check == ExponentialRowCheck::FirstFailure {
            break;
        }
    }
    Ok(report)
}
