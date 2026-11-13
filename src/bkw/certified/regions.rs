//! Exact membership in versioned Schmidt-Spitzer and circle-predicate regions.
//!
//! Only the two-smallest-root definition determines the limit-spectrum region.
//! Circle nonnegativity is an independent algebraic predicate. An explicit
//! circle certificate applies the unpublished manuscript theorem; limit
//! geometry alone uses a separate conditional hypothesis. We use rational invariants
//! c=b_2^2/|b_1|^3 and d=b_3/|b_1|^2, so no irrational square root is rounded.

use super::*;
use std::collections::BTreeMap;
use std::sync::OnceLock;

/// Highest nontrivial tail bandwidth supported by the bundled region charts.
pub const MAXIMUM_PRECOMPUTED_TOEPLITZ_BANDWIDTH: usize = 3;

/// Explicit all-row assumption for an exactly real Hessenberg limit spectrum.
pub const LIMIT_SPECTRUM_HYPOTHESIS: &str = "hessenberg_real_limit_spectrum: for a real Hessenberg Toeplitz symbol b, Lambda(b) subset R implies that every standard row is real-rooted";

/// Exact evaluation of one versioned integer-polynomial Boolean expression.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct RegionEvaluation {
    /// Versioned chart identifier.
    pub id: String,
    /// Variable order of the polynomial exponent vectors.
    pub variables: Vec<String>,
    /// Exact rational point, in that same variable order.
    pub point: Vec<String>,
    /// Whether this region's predicate holds at the point.
    pub predicate: bool,
    /// SHA-256 of the exact quantifier-elimination input text.
    pub input_sha256: String,
    /// SHA-256 of the raw quantifier-elimination output artifact.
    pub output_sha256: String,
}

/// Replayable exact geometry, separate from any all-row conjecture.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ToeplitzRegionCertificate {
    /// `polytool.toeplitz-region-certificate.v1`.
    pub schema: String,
    /// Rational b_1,...,b_s, with trailing zero tails removed.
    pub tails: Vec<String>,
    /// Effective s in z^-1+b_0+sum_(j=1)^s b_j z^j.
    pub bandwidth: usize,
    /// Sign chart of b_1: -1, 0, or 1.
    pub first_tail_chart: i8,
    /// Exact normalization and the rational-invariant convention.
    pub normalization: String,
    /// True/false for a completed exact limit-region chart; None if unresolved.
    pub limit_real: Option<bool>,
    /// Exact off-real smallest-root-tie predicate; absent for closed-form cases.
    pub off_real_tie: Option<RegionEvaluation>,
    /// Independent exact existence predicate for circle nonnegativity.
    pub circle_condition: Option<RegionEvaluation>,
}

impl ToeplitzRegionCertificate {
    /// Recompute the normalized point and every region sign.
    ///
    /// # Errors
    /// Invalid rational tails or a changed certificate.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        if certify_toeplitz_regions(&self.tails)?.as_ref() != Some(self) {
            return Err(RecurrenceError::VerificationFailed(
                "Toeplitz region certificate does not replay".into(),
            ));
        }
        Ok(())
    }
}

#[derive(Deserialize)]
struct RegionData {
    entries: BTreeMap<String, RegionEntry>,
}
#[derive(Deserialize)]
struct RegionEntry {
    variables: Vec<String>,
    formula: Expression,
    input_sha256: String,
    output_sha256: String,
}
#[derive(Deserialize)]
#[serde(tag = "op")]
enum Expression {
    #[serde(rename = "true")]
    True,
    #[serde(rename = "false")]
    False,
    #[serde(rename = "and")]
    All { args: Vec<Expression> },
    #[serde(rename = "or")]
    Any { args: Vec<Expression> },
    #[serde(rename = "not")]
    Not { arg: Box<Expression> },
    #[serde(rename = "sign")]
    Sign {
        relation: Relation,
        polynomial: Vec<Monomial>,
    },
}
#[derive(Deserialize)]
enum Relation {
    Equal,
    Unequal,
    Greater,
    GreaterEqual,
    Less,
    LessEqual,
}
#[derive(Deserialize)]
struct Monomial {
    coefficient: String,
    powers: Vec<usize>,
}

fn data() -> &'static RegionData {
    static DATA: OnceLock<RegionData> = OnceLock::new();
    DATA.get_or_init(|| {
        serde_json::from_str(include_str!("../../../data/bkw/regions-v1.json"))
            .expect("validated bundled BKW regions")
    })
}
fn evaluate(expression: &Expression, point: &[Q]) -> bool {
    match expression {
        Expression::True => true,
        Expression::False => false,
        Expression::All { args } => args.iter().all(|e| evaluate(e, point)),
        Expression::Any { args } => args.iter().any(|e| evaluate(e, point)),
        Expression::Not { arg } => !evaluate(arg, point),
        Expression::Sign {
            relation,
            polynomial,
        } => {
            let value = polynomial.iter().fold(Q::zero(), |total, term| {
                assert_eq!(
                    term.powers.len(),
                    point.len(),
                    "bundled polynomial dimension"
                );
                let coefficient = term
                    .coefficient
                    .parse::<num_bigint::BigInt>()
                    .expect("bundled integer coefficient");
                total
                    + term
                        .powers
                        .iter()
                        .zip(point)
                        .fold(Q::from_integer(coefficient), |value, (&power, x)| {
                            value * x.pow(power as i32)
                        })
            });
            match relation {
                Relation::Equal => value.is_zero(),
                Relation::Unequal => !value.is_zero(),
                Relation::Greater => value.is_positive(),
                Relation::GreaterEqual => !value.is_negative(),
                Relation::Less => value.is_negative(),
                Relation::LessEqual => !value.is_positive(),
            }
        }
    }
}
fn membership(id: &str, c: &Q, d: &Q) -> Option<RegionEvaluation> {
    let entry = data().entries.get(id)?;
    let point: Vec<_> = entry
        .variables
        .iter()
        .map(|variable| match variable.as_str() {
            "c" => c.clone(),
            "d" => d.clone(),
            _ => unreachable!("bundled region variable"),
        })
        .collect();
    Some(RegionEvaluation {
        id: id.into(),
        variables: entry.variables.clone(),
        point: point.iter().map(ToString::to_string).collect(),
        predicate: evaluate(&entry.formula, &point),
        input_sha256: entry.input_sha256.clone(),
        output_sha256: entry.output_sha256.clone(),
    })
}

/// Certify limit-spectrum and circle-region membership for rational tails.
///
/// A real shift b_0 and an invertible real affine spectral coordinate do not
/// change reality of the limit set. For b_1!=0 use rho=1/sqrt(|b_1|) and the
/// symbol rho*(b(rho*z)-b_0). Reflection makes the limit and circle regions even
/// in b_2. Their coordinates c=b_2^2/|b_1|^3 and d=b_3/|b_1|^2 are rational.
/// If b_1=0 use rho=1, c=b_2^2 and d=b_3. Trailing zero tails reduce bandwidth.
/// Outside the bundled range this returns None, never a conjectural verdict.
///
/// # Errors
/// Invalid rational coefficient text.
pub fn certify_toeplitz_regions(
    tails: &[String],
) -> Result<Option<ToeplitzRegionCertificate>, RecurrenceError> {
    let mut b = decode(tails)?.coeffs().to_vec();
    while b.last().is_some_and(Zero::is_zero) {
        b.pop();
    }
    let bandwidth = b.len();
    if bandwidth > MAXIMUM_PRECOMPUTED_TOEPLITZ_BANDWIDTH {
        return Ok(None);
    }
    let first = b.first().cloned().unwrap_or_else(Q::zero);
    let chart = if first.is_positive() {
        1
    } else if first.is_negative() {
        -1
    } else {
        0
    };
    let scale = if first.is_zero() {
        Q::one()
    } else {
        first.abs()
    };
    let c = b.get(1).cloned().unwrap_or_else(Q::zero).pow(2) / scale.pow(3);
    let d = b.get(2).cloned().unwrap_or_else(Q::zero) / scale.pow(2);
    let suffix = match chart {
        1 => "plus",
        -1 => "minus",
        _ => "zero",
    };
    let (off_real_tie, circle_condition, limit_real) = if bandwidth <= 1 {
        (None, None, Some(bandwidth == 0 || first.is_positive()))
    } else {
        let mode = if bandwidth == 2 {
            "s2"
        } else if c.is_zero() {
            "s3zero"
        } else {
            "s3"
        };
        let off = membership(&format!("{mode}{suffix}"), &c, &d);
        let circle = membership(&format!("circle{bandwidth}{suffix}"), &c, &d);
        let real = off.as_ref().map(|e| !e.predicate);
        (off, circle, real)
    };
    Ok(Some(ToeplitzRegionCertificate {
        schema: "polytool.toeplitz-region-certificate.v1".into(),
        tails: b.iter().map(ToString::to_string).collect(),
        bandwidth,
        first_tail_chart: chart,
        normalization: if first.is_zero() {
            "rho=1; b_0 removed; c=b_2^2, d=b_3"
        } else {
            "rho=1/sqrt(abs(b_1)); rho*(b(rho*z)-b_0); c=b_2^2/abs(b_1)^3, d=b_3/abs(b_1)^2"
        }
        .into(),
        limit_real,
        off_real_tie,
        circle_condition,
    }))
}
