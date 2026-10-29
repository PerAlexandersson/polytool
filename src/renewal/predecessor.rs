//! Symbolic variable-predecessor recurrences
//! `Q_0 = 1`, `Q_n = Q_(n-1) + w(n) t Q_(p(n))` for `n >= 1`.
//!
//! If `w(n) >= 0`, `0 <= p(n) < n` and `p` is nondecreasing, then every `Q_n`
//! has nonnegative coefficients and constant term 1, is real-rooted, and
//! `Q_i << Q_n` (oriented weak interlacing) for `p(n) <= i <= n`.  For the
//! fixed delay `p(n) = max(0, n-L)` this gives the window `Q_A << Q_B` for
//! `0 <= A <= B`, `B - A <= L`.  The paper proof is in
//! <https://github.com/PerAlexandersson/RealRooted/issues/2037>; the window
//! theorem (and its scaled form) is now formalized in the RealRooted Lean
//! library (`RealRooted/Interlacing/VariableDelay.lean`).  Individual matches
//! reported here are not replayed in Lean, which is what
//! [`PAPER_PROOF_STATUS`] records.
//!
//! [`VariablePredecessorRecurrence`] keeps the predecessor rule and the weight
//! separate from the fixed-offset [`crate::recurrence::RecurrenceTerm`].
//! Supported symbolic forms are [`PredecessorRule::FixedDelay`] together with
//! [`PredecessorWeight::Polynomial`] or [`PredecessorWeight::Periodic`].
//! Finite lists ([`PredecessorRule::Explicit`], [`PredecessorWeight::Explicit`])
//! certify only the supplied prefix; the certificate scope says so, and
//! [`VariablePredecessorRecurrence::certify_all_indices`] rejects them.
//!
//! Nonnegativity of a polynomial weight on the integers `n >= 1` is proven
//! exactly, never by sampling: the certificate stores an integer shift `c`
//! such that `w(c+m)` has nonnegative coefficients in `m` (so `w >= 0` on
//! `[c, oo)`), together with the exact values `w(1), ..., w(c-1)`.
//! [`VariablePredecessorRecurrence::verify_certificate`] replays it.
//!
//! [`VariablePredecessorRecurrence::rows`] is the canonical generator; rows
//! are `Vec<BigRational>` in ascending powers of `t`.
//! [`recognize_fixed_delay_predecessor`] recognizes
//! `P_n = P_(n-1) + w(n) t P_(n-L)` among fitted [`Recurrence`] values.

use super::RenewalError;
use crate::polynomial::Polynomial;
use crate::recurrence::{
    format_rational_coeff, parse_rational_coeff, BigRational, BivarPoly, Recurrence, RecurrenceSign,
};
use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Largest certificate shift `c` accepted for a polynomial weight; the values
/// `w(1), ..., w(c-1)` are checked one by one.  Larger shifts are rejected
/// conservatively ("not certified"), never accepted by sampling.
pub const MAX_POLYNOMIAL_WEIGHT_SHIFT: usize = 1 << 16;

/// Reference for the variable-predecessor and weighted renewal theorems.
pub const VARIABLE_PREDECESSOR_REFERENCE: &str =
    "https://github.com/PerAlexandersson/RealRooted/issues/2037";

/// Proof status reported with every certificate and theorem match.
pub const PAPER_PROOF_STATUS: &str = "paper_proof_not_lean_checked";

/// The predecessor function `p(n)`, defined for `n >= 1`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PredecessorRule {
    /// `p(n) = max(0, n - delay)`; symbolic, valid for every `n` when
    /// `delay >= 1`.
    FixedDelay {
        /// The delay `L`.
        delay: usize,
    },
    /// Finite list: `values[k] = p(k+1)`.  Certifies only `n <= values.len()`.
    Explicit {
        /// Predecessor values `p(1), p(2), ...`.
        values: Vec<usize>,
    },
}

/// The weight `w(n)`, defined for `n >= 1`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PredecessorWeight {
    /// `w(n) = sum_k coefficients[k] n^k` (ascending powers of `n`).
    Polynomial {
        /// Exact coefficients in ascending powers of `n`.
        coefficients: Vec<BigRational>,
    },
    /// `w(n) = values[(n-1) mod values.len()]`; the period must be positive.
    Periodic {
        /// One value per residue class, starting with `w(1)`.
        values: Vec<BigRational>,
    },
    /// Finite list: `values[k] = w(k+1)`.  Certifies only `n <= values.len()`.
    Explicit {
        /// Weight values `w(1), w(2), ...`.
        values: Vec<BigRational>,
    },
}

/// `Q_0 = 1`, `Q_n = Q_(n-1) + w(n) t Q_(p(n))`.
///
/// Construction does not validate the hypotheses; use
/// [`certify`](Self::certify).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VariablePredecessorRecurrence {
    /// The predecessor rule `p(n)`.
    pub predecessor: PredecessorRule,
    /// The weight `w(n)`.
    pub weight: PredecessorWeight,
}

/// Index range covered by a certificate.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum CertificateScope {
    /// Symbolic certificate: the hypotheses hold for every `n >= 1`.
    AllIndices,
    /// Finite data: the hypotheses were checked for `1 <= n <= last_index`
    /// only, so the conclusion is certified for `Q_0, ..., Q_last_index`.
    Prefix {
        /// Largest index `n` covered by the supplied finite data.
        last_index: usize,
    },
}

/// Exact certificate for the predecessor hypotheses `0 <= p(n) < n` and
/// monotonicity.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum PredecessorCertificate {
    /// `p(n) = max(0, n - delay)` with `delay >= 1`: symbolically
    /// nondecreasing with `0 <= p(n) <= n - 1`.
    FixedDelay {
        /// The delay `L >= 1`.
        delay: usize,
    },
    /// Every supplied value satisfies `p(n) < n`, and the list is
    /// nondecreasing.  Covers only the supplied prefix.
    ExplicitPrefix {
        /// Checked values `p(1), p(2), ...`.
        values: Vec<usize>,
    },
}

/// Exact certificate for `w(n) >= 0`.  Rationals are integer or `num/den`
/// strings.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum WeightCertificate {
    /// `w(n) = sum coefficients[k] n^k` is nonnegative on every integer
    /// `n >= 1`: `w(shift + m) = sum shifted_coefficients[k] m^k` has only
    /// nonnegative coefficients (so `w >= 0` on the real ray `[shift, oo)`),
    /// and `checked_values = [w(1), ..., w(shift-1)]` are all nonnegative.
    /// The zero weight has empty coefficient lists and `shift = 1`.
    Polynomial {
        /// Trimmed coefficients of `w` in ascending powers of `n`.
        coefficients: Vec<String>,
        /// The integer `c >= 1` of the Taylor-shift certificate.
        shift: usize,
        /// Coefficients of `w(c+m)` in ascending powers of `m`.
        shifted_coefficients: Vec<String>,
        /// Exact values `w(1), ..., w(c-1)`.
        checked_values: Vec<String>,
    },
    /// Every residue value of a periodic weight is nonnegative.
    Periodic {
        /// `w(1), ..., w(period)`.
        values: Vec<String>,
    },
    /// Every supplied value is nonnegative.  Covers only the supplied prefix.
    ExplicitPrefix {
        /// Checked values `w(1), w(2), ...`.
        values: Vec<String>,
    },
}

/// Exact certificate of the variable-predecessor hypotheses and the resulting
/// conditional conclusion.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct VariablePredecessorCertificate {
    /// Indices covered: all `n` for symbolic data, a prefix for finite lists.
    pub scope: CertificateScope,
    /// Predecessor hypotheses.
    pub predecessor: PredecessorCertificate,
    /// Weight nonnegativity.
    pub weight: WeightCertificate,
    /// `L` for a fixed delay: `Q_A << Q_B` whenever `0 <= A <= B` and
    /// `B - A <= L` (within the scope).  `None` for explicit predecessors.
    pub interlacing_window: Option<usize>,
    /// Plain-text statement of the certified conclusion.
    pub conclusion: String,
    /// Always `paper_proof_not_lean_checked`.
    pub proof_status: String,
    /// Source of the paper proof.
    pub reference: String,
}

fn trim(mut coefficients: Vec<BigRational>) -> Vec<BigRational> {
    while coefficients.last().is_some_and(Zero::is_zero) {
        coefficients.pop();
    }
    coefficients
}

fn rational(n: usize) -> BigRational {
    BigRational::from_integer(BigInt::from(n))
}

fn strings(values: &[BigRational]) -> Vec<String> {
    values.iter().map(format_rational_coeff).collect()
}

fn taylor_shift(coefficients: &[BigRational], shift: usize) -> Vec<BigRational> {
    Polynomial::new(coefficients.to_vec())
        .shift(&rational(shift))
        .coeffs()
        .to_vec()
}

fn evaluate(coefficients: &[BigRational], n: usize) -> BigRational {
    crate::vec_poly::evaluate(coefficients, &rational(n))
}

/// Find the certificate shift for a nonzero polynomial weight with positive
/// leading coefficient, or explain why none was found.
fn polynomial_weight_shift(coefficients: &[BigRational]) -> Result<usize, RenewalError> {
    let ok = |c: usize| {
        taylor_shift(coefficients, c)
            .iter()
            .all(|x| !x.is_negative())
    };
    let mut hi = 1usize;
    while !ok(hi) {
        if hi > MAX_POLYNOMIAL_WEIGHT_SHIFT {
            return Err(RenewalError::WeightShiftNotFound {
                max_shift: MAX_POLYNOMIAL_WEIGHT_SHIFT,
            });
        }
        hi *= 2;
    }
    // ok is monotone in c: substituting m -> m + 1 keeps coefficients nonnegative.
    let mut lo = hi / 2; // not ok (or zero)
    while hi - lo > 1 {
        let mid = lo + (hi - lo) / 2;
        if ok(mid) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    if hi > MAX_POLYNOMIAL_WEIGHT_SHIFT {
        return Err(RenewalError::WeightShiftTooLarge {
            shift: hi,
            max_shift: MAX_POLYNOMIAL_WEIGHT_SHIFT,
        });
    }
    Ok(hi)
}

fn certify_polynomial_weight(
    coefficients: &[BigRational],
) -> Result<WeightCertificate, RenewalError> {
    let coefficients = trim(coefficients.to_vec());
    let Some(leading) = coefficients.last() else {
        return Ok(WeightCertificate::Polynomial {
            coefficients: Vec::new(),
            shift: 1,
            shifted_coefficients: Vec::new(),
            checked_values: Vec::new(),
        });
    };
    if leading.is_negative() {
        // Cauchy bound: every real zero is smaller than 1 + max |a_i / a_d|,
        // so w has the sign of its leading coefficient from there on.
        let bound = coefficients
            .iter()
            .map(|c| (c / leading).abs())
            .fold(BigRational::zero(), |a, b| a.max(b))
            + BigRational::one();
        let witness = bound.ceil().to_integer().max(BigInt::one());
        let value =
            crate::vec_poly::evaluate(&coefficients, &BigRational::from_integer(witness.clone()));
        debug_assert!(value.is_negative());
        return Err(RenewalError::NegativeLeadingCoefficient {
            index: witness,
            value: format_rational_coeff(&value),
        });
    }
    let shift = polynomial_weight_shift(&coefficients)?;
    let mut checked_values = Vec::with_capacity(shift - 1);
    for n in 1..shift {
        let value = evaluate(&coefficients, n);
        if value.is_negative() {
            return Err(RenewalError::NegativeWeight {
                kind: None,
                index: n,
                value: format_rational_coeff(&value),
            });
        }
        checked_values.push(value);
    }
    Ok(WeightCertificate::Polynomial {
        shifted_coefficients: strings(&taylor_shift(&coefficients, shift)),
        coefficients: strings(&coefficients),
        shift,
        checked_values: strings(&checked_values),
    })
}

fn certify_values(values: &[BigRational], kind: &'static str) -> Result<(), RenewalError> {
    match values.iter().position(Signed::is_negative) {
        Some(k) => Err(RenewalError::NegativeWeight {
            kind: Some(kind),
            index: k + 1,
            value: format_rational_coeff(&values[k]),
        }),
        None => Ok(()),
    }
}

impl VariablePredecessorRecurrence {
    /// The fixed-delay family `p(n) = max(0, n - delay)` with weight `weight`.
    pub fn fixed_delay(delay: usize, weight: PredecessorWeight) -> Self {
        Self {
            predecessor: PredecessorRule::FixedDelay { delay },
            weight,
        }
    }

    /// Largest index `n` for which `p(n)` and `w(n)` are defined; `None` when
    /// both are symbolic (every `n`).
    pub fn defined_through(&self) -> Option<usize> {
        let p = match &self.predecessor {
            PredecessorRule::FixedDelay { .. } => None,
            PredecessorRule::Explicit { values } => Some(values.len()),
        };
        let w = match &self.weight {
            PredecessorWeight::Explicit { values } => Some(values.len()),
            _ => None,
        };
        match (p, w) {
            (Some(a), Some(b)) => Some(a.min(b)),
            (a, b) => a.or(b),
        }
    }

    /// `p(n)` for `n >= 1`; `None` for `n = 0`, outside explicit data, or an
    /// empty period.
    pub fn predecessor_at(&self, n: usize) -> Option<usize> {
        if n == 0 {
            return None;
        }
        match &self.predecessor {
            PredecessorRule::FixedDelay { delay } => Some(n.saturating_sub(*delay)),
            PredecessorRule::Explicit { values } => values.get(n - 1).copied(),
        }
    }

    /// `w(n)` for `n >= 1`; `None` for `n = 0`, outside explicit data, or an
    /// empty period.
    pub fn weight_at(&self, n: usize) -> Option<BigRational> {
        if n == 0 {
            return None;
        }
        match &self.weight {
            PredecessorWeight::Polynomial { coefficients } => Some(evaluate(coefficients, n)),
            PredecessorWeight::Periodic { values } => {
                (!values.is_empty()).then(|| values[(n - 1) % values.len()].clone())
            }
            PredecessorWeight::Explicit { values } => values.get(n - 1).cloned(),
        }
    }

    /// Certify the hypotheses exactly.  Errors name a violating index or
    /// explain why no exact certificate was found; failure is never decided by
    /// sampling.
    pub fn certify(&self) -> Result<VariablePredecessorCertificate, RenewalError> {
        let predecessor = match &self.predecessor {
            PredecessorRule::FixedDelay { delay: 0 } => return Err(RenewalError::ZeroDelay),
            PredecessorRule::FixedDelay { delay } => {
                PredecessorCertificate::FixedDelay { delay: *delay }
            }
            PredecessorRule::Explicit { values } => {
                for (k, &p) in values.iter().enumerate() {
                    let n = k + 1;
                    if p >= n {
                        return Err(RenewalError::PredecessorNotBelowIndex {
                            index: n,
                            predecessor: p,
                        });
                    }
                    if k > 0 && p < values[k - 1] {
                        return Err(RenewalError::PredecessorDecreases {
                            index: n,
                            previous: values[k - 1],
                            predecessor: p,
                        });
                    }
                }
                PredecessorCertificate::ExplicitPrefix {
                    values: values.clone(),
                }
            }
        };
        let weight = match &self.weight {
            PredecessorWeight::Polynomial { coefficients } => {
                certify_polynomial_weight(coefficients)?
            }
            PredecessorWeight::Periodic { values } => {
                if values.is_empty() {
                    return Err(RenewalError::EmptyPeriod);
                }
                certify_values(values, "periodic")?;
                WeightCertificate::Periodic {
                    values: strings(values),
                }
            }
            PredecessorWeight::Explicit { values } => {
                certify_values(values, "explicit")?;
                WeightCertificate::ExplicitPrefix {
                    values: strings(values),
                }
            }
        };
        let scope = match self.defined_through() {
            None => CertificateScope::AllIndices,
            Some(last_index) => CertificateScope::Prefix { last_index },
        };
        let interlacing_window = match self.predecessor {
            PredecessorRule::FixedDelay { delay } => Some(delay),
            PredecessorRule::Explicit { .. } => None,
        };
        let range = match &scope {
            CertificateScope::AllIndices => "every n >= 1".to_string(),
            CertificateScope::Prefix { last_index } => {
                format!("1 <= n <= {last_index} only (finite data)")
            }
        };
        let mut conclusion = format!(
            "For {range}: Q_n has nonnegative coefficients and constant term 1, is \
             real-rooted, and Q_i << Q_n (oriented weak interlacing) for p(n) <= i <= n."
        );
        if let Some(delay) = interlacing_window {
            let upper = match &scope {
                CertificateScope::AllIndices => String::new(),
                CertificateScope::Prefix { last_index } => format!(" <= {last_index}"),
            };
            conclusion.push_str(&format!(
                " Hence Q_A << Q_B whenever 0 <= A <= B{upper} and B - A <= {delay}."
            ));
        }
        Ok(VariablePredecessorCertificate {
            scope,
            predecessor,
            weight,
            interlacing_window,
            conclusion,
            proof_status: PAPER_PROOF_STATUS.into(),
            reference: VARIABLE_PREDECESSOR_REFERENCE.into(),
        })
    }

    /// Certify the hypotheses for every `n >= 1`.  Finite predecessor or
    /// weight lists are rejected: they certify only a prefix.
    pub fn certify_all_indices(&self) -> Result<VariablePredecessorCertificate, RenewalError> {
        let certificate = self.certify()?;
        match certificate.scope {
            CertificateScope::AllIndices => Ok(certificate),
            CertificateScope::Prefix { last_index } => {
                Err(RenewalError::FiniteDataOnly { last_index })
            }
        }
    }

    /// Replay a certificate against this recurrence without trusting it.
    ///
    /// For a polynomial weight any valid shift is accepted: the stored
    /// coefficients must equal the weight, the Taylor shift is recomputed and
    /// all stored values are re-evaluated.  The scope, window and other
    /// components must agree with this recurrence.
    pub fn verify_certificate(
        &self,
        certificate: &VariablePredecessorCertificate,
    ) -> Result<(), RenewalError> {
        let parse = |values: &[String]| -> Result<Vec<BigRational>, RenewalError> {
            values
                .iter()
                .map(|v| parse_rational_coeff(v).map_err(RenewalError::from))
                .collect()
        };
        let derived = self.certify()?;
        if certificate.scope != derived.scope
            || certificate.predecessor != derived.predecessor
            || certificate.interlacing_window != derived.interlacing_window
        {
            return Err(RenewalError::CertificateScopeMismatch);
        }
        match (&self.weight, &certificate.weight) {
            (
                PredecessorWeight::Polynomial { coefficients },
                WeightCertificate::Polynomial {
                    coefficients: stored,
                    shift,
                    shifted_coefficients,
                    checked_values,
                },
            ) => {
                let coefficients = trim(coefficients.clone());
                if parse(stored)? != coefficients || *shift == 0 {
                    return Err(RenewalError::CertificateWeightPolynomialMismatch);
                }
                let shifted = parse(shifted_coefficients)?;
                if shifted != taylor_shift(&coefficients, *shift)
                    || shifted.iter().any(Signed::is_negative)
                {
                    return Err(RenewalError::CertificateShiftedCoefficientsInvalid);
                }
                let values = parse(checked_values)?;
                if values.len() != shift - 1
                    || values.iter().any(Signed::is_negative)
                    || values
                        .iter()
                        .enumerate()
                        .any(|(k, v)| *v != evaluate(&coefficients, k + 1))
                {
                    return Err(RenewalError::CertificateCheckedValuesInvalid);
                }
                Ok(())
            }
            (_, stored) if *stored == derived.weight => Ok(()),
            _ => Err(RenewalError::CertificateWeightMismatch),
        }
    }

    /// `L` when `p(n) = max(0, n-L)`: then `Q_A << Q_B` for `B - A <= L`.
    pub fn interlacing_window(&self) -> Option<usize> {
        match self.predecessor {
            PredecessorRule::FixedDelay { delay } => Some(delay),
            PredecessorRule::Explicit { .. } => None,
        }
    }

    /// Canonical rows `Q_0, ..., Q_(total_rows-1)` (ascending powers of `t`).
    /// This is the library generator for the family.  Errors when an index is
    /// outside the supplied finite data or when `p(n) >= n`; the weight sign is
    /// not checked here (see [`certify`](Self::certify)).
    pub fn rows(&self, total_rows: usize) -> Result<Vec<Vec<BigRational>>, RenewalError> {
        let mut rows: Vec<Vec<BigRational>> = Vec::with_capacity(total_rows);
        for n in 0..total_rows {
            if n == 0 {
                rows.push(vec![BigRational::one()]);
                continue;
            }
            let (Some(p), Some(w)) = (self.predecessor_at(n), self.weight_at(n)) else {
                return Err(RenewalError::UndefinedIndex { index: n });
            };
            if p >= n {
                return Err(RenewalError::PredecessorNotBelowIndex {
                    index: n,
                    predecessor: p,
                });
            }
            let mut row = rows[n - 1].clone();
            let previous = &rows[p];
            row.resize(row.len().max(previous.len() + 1), BigRational::zero());
            for (k, c) in previous.iter().enumerate() {
                row[k + 1] += &w * c;
            }
            while row.len() > 1 && row.last().is_some_and(Zero::is_zero) {
                row.pop();
            }
            rows.push(row);
        }
        Ok(rows)
    }

    /// The same family after relabeling `n -> n + shift`: weight
    /// `w(n + shift)` with the same fixed delay.  Only fixed-delay rules are
    /// supported.  Nonnegativity for `n >= 1` is inherited.
    pub(crate) fn with_index_shift(&self, shift: usize) -> Option<Self> {
        let PredecessorRule::FixedDelay { delay } = self.predecessor else {
            return None;
        };
        let weight = match &self.weight {
            PredecessorWeight::Polynomial { coefficients } => PredecessorWeight::Polynomial {
                coefficients: taylor_shift(&trim(coefficients.clone()), shift),
            },
            PredecessorWeight::Periodic { values } if !values.is_empty() => {
                let mut values = values.clone();
                let len = values.len();
                values.rotate_left(shift % len);
                PredecessorWeight::Periodic { values }
            }
            PredecessorWeight::Explicit { values } => PredecessorWeight::Explicit {
                values: values.iter().skip(shift).cloned().collect(),
            },
            PredecessorWeight::Periodic { .. } => return None,
        };
        Some(Self::fixed_delay(delay, weight))
    }

    /// Plain-text description, e.g.
    /// `Q_0 = 1, Q_n = Q_(n-1) + w(n) t Q_(max(0, n-2)), w(n) = 1 + n`.
    pub fn describe(&self) -> String {
        let p = match &self.predecessor {
            PredecessorRule::FixedDelay { delay } => format!("max(0, n-{delay})"),
            PredecessorRule::Explicit { values } => {
                format!("p(n), p(1..={}) = {values:?}", values.len())
            }
        };
        let w = match &self.weight {
            PredecessorWeight::Polynomial { coefficients } => {
                let coeffs = trim(coefficients.clone())
                    .into_iter()
                    .map(|c| vec![c])
                    .collect::<Vec<_>>();
                if coeffs.is_empty() {
                    "0".into()
                } else {
                    BivarPoly { coeffs }.to_string()
                }
            }
            PredecessorWeight::Periodic { values } => format!(
                "w(1..={}) = [{}] repeated periodically",
                values.len(),
                strings(values).join(", ")
            ),
            PredecessorWeight::Explicit { values } => format!(
                "w(1..={}) = [{}] (finite list)",
                values.len(),
                strings(values).join(", ")
            ),
        };
        format!("Q_0 = 1, Q_n = Q_(n-1) + w(n) t Q_({p}), w(n) = {w}")
    }
}

/// Recognize `P_n = P_(n-1) + w(n) t P_(n-L)` with `L >= 1` and `w`
/// polynomial in `n`, as a fixed-delay [`VariablePredecessorRecurrence`].
///
/// A constant nonzero LHS factor and split/duplicate terms are normalized;
/// derivatives, alternating signs, inhomogeneous terms, `t`-powers other than
/// 0 and 1, more than one weighted lag, and a coefficient of `P_(n-1)` other
/// than `1 (+ w(n) t if L = 1)` are rejected.  The weight sign is not checked
/// here (see [`VariablePredecessorRecurrence::certify`]).
///
/// The canonical family uses `Q_(max(0, n-L))`, so for `n < L` it refers to
/// `Q_0 = 1`, whereas a fixed-offset recurrence uses zero rows at negative
/// indices; compare initial rows before applying the theorem.
pub fn recognize_fixed_delay_predecessor(
    recurrence: &Recurrence,
) -> Result<VariablePredecessorRecurrence, RenewalError> {
    let scalar = |p: &BivarPoly| -> Option<BigRational> {
        let mut value = BigRational::zero();
        for (i, row) in p.coeffs.iter().enumerate() {
            for (j, c) in row.iter().enumerate() {
                if i == 0 && j == 0 {
                    value = c.clone();
                } else if !c.is_zero() {
                    return None;
                }
            }
        }
        Some(value)
    };
    let denominator = match &recurrence.denominator {
        None => BigRational::one(),
        Some(p) => scalar(p)
            .filter(|c| !c.is_zero())
            .ok_or(RenewalError::NonconstantLeftFactor)?,
    };
    if recurrence
        .inhomogeneous
        .as_ref()
        .is_some_and(|p| scalar(p) != Some(BigRational::zero()))
    {
        return Err(RenewalError::NotHomogeneous);
    }
    // (offset, n-power, t-power) -> coefficient
    let mut terms = BTreeMap::<(usize, usize, usize), BigRational>::new();
    for term in &recurrence.terms {
        if term.deriv_order != 0 || term.sign != RecurrenceSign::None {
            return Err(RenewalError::UnsupportedTermKind);
        }
        for (i, row) in term.coeff.coeffs.iter().enumerate() {
            for (j, c) in row.iter().enumerate() {
                if !c.is_zero() {
                    *terms.entry((term.offset, i, j)).or_default() += c / &denominator;
                }
            }
        }
    }
    terms.retain(|_, c| !c.is_zero());
    if terms.remove(&(1, 0, 0)) != Some(BigRational::one()) {
        return Err(RenewalError::PreviousCoefficientNotOneWithWeight);
    }
    let mut delay = None;
    let mut weight = Vec::new();
    for ((offset, i, j), c) in terms {
        if j != 1 || offset == 0 {
            return Err(RenewalError::UnsupportedWeightedTerm);
        }
        if delay.is_some_and(|d| d != offset) {
            return Err(RenewalError::MultipleWeightedLags);
        }
        delay = Some(offset);
        weight.resize(weight.len().max(i + 1), BigRational::zero());
        weight[i] = c;
    }
    let delay = delay.ok_or(RenewalError::NoWeightedTerm)?;
    Ok(VariablePredecessorRecurrence::fixed_delay(
        delay,
        PredecessorWeight::Polynomial {
            coefficients: weight,
        },
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::check_weak_interlacing_bigint_coeffs;
    use crate::recurrence::RecurrenceTerm;
    use num_integer::Integer;

    fn q(n: i64) -> BigRational {
        BigRational::from_integer(n.into())
    }
    fn qq(n: i64, d: i64) -> BigRational {
        BigRational::new(n.into(), d.into())
    }
    fn poly(coefficients: &[BigRational]) -> PredecessorWeight {
        PredecessorWeight::Polynomial {
            coefficients: coefficients.to_vec(),
        }
    }
    fn integer_row(row: &[BigRational]) -> Vec<BigInt> {
        let scale = row.iter().fold(BigInt::one(), |a, b| a.lcm(b.denom()));
        row.iter()
            .map(|c| c.numer() * (&scale / c.denom()))
            .collect()
    }
    /// Finite sanity check of the proven window on generated rows.
    fn assert_window(rows: &[Vec<BigRational>], window: usize) {
        for b in 0..rows.len() {
            for a in b.saturating_sub(window)..=b {
                assert_eq!(
                    check_weak_interlacing_bigint_coeffs(
                        &integer_row(&rows[a]),
                        &integer_row(&rows[b])
                    ),
                    Some(true),
                    "Q_{a} << Q_{b}"
                );
            }
        }
    }

    #[test]
    fn fixed_delay_polynomial_weight_is_certified_symbolically() {
        // w(n) = (n - 5/2)^2 - 1/16 is negative on a real interval inside (2, 3)
        // but nonnegative at every integer n >= 1.
        let w = vec![qq(99, 16), q(-5), q(1)];
        let family = VariablePredecessorRecurrence::fixed_delay(2, poly(&w));
        let certificate = family.certify_all_indices().unwrap();
        assert_eq!(certificate.scope, CertificateScope::AllIndices);
        assert_eq!(certificate.interlacing_window, Some(2));
        let WeightCertificate::Polynomial {
            shift,
            ref shifted_coefficients,
            ref checked_values,
            ..
        } = certificate.weight
        else {
            panic!("polynomial certificate expected")
        };
        assert_eq!(shift, 3);
        assert_eq!(shifted_coefficients, &["3/16", "1", "1"]);
        assert_eq!(checked_values, &["35/16", "3/16"]);
        family.verify_certificate(&certificate).unwrap();
        let rows = family.rows(9).unwrap();
        assert_eq!(rows[1], vec![q(1), qq(35, 16)]);
        assert_window(&rows, 2);
    }

    #[test]
    fn periodic_and_zero_weights() {
        let family = VariablePredecessorRecurrence::fixed_delay(
            3,
            PredecessorWeight::Periodic {
                values: vec![q(2), q(0), qq(1, 3)],
            },
        );
        let certificate = family.certify_all_indices().unwrap();
        family.verify_certificate(&certificate).unwrap();
        assert_window(&family.rows(12).unwrap(), 3);
        // p(n) = n - 1 with w = 1 gives (1+t)^n.
        let binomial = VariablePredecessorRecurrence::fixed_delay(1, poly(&[q(1)]));
        assert_eq!(binomial.rows(4).unwrap()[3], vec![q(1), q(3), q(3), q(1)]);
        let zero = VariablePredecessorRecurrence::fixed_delay(1, poly(&[q(0)]));
        zero.verify_certificate(&zero.certify_all_indices().unwrap())
            .unwrap();
        assert_eq!(zero.rows(3).unwrap(), vec![vec![q(1)]; 3]);
    }

    #[test]
    fn explicit_data_certifies_only_a_prefix() {
        let family = VariablePredecessorRecurrence {
            predecessor: PredecessorRule::Explicit {
                values: vec![0, 0, 1, 1, 3],
            },
            weight: poly(&[q(1), q(1)]),
        };
        let certificate = family.certify().unwrap();
        assert_eq!(
            certificate.scope,
            CertificateScope::Prefix { last_index: 5 }
        );
        assert!(certificate.conclusion.contains("1 <= n <= 5 only"));
        assert!(family.certify_all_indices().is_err());
        let rows = family.rows(6).unwrap();
        assert!(family.rows(7).is_err());
        for n in 1..rows.len() {
            for i in family.predecessor_at(n).unwrap()..=n {
                assert_eq!(
                    check_weak_interlacing_bigint_coeffs(
                        &integer_row(&rows[i]),
                        &integer_row(&rows[n])
                    ),
                    Some(true)
                );
            }
        }
        // Sampled weight values with a symbolic delay are still a prefix.
        let sampled = VariablePredecessorRecurrence::fixed_delay(
            2,
            PredecessorWeight::Explicit {
                values: vec![q(1); 10],
            },
        );
        assert_eq!(
            sampled.certify().unwrap().scope,
            CertificateScope::Prefix { last_index: 10 }
        );
        assert!(sampled.certify_all_indices().is_err());
    }

    #[test]
    fn rejected_hypotheses() {
        // Positive at n = 1..9 (so sampling would accept it), negative at n = 10.
        let late = VariablePredecessorRecurrence::fixed_delay(2, poly(&[qq(199, 2), q(-20), q(1)]));
        let error = late.certify().unwrap_err();
        assert!(error.to_string().contains("w(10) = -1/2"), "{error}");
        // Negative leading coefficient: witness beyond the Cauchy bound.
        let falling = VariablePredecessorRecurrence::fixed_delay(1, poly(&[q(100), q(-1)]));
        assert!(falling
            .certify()
            .unwrap_err()
            .to_string()
            .contains("w(101)"));
        let periodic = VariablePredecessorRecurrence::fixed_delay(
            2,
            PredecessorWeight::Periodic {
                values: vec![q(1), q(-1)],
            },
        );
        assert!(periodic.certify().unwrap_err().to_string().contains("w(2)"));
        let empty = VariablePredecessorRecurrence::fixed_delay(
            2,
            PredecessorWeight::Periodic { values: vec![] },
        );
        assert!(empty.certify().is_err());
        let zero_delay = VariablePredecessorRecurrence::fixed_delay(0, poly(&[q(1)]));
        assert!(zero_delay.certify().is_err());
        let too_large = VariablePredecessorRecurrence {
            predecessor: PredecessorRule::Explicit { values: vec![0, 2] },
            weight: poly(&[q(1)]),
        };
        assert!(too_large
            .certify()
            .unwrap_err()
            .to_string()
            .contains("p(2) = 2"));
        assert!(too_large.rows(3).is_err());
        let decreasing = VariablePredecessorRecurrence {
            predecessor: PredecessorRule::Explicit {
                values: vec![0, 1, 0],
            },
            weight: poly(&[q(1)]),
        };
        assert!(decreasing
            .certify()
            .unwrap_err()
            .to_string()
            .contains("decreases"));
    }

    #[test]
    fn tampered_certificates_are_rejected() {
        let family =
            VariablePredecessorRecurrence::fixed_delay(2, poly(&[qq(99, 16), q(-5), q(1)]));
        let mut certificate = family.certify().unwrap();
        // A larger shift is still a valid certificate.
        if let WeightCertificate::Polynomial {
            shift,
            shifted_coefficients,
            checked_values,
            ..
        } = &mut certificate.weight
        {
            *shift = 4;
            *shifted_coefficients = vec!["35/16".into(), "3".into(), "1".into()];
            checked_values.push("3/16".into());
        }
        family.verify_certificate(&certificate).unwrap();
        if let WeightCertificate::Polynomial { shift, .. } = &mut certificate.weight {
            *shift = 2;
        }
        assert!(family.verify_certificate(&certificate).is_err());
        let mut certificate = family.certify().unwrap();
        certificate.scope = CertificateScope::Prefix { last_index: 3 };
        assert!(family.verify_certificate(&certificate).is_err());
    }

    #[test]
    fn recognize_index_dependent_fixed_delay_recurrence() {
        // P_n = P_(n-1) + (1+n) t P_(n-2), written with LHS factor 2.
        let rec = Recurrence {
            denominator: Some(BivarPoly {
                coeffs: vec![vec![q(2)]],
            }),
            inhomogeneous: None,
            terms: vec![
                RecurrenceTerm::new(
                    1,
                    0,
                    BivarPoly {
                        coeffs: vec![vec![q(2)]],
                    },
                ),
                RecurrenceTerm::new(
                    2,
                    0,
                    BivarPoly {
                        coeffs: vec![vec![q(0), q(2)], vec![q(0), q(2)]],
                    },
                ),
            ],
        };
        let family = recognize_fixed_delay_predecessor(&rec).unwrap();
        assert_eq!(family.interlacing_window(), Some(2));
        assert_eq!(family.weight_at(3), Some(q(4)));
        assert!(family.describe().contains("max(0, n-2)"));
        // L = 1: P_n = (1 + n t) P_(n-1).
        let rec = Recurrence {
            denominator: None,
            inhomogeneous: None,
            terms: vec![RecurrenceTerm::new(
                1,
                0,
                BivarPoly {
                    coeffs: vec![vec![q(1), q(0)], vec![q(0), q(1)]],
                },
            )],
        };
        let family = recognize_fixed_delay_predecessor(&rec).unwrap();
        assert_eq!(family.rows(3).unwrap()[2], vec![q(1), q(3), q(2)]);
        // t^2 terms, missing P_(n-1) and two weighted lags are rejected.
        let mut bad = rec.clone();
        bad.terms[0].coeff.coeffs[0].push(q(1));
        assert!(recognize_fixed_delay_predecessor(&bad).is_err());
        let mut bad = rec.clone();
        bad.terms[0].coeff.coeffs[0][0] = q(2);
        assert!(recognize_fixed_delay_predecessor(&bad).is_err());
        let mut bad = rec;
        bad.terms.push(RecurrenceTerm::new(
            3,
            0,
            BivarPoly {
                coeffs: vec![vec![q(0), q(1)]],
            },
        ));
        assert!(recognize_fixed_delay_predecessor(&bad)
            .unwrap_err()
            .to_string()
            .contains("more than one"));
    }
}
