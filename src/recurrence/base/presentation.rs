//! Exact validity boundaries and primitive integral presentation of fitted rules.

use super::*;
use crate::Polynomial;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use serde::{Deserialize, Serialize};

/// Exact seed requirements for a rule on a finite supplied prefix.
///
/// Replay is finite evidence, not a proof of the identity at later indices.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RecurrencePresentation {
    /// Earliest absolute index after every failed supplied identity and available lag.
    pub earliest_valid_index: usize,
    /// All nonnegative absolute indices where the leader is identically zero in t.
    /// Exact integers are strings, including those beyond machine index bounds.
    /// Empty for an everywhere-zero leader; see `left_hand_identically_zero`.
    pub left_hand_zeros: Vec<String>,
    /// True precisely when the leader is identically zero for every index.
    pub left_hand_identically_zero: bool,
    /// Start after replay failures and the last zero; absent if the leader is zero
    /// everywhere or the required seed rows have not all been supplied.
    pub start_index: Option<usize>,
    /// `ready`, `insufficient_initial_rows`, or `zero_left_hand_factor`.
    pub status: String,
    /// Explanation when no defining range can be stated.
    pub reason: Option<String>,
}

impl Recurrence {
    /// Clear denominators, divide the content of the entire identity, and make
    /// the leading coefficient of the leader positive (n first, then t).
    /// The original fitted scaling is not modified. No polynomial factors are
    /// cancelled: their zeros carry seed information.
    pub fn normalized_presentation(&self) -> Result<Self, RecurrenceError> {
        let one = BivarPoly {
            coeffs: vec![vec![BigRational::one()]],
        };
        let leader = self.denominator.as_ref().unwrap_or(&one);
        let leading = leader
            .coeffs
            .iter()
            .rev()
            .flat_map(|r| r.iter().rev())
            .find(|c| !c.is_zero())
            .ok_or_else(|| {
                RecurrenceError::InvalidInput(
                    "cannot normalize an identically zero left-hand factor".into(),
                )
            })?;
        let polynomials = std::iter::once(leader)
            .chain(self.terms.iter().map(|term| &term.coeff))
            .chain(self.inhomogeneous.iter());
        let coefficients: Vec<_> = polynomials
            .flat_map(|p| p.coeffs.iter().flatten())
            .collect();
        let denominator = coefficients
            .iter()
            .fold(BigInt::one(), |a, c| a.lcm(c.denom()));
        let content = coefficients.iter().fold(BigInt::zero(), |a, c| {
            a.gcd(&(c.numer() * (&denominator / c.denom())))
        });
        let mut scale = BigRational::new(denominator, content);
        if leading.is_negative() {
            scale = -scale;
        }
        let scaled = |p: &BivarPoly| {
            BivarPoly::from_coefficients_normalized(
                p.coeffs
                    .iter()
                    .map(|row| row.iter().map(|c| c * &scale).collect())
                    .collect(),
            )
        };
        let mut result = self.clone();
        let normalized_leader = scaled(leader);
        result.denominator = (!normalized_leader.is_one()).then_some(normalized_leader);
        for term in &mut result.terms {
            term.coeff = scaled(&term.coeff);
        }
        result.inhomogeneous = self.inhomogeneous.as_ref().map(scaled);
        Ok(result)
    }

    /// Find all nonnegative integers n where the left-hand factor vanishes
    /// identically as a polynomial in t, using the gcd of its t-coefficients.
    /// A t-dependent factor has no scalar "zero at n" unless every t-coefficient
    /// vanishes there. An identically zero leader is an error (infinitely many zeros).
    pub fn left_hand_nonnegative_integer_zeros(&self) -> Result<Vec<BigInt>, RecurrenceError> {
        let Some(leader) = &self.denominator else {
            return Ok(Vec::new());
        };
        let mut gcd = Polynomial::<BigRational>::zero();
        for j in 0..leader.coefficient_width() {
            gcd = gcd.gcd(&Polynomial::new(
                leader
                    .coeffs
                    .iter()
                    .map(|row| row.get(j).cloned().unwrap_or_else(BigRational::zero))
                    .collect(),
            ));
        }
        if gcd.is_zero() {
            return Err(RecurrenceError::InvalidInput(
                "left-hand factor is identically zero at every index".into(),
            ));
        }
        Ok(
            crate::recurrence::implicit::algebra::nonnegative_integer_roots(gcd.coeffs())
                .into_iter()
                .collect(),
        )
    }

    /// Replay every supplied identity with available lags, then require seeds
    /// through the last nonnegative integer zero of the leader, including future
    /// zeros beyond the supplied data. Indices use `first_index` throughout.
    pub fn presentation_for_rows(
        &self,
        rows: &[Vec<BigRational>],
        first_index: usize,
    ) -> Result<RecurrencePresentation, RecurrenceError> {
        let end = first_index.checked_add(rows.len()).ok_or_else(|| {
            RecurrenceError::InvalidInput("recurrence presentation index overflow".into())
        })?;
        let shift = i64::try_from(first_index)
            .ok()
            .and_then(|n| n.checked_sub(1))
            .ok_or_else(|| RecurrenceError::InvalidInput("recurrence index exceeds i64".into()))?;
        let local = self.shift_index(shift);
        let mut count = self.max_offset();
        for offset in self.max_offset()..rows.len() {
            if !recurrence_holds_at_rational(rows, &local, offset + 1) {
                count = offset + 1;
            }
        }
        let earliest = first_index.checked_add(count).ok_or_else(|| {
            RecurrenceError::InvalidInput("recurrence seed count overflow".into())
        })?;
        let zeros = match self.left_hand_nonnegative_integer_zeros() {
            Ok(zeros) => zeros,
            Err(error) => {
                return Ok(RecurrencePresentation {
                    earliest_valid_index: earliest,
                    left_hand_zeros: Vec::new(),
                    left_hand_identically_zero: true,
                    start_index: None,
                    status: "zero_left_hand_factor".into(),
                    reason: Some(error.to_string()),
                })
            }
        };
        let required = zeros.last().map_or_else(
            || BigInt::from(earliest),
            |last| (last + BigInt::one()).max(BigInt::from(earliest)),
        );
        let enough = required <= BigInt::from(end);
        Ok(RecurrencePresentation {
            earliest_valid_index: earliest,
            left_hand_zeros: zeros.iter().map(ToString::to_string).collect(),
            left_hand_identically_zero: false,
            start_index: if enough {
                required.to_string().parse().ok()
            } else {
                None
            },
            status: if enough {
                "ready"
            } else {
                "insufficient_initial_rows"
            }
            .into(),
            reason: (!enough).then(|| {
                format!(
                "supplied data ends before all required initial rows: need rows through index {}",
                required - BigInt::one())
            }),
        })
    }
}
