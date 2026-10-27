//! Exact recognition of `P_n = sum_j r_j t^j P_(n-1-ell*j)`.
//!
//! With canonical initial data, `R(z)=sum r_j z^j`, `R(0)=1`, nonnegative
//! coefficients and real zeros imply `P_n << P_(n+s)` for `degree(R)*s <=
//! degree(R)*ell+1`. This follows by periodic weighted gap subsets and the
//! nondecreasing-predecessor interlacing induction. See Fisk, arXiv:math/0612833v2,
//! Example 16.10, for the unweighted fixed-delay antecedent.
//!
//! Recognition certifies the structure of an input recurrence, not that a
//! recurrence fitted to finitely many rows holds for the intended sequence.
//!
//! Entry points: [`recognize_weighted_renewal`] turns a
//! [`Recurrence`] (see [`crate::recurrence`]) into a [`WeightedRenewal`], whose
//! [`WeightedRenewal::rows`] generates the canonical rows exactly and whose
//! [`WeightedRenewal::max_interlacing_step`] gives the proven interlacing
//! distance.  Rows are `Vec<BigRational>` in ascending powers of `t`; recognition
//! and certification failures are returned as [`RenewalError`].
//!
//! The symbolic variable-predecessor family `Q_0 = 1`,
//! `Q_n = Q_(n-1) + w(n) t Q_(p(n))` is [`VariablePredecessorRecurrence`], with
//! exact hypotheses certificates ([`VariablePredecessorCertificate`]) and the
//! canonical generator [`VariablePredecessorRecurrence::rows`].
//! [`match_recurrence_theorems`] and [`RecurrenceTheoremReport`] turn a found or
//! supplied recurrence into explicit conditional theorem matches, kept separate
//! from finite fit evidence.

use crate::is_real_rooted_bigint_coeffs;
use crate::recurrence::{
    BivarPoly, ParseRationalError, Recurrence, RecurrenceError, RecurrenceSign,
};
use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Signed, Zero};
use std::collections::BTreeMap;
use std::fmt;

mod predecessor;
pub use predecessor::*;
mod theorem;
pub use theorem::*;

/// Why a renewal recurrence was not recognized, or a variable-predecessor
/// hypothesis was not certified, rows were not generated or a certificate
/// failed its replay.
///
/// The [`Display`](fmt::Display) text of each variant is the message of the
/// former `String` error.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RenewalError {
    /// The delay `ell` of [`WeightedRenewal::new`] is zero.
    NonPositiveDelay,
    /// The seed does not start with 1 or has a negative coefficient.
    InvalidSeed,
    /// The seed polynomial is not real-rooted.
    SeedNotRealRooted,
    /// `1 + ell * degree(seed)` overflows `usize`.
    OffsetOverflow,
    /// The left-hand-side factor of the recurrence is not a nonzero constant.
    NonconstantLeftFactor,
    /// The recurrence has a nonzero inhomogeneous term.
    NotHomogeneous,
    /// A term has a derivative or an alternating sign.
    UnsupportedTermKind,
    /// A coefficient depends on the index `n`.
    IndexDependentCoefficient,
    /// After normalization the coefficient of `P_(n-1)` is not 1.
    PreviousCoefficientNotOne,
    /// A term is not of the form `t^j P_(n-1-ell*j)`.
    UnsupportedSupport,
    /// Two terms imply different delays `ell`.
    InconsistentDelays,
    /// After normalization the coefficient of `P_(n-1)` is not `1 (+ w(n) t)`.
    PreviousCoefficientNotOneWithWeight,
    /// A term is neither `P_(n-1)` nor the single weighted term `w(n) t P_(n-L)`.
    UnsupportedWeightedTerm,
    /// Weighted terms have more than one lag.
    MultipleWeightedLags,
    /// There is no weighted term `w(n) t P_(n-L)`.
    NoWeightedTerm,
    /// A fixed delay of zero gives `p(n) = n`.
    ZeroDelay,
    /// The predecessor of index `n` is not smaller than `n`.
    PredecessorNotBelowIndex {
        /// The index `n`.
        index: usize,
        /// The predecessor `p(n)`.
        predecessor: usize,
    },
    /// An explicit predecessor list decreases at index `n`.
    PredecessorDecreases {
        /// The index `n` where `p(n) < p(n - 1)`.
        index: usize,
        /// The previous predecessor `p(n - 1)`.
        previous: usize,
        /// The predecessor `p(n)`.
        predecessor: usize,
    },
    /// A periodic weight has no values.
    EmptyPeriod,
    /// The leading coefficient of a polynomial weight is negative, so the weight
    /// is negative from `index` on.
    NegativeLeadingCoefficient {
        /// An index at which the weight is negative.
        index: BigInt,
        /// The exact value `w(index)`.
        value: String,
    },
    /// A weight value is negative.
    NegativeWeight {
        /// `"periodic"` or `"explicit"` for a value list, `None` for a polynomial.
        kind: Option<&'static str>,
        /// The index `n >= 1`.
        index: usize,
        /// The exact value `w(index)`.
        value: String,
    },
    /// No certificate shift `c <= max_shift` makes `w(c + m)` coefficientwise
    /// nonnegative.
    WeightShiftNotFound {
        /// The search bound.
        max_shift: usize,
    },
    /// The minimal certificate shift exceeds the supported bound.
    WeightShiftTooLarge {
        /// The minimal shift.
        shift: usize,
        /// The supported bound.
        max_shift: usize,
    },
    /// Finite predecessor or weight data certifies only a prefix, not all indices.
    FiniteDataOnly {
        /// The last certified index.
        last_index: usize,
    },
    /// The rows need `p(n)` or `w(n)` outside the supplied finite data.
    UndefinedIndex {
        /// The index `n`.
        index: usize,
    },
    /// A certificate's scope, predecessor or window differs from the recurrence.
    CertificateScopeMismatch,
    /// A certificate's stored weight polynomial differs from the weight.
    CertificateWeightPolynomialMismatch,
    /// A certificate's shifted coefficients are wrong or negative.
    CertificateShiftedCoefficientsInvalid,
    /// A certificate's checked weight values are wrong or negative.
    CertificateCheckedValuesInvalid,
    /// A certificate's weight differs from the recurrence's weight.
    CertificateWeightMismatch,
    /// A stored coefficient string is not an exact rational.
    Rational(ParseRationalError),
    /// The supplied recurrence record could not be decoded.
    Recurrence(RecurrenceError),
}

impl fmt::Display for RenewalError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonPositiveDelay => write!(f, "ell must be positive"),
            Self::InvalidSeed => {
                write!(f, "seed must have constant 1 and nonnegative coefficients")
            }
            Self::SeedNotRealRooted => write!(f, "seed is not real-rooted"),
            Self::OffsetOverflow => write!(f, "recurrence offset overflows usize"),
            Self::NonconstantLeftFactor => write!(f, "LHS must be a nonzero constant"),
            Self::NotHomogeneous => write!(f, "recurrence must be homogeneous"),
            Self::UnsupportedTermKind => {
                write!(f, "derivatives and alternating signs are unsupported")
            }
            Self::IndexDependentCoefficient => {
                write!(f, "coefficients must be index-independent")
            }
            Self::PreviousCoefficientNotOne => {
                write!(f, "coefficient of P_(n-1) must be 1 after normalization")
            }
            Self::UnsupportedSupport => {
                write!(f, "support must satisfy offset = 1 + ell * t_degree")
            }
            Self::InconsistentDelays => write!(f, "inconsistent delays"),
            Self::PreviousCoefficientNotOneWithWeight => write!(
                f,
                "coefficient of P_(n-1) must be 1 (+ w(n) t) after normalization"
            ),
            Self::UnsupportedWeightedTerm => {
                write!(f, "only P_(n-1) and one term w(n) t P_(n-L) are supported")
            }
            Self::MultipleWeightedLags => write!(f, "more than one weighted lag"),
            Self::NoWeightedTerm => write!(f, "no weighted term w(n) t P_(n-L)"),
            Self::ZeroDelay => write!(
                f,
                "predecessor p(n) = n violates p(n) < n (delay must be >= 1)"
            ),
            Self::PredecessorNotBelowIndex { index, predecessor } => {
                write!(
                    f,
                    "predecessor p({index}) = {predecessor} violates p(n) < n"
                )
            }
            Self::PredecessorDecreases {
                index,
                previous,
                predecessor,
            } => write!(
                f,
                "predecessor decreases: p({}) = {previous} > p({index}) = {predecessor}",
                index.saturating_sub(1)
            ),
            Self::EmptyPeriod => write!(f, "periodic weight needs a positive period"),
            Self::NegativeLeadingCoefficient { index, value } => write!(
                f,
                "weight is negative: w({index}) = {value} (negative leading coefficient)"
            ),
            Self::NegativeWeight {
                kind: None,
                index,
                value,
            } => write!(f, "weight is negative: w({index}) = {value}"),
            Self::NegativeWeight {
                kind: Some(kind),
                index,
                value,
            } => write!(f, "weight is negative: {kind} w({index}) = {value}"),
            Self::WeightShiftNotFound { max_shift } => write!(
                f,
                "weight nonnegativity not certified: no shift c <= {max_shift} \
                 makes w(c+m) coefficientwise nonnegative"
            ),
            Self::WeightShiftTooLarge { shift, max_shift } => write!(
                f,
                "weight nonnegativity not certified: the certificate shift {shift} exceeds \
                 {max_shift}"
            ),
            Self::FiniteDataOnly { last_index } => write!(
                f,
                "finite predecessor/weight data certifies only n <= {last_index}; a global \
                 claim needs a symbolic rule"
            ),
            Self::UndefinedIndex { index } => {
                write!(f, "p({index}) or w({index}) is not defined")
            }
            Self::CertificateScopeMismatch => {
                write!(f, "certificate scope, predecessor or window does not match")
            }
            Self::CertificateWeightPolynomialMismatch => {
                write!(f, "certificate weight polynomial does not match")
            }
            Self::CertificateShiftedCoefficientsInvalid => {
                write!(f, "shifted weight coefficients are wrong or negative")
            }
            Self::CertificateCheckedValuesInvalid => {
                write!(f, "checked weight values are wrong or negative")
            }
            Self::CertificateWeightMismatch => write!(f, "certificate weight does not match"),
            Self::Rational(error) => error.fmt(f),
            Self::Recurrence(error) => error.fmt(f),
        }
    }
}

impl std::error::Error for RenewalError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Rational(error) => Some(error),
            Self::Recurrence(error) => Some(error),
            _ => None,
        }
    }
}

impl From<ParseRationalError> for RenewalError {
    fn from(error: ParseRationalError) -> Self {
        Self::Rational(error)
    }
}

impl From<RecurrenceError> for RenewalError {
    fn from(error: RecurrenceError) -> Self {
        Self::Recurrence(error)
    }
}

/// A PF seed and positive delay, checked exactly at construction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WeightedRenewal {
    ell: usize,
    seed: Vec<BigRational>,
}

fn constant(p: &BivarPoly) -> Option<BigRational> {
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
}

impl WeightedRenewal {
    /// Validate a seed in ascending degree order. Trailing zeros are removed.
    pub fn new(ell: usize, mut seed: Vec<BigRational>) -> Result<Self, RenewalError> {
        if ell == 0 {
            return Err(RenewalError::NonPositiveDelay);
        }
        while seed.last().is_some_and(Zero::is_zero) {
            seed.pop();
        }
        if seed.first() != Some(&BigRational::one()) || seed.iter().any(Signed::is_negative) {
            return Err(RenewalError::InvalidSeed);
        }
        let scale = seed.iter().fold(BigInt::one(), |a, b| a.lcm(b.denom()));
        let integers: Vec<_> = seed
            .iter()
            .map(|c| c.numer() * (&scale / c.denom()))
            .collect();
        if !is_real_rooted_bigint_coeffs(&integers) {
            return Err(RenewalError::SeedNotRealRooted);
        }
        ell.checked_mul(seed.len() - 1)
            .and_then(|n| n.checked_add(1))
            .ok_or(RenewalError::OffsetOverflow)?;
        Ok(Self { ell, seed })
    }

    /// Delay parameter in the generating function `1/(1-x R(t x^ell))`.
    pub fn ell(&self) -> usize {
        self.ell
    }

    /// The certified PF seed, in ascending degree order.
    pub fn seed(&self) -> &[BigRational] {
        &self.seed
    }

    /// Largest guaranteed index step; `None` means every step (constant seed).
    pub fn max_interlacing_step(&self) -> Option<usize> {
        let d = self.seed.len() - 1;
        (d * self.ell + 1).checked_div(d)
    }

    /// Canonical rows starting at `P_0=1`, with negative-index rows zero.
    /// This is the library generator for the recognized family.
    pub fn rows(&self, total_rows: usize) -> Vec<Vec<BigRational>> {
        let mut rows: Vec<Vec<BigRational>> = Vec::with_capacity(total_rows);
        for n in 0..total_rows {
            if n == 0 {
                rows.push(vec![BigRational::one()]);
                continue;
            }
            let mut row = rows[n - 1].clone();
            for (j, coefficient) in self.seed.iter().enumerate().skip(1) {
                let offset = 1 + self.ell * j;
                if offset > n {
                    break;
                }
                let previous = &rows[n - offset];
                row.resize(row.len().max(previous.len() + j), BigRational::zero());
                for (k, c) in previous.iter().enumerate() {
                    row[k + j] += coefficient * c;
                }
            }
            while row.len() > 1 && row.last().is_some_and(Zero::is_zero) {
                row.pop();
            }
            rows.push(row);
        }
        rows
    }

    /// Check enough initial rows to determine the recurrence, starting at zero.
    /// Extra supplied rows are checked too. This does not prove a fitted identity.
    pub fn matches_initial_rows(&self, first_index: usize, rows: &[Vec<BigRational>]) -> bool {
        let required = 1 + self.ell * (self.seed.len() - 1);
        if first_index != 0 || rows.len() < required {
            return false;
        }
        self.rows(rows.len())
            .iter()
            .zip(rows)
            .all(|(expected, actual)| {
                let mut end = actual.len();
                while end > 0 && actual[end - 1].is_zero() {
                    end -= 1;
                }
                expected.as_slice() == &actual[..end]
            })
    }
}

/// Recognize a stationary homogeneous renewal recurrence and check its seed.
/// Constant nonzero LHS factors and split/duplicate RHS terms are normalized.
/// Derivatives, alternating signs, index dependence and other lag supports are
/// rejected conservatively. For `R=1`, the immaterial delay is returned as 1.
pub fn recognize_weighted_renewal(
    recurrence: &Recurrence,
) -> Result<WeightedRenewal, RenewalError> {
    let denominator = match &recurrence.denominator {
        None => BigRational::one(),
        Some(p) => constant(p)
            .filter(|c| !c.is_zero())
            .ok_or(RenewalError::NonconstantLeftFactor)?,
    };
    if recurrence
        .inhomogeneous
        .as_ref()
        .is_some_and(|p| constant(p) != Some(BigRational::zero()))
    {
        return Err(RenewalError::NotHomogeneous);
    }
    let mut terms = BTreeMap::<(usize, usize), BigRational>::new();
    for term in &recurrence.terms {
        if term.deriv_order != 0 || term.sign != RecurrenceSign::None {
            return Err(RenewalError::UnsupportedTermKind);
        }
        for (i, row) in term.coeff.coeffs.iter().enumerate() {
            for (j, c) in row.iter().enumerate() {
                if c.is_zero() {
                    continue;
                }
                if i != 0 {
                    return Err(RenewalError::IndexDependentCoefficient);
                }
                *terms.entry((term.offset, j)).or_default() += c / &denominator;
            }
        }
    }
    terms.retain(|_, c| !c.is_zero());
    if terms.remove(&(1, 0)) != Some(BigRational::one()) {
        return Err(RenewalError::PreviousCoefficientNotOne);
    }
    let mut ell = None;
    let mut seed = vec![BigRational::one()];
    for ((offset, j), c) in terms {
        if j == 0 || offset <= 1 || (offset - 1) % j != 0 {
            return Err(RenewalError::UnsupportedSupport);
        }
        let delay = (offset - 1) / j;
        if ell.is_some_and(|old| old != delay) {
            return Err(RenewalError::InconsistentDelays);
        }
        ell = Some(delay);
        seed.resize(seed.len().max(j + 1), BigRational::zero());
        seed[j] += c;
    }
    WeightedRenewal::new(ell.unwrap_or(1), seed)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::recurrence::RecurrenceTerm;
    fn q(n: i64) -> BigRational {
        BigRational::from_integer(n.into())
    }
    fn recurrence(ell: usize, seed: &[i64]) -> Recurrence {
        Recurrence {
            denominator: None,
            inhomogeneous: None,
            terms: seed
                .iter()
                .enumerate()
                .map(|(j, c)| {
                    let mut row = vec![q(0); j + 1];
                    row[j] = q(*c);
                    RecurrenceTerm::new(1 + ell * j, 0, BivarPoly { coeffs: vec![row] })
                })
                .collect(),
        }
    }
    #[test]
    fn renewal_errors_are_typed_and_keep_their_messages() {
        assert_eq!(
            WeightedRenewal::new(0, vec![q(1), q(1)]).unwrap_err(),
            RenewalError::NonPositiveDelay
        );
        assert_eq!(
            WeightedRenewal::new(1, vec![q(2)]).unwrap_err().to_string(),
            "seed must have constant 1 and nonnegative coefficients"
        );
        assert_eq!(
            WeightedRenewal::new(1, vec![q(1), q(0), q(1)]).unwrap_err(),
            RenewalError::SeedNotRealRooted
        );
        let mut inhomogeneous = recurrence(1, &[1, 1]);
        inhomogeneous.inhomogeneous = Some(BivarPoly {
            coeffs: vec![vec![q(1)]],
        });
        assert_eq!(
            recognize_weighted_renewal(&inhomogeneous).unwrap_err(),
            RenewalError::NotHomogeneous
        );
        let source = RenewalError::from(crate::recurrence::ParseRationalError::Empty);
        assert_eq!(source.to_string(), "empty rational coefficient");
        assert!(std::error::Error::source(&source).is_some());
    }

    #[test]
    fn pf_seed_and_independent_coefficient_formula() {
        let family = recognize_weighted_renewal(&recurrence(1, &[1, 3, 2])).unwrap();
        // [t^j] P_n = [z^j] (1+3z+2z^2)^(n-j).
        assert_eq!(
            family.rows(6),
            vec![
                vec![q(1)],
                vec![q(1)],
                vec![q(1), q(3)],
                vec![q(1), q(6), q(2)],
                vec![q(1), q(9), q(13)],
                vec![q(1), q(12), q(33), q(12)]
            ]
        );
        assert_eq!(family.max_interlacing_step(), Some(1));
        assert!(family.matches_initial_rows(0, &family.rows(5)));
        let mut wrong = family.rows(5);
        wrong[1] = vec![q(1), q(1)];
        assert!(!family.matches_initial_rows(0, &wrong));
        assert!(!family.matches_initial_rows(1, &family.rows(5)));
        assert!(!family.matches_initial_rows(0, &family.rows(2)));
    }
    #[test]
    fn normalize_duplicates_and_rational_seed() {
        let mut rec = recurrence(2, &[1, 2, 1]);
        rec.terms.extend(rec.terms.clone());
        rec.denominator = Some(BivarPoly {
            coeffs: vec![vec![q(2)]],
        });
        assert_eq!(
            recognize_weighted_renewal(&rec)
                .unwrap()
                .max_interlacing_step(),
            Some(2)
        );
        assert!(WeightedRenewal::new(1, vec![q(1), q(1), q(1) / q(4)]).is_ok());
        assert!(WeightedRenewal::new(1, vec![q(1), q(1), q(1) / q(3)]).is_err());
    }
    #[test]
    fn reject_false_positive_structures() {
        assert!(recognize_weighted_renewal(&recurrence(1, &[1, 1, 1])).is_err());
        assert!(recognize_weighted_renewal(&recurrence(1, &[1, -1])).is_err());
        assert!(recognize_weighted_renewal(&recurrence(0, &[1, 1])).is_err());
        let mut rec = recurrence(1, &[1, 2, 1]);
        rec.terms[2].offset += 1;
        assert!(recognize_weighted_renewal(&rec).is_err());
        rec = recurrence(1, &[1, 1]);
        rec.terms[1].deriv_order = 1;
        assert!(recognize_weighted_renewal(&rec).is_err());
        rec.terms[1].deriv_order = 0;
        rec.terms[1].sign = RecurrenceSign::AlternatingN;
        assert!(recognize_weighted_renewal(&rec).is_err());
        rec = recurrence(1, &[1, 1]);
        rec.terms[1].coeff.coeffs.push(vec![q(1)]);
        assert!(recognize_weighted_renewal(&rec).is_err());
        rec = recurrence(1, &[1, 1]);
        rec.inhomogeneous = Some(BivarPoly {
            coeffs: vec![vec![q(1)]],
        });
        assert!(recognize_weighted_renewal(&rec).is_err());
    }
    #[test]
    fn constant_and_long_delay_boundaries() {
        let family = recognize_weighted_renewal(&recurrence(1, &[1])).unwrap();
        assert_eq!(family.max_interlacing_step(), None);
        assert_eq!(family.rows(3), vec![vec![q(1)]; 3]);
        let family = recognize_weighted_renewal(&recurrence(3, &[1, 2])).unwrap();
        assert_eq!(family.max_interlacing_step(), Some(4));
        assert_eq!(family.rows(5)[4], vec![q(1), q(2)]);
        assert!(WeightedRenewal::new(usize::MAX, vec![q(1), q(1)]).is_err());
    }
}
