//! Typed errors of the recurrence reports, searches, certificates and replays.
//!
//! The core types of this module keep their finer-grained errors
//! ([`RecurrenceEvaluationError`], [`VectorRecurrenceFitError`],
//! [`TransformError`], ...).  [`RecurrenceError`] covers the report, search and
//! verification layers, whose failures are descriptive messages naming the
//! offending field, index or row; [`ParseRationalError`] covers the exact
//! coefficient parser they share.

use super::TransformError;
use crate::oeis::OeisSequenceError;
use std::fmt;

/// Why an exact rational coefficient string could not be parsed
/// ([`parse_rational_coeff`](crate::recurrence::parse_rational_coeff)).
///
/// The [`Display`](fmt::Display) text of each variant is the message of the
/// former `String` error.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ParseRationalError {
    /// The string is empty or only whitespace.
    Empty,
    /// The numerator of `num/den` is not an integer.
    InvalidNumerator {
        /// The trimmed numerator text.
        text: String,
        /// The integer parser's message.
        message: String,
    },
    /// The denominator of `num/den` is not an integer.
    InvalidDenominator {
        /// The trimmed denominator text.
        text: String,
        /// The integer parser's message.
        message: String,
    },
    /// The denominator of `num/den` is zero.
    ZeroDenominator {
        /// The trimmed coefficient text.
        token: String,
    },
    /// A coefficient without `/` is not an integer.
    InvalidInteger {
        /// The trimmed coefficient text.
        token: String,
        /// The integer parser's message.
        message: String,
    },
}

impl fmt::Display for ParseRationalError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Empty => write!(f, "empty rational coefficient"),
            Self::InvalidNumerator { text, message } => {
                write!(f, "invalid numerator `{text}`: {message}")
            }
            Self::InvalidDenominator { text, message } => {
                write!(f, "invalid denominator `{text}`: {message}")
            }
            Self::ZeroDenominator { token } => {
                write!(f, "zero denominator in coefficient `{token}`")
            }
            Self::InvalidInteger { token, message } => {
                write!(f, "invalid integer `{token}`: {message}")
            }
        }
    }
}

impl std::error::Error for ParseRationalError {}

/// Why a recurrence report, search, transport or replay failed.
///
/// Each message variant carries the human-readable explanation (it names the
/// offending field, index or row), and its [`Display`](fmt::Display) text is
/// that explanation alone, as in the former `String` errors.  Failures of an
/// exact recomputation (`verify_*`, replay) are [`VerificationFailed`]; a
/// malformed or inconsistent request, option or JSON record is
/// [`InvalidInput`]; a documented ceiling on the request that was exceeded is
/// [`LimitExceeded`]; an exhausted computation budget is [`BudgetExhausted`]; and a construction that does not apply to the supplied
/// recurrence is [`Unsupported`].
///
/// [`VerificationFailed`]: RecurrenceError::VerificationFailed
/// [`InvalidInput`]: RecurrenceError::InvalidInput
/// [`LimitExceeded`]: RecurrenceError::LimitExceeded
/// [`BudgetExhausted`]: RecurrenceError::BudgetExhausted
/// [`Unsupported`]: RecurrenceError::Unsupported
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RecurrenceError {
    /// A request, option, report or JSON record is malformed, inconsistent or
    /// out of range.
    InvalidInput(String),
    /// The construction does not apply to the supplied recurrence or data.
    Unsupported(String),
    /// A documented limit or service ceiling on the request was exceeded.
    LimitExceeded(String),
    /// A resource budget of the exact computation (entries, degrees, coefficient
    /// bits, rows) ran out, so the result is incomplete rather than invalid.
    BudgetExhausted(String),
    /// A certificate, report or replay disagrees with the exact recomputation.
    VerificationFailed(String),
    /// An exact rational coefficient string could not be parsed.
    Rational(ParseRationalError),
    /// A coordinate transformation, inverse or transport was rejected.
    Transform(TransformError),
    /// A bundled OEIS family could not supply its recurrence.
    Oeis(OeisSequenceError),
}

impl RecurrenceError {
    /// Prefix the message with `context` (`"{context}: {message}"`), keeping the
    /// variant, so a caller that distinguishes an exhausted budget from invalid
    /// input still can.  The wrapping variants ([`Transform`](Self::Transform),
    /// [`Oeis`](Self::Oeis), [`Rational`](Self::Rational)) become
    /// [`InvalidInput`](Self::InvalidInput).
    pub(crate) fn with_context(self, context: impl fmt::Display) -> Self {
        match self {
            Self::InvalidInput(message) => Self::InvalidInput(format!("{context}: {message}")),
            Self::Unsupported(message) => Self::Unsupported(format!("{context}: {message}")),
            Self::LimitExceeded(message) => Self::LimitExceeded(format!("{context}: {message}")),
            Self::BudgetExhausted(message) => {
                Self::BudgetExhausted(format!("{context}: {message}"))
            }
            Self::VerificationFailed(message) => {
                Self::VerificationFailed(format!("{context}: {message}"))
            }
            other => Self::InvalidInput(format!("{context}: {other}")),
        }
    }
}

impl fmt::Display for RecurrenceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidInput(message)
            | Self::Unsupported(message)
            | Self::LimitExceeded(message)
            | Self::BudgetExhausted(message)
            | Self::VerificationFailed(message) => f.write_str(message),
            Self::Rational(error) => error.fmt(f),
            Self::Transform(error) => error.fmt(f),
            Self::Oeis(error) => error.fmt(f),
        }
    }
}

impl std::error::Error for RecurrenceError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Rational(error) => Some(error),
            Self::Transform(error) => Some(error),
            Self::Oeis(error) => Some(error),
            _ => None,
        }
    }
}

impl From<ParseRationalError> for RecurrenceError {
    fn from(error: ParseRationalError) -> Self {
        Self::Rational(error)
    }
}

impl From<OeisSequenceError> for RecurrenceError {
    fn from(error: OeisSequenceError) -> Self {
        Self::Oeis(error)
    }
}

impl From<TransformError> for RecurrenceError {
    fn from(error: TransformError) -> Self {
        Self::Transform(error)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::recurrence::parse_rational_coeff;
    use std::error::Error;

    #[test]
    fn rational_parse_errors_keep_their_messages() {
        assert_eq!(
            parse_rational_coeff("  ").unwrap_err(),
            ParseRationalError::Empty
        );
        assert_eq!(
            parse_rational_coeff("1/0").unwrap_err().to_string(),
            "zero denominator in coefficient `1/0`"
        );
        assert!(matches!(
            parse_rational_coeff("x/2").unwrap_err(),
            ParseRationalError::InvalidNumerator { .. }
        ));
        assert!(parse_rational_coeff("1/y")
            .unwrap_err()
            .to_string()
            .starts_with("invalid denominator `y`: "));
        assert!(parse_rational_coeff("q")
            .unwrap_err()
            .to_string()
            .starts_with("invalid integer `q`: "));
    }

    #[test]
    fn display_is_the_bare_message_and_sources_are_exposed() {
        let message = RecurrenceError::InvalidInput("row 3 is malformed".into());
        assert_eq!(message.to_string(), "row 3 is malformed");
        assert!(message.source().is_none());
        let wrapped = RecurrenceError::from(ParseRationalError::Empty);
        assert_eq!(wrapped.to_string(), "empty rational coefficient");
        assert!(wrapped.source().is_some());
        let transform = RecurrenceError::from(TransformError::SingularMobius);
        assert_eq!(
            transform.to_string(),
            TransformError::SingularMobius.to_string()
        );
        assert!(transform.source().is_some());
    }
}
