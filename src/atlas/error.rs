//! The error type of the atlas engine.

use crate::recurrence::RecurrenceError;
use std::fmt;

/// Errors of the atlas engine: storage, malformed records, certificates that
/// do not replay and claims that do not follow from their certificates.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AtlasError {
    /// Reading or writing a data file failed.
    Io(String),
    /// A record does not parse or violates the format.
    InvalidRecord(String),
    /// A certificate does not replay, or a claim is not backed.
    VerificationFailed(String),
    /// The input is outside what the engine supports (not an error of the
    /// data).
    Unsupported(String),
    /// An error of the recurrence machinery.
    Recurrence(RecurrenceError),
}

impl fmt::Display for AtlasError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(message) => write!(f, "atlas I/O error: {message}"),
            Self::InvalidRecord(message) => write!(f, "invalid atlas record: {message}"),
            Self::VerificationFailed(message) => {
                write!(f, "atlas verification failed: {message}")
            }
            Self::Unsupported(message) => write!(f, "unsupported by the atlas engine: {message}"),
            Self::Recurrence(error) => write!(f, "{error}"),
        }
    }
}

impl std::error::Error for AtlasError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Recurrence(error) => Some(error),
            _ => None,
        }
    }
}

impl From<RecurrenceError> for AtlasError {
    fn from(error: RecurrenceError) -> Self {
        Self::Recurrence(error)
    }
}

impl From<std::io::Error> for AtlasError {
    fn from(error: std::io::Error) -> Self {
        Self::Io(error.to_string())
    }
}

impl From<serde_json::Error> for AtlasError {
    fn from(error: serde_json::Error) -> Self {
        Self::InvalidRecord(error.to_string())
    }
}
