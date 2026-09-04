//! Recurrence-backed OEIS polynomial families.
//!
//! Every catalog entry is a small sparse recurrence description plus the
//! initial rows needed to replay it.  The generated ID functions, such as
//! [`A008292`], make individual families convenient to use while [`by_id`] and
//! [`catalog`] support dynamic discovery.

use crate::recurrence::{
    parse_rational_coeff, BigRational, BivarPoly, Recurrence, RecurrenceEvaluationError,
    RecurrenceSign, RecurrenceTerm,
};
use num_bigint::BigInt;
use num_traits::One;
use std::collections::BTreeSet;
use std::fmt;

/// Publication status of one bundled recurrence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OeisSequenceStatus {
    /// The recurrence was checked on rows excluded from fitting.
    Verified,
    /// The recurrence matches a current OEIS prefix, but fitting provenance is unavailable.
    Validated,
    /// The recurrence generates known rows but still lacks a genuine holdout.
    Experimental,
}

impl OeisSequenceStatus {
    /// Stable lowercase name used by the CLI and MCP server.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Verified => "verified",
            Self::Validated => "validated",
            Self::Experimental => "experimental",
        }
    }
}

/// Shape of the OEIS data represented by polynomial coefficient rows.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OeisLayout {
    /// OEIS `tabl`: successive rows have lengths `1, 2, 3, ...`.
    RegularTriangle,
    /// OEIS `tabf` or another explicitly row-grouped table.
    Table,
}

impl OeisLayout {
    /// Stable lowercase name used by structured output.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::RegularTriangle => "regular_triangle",
            Self::Table => "table",
        }
    }
}

/// One nonzero coefficient `value * n^n_degree * t^variable_degree`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SparseCoefficient {
    pub n_degree: usize,
    pub variable_degree: usize,
    pub value: &'static str,
}

impl SparseCoefficient {
    /// Construct one sparse bivariate monomial coefficient.
    pub const fn new(n_degree: usize, variable_degree: usize, value: &'static str) -> Self {
        Self {
            n_degree,
            variable_degree,
            value,
        }
    }
}

/// One differential-polynomial term in a scalar recurrence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SparseRecurrenceTerm {
    pub offset: usize,
    pub derivative_order: usize,
    pub alternating_sign: bool,
    pub coefficient: &'static [SparseCoefficient],
}

/// Sparse executable form of a recurrence found by polytool.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SparseRecurrence {
    pub terms: &'static [SparseRecurrenceTerm],
    pub denominator: Option<&'static [SparseCoefficient]>,
    pub inhomogeneous: Option<&'static [SparseCoefficient]>,
}

/// One bundled OEIS polynomial family.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OeisSequenceDefinition {
    pub id: &'static str,
    pub name: &'static str,
    pub status: OeisSequenceStatus,
    pub layout: OeisLayout,
    /// Displayed index of the first bundled row.
    pub first_row: i64,
    /// Flattened OEIS index of the first bundled coefficient.
    pub flattened_offset: i64,
    /// Whether the bundled flattened prefix agrees with local OEIS data.
    pub bfile_prefix_verified: bool,
    /// Recurrence benchmark slug used for exact regression testing.
    pub fixture_slug: &'static str,
    /// Number of canonical rows preceding the first fixture row.
    pub fixture_row_offset: usize,
    /// Complete rows before the recurrence's first row.
    pub prefix_rows: &'static [&'static [&'static str]],
    /// Index at which recurrence coefficients evaluate their `n` variable.
    pub recurrence_first_index: usize,
    /// Structural width of the first recurrence-generated row.
    pub recurrence_first_width: usize,
    /// Initial rows supplied to the recurrence evaluator.
    pub initial_rows: &'static [&'static [&'static str]],
    pub recurrence: SparseRecurrence,
    pub source_rows: usize,
    pub verification_rows: usize,
    pub rows_sha256: &'static str,
}

/// Failure while decoding or evaluating a bundled sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OeisSequenceError {
    InvalidDefinition { id: &'static str, message: String },
    Evaluation(RecurrenceEvaluationError),
    NonIntegralCoefficient { id: &'static str, row: i64 },
    RowBeforeStart { requested: i64, first: i64 },
    BFilePrefixNotVerified { id: &'static str },
}

impl fmt::Display for OeisSequenceError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidDefinition { id, message } => {
                write!(formatter, "invalid OEIS catalog entry {id}: {message}")
            }
            Self::Evaluation(error) => error.fmt(formatter),
            Self::NonIntegralCoefficient { id, row } => {
                write!(formatter, "{id} generated a nonintegral coefficient in row {row}")
            }
            Self::RowBeforeStart { requested, first } => {
                write!(formatter, "requested row {requested} is before first row {first}")
            }
            Self::BFilePrefixNotVerified { id } => write!(
                formatter,
                "strict b-file output for {id} is disabled because its OEIS prefix mapping is not verified"
            ),
        }
    }
}

impl std::error::Error for OeisSequenceError {}

impl From<RecurrenceEvaluationError> for OeisSequenceError {
    fn from(error: RecurrenceEvaluationError) -> Self {
        Self::Evaluation(error)
    }
}

fn decode_sparse_polynomial(
    id: &'static str,
    terms: &[SparseCoefficient],
) -> Result<BivarPoly, OeisSequenceError> {
    let n_degree = terms.iter().map(|term| term.n_degree).max().unwrap_or(0);
    let variable_degree = terms
        .iter()
        .map(|term| term.variable_degree)
        .max()
        .unwrap_or(0);
    let mut coeffs =
        vec![vec![BigRational::from_integer(BigInt::from(0)); variable_degree + 1]; n_degree + 1];
    let mut seen = BTreeSet::new();
    for term in terms {
        if !seen.insert((term.n_degree, term.variable_degree)) {
            return Err(OeisSequenceError::InvalidDefinition {
                id,
                message: format!(
                    "duplicate sparse monomial n^{} t^{}",
                    term.n_degree, term.variable_degree
                ),
            });
        }
        coeffs[term.n_degree][term.variable_degree] =
            parse_rational_coeff(term.value).map_err(|message| {
                OeisSequenceError::InvalidDefinition {
                    id,
                    message: format!("invalid recurrence coefficient `{}`: {message}", term.value),
                }
            })?;
    }
    Ok(BivarPoly { coeffs })
}

fn decode_rows(
    id: &'static str,
    rows: &[&[&str]],
) -> Result<Vec<Vec<BigRational>>, OeisSequenceError> {
    rows.iter()
        .enumerate()
        .map(|(row_index, row)| {
            row.iter()
                .enumerate()
                .map(|(coefficient_index, value)| {
                    parse_rational_coeff(value).map_err(|message| {
                        OeisSequenceError::InvalidDefinition {
                            id,
                            message: format!(
                                "invalid initial coefficient at row {row_index}, coefficient \
                                 {coefficient_index}: {message}"
                            ),
                        }
                    })
                })
                .collect()
        })
        .collect()
}

impl OeisSequenceDefinition {
    /// Decode the sparse static description into the recurrence evaluator.
    pub fn recurrence_parts(
        &self,
    ) -> Result<(Recurrence, Vec<Vec<BigRational>>), OeisSequenceError> {
        let terms = self
            .recurrence
            .terms
            .iter()
            .map(|term| {
                Ok(RecurrenceTerm {
                    offset: term.offset,
                    deriv_order: term.derivative_order,
                    sign: if term.alternating_sign {
                        RecurrenceSign::AlternatingN
                    } else {
                        RecurrenceSign::None
                    },
                    coeff: decode_sparse_polynomial(self.id, term.coefficient)?,
                })
            })
            .collect::<Result<Vec<_>, OeisSequenceError>>()?;
        let denominator = self
            .recurrence
            .denominator
            .map(|terms| decode_sparse_polynomial(self.id, terms))
            .transpose()?;
        let inhomogeneous = self
            .recurrence
            .inhomogeneous
            .map(|terms| decode_sparse_polynomial(self.id, terms))
            .transpose()?;
        let initial_rows = decode_rows(self.id, self.initial_rows)?;
        Ok((
            Recurrence {
                terms,
                denominator,
                inhomogeneous,
            },
            initial_rows,
        ))
    }

    /// Generate the first `row_count` complete rows in catalog order.
    pub fn generate_rows(&self, row_count: usize) -> Result<Vec<Vec<BigInt>>, OeisSequenceError> {
        let prefix = decode_rows(self.id, self.prefix_rows)?;
        let needed_recurrence_rows = row_count.saturating_sub(prefix.len());
        let recurrence_rows = if needed_recurrence_rows == 0 {
            Vec::new()
        } else {
            let (recurrence, initial_rows) = self.recurrence_parts()?;
            recurrence.generate_rows_rational(
                &initial_rows,
                self.recurrence_first_index,
                needed_recurrence_rows,
            )?
        };
        let prefix_count = prefix.len();
        prefix
            .into_iter()
            .take(row_count)
            .chain(recurrence_rows)
            .enumerate()
            .map(|(index, row)| {
                let displayed_row = self.first_row + index as i64;
                let mut integers = row
                    .into_iter()
                    .map(|coefficient| {
                        if coefficient.denom() != &BigInt::one() {
                            Err(OeisSequenceError::NonIntegralCoefficient {
                                id: self.id,
                                row: displayed_row,
                            })
                        } else {
                            Ok(coefficient.numer().clone())
                        }
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                if self.bfile_prefix_verified
                    && self.layout == OeisLayout::RegularTriangle
                    && index >= prefix_count
                {
                    let recurrence_index = index - prefix_count;
                    let expected_width = self
                        .recurrence_first_width
                        .checked_add(recurrence_index)
                        .ok_or_else(|| OeisSequenceError::InvalidDefinition {
                            id: self.id,
                            message: "row width overflow".to_string(),
                        })?;
                    integers.resize(expected_width, BigInt::from(0));
                }
                Ok(integers)
            })
            .collect()
    }

    /// Generate `row_count` rows beginning at a displayed OEIS row index.
    pub fn generate_rows_from(
        &self,
        first_row: i64,
        row_count: usize,
    ) -> Result<Vec<Vec<BigInt>>, OeisSequenceError> {
        if first_row < self.first_row {
            return Err(OeisSequenceError::RowBeforeStart {
                requested: first_row,
                first: self.first_row,
            });
        }
        let skipped = usize::try_from(first_row - self.first_row).unwrap_or(usize::MAX);
        let rows = self.generate_rows(skipped.saturating_add(row_count))?;
        Ok(rows.into_iter().skip(skipped).collect())
    }

    /// Flattened index corresponding to the beginning of `displayed_row`.
    pub fn flattened_offset_for_row(&self, displayed_row: i64) -> Result<i64, OeisSequenceError> {
        if !self.bfile_prefix_verified {
            return Err(OeisSequenceError::BFilePrefixNotVerified { id: self.id });
        }
        if displayed_row < self.first_row {
            return Err(OeisSequenceError::RowBeforeStart {
                requested: displayed_row,
                first: self.first_row,
            });
        }
        let skipped = usize::try_from(displayed_row - self.first_row).unwrap_or(usize::MAX);
        let rows = self.generate_rows(skipped)?;
        let skipped_terms = rows.iter().try_fold(0i64, |sum, row| {
            i64::try_from(row.len())
                .ok()
                .and_then(|length| sum.checked_add(length))
        });
        skipped_terms
            .and_then(|terms| self.flattened_offset.checked_add(terms))
            .ok_or_else(|| OeisSequenceError::InvalidDefinition {
                id: self.id,
                message: "flattened offset overflow".to_string(),
            })
    }
}

include!("oeis_catalog_generated.rs");

/// All bundled OEIS sequence definitions, sorted by A-number.
pub fn catalog() -> &'static [OeisSequenceDefinition] {
    OEIS_CATALOG
}

/// Look up an A-number case-insensitively.
pub fn by_id(id: &str) -> Option<&'static OeisSequenceDefinition> {
    let normalized = id.to_ascii_uppercase();
    OEIS_CATALOG
        .binary_search_by_key(&normalized.as_str(), |entry| entry.id)
        .ok()
        .map(|index| &OEIS_CATALOG[index])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn catalog_is_sorted_unique_and_decodable() {
        assert_eq!(catalog().len(), 785);
        assert_eq!(
            catalog()
                .iter()
                .filter(|entry| entry.status == OeisSequenceStatus::Verified)
                .count(),
            125
        );
        assert_eq!(
            catalog()
                .iter()
                .filter(|entry| entry.status == OeisSequenceStatus::Validated)
                .count(),
            630
        );
        assert_eq!(
            catalog()
                .iter()
                .filter(|entry| entry.status == OeisSequenceStatus::Experimental)
                .count(),
            30
        );
        for pair in catalog().windows(2) {
            assert!(pair[0].id < pair[1].id);
        }
        for entry in catalog() {
            entry.recurrence_parts().expect(entry.id);
            if entry.status == OeisSequenceStatus::Verified {
                assert!(
                    entry.verification_rows > 0,
                    "{} lacks holdout rows",
                    entry.id
                );
            }
        }
    }

    #[test]
    fn direct_id_function_and_lookup_agree() {
        assert_eq!(A008292(), by_id("a008292").unwrap());
    }

    #[test]
    fn a008292_begins_with_eulerian_triangle_rows() {
        let rows = A008292().generate_rows(4).unwrap();
        assert_eq!(
            rows,
            vec![
                vec![1.into()],
                vec![1.into(), 1.into()],
                vec![1.into(), 4.into(), 1.into()],
                vec![1.into(), 11.into(), 11.into(), 1.into()],
            ]
        );
    }

    #[test]
    fn bfile_offsets_account_for_complete_rows() {
        let sequence = A008292();
        assert!(sequence.bfile_prefix_verified);
        assert_eq!(sequence.flattened_offset_for_row(1).unwrap(), 1);
        assert_eq!(sequence.flattened_offset_for_row(4).unwrap(), 7);
    }

    #[test]
    fn every_sparse_definition_reproduces_its_fixture_rows() {
        fn trim(mut row: Vec<BigInt>) -> Vec<BigInt> {
            while row.last().is_some_and(|value| value == &BigInt::from(0)) {
                row.pop();
            }
            if row.is_empty() {
                row.push(BigInt::from(0));
            }
            row
        }

        let fixtures = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("fixtures/recurrence-benchmarks/rows");
        for sequence in catalog() {
            if sequence.fixture_slug.is_empty() {
                continue;
            }
            let path = fixtures.join(format!("{}.txt", sequence.fixture_slug));
            let expected = std::fs::read_to_string(&path)
                .unwrap_or_else(|error| panic!("failed to read {}: {error}", path.display()))
                .lines()
                .filter(|line| !line.trim().is_empty())
                .map(|line| {
                    line.split(',')
                        .map(|value| value.trim().parse::<BigInt>().unwrap())
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();
            let generated = sequence
                .generate_rows(sequence.fixture_row_offset + expected.len())
                .unwrap_or_else(|error| panic!("{}: {error}", sequence.id));
            let actual = generated[sequence.fixture_row_offset..]
                .iter()
                .cloned()
                .map(trim)
                .collect::<Vec<_>>();
            let expected = expected.into_iter().map(trim).collect::<Vec<_>>();
            assert_eq!(actual, expected, "{}", sequence.id);
        }
    }

    #[test]
    fn every_imported_lean_definition_reproduces_its_validation_row() {
        for (id, row_index, expected) in OEIS_IMPORTED_VALIDATION_ROWS {
            let generated = by_id(id)
                .unwrap()
                .generate_rows(row_index + 1)
                .unwrap_or_else(|error| panic!("{id}: {error}"));
            let actual = &generated[*row_index];
            let expected = expected
                .iter()
                .map(|value| value.parse::<BigInt>().unwrap())
                .collect::<Vec<_>>();
            assert_eq!(actual, &expected, "{id} row {row_index}");
        }
    }
}
