//! Infer coefficient rows from an OEIS internal-format record.

use super::{parse_bfile, parse_oeis_internal_format, OeisSourceError, OeisTermsSource};
use num_bigint::BigInt;

/// A regular rule for successive row lengths.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OeisRowLengthRule {
    /// Each successive length increases by this amount.
    ConstantIncrease(usize),
    /// Length increases by one every `period` rows, at the given zero-based phase.
    PeriodicIncrease {
        /// Number of transitions per increase.
        period: usize,
        /// Transition index (from zero) modulo `period` where increases occur.
        phase: usize,
    },
}

/// Whether complete rows were inferred or inference was refused.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OeisRowsStatus {
    /// At least one complete row was inferred.
    Ok,
    /// The record did not support a row inference.
    Refused,
}

/// Reason that row inference was refused.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum OeisRowsRefusalReason {
    /// Neither triangle keyword was present.
    NotTriangleKeyword,
    /// No numeric example table was found.
    NoExampleTable,
    /// An example table starts after the beginning of the terms.
    TableDoesNotStartAtFirstTerm,
    /// There is not enough evidence to establish a rule.
    TooFewExampleRows,
    /// The example row lengths do not follow a supported regular rule.
    IrregularLengths,
    /// An example table disagrees with the terms.
    TermsDisagreeWithTable,
    /// No complete row is available.
    NoCompleteRows,
}

impl OeisRowsRefusalReason {
    /// Stable machine-readable identifier.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::NotTriangleKeyword => "not_triangle_keyword",
            Self::NoExampleTable => "no_example_table",
            Self::TableDoesNotStartAtFirstTerm => "table_does_not_start_at_first_term",
            Self::TooFewExampleRows => "too_few_example_rows",
            Self::IrregularLengths => "irregular_lengths",
            Self::TermsDisagreeWithTable => "terms_disagree_with_table",
            Self::NoCompleteRows => "no_complete_rows",
        }
    }
    /// Human-readable explanation.
    #[must_use]
    pub fn message(self) -> &'static str {
        match self {
            Self::NotTriangleKeyword => "not tabl or tabf",
            Self::NoExampleTable => "no table in %e",
            Self::TableDoesNotStartAtFirstTerm => "table in %e does not start at the first term",
            Self::TooFewExampleRows => "too few example rows to establish a rule",
            Self::IrregularLengths => "irregular lengths",
            Self::TermsDisagreeWithTable => "terms disagree with the table",
            Self::NoCompleteRows => "no complete rows",
        }
    }
}

impl std::fmt::Display for OeisRowsRefusalReason {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.message())
    }
}

/// Why an available b-file could not be used.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum OeisBFileRejection {
    /// The file could not be read as UTF-8.
    Unreadable,
    /// The file does not have valid b-file syntax.
    InvalidBFile,
    /// Indices are not consecutive.
    NonConsecutiveIndices,
    /// The values do not begin with the internal-format terms.
    PrefixMismatch,
}
impl OeisBFileRejection {
    /// Stable machine-readable identifier.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Unreadable => "unreadable",
            Self::InvalidBFile => "invalid_bfile",
            Self::NonConsecutiveIndices => "non_consecutive_indices",
            Self::PrefixMismatch => "prefix_mismatch",
        }
    }
    /// Human-readable explanation.
    #[must_use]
    pub fn message(self) -> &'static str {
        match self {
            Self::Unreadable => "b-file is unreadable or not UTF-8",
            Self::InvalidBFile => "not a valid b-file",
            Self::NonConsecutiveIndices => "b-file indices are not consecutive",
            Self::PrefixMismatch => "b-file values do not match the term prefix",
        }
    }
}

/// Outcome of inferring complete rows from a local OEIS record.
#[derive(Debug, Clone)]
pub struct OeisInferredRows {
    /// Sequence title.
    pub name: String,
    /// Keywords from the record.
    pub keywords: Vec<String>,
    /// Accepted source of the terms.
    pub term_source: OeisTermsSource,
    /// Inferred rule, absent when inference was refused.
    pub row_length_rule: Option<OeisRowLengthRule>,
    /// Complete coefficient rows in their given order.
    pub rows: Vec<Vec<BigInt>>,
    /// Refusal reason, or `None` when rows were inferred.
    pub reason: Option<OeisRowsRefusalReason>,
    /// Why a supplied b-file was not used, if any.
    pub bfile_rejection: Option<OeisBFileRejection>,
}

impl OeisInferredRows {
    /// Status of the inference.
    #[must_use]
    pub fn status(&self) -> OeisRowsStatus {
        if self.reason.is_some() {
            OeisRowsStatus::Refused
        } else {
            OeisRowsStatus::Ok
        }
    }
}

fn table_row(line: &str) -> Option<Vec<BigInt>> {
    let mut parts = line.splitn(3, ' ');
    if parts.next()? != "%e" {
        return None;
    }
    parts.next()?;
    let mut text = parts.next()?.trim();
    if let Some(rest) = text.strip_prefix('[') {
        let (label, rest) = rest.split_once(']')?;
        label.trim().parse::<i64>().ok()?;
        text = rest.trim_start();
    } else if let Some((label, rest)) = text.split_once([':', '|']) {
        let label = label.trim().strip_prefix("n=").unwrap_or(label.trim());
        label.trim().parse::<i64>().ok()?;
        text = rest.trim_start();
    } else if let Some(rest) = text.strip_prefix("n=") {
        let (label, rest) = rest.split_once(char::is_whitespace)?;
        label.trim().parse::<i64>().ok()?;
        text = rest.trim_start();
    } else if let Some((label, rest)) = text.split_once('=') {
        label.trim().parse::<i64>().ok()?;
        text = rest.trim_start();
    }
    let text = text.trim_end_matches(['.', ';', ',', ' ']);
    let tokens: Vec<_> = text
        .split(|c: char| c == ',' || c.is_ascii_whitespace())
        .filter(|s| !s.is_empty())
        .collect();
    if tokens.is_empty() {
        return None;
    }
    tokens
        .into_iter()
        .map(|s| s.parse::<BigInt>().ok())
        .collect()
}

fn regular_rule(lengths: &[usize]) -> Option<OeisRowLengthRule> {
    if lengths.len() < 3 {
        return None;
    }
    let differences: Vec<_> = lengths.windows(2).map(|w| w[1].checked_sub(w[0])).collect();
    if let Some(Some(step)) = differences.first() {
        if *step > 0 && differences.iter().all(|d| *d == Some(*step)) {
            return Some(OeisRowLengthRule::ConstantIncrease(*step));
        }
    }
    if differences.iter().any(|d| !matches!(d, Some(0 | 1))) {
        return None;
    }
    // Require two observed increases to identify a period and its phase.
    let ones: Vec<_> = differences
        .iter()
        .enumerate()
        .filter_map(|(i, d)| (*d == Some(1)).then_some(i))
        .collect();
    if ones.len() < 2 {
        return None;
    }
    let period = ones[1] - ones[0];
    if period < 2
        || differences
            .iter()
            .enumerate()
            .any(|(i, d)| (*d == Some(1)) != (i % period == ones[0] % period))
    {
        return None;
    }
    Some(OeisRowLengthRule::PeriodicIncrease {
        period,
        phase: ones[0] % period,
    })
}

/// Infer complete rows from internal-format text and an optional b-file.
///
/// For `tabf`, only a run of consecutive example lines matching the beginning
/// of the terms is used; at least three rows and two periodic increases are
/// required to infer a rule. Incomplete final rows are discarded.
///
/// # Errors
/// Returns the internal-format parser error if the record is invalid.
pub fn infer_oeis_rows(
    text: &str,
    id: &str,
    bfile_text: Option<&str>,
) -> Result<OeisInferredRows, OeisSourceError> {
    let mut record = parse_oeis_internal_format(text, id)?;
    let mut bfile_rejection = None;
    if let Some(text) = bfile_text {
        let valid_shape = text.lines().all(|line| {
            let line = line.trim();
            line.is_empty() || line.starts_with('#') || line.split_whitespace().count() == 2
        });
        let parsed = if valid_shape {
            parse_bfile(text)
        } else {
            Err(OeisSourceError::EmptyBFile)
        };
        match parsed {
            Ok(bfile) if bfile.terms.starts_with(&record.terms) => {
                record.terms = bfile.terms;
                record.terms_source = OeisTermsSource::BFile;
            }
            Ok(_) => bfile_rejection = Some(OeisBFileRejection::PrefixMismatch),
            Err(OeisSourceError::NonConsecutiveBFileIndex { .. }) => {
                bfile_rejection = Some(OeisBFileRejection::NonConsecutiveIndices)
            }
            Err(_) => bfile_rejection = Some(OeisBFileRejection::InvalidBFile),
        }
    }
    let mut result = OeisInferredRows {
        name: record.name.clone(),
        keywords: record.keywords.clone(),
        term_source: record.terms_source,
        row_length_rule: None,
        rows: Vec::new(),
        reason: None,
        bfile_rejection,
    };
    let (first, rule) = if record.has_keyword("tabl") {
        (1, OeisRowLengthRule::ConstantIncrease(1))
    } else if record.has_keyword("tabf") {
        let lines: Vec<_> = text.lines().collect();
        let mut best: Vec<Vec<BigInt>> = Vec::new();
        let mut saw_table = false;
        let mut starts_later = false;
        let mut disagreement = false;
        let mut best_mismatch = false;
        for start in 0..lines.len() {
            let Some(first_row) = table_row(lines[start]) else {
                continue;
            };
            saw_table = true;
            if !record.terms.starts_with(&first_row) {
                if record
                    .terms
                    .windows(first_row.len())
                    .skip(1)
                    .any(|window| window == first_row)
                {
                    starts_later = true;
                }
                continue;
            }
            let mut candidate = Vec::new();
            let mut offset = 0;
            let mut mismatch = false;
            for line in &lines[start..] {
                let Some(row) = table_row(line) else {
                    break;
                };
                let tail = &record.terms[offset..];
                if !tail.starts_with(&row) {
                    if !row.starts_with(tail) {
                        mismatch = true;
                    }
                    break;
                }
                offset += row.len();
                candidate.push(row);
                if offset == record.terms.len() {
                    break;
                }
            }
            if mismatch {
                disagreement = true;
            }
            if candidate.len() > best.len()
                || (candidate.len() == best.len() && best_mismatch && !mismatch)
            {
                best = candidate;
                best_mismatch = mismatch;
            }
        }
        if best_mismatch {
            result.reason = Some(OeisRowsRefusalReason::TermsDisagreeWithTable);
            return Ok(result);
        }
        if best.is_empty() {
            result.reason = Some(if disagreement {
                OeisRowsRefusalReason::TermsDisagreeWithTable
            } else if starts_later {
                OeisRowsRefusalReason::TableDoesNotStartAtFirstTerm
            } else if saw_table {
                OeisRowsRefusalReason::TermsDisagreeWithTable
            } else {
                OeisRowsRefusalReason::NoExampleTable
            });
            return Ok(result);
        }
        let lengths: Vec<_> = best.iter().map(Vec::len).collect();
        if lengths.len() < 3 {
            result.reason = Some(OeisRowsRefusalReason::TooFewExampleRows);
            return Ok(result);
        }
        let Some(rule) = regular_rule(&lengths) else {
            let differences: Vec<_> = lengths.windows(2).map(|w| w[1].checked_sub(w[0])).collect();
            let periodic_candidate = differences.iter().all(|d| matches!(d, Some(0 | 1)))
                && differences.iter().filter(|d| **d == Some(1)).count() < 2;
            result.reason = Some(if periodic_candidate {
                OeisRowsRefusalReason::TooFewExampleRows
            } else {
                OeisRowsRefusalReason::IrregularLengths
            });
            return Ok(result);
        };
        (lengths[0], rule)
    } else {
        result.reason = Some(OeisRowsRefusalReason::NotTriangleKeyword);
        return Ok(result);
    };
    let mut width = first;
    let mut offset = 0usize;
    let mut transition = 0usize;
    while width <= record.terms.len().saturating_sub(offset) {
        result
            .rows
            .push(record.terms[offset..offset + width].to_vec());
        offset += width;
        let increase = match rule {
            OeisRowLengthRule::ConstantIncrease(c) => c,
            OeisRowLengthRule::PeriodicIncrease { period, phase } => {
                usize::from(transition % period == phase)
            }
        };
        let Some(next) = width.checked_add(increase) else {
            break;
        };
        width = next;
        transition += 1;
    }
    result.row_length_rule = Some(rule);
    if result.rows.is_empty() {
        result.reason = Some(OeisRowsRefusalReason::NoCompleteRows);
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn record(keyword: &str, terms: &str, examples: &str) -> String {
        format!("%N A000001 Example\n%O A000001 0,1\n%K A000001 {keyword}\n%S A000001 {terms}\n{examples}")
    }
    #[test]
    fn sources_and_triangle() {
        let text = record("tabl", "1,2,3,4,5,6", "");
        let run = |b| infer_oeis_rows(&text, "A000001", b).unwrap();
        assert_eq!(run(None).rows.len(), 3);
        assert_eq!(
            run(Some("0 1\n1 2\n2 3\n3 4\n4 5\n5 6\n6 7\n")).term_source,
            OeisTermsSource::BFile
        );
        assert_eq!(
            run(Some("0 1\n2 2\n")).term_source,
            OeisTermsSource::InternalFormat
        );
        assert_eq!(
            run(Some("0 9\n1 2\n")).term_source,
            OeisTermsSource::InternalFormat
        );
        assert_eq!(
            run(Some("0 1 extra\n")).term_source,
            OeisTermsSource::InternalFormat
        );
        assert_eq!(
            infer_oeis_rows(&record("nonn", "1,2", ""), "A000001", None)
                .unwrap()
                .reason,
            Some(OeisRowsRefusalReason::NotTriangleKeyword)
        );
    }
    #[test]
    fn table_growth_and_refusals() {
        let constant = record(
            "tabf",
            "1,2,3,4,5,6,7,8,9",
            "%e A000001 1\n%e A000001 2,3,4\n%e A000001 5,6,7,8,9\n",
        );
        assert_eq!(
            infer_oeis_rows(&constant, "A000001", None)
                .unwrap()
                .row_length_rule,
            Some(OeisRowLengthRule::ConstantIncrease(2))
        );
        let periodic = record(
            "tabf",
            "1,2,3,4,5,6,7,8,9",
            "%e A000001 1\n%e A000001 2\n%e A000001 3,4\n%e A000001 5,6\n%e A000001 7,8,9\n",
        );
        assert_eq!(
            infer_oeis_rows(&periodic, "A000001", None)
                .unwrap()
                .row_length_rule,
            Some(OeisRowLengthRule::PeriodicIncrease {
                period: 2,
                phase: 1
            })
        );
        let irregular = record(
            "tabf",
            "1,2,3,4,5,6,7",
            "%e A000001 1\n%e A000001 2,3\n%e A000001 4,5,6,7\n",
        );
        assert_eq!(
            infer_oeis_rows(&irregular, "A000001", None).unwrap().reason,
            Some(OeisRowsRefusalReason::IrregularLengths)
        );
        assert_eq!(
            infer_oeis_rows(&record("tabf", "1,2", ""), "A000001", None)
                .unwrap()
                .reason,
            Some(OeisRowsRefusalReason::NoExampleTable)
        );
    }
    #[test]
    fn examples_accept_labels_and_terminal_punctuation_only() {
        for row in [
            "1: 4, 5, 6.",
            "n=1: 4, 5, 6;,.",
            "n=1 4, 5, 6.",
            "1= 4, 5, 6.",
            "1 | 4, 5, 6;",
            "[1] 4, 5, 6,",
        ] {
            assert_eq!(
                table_row(&format!("%e A000001 {row}")),
                Some(vec![4.into(), 5.into(), 6.into()])
            );
        }
        assert!(table_row("%e A000001 [1] 4, 5, 6. trailing text").is_none());
        assert!(table_row("%e A000001 4, 5x, 6.").is_none());
    }

    #[test]
    fn refusals_are_distinct_and_periods_need_evidence() {
        assert_eq!(
            serde_json::to_string(&OeisRowsRefusalReason::TooFewExampleRows).unwrap(),
            "\"too_few_example_rows\""
        );
        let run = |terms, examples| {
            infer_oeis_rows(&record("tabf", terms, examples), "A000001", None)
                .unwrap()
                .reason
        };
        assert_eq!(
            run("1,2,3", ""),
            Some(OeisRowsRefusalReason::NoExampleTable)
        );
        assert_eq!(
            run("1,2,3", "%e A000001 2,3\n"),
            Some(OeisRowsRefusalReason::TableDoesNotStartAtFirstTerm)
        );
        assert_eq!(
            run("1,2,3", "%e A000001 1\n%e A000001 2,3\n"),
            Some(OeisRowsRefusalReason::TooFewExampleRows)
        );
        assert_eq!(
            run(
                "1,2,3,4,5,6",
                "%e A000001 1\n%e A000001 2\n%e A000001 3,4\n%e A000001 5,6\n"
            ),
            Some(OeisRowsRefusalReason::TooFewExampleRows)
        );
        assert_eq!(
            run(
                "1,2,3,4,5,6",
                "%e A000001 1\n%e A000001 2,3\n%e A000001 4,9,6\n"
            ),
            Some(OeisRowsRefusalReason::TermsDisagreeWithTable)
        );
    }

    #[test]
    fn bfile_diagnostics_and_blank_crlf_lines() {
        let text = record("tabl", "1,2,3", "");
        let run = |b| infer_oeis_rows(&text, "A000001", Some(b)).unwrap();
        assert_eq!(
            run("# comment\r\n \r\n0 1\r\n1 2\r\n2 3\r\n3 4\r\n").term_source,
            OeisTermsSource::BFile
        );
        assert_eq!(
            run("version https://git-lfs.github.com/spec/v1\n").bfile_rejection,
            Some(OeisBFileRejection::InvalidBFile)
        );
        assert_eq!(
            run("0 1\n2 2\n").bfile_rejection,
            Some(OeisBFileRejection::NonConsecutiveIndices)
        );
        assert_eq!(
            run("0 9\n1 2\n2 3\n").bfile_rejection,
            Some(OeisBFileRejection::PrefixMismatch)
        );
        assert_eq!(
            run("0 1 extra\n").bfile_rejection,
            Some(OeisBFileRejection::InvalidBFile)
        );
    }
}
