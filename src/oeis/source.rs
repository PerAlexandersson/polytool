//! Reading OEIS terms from the plain-text formats of the OEIS.
//!
//! [`parse_oeis_internal_format`] reads one record of the OEIS internal format
//! (the `%N`, `%O`, `%K` and `%S/%T/%U` lines of an `A123456.seq` file) and
//! [`parse_bfile`] reads a b-file (`n a(n)` lines).  Both are independent of
//! the bundled catalog and of the file system; they back
//! `polytool oeis add` and `polytool oeis check`.
//! [`parse_oeis_internal_terms`] preserves term spelling for sequence indexes
//! and uses the separate strict-term, optional-metadata policy documented there.
//!
//! ```
//! use polytool::oeis::{parse_bfile, parse_oeis_internal_format};
//!
//! let record = parse_oeis_internal_format(
//!     "%S A000001 1,1,1,1,4,1\n%N A000001 Example triangle.\n\
//!      %K A000001 nonn,tabl\n%O A000001 1,5\n",
//!     "A000001",
//! )
//! .unwrap();
//! assert_eq!(record.name, "Example triangle.");
//! assert_eq!(record.offset, 1);
//! assert!(record.has_keyword("tabl"));
//! assert_eq!(record.terms.len(), 6);
//!
//! let bfile = parse_bfile("# comment\n0 1\n1 -2\n").unwrap();
//! assert_eq!(bfile.first_index, 0);
//! assert_eq!(bfile.terms.len(), 2);
//! assert!(parse_bfile("0 1\n2 1\n").is_err());
//! ```

use num_bigint::BigInt;
use std::fmt;
use std::path::Path;

/// Read a local b-file, optionally preferring a separate b-file directory.
///
/// For an A-number such as `A123456`, look first in
/// `bfile_dir/A123/b123456.txt`, then in `data_dir/files/A123/b123456.txt`.
/// Only an absent file falls back; an unreadable override returns its I/O
/// error. `Ok(None)` means neither file exists. Contents are returned unchanged,
/// including Git LFS pointers; callers still apply their usual parsing and
/// term-agreement rules (for example, [`super::infer_oeis_rows`]).
///
/// # Errors
///
/// Returns `InvalidInput` for an identifier other than `A` and six ASCII digits,
/// or the I/O error if an existing file cannot be read as UTF-8.
pub fn read_oeis_bfile(
    data_dir: &Path,
    id: &str,
    bfile_dir: Option<&Path>,
) -> std::io::Result<Option<String>> {
    if id.len() != 7 || !id.starts_with('A') || !id[1..].bytes().all(|b| b.is_ascii_digit()) {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "expected an A-number with six digits",
        ));
    }
    let relative = Path::new(&id[..4]).join(format!("b{}.txt", &id[1..]));
    let fallback = data_dir.join("files");
    for directory in bfile_dir
        .into_iter()
        .chain(std::iter::once(fallback.as_path()))
    {
        match std::fs::read_to_string(directory.join(&relative)) {
            Ok(text) => return Ok(Some(text)),
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
            Err(error) => return Err(error),
        }
    }
    Ok(None)
}

/// Where the terms of an [`OeisSourceRecord`] come from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OeisTermsSource {
    /// The `%S/%T/%U` lines of the internal format.
    InternalFormat,
    /// A b-file.
    BFile,
}

impl OeisTermsSource {
    /// The label recorded for the source: `oeis_internal_format` or `bfile`.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::InternalFormat => "oeis_internal_format",
            Self::BFile => "bfile",
        }
    }
}

impl fmt::Display for OeisTermsSource {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Terms and metadata of one OEIS sequence.
#[derive(Debug, Clone)]
pub struct OeisSourceRecord {
    /// OEIS title (`%N`).
    pub name: String,
    /// First index of the flattened sequence (first value of `%O`).
    pub offset: i64,
    /// OEIS keywords (`%K`), such as `tabl` or `tabf`.
    pub keywords: Vec<String>,
    /// Flattened terms starting at `offset`.
    pub terms: Vec<BigInt>,
    /// Where the terms were read from.
    pub terms_source: OeisTermsSource,
}

impl OeisSourceRecord {
    /// Whether `keyword` is among the OEIS keywords of the record.
    #[must_use]
    pub fn has_keyword(&self, keyword: &str) -> bool {
        self.keywords.iter().any(|candidate| candidate == keyword)
    }
}

/// The terms of a b-file with the index of the first one.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParsedBFile {
    /// Index of the first term.
    pub first_index: i64,
    /// The terms, with consecutive indices from `first_index`.
    pub terms: Vec<BigInt>,
}

/// Internal-format terms and metadata for a text-preserving sequence index.
///
/// Unlike [`OeisSourceRecord`], this needs neither a title nor an offset.
/// Terms retain their spelling (including leading zeros and `-0`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OeisInternalTerms {
    /// First value of the last `%O` line, or zero if missing or invalid.
    pub offset: i64,
    /// Comma-separated keywords from the last `%K` line, trimmed.
    pub keywords: Vec<String>,
    /// Nonempty, trimmed `%S/%T/%U` tokens in their input order.
    pub terms: Vec<String>,
}

/// Error of the internal-format and b-file parsers in this module.
///
/// Line numbers are 1-based.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OeisSourceError {
    /// An internal-format term is not ASCII digits with an optional minus.
    InvalidInternalTerm {
        /// The trimmed token.
        term: String,
    },
    /// The internal-format record contains no nonempty terms.
    EmptyInternalTerms,
    /// The first value of a `%O` line is not an integer.
    InvalidOffsetLine {
        /// The sequence identifier given to the parser.
        id: String,
        /// The offending line.
        line: String,
    },
    /// The record has no `%O` line.
    MissingOffset {
        /// The sequence identifier given to the parser.
        id: String,
    },
    /// The record has no (or an empty) `%N` line.
    MissingName {
        /// The sequence identifier given to the parser.
        id: String,
    },
    /// The text is a Git LFS pointer, not b-file data.
    GitLfsPointer,
    /// A b-file line does not have the form `n a(n)`.
    MalformedBFileLine {
        /// Line number.
        line: usize,
    },
    /// The index of a b-file line is not an integer.
    InvalidBFileIndex {
        /// Line number.
        line: usize,
    },
    /// The term of a b-file line is not an integer.
    InvalidBFileTerm {
        /// Line number.
        line: usize,
    },
    /// The index of a b-file line is not the successor of the previous one.
    NonConsecutiveBFileIndex {
        /// Line number.
        line: usize,
        /// The index found.
        index: i64,
        /// The index of the previous line.
        previous: i64,
    },
    /// The b-file has no terms.
    EmptyBFile,
}

impl fmt::Display for OeisSourceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidInternalTerm { term } => {
                write!(f, "invalid internal-format term `{term}`")
            }
            Self::EmptyInternalTerms => f.write_str("internal-format record has no terms"),
            Self::InvalidOffsetLine { id, line } => {
                write!(f, "{id}: invalid offset line `{line}`")
            }
            Self::MissingOffset { id } => write!(f, "{id}: missing %O offset"),
            Self::MissingName { id } => write!(f, "{id}: missing %N name"),
            Self::GitLfsPointer => write!(f, "the b-file is a Git LFS pointer, not data"),
            Self::MalformedBFileLine { line } => {
                write!(f, "b-file line {line}: expected `n a(n)`")
            }
            Self::InvalidBFileIndex { line } => write!(f, "b-file line {line}: invalid index"),
            Self::InvalidBFileTerm { line } => write!(f, "b-file line {line}: invalid term"),
            Self::NonConsecutiveBFileIndex {
                line,
                index,
                previous,
            } => write!(
                f,
                "b-file line {line}: index {index} does not follow {previous}"
            ),
            Self::EmptyBFile => write!(f, "b-file has no terms"),
        }
    }
}

impl std::error::Error for OeisSourceError {}

fn parse_term(token: &str) -> Option<BigInt> {
    let token = token.trim();
    let digits = token.strip_prefix(['-', '+']).unwrap_or(token);
    if digits.is_empty() || !digits.bytes().all(|byte| byte.is_ascii_digit()) {
        return None;
    }
    token.parse().ok()
}

fn internal_format_fields(text: &str) -> impl Iterator<Item = (&str, &str, &str)> {
    text.lines().filter_map(|line| {
        let mut parts = line.splitn(3, ' ');
        let (Some(tag), Some(_), Some(rest)) = (parts.next(), parts.next(), parts.next()) else {
            return None;
        };
        Some((line, tag, rest))
    })
}

/// Parse internal-format terms for a text-preserving sequence index.
///
/// This stricter term policy rejects the entire record if a nonempty term is
/// not ASCII digits with an optional leading minus; a leading plus is rejected.
/// Empty tokens are skipped. Missing or invalid offsets become zero, the last
/// `%K` line supplies the keywords, and a `%N` title is not required.
/// Use [`parse_oeis_internal_format`] for a titled record with a required offset,
/// numeric terms, and invalid term tokens skipped instead.
///
/// # Errors
///
/// [`OeisSourceError::InvalidInternalTerm`] for an invalid nonempty term, or
/// [`OeisSourceError::EmptyInternalTerms`] if there are no terms.
pub fn parse_oeis_internal_terms(text: &str) -> Result<OeisInternalTerms, OeisSourceError> {
    let mut offset = None;
    let mut keywords = Vec::new();
    let mut terms = Vec::new();
    for (_, tag, rest) in internal_format_fields(text) {
        match tag {
            "%S" | "%T" | "%U" => terms.extend(
                rest.split(',')
                    .map(str::trim)
                    .filter(|term| !term.is_empty())
                    .map(str::to_string),
            ),
            "%O" => offset = rest.split(',').next().and_then(|v| v.trim().parse().ok()),
            "%K" => keywords = rest.split(',').map(|k| k.trim().to_string()).collect(),
            _ => {}
        }
    }
    for term in &terms {
        let digits = term.strip_prefix('-').unwrap_or(term);
        if digits.is_empty() || !digits.bytes().all(|byte| byte.is_ascii_digit()) {
            return Err(OeisSourceError::InvalidInternalTerm { term: term.clone() });
        }
    }
    if terms.is_empty() {
        return Err(OeisSourceError::EmptyInternalTerms);
    }
    Ok(OeisInternalTerms {
        offset: offset.unwrap_or(0),
        keywords,
        terms,
    })
}

/// Parse one record in the OEIS internal format.
///
/// The `%N` line gives the name, the first value of `%O` the offset, `%K` the
/// keywords and the `%S`, `%T` and `%U` lines the terms (tokens that are not
/// integers are skipped); other lines are ignored.  `id` only labels the
/// errors.
///
/// # Errors
///
/// [`OeisSourceError::InvalidOffsetLine`], [`OeisSourceError::MissingOffset`]
/// and [`OeisSourceError::MissingName`].
pub fn parse_oeis_internal_format(
    text: &str,
    id: &str,
) -> Result<OeisSourceRecord, OeisSourceError> {
    let mut name = String::new();
    let mut offset = None;
    let mut keywords = Vec::new();
    let mut terms = Vec::new();
    for (line, tag, rest) in internal_format_fields(text) {
        match tag {
            "%N" => name = rest.trim().to_string(),
            "%O" => {
                let first = rest.split(',').next().unwrap_or("").trim();
                offset =
                    Some(
                        first
                            .parse::<i64>()
                            .map_err(|_| OeisSourceError::InvalidOffsetLine {
                                id: id.to_string(),
                                line: line.to_string(),
                            })?,
                    );
            }
            "%K" => keywords.extend(rest.split(',').map(|keyword| keyword.trim().to_string())),
            "%S" | "%T" | "%U" => terms.extend(rest.split(',').filter_map(parse_term)),
            _ => {}
        }
    }
    let Some(offset) = offset else {
        return Err(OeisSourceError::MissingOffset { id: id.to_string() });
    };
    if name.is_empty() {
        return Err(OeisSourceError::MissingName { id: id.to_string() });
    }
    Ok(OeisSourceRecord {
        name,
        offset,
        keywords,
        terms,
        terms_source: OeisTermsSource::InternalFormat,
    })
}

/// Parse a b-file: `n a(n)` lines with consecutive indices; `#` comments and
/// blank lines are ignored.
///
/// # Errors
///
/// [`OeisSourceError::GitLfsPointer`] for a Git LFS pointer file, a line error
/// ([`OeisSourceError::MalformedBFileLine`],
/// [`OeisSourceError::InvalidBFileIndex`], [`OeisSourceError::InvalidBFileTerm`],
/// [`OeisSourceError::NonConsecutiveBFileIndex`]) and
/// [`OeisSourceError::EmptyBFile`].
pub fn parse_bfile(text: &str) -> Result<ParsedBFile, OeisSourceError> {
    if text.starts_with("version https://git-lfs") {
        return Err(OeisSourceError::GitLfsPointer);
    }
    let mut first = None;
    let mut terms = Vec::new();
    for (line_index, line) in text.lines().enumerate() {
        let line_number = line_index + 1;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut parts = line.split_whitespace();
        let (Some(index), Some(value)) = (parts.next(), parts.next()) else {
            return Err(OeisSourceError::MalformedBFileLine { line: line_number });
        };
        let index = index
            .parse::<i64>()
            .map_err(|_| OeisSourceError::InvalidBFileIndex { line: line_number })?;
        let value =
            parse_term(value).ok_or(OeisSourceError::InvalidBFileTerm { line: line_number })?;
        let expected = first.map(|first: i64| first + terms.len() as i64);
        if let Some(expected) = expected {
            if index != expected {
                return Err(OeisSourceError::NonConsecutiveBFileIndex {
                    line: line_number,
                    index,
                    previous: expected - 1,
                });
            }
        } else {
            first = Some(index);
        }
        terms.push(value);
    }
    let first_index = first.ok_or(OeisSourceError::EmptyBFile)?;
    Ok(ParsedBFile { first_index, terms })
}

#[cfg(test)]
mod tests {
    use super::*;

    const RECORD: &str = "%I A000001 #1 Jan 01 2000\n\
        %S A000001 1,1,1,1,4,1,\n\
        %T A000001 1,11,11,1\n\
        %N A000001 Example triangle.\n\
        %K A000001 nonn,tabl\n\
        %O A000001 1,5\n";

    #[test]
    fn indexed_terms_preserve_spelling_and_optional_metadata() {
        let record = parse_oeis_internal_terms(
            "%S A1 001, -0,, -0002,\n%T A1 3,\n%U A1 4\n\
             %O A1 -2,5\n%O A1 bad\n%K A1 old\n%K A1 nonn, tabl,\n",
        )
        .unwrap();
        assert_eq!(record.offset, 0);
        assert_eq!(record.keywords, ["nonn", "tabl", ""]);
        assert_eq!(record.terms, ["001", "-0", "-0002", "3", "4"]);
        assert_eq!(parse_oeis_internal_terms("%S A1 1\n").unwrap().offset, 0);
        assert_eq!(parse_oeis_internal_terms(RECORD).unwrap().offset, 1);
        for invalid in ["+", "+1", "-", "--1", "1 2", "x", "١"] {
            assert_eq!(
                parse_oeis_internal_terms(&format!("%S A1 1,{invalid},2\n")),
                Err(OeisSourceError::InvalidInternalTerm {
                    term: invalid.into()
                }),
            );
        }
        assert_eq!(
            parse_oeis_internal_terms("%S A1 , ,\n%N A1 title\n"),
            Err(OeisSourceError::EmptyInternalTerms),
        );
    }

    #[test]
    fn titled_records_keep_their_distinct_term_and_metadata_policy() {
        let record = parse_oeis_internal_format(
            "%N A1 title\n%O A1 2\n%K A1 one\n%K A1 two\n%S A1 001,+2,x,-0\n",
            "A1",
        )
        .unwrap();
        assert_eq!(record.offset, 2);
        assert_eq!(record.keywords, ["one", "two"]);
        assert_eq!(record.terms, [1, 2, 0].map(BigInt::from));
    }

    #[test]
    fn internal_format_record_is_parsed() {
        let record = parse_oeis_internal_format(RECORD, "A000001").unwrap();
        assert_eq!(record.name, "Example triangle.");
        assert_eq!(record.offset, 1);
        assert!(record.has_keyword("tabl"));
        assert_eq!(record.terms.len(), 10);
        assert_eq!(record.terms_source.as_str(), "oeis_internal_format");
    }

    #[test]
    fn internal_format_errors_keep_their_messages() {
        let error = |text: &str| {
            parse_oeis_internal_format(text, "A1")
                .unwrap_err()
                .to_string()
        };
        assert_eq!(error("%N A1 x\n"), "A1: missing %O offset");
        assert_eq!(error("%O A1 0\n"), "A1: missing %N name");
        assert_eq!(
            error("%N A1 x\n%O A1 z,1\n"),
            "A1: invalid offset line `%O A1 z,1`"
        );
    }

    #[test]
    fn bfiles_need_consecutive_indices() {
        assert_eq!(
            parse_bfile("# c\n0 1\n1 -2\n").unwrap(),
            ParsedBFile {
                first_index: 0,
                terms: vec![BigInt::from(1), BigInt::from(-2)]
            }
        );
        let error = |text: &str| parse_bfile(text).unwrap_err();
        assert_eq!(
            error("0 1\n2 1\n"),
            OeisSourceError::NonConsecutiveBFileIndex {
                line: 2,
                index: 2,
                previous: 0
            }
        );
        assert_eq!(
            error("0 1\n2 1\n").to_string(),
            "b-file line 2: index 2 does not follow 0"
        );
        assert_eq!(
            error("version https://git-lfs.github.com/spec/v1\n").to_string(),
            "the b-file is a Git LFS pointer, not data"
        );
        assert_eq!(error("0\n").to_string(), "b-file line 1: expected `n a(n)`");
        assert_eq!(error("x 1\n").to_string(), "b-file line 1: invalid index");
        assert_eq!(error("0 y\n").to_string(), "b-file line 1: invalid term");
        assert_eq!(error("# nothing\n").to_string(), "b-file has no terms");
    }
}
