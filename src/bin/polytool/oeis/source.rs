//! Reading OEIS terms for `polytool oeis add` and `polytool oeis check`.
//!
//! Two inputs are supported:
//!
//! - a local copy of the OEIS in its internal format (`--oeis-data DIR` or the
//!   `OEIS_DATA_DIR` environment variable): one `A123456.seq` file per
//!   sequence under `DIR/seq/A123/` (or directly under `DIR/A123/`), whose
//!   `%N`, `%O`, `%K` and `%S/%T/%U` lines give the name, offset, keywords and
//!   terms.  When the copy also holds the b-file (`DIR/files/A123/b123456.txt`)
//!   and the b-file agrees with the `%S/%T/%U` terms, its longer term list is
//!   used instead;
//! - an explicit b-file (`--bfile FILE`, lines `n a(n)`), with the name,
//!   offset and layout given on the command line.

use num_bigint::BigInt;
use std::path::{Path, PathBuf};

/// Terms and metadata of one OEIS sequence.
#[derive(Debug, Clone)]
pub(crate) struct OeisSourceRecord {
    /// OEIS title (`%N`).
    pub(crate) name: String,
    /// First index of the flattened sequence (first value of `%O`).
    pub(crate) offset: i64,
    /// OEIS keywords (`%K`), such as `tabl` or `tabf`.
    pub(crate) keywords: Vec<String>,
    /// Flattened terms starting at `offset`.
    pub(crate) terms: Vec<BigInt>,
    /// `oeis_internal_format` or `bfile`.
    pub(crate) terms_source: &'static str,
}

impl OeisSourceRecord {
    pub(crate) fn has_keyword(&self, keyword: &str) -> bool {
        self.keywords.iter().any(|value| value == keyword)
    }
}

fn sequence_directory(id: &str) -> &str {
    &id[..4]
}

/// Path of the internal-format file of `id`, if it exists.
pub(crate) fn sequence_file(root: &Path, id: &str) -> Option<PathBuf> {
    let file = format!("{id}.seq");
    [
        root.join("seq").join(sequence_directory(id)).join(&file),
        root.join(sequence_directory(id)).join(&file),
    ]
    .into_iter()
    .find(|path| path.is_file())
}

fn local_bfile(root: &Path, id: &str) -> Option<PathBuf> {
    let file = format!("b{}.txt", &id[1..]);
    [
        root.join("files").join(sequence_directory(id)).join(&file),
        root.join(sequence_directory(id)).join(&file),
    ]
    .into_iter()
    .find(|path| path.is_file())
}

fn parse_term(token: &str) -> Option<BigInt> {
    let token = token.trim();
    let digits = token.strip_prefix(['-', '+']).unwrap_or(token);
    if digits.is_empty() || !digits.bytes().all(|byte| byte.is_ascii_digit()) {
        return None;
    }
    token.parse().ok()
}

/// Parse one record in the OEIS internal format.
pub(crate) fn parse_internal_format(text: &str, id: &str) -> Result<OeisSourceRecord, String> {
    let mut name = String::new();
    let mut offset = None;
    let mut keywords = Vec::new();
    let mut terms = Vec::new();
    for line in text.lines() {
        let mut parts = line.splitn(3, ' ');
        let (Some(tag), Some(_), Some(rest)) = (parts.next(), parts.next(), parts.next()) else {
            continue;
        };
        match tag {
            "%N" => name = rest.trim().to_string(),
            "%O" => {
                let first = rest.split(',').next().unwrap_or("").trim();
                offset = Some(
                    first
                        .parse::<i64>()
                        .map_err(|_| format!("{id}: invalid offset line `{line}`"))?,
                );
            }
            "%K" => keywords.extend(rest.split(',').map(|keyword| keyword.trim().to_string())),
            "%S" | "%T" | "%U" => terms.extend(rest.split(',').filter_map(parse_term)),
            _ => {}
        }
    }
    let Some(offset) = offset else {
        return Err(format!("{id}: missing %O offset"));
    };
    if name.is_empty() {
        return Err(format!("{id}: missing %N name"));
    }
    Ok(OeisSourceRecord {
        name,
        offset,
        keywords,
        terms,
        terms_source: "oeis_internal_format",
    })
}

/// Parse a b-file: `n a(n)` lines with consecutive indices; `#` comments and
/// blank lines are ignored.  Returns the first index and the terms.
pub(crate) fn parse_bfile(text: &str) -> Result<(i64, Vec<BigInt>), String> {
    if text.starts_with("version https://git-lfs") {
        return Err("the b-file is a Git LFS pointer, not data".to_string());
    }
    let mut first = None;
    let mut terms = Vec::new();
    for (line_index, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut parts = line.split_whitespace();
        let (Some(index), Some(value)) = (parts.next(), parts.next()) else {
            return Err(format!("b-file line {}: expected `n a(n)`", line_index + 1));
        };
        let index = index
            .parse::<i64>()
            .map_err(|_| format!("b-file line {}: invalid index", line_index + 1))?;
        let value = parse_term(value)
            .ok_or_else(|| format!("b-file line {}: invalid term", line_index + 1))?;
        let expected = first.map(|first: i64| first + terms.len() as i64);
        if let Some(expected) = expected {
            if index != expected {
                return Err(format!(
                    "b-file line {}: index {index} does not follow {}",
                    line_index + 1,
                    expected - 1
                ));
            }
        } else {
            first = Some(index);
        }
        terms.push(value);
    }
    let first = first.ok_or_else(|| "b-file has no terms".to_string())?;
    Ok((first, terms))
}

/// Read `id` from a local internal-format copy, preferring a usable b-file.
pub(crate) fn read_local_record(root: &Path, id: &str) -> Result<OeisSourceRecord, String> {
    let path = sequence_file(root, id)
        .ok_or_else(|| format!("{id}: no local OEIS record under {}", root.display()))?;
    let text = std::fs::read(&path)
        .map_err(|error| format!("{id}: failed to read {}: {error}", path.display()))?;
    let mut record = parse_internal_format(&String::from_utf8_lossy(&text), id)?;
    if let Some(bfile) = local_bfile(root, id) {
        if let Ok(text) = std::fs::read_to_string(&bfile) {
            if let Ok((first, terms)) = parse_bfile(&text) {
                if first == record.offset
                    && terms.len() > record.terms.len()
                    && terms.starts_with(&record.terms)
                {
                    record.terms = terms;
                    record.terms_source = "bfile";
                }
            }
        }
    }
    Ok(record)
}

/// Row-width rule `width(i) = first + step * floor(i / every)` for row
/// `i = 0, 1, ...`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct RowWidths {
    pub(crate) first: usize,
    pub(crate) step: usize,
    pub(crate) every: usize,
}

impl RowWidths {
    /// The regular triangle (`tabl`): widths 1, 2, 3, ...
    pub(crate) const TRIANGLE: Self = Self {
        first: 1,
        step: 1,
        every: 1,
    };

    /// Parse `W`, `W+S` or `W+S/R`: for example `1+2` for rows of lengths
    /// 1, 3, 5, ... and `1+1/2` for 1, 1, 2, 2, 3, 3, ...
    pub(crate) fn parse(text: &str) -> Result<Self, String> {
        let invalid = || format!("invalid row-width rule `{text}` (expected W, W+S or W+S/R)");
        let (first, rest) = text.split_once('+').unwrap_or((text, "0"));
        let (step, every) = rest.split_once('/').unwrap_or((rest, "1"));
        let number = |part: &str| part.trim().parse::<usize>().map_err(|_| invalid());
        let (first, step, every) = (number(first)?, number(step)?, number(every)?);
        if first == 0 || every == 0 {
            return Err(invalid());
        }
        Ok(Self { first, step, every })
    }

    /// Width of row `index`.
    pub(crate) fn width(&self, index: usize) -> usize {
        self.first + self.step * (index / self.every)
    }

    /// Split terms into complete rows; returns the rows and the number of
    /// trailing terms that do not fill a row.
    pub(crate) fn split<'a>(&self, terms: &'a [BigInt]) -> (Vec<&'a [BigInt]>, usize) {
        let mut rows = Vec::new();
        let mut cursor = 0;
        loop {
            let width = self.width(rows.len());
            if cursor + width > terms.len() {
                return (rows, terms.len() - cursor);
            }
            rows.push(&terms[cursor..cursor + width]);
            cursor += width;
        }
    }
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
    fn internal_format_record_is_parsed() {
        let record = parse_internal_format(RECORD, "A000001").unwrap();
        assert_eq!(record.name, "Example triangle.");
        assert_eq!(record.offset, 1);
        assert!(record.has_keyword("tabl"));
        assert_eq!(record.terms.len(), 10);
        let (rows, rest) = RowWidths::TRIANGLE.split(&record.terms);
        assert_eq!(rows.len(), 4);
        assert_eq!(rest, 0);
    }

    #[test]
    fn bfiles_need_consecutive_indices() {
        assert_eq!(
            parse_bfile("# c\n0 1\n1 -2\n").unwrap(),
            (0, vec![BigInt::from(1), BigInt::from(-2)])
        );
        assert!(parse_bfile("0 1\n2 1\n").is_err());
        assert!(parse_bfile("version https://git-lfs.github.com/spec/v1\n").is_err());
    }

    #[test]
    fn row_width_rules_parse() {
        let widths = |rule: &str| {
            let rule = RowWidths::parse(rule).unwrap();
            (0..6).map(|index| rule.width(index)).collect::<Vec<_>>()
        };
        assert_eq!(widths("1+2"), [1, 3, 5, 7, 9, 11]);
        assert_eq!(widths("3"), [3; 6]);
        assert_eq!(widths("1+1/2"), [1, 1, 2, 2, 3, 3]);
        assert_eq!(RowWidths::parse("1+1").unwrap(), RowWidths::TRIANGLE);
        assert!(RowWidths::parse("0+1").is_err());
        assert!(RowWidths::parse("1+1/0").is_err());
    }
}
