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
pub(crate) use polytool::oeis::OeisSourceRecord;
use polytool::oeis::{parse_bfile, parse_oeis_internal_format, OeisTermsSource};
use std::path::{Path, PathBuf};

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

/// Read `id` from a local internal-format copy, preferring a usable b-file.
pub(crate) fn read_local_record(root: &Path, id: &str) -> Result<OeisSourceRecord, String> {
    let path = sequence_file(root, id)
        .ok_or_else(|| format!("{id}: no local OEIS record under {}", root.display()))?;
    let text = std::fs::read(&path)
        .map_err(|error| format!("{id}: failed to read {}: {error}", path.display()))?;
    let mut record = parse_oeis_internal_format(&String::from_utf8_lossy(&text), id)
        .map_err(|error| error.to_string())?;
    if let Some(bfile) = local_bfile(root, id) {
        if let Ok(text) = std::fs::read_to_string(&bfile) {
            if let Ok(bfile) = parse_bfile(&text) {
                if bfile.first_index == record.offset
                    && bfile.terms.len() > record.terms.len()
                    && bfile.terms.starts_with(&record.terms)
                {
                    record.terms = bfile.terms;
                    record.terms_source = OeisTermsSource::BFile;
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
    fn internal_format_record_splits_into_rows() {
        let record = parse_oeis_internal_format(RECORD, "A000001").unwrap();
        let (rows, rest) = RowWidths::TRIANGLE.split(&record.terms);
        assert_eq!(rows.len(), 4);
        assert_eq!(rest, 0);
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
