//! OEIS b-file text from rows of integer coefficients.
//!
//! [`flatten_rows_to_bfile`] flattens a list of coefficient rows (a triangle
//! read by rows, or any sequence of polynomial coefficient vectors) into the
//! plain `index value` lines of an OEIS b-file.  It is independent of the
//! bundled catalog, so it works for the rows of any recurrence, and it is the
//! formatter behind `polytool oeis generate --format bfile` and
//! `polytool recurrence-generate --format bfile`.
//!
//! Conventions: rows are coefficient vectors in ascending powers (entry `k` is
//! the coefficient of `t^k`) exactly as stored; [`BFileRowOrientation`] can
//! read each row from the other end, and [`BFileZeroPolicy`] can drop zero
//! coefficients at either end of each row (in stored order, before the
//! orientation is applied).  An optional term cap ([`BFileOptions::max_terms`])
//! never splits a row: a row is either written completely or not at all, and
//! [`BFileStoppingRule`] says whether the cap is a hard upper bound or a target
//! for the nearest whole row.  Indices are consecutive, beginning at
//! [`BFileOptions::first_index`] (the OEIS offset).

use num_bigint::BigInt;
use num_traits::Zero;
use std::fmt;
use std::fmt::Write as _;

/// The order in which the coefficients of each row are written.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BFileRowOrientation {
    /// Ascending powers, as the rows are stored (the default).
    #[default]
    AsStored,
    /// Each row read from its highest stored coefficient down to its lowest.
    Reversed,
}

/// Which zero coefficients are dropped from the ends of each row.
///
/// "Low" is the start of the stored row (the constant coefficient), "high" its
/// end (the leading coefficient and any structural padding).  The policy is
/// applied to the stored row, before [`BFileRowOrientation`].  A row that
/// consists only of zeros becomes empty under any trimming policy and then
/// contributes no terms.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BFileZeroPolicy {
    /// Keep every coefficient, including zeros at both ends (the default; what
    /// `polytool oeis generate --format bfile` writes).
    #[default]
    KeepAll,
    /// Drop zeros at the low (constant-term) end of each row.
    TrimLowOrder,
    /// Drop zeros at the high (leading-coefficient) end of each row.
    TrimHighOrder,
    /// Drop zeros at both ends of each row.
    TrimBothEnds,
}

/// How a term cap ([`BFileOptions::max_terms`]) selects the last row.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BFileStoppingRule {
    /// Write the longest prefix of whole rows with at most the cap many terms
    /// (the default; the behaviour of `polytool oeis generate`).
    #[default]
    AtMostCap,
    /// Write the prefix of whole rows whose term count is nearest to the cap,
    /// which may exceed it; a tie is resolved towards the shorter prefix.
    NearestWholeRow,
}

/// Options of [`flatten_rows_to_bfile`].
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct BFileOptions {
    /// Index (the OEIS offset) of the first written term.
    pub first_index: i64,
    /// Order of the coefficients within each row.
    pub orientation: BFileRowOrientation,
    /// Zero coefficients dropped from the ends of each row.
    pub zero_policy: BFileZeroPolicy,
    /// Cap on the number of terms, counted after the zero policy; `None` writes
    /// every row.  A row is never split.
    pub max_terms: Option<usize>,
    /// How the cap picks the last row; ignored when `max_terms` is `None`.
    pub stopping_rule: BFileStoppingRule,
}

/// Result of [`flatten_rows_to_bfile`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BFileText {
    /// The b-file lines, each `index value` and terminated by a newline.
    pub text: String,
    /// Number of leading rows written completely (rows that the zero policy
    /// empties are counted).
    pub rows_included: usize,
    /// Number of terms (lines) written.
    pub terms_included: usize,
}

/// Error of [`flatten_rows_to_bfile`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BFileError {
    /// A term index does not fit in `i64`.
    IndexOverflow,
}

impl fmt::Display for BFileError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::IndexOverflow => write!(f, "b-file index overflows i64"),
        }
    }
}

impl std::error::Error for BFileError {}

/// The part of `row` that the zero policy keeps (in stored order).
///
/// # Examples
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::oeis::{bfile_row_terms, BFileZeroPolicy};
///
/// let row: Vec<BigInt> = [0, 2, 0, 3, 0].iter().map(|&x| BigInt::from(x)).collect();
/// assert_eq!(bfile_row_terms(&row, BFileZeroPolicy::KeepAll).len(), 5);
/// assert_eq!(bfile_row_terms(&row, BFileZeroPolicy::TrimLowOrder).len(), 4);
/// assert_eq!(bfile_row_terms(&row, BFileZeroPolicy::TrimBothEnds).len(), 3);
/// ```
pub fn bfile_row_terms(row: &[BigInt], zero_policy: BFileZeroPolicy) -> &[BigInt] {
    let trim_low = matches!(
        zero_policy,
        BFileZeroPolicy::TrimLowOrder | BFileZeroPolicy::TrimBothEnds
    );
    let trim_high = matches!(
        zero_policy,
        BFileZeroPolicy::TrimHighOrder | BFileZeroPolicy::TrimBothEnds
    );
    let mut start = 0;
    let mut end = row.len();
    if trim_low {
        while start < end && row[start].is_zero() {
            start += 1;
        }
    }
    if trim_high {
        while end > start && row[end - 1].is_zero() {
            end -= 1;
        }
    }
    &row[start..end]
}

/// Flatten coefficient rows into OEIS b-file text.
///
/// Row `r` contributes the coefficients kept by [`BFileOptions::zero_policy`],
/// read in [`BFileOptions::orientation`], one `index value` line each, with
/// consecutive indices from [`BFileOptions::first_index`].  If
/// [`BFileOptions::max_terms`] is set, the output ends after a whole row chosen
/// by [`BFileOptions::stopping_rule`]; a row is never split.  With
/// [`BFileStoppingRule::AtMostCap`], a first row longer than the cap yields no
/// output (`rows_included == 0`); with [`BFileStoppingRule::NearestWholeRow`] it
/// is written if that is nearer to the cap than writing nothing.
///
/// # Errors
///
/// [`BFileError::IndexOverflow`] if an index does not fit in `i64`.
///
/// # Examples
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::oeis::{flatten_rows_to_bfile, BFileOptions};
///
/// let rows: Vec<Vec<BigInt>> = vec![
///     vec![BigInt::from(1)],
///     vec![BigInt::from(1), BigInt::from(1)],
///     vec![BigInt::from(1), BigInt::from(4), BigInt::from(1)],
/// ];
/// let options = BFileOptions {
///     first_index: 1,
///     max_terms: Some(5),
///     ..BFileOptions::default()
/// };
/// let file = flatten_rows_to_bfile(&rows, &options).unwrap();
/// // Three terms fit; the third row would cross the cap and is dropped whole.
/// assert_eq!(file.text, "1 1\n2 1\n3 1\n");
/// assert_eq!(file.rows_included, 2);
/// ```
pub fn flatten_rows_to_bfile(
    rows: &[Vec<BigInt>],
    options: &BFileOptions,
) -> Result<BFileText, BFileError> {
    let mut text = String::new();
    let mut index = Some(options.first_index);
    let mut rows_included = 0usize;
    let mut terms_included = 0usize;
    for row in rows {
        let kept = bfile_row_terms(row, options.zero_policy);
        let mut last_row = false;
        if let Some(cap) = options.max_terms {
            if terms_included.saturating_add(kept.len()) > cap {
                let overshoot = terms_included.saturating_add(kept.len()) - cap;
                let undershoot = cap.saturating_sub(terms_included);
                if options.stopping_rule == BFileStoppingRule::NearestWholeRow
                    && overshoot < undershoot
                {
                    last_row = true;
                } else {
                    break;
                }
            }
        }
        let mut write_term = |value: &BigInt| -> Result<(), BFileError> {
            let current = index.ok_or(BFileError::IndexOverflow)?;
            writeln!(text, "{current} {value}").expect("writing to a String cannot fail");
            index = current.checked_add(1);
            Ok(())
        };
        match options.orientation {
            BFileRowOrientation::AsStored => kept.iter().try_for_each(&mut write_term)?,
            BFileRowOrientation::Reversed => kept.iter().rev().try_for_each(&mut write_term)?,
        }
        terms_included += kept.len();
        rows_included += 1;
        if last_row {
            break;
        }
    }
    Ok(BFileText {
        text,
        rows_included,
        terms_included,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rows(data: &[&[i64]]) -> Vec<Vec<BigInt>> {
        data.iter()
            .map(|row| row.iter().map(|&x| BigInt::from(x)).collect())
            .collect()
    }

    fn lines(file: &BFileText) -> Vec<&str> {
        file.text.lines().collect()
    }

    #[test]
    fn writes_consecutive_indices_from_the_offset() {
        let data = rows(&[&[1], &[1, 1], &[1, 4, 1]]);
        let options = BFileOptions {
            first_index: 7,
            ..BFileOptions::default()
        };
        let file = flatten_rows_to_bfile(&data, &options).unwrap();
        assert_eq!(lines(&file), ["7 1", "8 1", "9 1", "10 1", "11 4", "12 1"]);
        assert_eq!((file.rows_included, file.terms_included), (3, 6));
    }

    #[test]
    fn negative_offset_and_values() {
        let data = rows(&[&[-3, 0]]);
        let options = BFileOptions {
            first_index: -1,
            ..BFileOptions::default()
        };
        assert_eq!(
            flatten_rows_to_bfile(&data, &options).unwrap().text,
            "-1 -3\n0 0\n"
        );
    }

    #[test]
    fn reversed_orientation_reverses_each_row_only() {
        let data = rows(&[&[1, 2, 3], &[4, 5]]);
        let options = BFileOptions {
            first_index: 0,
            orientation: BFileRowOrientation::Reversed,
            ..BFileOptions::default()
        };
        let file = flatten_rows_to_bfile(&data, &options).unwrap();
        assert_eq!(lines(&file), ["0 3", "1 2", "2 1", "3 5", "4 4"]);
    }

    #[test]
    fn at_most_cap_never_splits_a_row() {
        let data = rows(&[&[1], &[1, 1], &[1, 2, 1], &[1, 3, 3, 1]]);
        for (cap, expected_rows, expected_terms) in [
            (0, 0, 0),
            (1, 1, 1),
            (2, 1, 1),
            (3, 2, 3),
            (5, 2, 3),
            (6, 3, 6),
            (9, 3, 6),
            (10, 4, 10),
            (100, 4, 10),
        ] {
            let options = BFileOptions {
                first_index: 1,
                max_terms: Some(cap),
                ..BFileOptions::default()
            };
            let file = flatten_rows_to_bfile(&data, &options).unwrap();
            assert_eq!(file.rows_included, expected_rows, "cap {cap}");
            assert_eq!(file.terms_included, expected_terms, "cap {cap}");
            assert_eq!(file.text.lines().count(), expected_terms, "cap {cap}");
        }
    }

    #[test]
    fn nearest_whole_row_may_exceed_the_cap_and_breaks_ties_down() {
        let data = rows(&[&[1], &[1, 1], &[1, 2, 1], &[1, 3, 3, 1]]);
        // Prefix term counts are 0, 1, 3, 6, 10.
        for (cap, expected_terms) in [
            (0, 0),
            (1, 1),
            (2, 1), // 1 and 3 are equally near: the shorter prefix wins
            (3, 3),
            (4, 3),
            (5, 6), // 6 is nearer than 3
            (6, 6),
            (7, 6),
            (8, 6), // 6 and 10 are equally near: the shorter prefix wins
            (9, 10),
            (100, 10),
        ] {
            let options = BFileOptions {
                max_terms: Some(cap),
                stopping_rule: BFileStoppingRule::NearestWholeRow,
                ..BFileOptions::default()
            };
            let file = flatten_rows_to_bfile(&data, &options).unwrap();
            assert_eq!(file.terms_included, expected_terms, "cap {cap}");
            assert_eq!(file.text.lines().count(), expected_terms, "cap {cap}");
        }
    }

    #[test]
    fn nearest_rule_can_write_a_first_row_longer_than_the_cap() {
        let data = rows(&[&[1, 2, 3, 4]]);
        let nearest = BFileOptions {
            max_terms: Some(3),
            stopping_rule: BFileStoppingRule::NearestWholeRow,
            ..BFileOptions::default()
        };
        assert_eq!(
            flatten_rows_to_bfile(&data, &nearest)
                .unwrap()
                .rows_included,
            1
        );
        let at_most = BFileOptions {
            max_terms: Some(3),
            ..BFileOptions::default()
        };
        assert_eq!(
            flatten_rows_to_bfile(&data, &at_most)
                .unwrap()
                .rows_included,
            0
        );
    }

    #[test]
    fn zero_policy_trims_stored_ends_before_the_orientation() {
        let data = rows(&[&[0, 0, 5, 0, 6, 0], &[0, 0]]);
        let write = |zero_policy, orientation| {
            let options = BFileOptions {
                orientation,
                zero_policy,
                ..BFileOptions::default()
            };
            flatten_rows_to_bfile(&data, &options)
                .unwrap()
                .text
                .lines()
                .map(|line| line.split(' ').nth(1).unwrap().to_string())
                .collect::<Vec<_>>()
                .join(",")
        };
        use BFileRowOrientation::{AsStored, Reversed};
        assert_eq!(write(BFileZeroPolicy::KeepAll, AsStored), "0,0,5,0,6,0,0,0");
        assert_eq!(write(BFileZeroPolicy::TrimLowOrder, AsStored), "5,0,6,0");
        assert_eq!(write(BFileZeroPolicy::TrimHighOrder, AsStored), "0,0,5,0,6");
        assert_eq!(write(BFileZeroPolicy::TrimBothEnds, AsStored), "5,0,6");
        assert_eq!(write(BFileZeroPolicy::TrimHighOrder, Reversed), "6,0,5,0,0");
    }

    #[test]
    fn cap_counts_terms_after_the_zero_policy() {
        let data = rows(&[&[1, 0, 0], &[2, 3, 0]]);
        let options = BFileOptions {
            zero_policy: BFileZeroPolicy::TrimHighOrder,
            max_terms: Some(3),
            ..BFileOptions::default()
        };
        let file = flatten_rows_to_bfile(&data, &options).unwrap();
        assert_eq!((file.rows_included, file.terms_included), (2, 3));
    }

    #[test]
    fn index_overflow_is_an_error() {
        let single = rows(&[&[1]]);
        let at_limit = BFileOptions {
            first_index: i64::MAX,
            ..BFileOptions::default()
        };
        assert!(flatten_rows_to_bfile(&single, &at_limit).is_ok());
        let data = rows(&[&[1, 2]]);
        let options = BFileOptions {
            first_index: i64::MAX,
            ..BFileOptions::default()
        };
        assert_eq!(
            flatten_rows_to_bfile(&data, &options),
            Err(BFileError::IndexOverflow)
        );
    }
}
