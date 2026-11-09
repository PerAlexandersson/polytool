//! OEIS-ready text for the rows of any recurrence: the `%e` "Triangle
//! begins:" block, the `%F` recurrence sentence, and `%t` (Mathematica) and
//! `%o` (PARI/GP) programs that reproduce the rows exactly.
//!
//! Entry points: [`OeisRecurrencePresentation`] re-indexes a recurrence (a
//! [`Recurrence`] with its seed rows, or a `polytool.recurrence.v1` record,
//! [`RecurrenceJson`]) to the row numbers of an OEIS entry and formats it
//! ([`OeisRecurrencePresentation::oeis_text`] gives every field at once as an
//! [`OeisText`]); [`OeisTriangleReading`] says how a coefficient row is read
//! as a row of the OEIS triangle; [`format_triangle_example`] lays out rows
//! that were already read.
//!
//! # Conventions
//!
//! * The polynomials are `P(n,x)`, coefficients in ascending powers of `x`
//!   exactly as the recurrence produces them.  Row `n` of the OEIS triangle is
//!   the coefficient list of `P(n,x)`, read as [`OeisTriangleReading`] says:
//!   zero coefficients dropped at either end by the [`BFileZeroPolicy`]
//!   (applied to the ascending list), optionally padded with zeros at the
//!   high end to the lengths `L0, L0 + 1, ...` of a regular triangle, and
//!   then read in ascending or descending powers ([`BFileRowOrientation`]).
//!   These are the same options as the b-file formatter
//!   [`crate::oeis::flatten_rows_to_bfile`], so a b-file of the read rows
//!   with the plain options agrees term for term with the text here.
//! * The OEIS row number of the first row is
//!   [`OeisTriangleReading::first_row_index`] (for a triangle read by rows, the
//!   offset `%O` of the entry).  The recurrence index `n` in the formula and
//!   the programs is this OEIS row number: the recurrence is shifted with
//!   [`Recurrence::shift_index`], so the coefficients are exactly those of the
//!   source recurrence at the corresponding row.
//! * Text is plain ASCII in OEIS style: `*` for products, `^` for powers,
//!   `d/dx P(n-1,x)` for derivatives, exact rationals `a/b`.  Coefficients of
//!   the recurrence are grouped by powers of `x`, with polynomial coefficients
//!   in `n` in descending powers, e.g. `(1+(2*n-2)*x)*P(n-1,x)`.
//! * The Mathematica program is a memoized `P[n_, x_] := P[n, x] =
//!   Expand[...]` with the seed rows as separate definitions; a left-hand
//!   factor `f(n,x)` is divided out exactly with `Together`.  The PARI/GP
//!   program fills a vector of polynomials (PARI reduces the exact quotient
//!   by `f(n,x)` to a polynomial).  Both end with the flattened rows
//!   `T(o..N)`.
//! * A zero row is written as the single term `0` by the b-file and the
//!   example block; Mathematica's `CoefficientList[0, x]` is the empty list.
//!
//! What is claimed: the text describes the supplied recurrence and seeds
//! exactly.  Whether they define the OEIS entry for every `n` is the
//! caller's identification (usually a finite comparison of terms).

use super::bfile::{bfile_row_terms, BFileRowOrientation, BFileZeroPolicy};
use crate::recurrence::{
    BigRational, BivarPoly, Recurrence, RecurrenceError, RecurrenceEvaluationError, RecurrenceJson,
    RecurrenceRowIterator, RecurrenceSign, RecurrenceTerm,
};
use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};
use serde::{Deserialize, Serialize};
use std::fmt;
use std::fmt::Write as _;

/// JSON schema tag of [`OeisText`].
pub const OEIS_TEXT_SCHEMA: &str = "polytool.oeis-text.v1";

/// How the coefficient rows of a polynomial sequence are read as the rows of
/// an OEIS triangle (see the conventions of [`crate::oeis`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct OeisTriangleReading {
    /// OEIS row number `n` of the first row (the offset `%O` of a triangle
    /// read by rows).
    pub first_row_index: i64,
    /// Ascending powers (as stored) or descending powers within each row.
    pub orientation: BFileRowOrientation,
    /// Zero coefficients dropped from the ends of each ascending row, before
    /// any padding and before the orientation is applied.
    pub zero_policy: BFileZeroPolicy,
    /// Read the rows as a regular triangle: row `i` (counted from the first
    /// row) has exactly `L0 + i` entries, `L0` the length of the first read
    /// row.  A shorter row is padded with zeros at its high end (when no low
    /// zeros are dropped); reading stops before the first row that is longer,
    /// or shorter while low zeros are dropped.
    pub regular_triangle: bool,
}

impl OeisTriangleReading {
    fn trims_low_zeros(&self) -> bool {
        matches!(
            self.zero_policy,
            BFileZeroPolicy::TrimLowOrder | BFileZeroPolicy::TrimBothEnds
        )
    }

    fn reversed(&self) -> bool {
        self.orientation == BFileRowOrientation::Reversed
    }

    /// The OEIS rows of `rows` (ascending coefficient rows, the first one being
    /// OEIS row [`first_row_index`](Self::first_row_index)).
    ///
    /// The result is the longest prefix the reading accepts: with
    /// [`regular_triangle`](Self::regular_triangle) it ends before the first
    /// row that breaks the pattern, otherwise every row is read.
    ///
    /// ```
    /// use num_bigint::BigInt;
    /// use polytool::oeis::{BFileRowOrientation, OeisTriangleReading};
    ///
    /// let rows: Vec<Vec<BigInt>> = [vec![1], vec![1, 2], vec![1, 10, 4]]
    ///     .iter()
    ///     .map(|row| row.iter().map(|&c| BigInt::from(c)).collect())
    ///     .collect();
    /// let reading = OeisTriangleReading {
    ///     orientation: BFileRowOrientation::Reversed,
    ///     ..OeisTriangleReading::default()
    /// };
    /// let read = reading.read_rows(&rows);
    /// assert_eq!(read[2], vec![BigInt::from(4), BigInt::from(10), BigInt::from(1)]);
    /// ```
    pub fn read_rows(&self, rows: &[Vec<BigInt>]) -> Vec<Vec<BigInt>> {
        let mut read = Vec::with_capacity(rows.len());
        let mut first_length = None;
        for (position, row) in rows.iter().enumerate() {
            let mut kept = bfile_row_terms(row, self.zero_policy).to_vec();
            if self.regular_triangle {
                let expected = *first_length.get_or_insert(kept.len()) + position;
                if kept.len() > expected || (kept.len() < expected && self.trims_low_zeros()) {
                    break;
                }
                kept.resize(expected, BigInt::zero());
            }
            if self.reversed() {
                kept.reverse();
            }
            read.push(kept);
        }
        read
    }
}

/// Error of [`OeisRecurrencePresentation`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OeisTextError {
    /// The recurrence record is malformed.
    Recurrence(RecurrenceError),
    /// Replaying the recurrence failed.
    Evaluation(RecurrenceEvaluationError),
    /// The first OEIS row lies before the first row of the recurrence.
    FirstRowBeforeRecurrence {
        /// Requested absolute index of the first OEIS row.
        first_row: usize,
        /// Absolute index of the first row of the recurrence.
        first_index: usize,
    },
    /// A generated row has a non-integral coefficient (OEIS rows are integers).
    NonIntegralCoefficient {
        /// OEIS row number of the row.
        row_index: i64,
    },
    /// An index does not fit in `i64`.
    IndexOverflow,
}

impl fmt::Display for OeisTextError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Recurrence(error) => write!(f, "{error}"),
            Self::Evaluation(error) => write!(f, "{error}"),
            Self::FirstRowBeforeRecurrence {
                first_row,
                first_index,
            } => write!(
                f,
                "the first OEIS row {first_row} lies before the first recurrence row {first_index}"
            ),
            Self::NonIntegralCoefficient { row_index } => {
                write!(f, "row {row_index} has a non-integral coefficient")
            }
            Self::IndexOverflow => write!(f, "row index overflows i64"),
        }
    }
}

impl std::error::Error for OeisTextError {}

impl From<RecurrenceError> for OeisTextError {
    fn from(error: RecurrenceError) -> Self {
        Self::Recurrence(error)
    }
}

impl From<RecurrenceEvaluationError> for OeisTextError {
    fn from(error: RecurrenceEvaluationError) -> Self {
        Self::Evaluation(error)
    }
}

/// Every OEIS text field of a recurrence, as returned by
/// [`OeisRecurrencePresentation::oeis_text`] (serialized by the CLI and the
/// web interface with schema [`OEIS_TEXT_SCHEMA`]).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct OeisText {
    /// Schema identifier, [`OEIS_TEXT_SCHEMA`].
    pub schema: String,
    /// OEIS row number of the first row.
    pub first_row_index: i64,
    /// Number of rows in [`triangle_example`](Self::triangle_example) and
    /// produced by the programs.
    pub rows_shown: usize,
    /// The `%e` block, "Triangle begins:" followed by the aligned rows.
    pub triangle_example: String,
    /// The `%F` sentence stating the recurrence and its initial rows.
    pub recurrence_sentence: String,
    /// The `%t` Mathematica program.
    pub mathematica: String,
    /// The `%o` PARI/GP program.
    pub pari: String,
}

/// A recurrence re-indexed to the row numbers of an OEIS triangle, with the
/// seed rows it needs from the first OEIS row on (see the
/// conventions of [`crate::oeis`]).
#[derive(Debug, Clone)]
pub struct OeisRecurrencePresentation {
    recurrence: Recurrence,
    first_row: usize,
    seeds: Vec<Vec<BigRational>>,
    reading: OeisTriangleReading,
    displayed: Recurrence,
}

impl OeisRecurrencePresentation {
    /// Present `recurrence`, whose rows start at absolute index `first_index`
    /// with `initial_rows`, from the absolute row `first_row` on, which becomes
    /// OEIS row `reading.first_row_index`.
    ///
    /// The seeds are the rows `first_row, first_row + 1, ...`: as many as the
    /// largest lag of the recurrence, and every supplied initial row from
    /// `first_row` on (a record keeps rows where its left-hand factor vanishes
    /// as initial rows).
    ///
    /// # Errors
    ///
    /// [`OeisTextError::FirstRowBeforeRecurrence`] when `first_row <
    /// first_index`, [`OeisTextError::Evaluation`] when the seeds cannot be
    /// generated, [`OeisTextError::IndexOverflow`] when an index does not fit
    /// in `i64`.
    pub fn new(
        recurrence: &Recurrence,
        first_index: usize,
        initial_rows: &[Vec<BigRational>],
        first_row: usize,
        reading: OeisTriangleReading,
    ) -> Result<Self, OeisTextError> {
        let skipped =
            first_row
                .checked_sub(first_index)
                .ok_or(OeisTextError::FirstRowBeforeRecurrence {
                    first_row,
                    first_index,
                })?;
        let seed_count = recurrence
            .max_offset()
            .max(initial_rows.len().saturating_sub(skipped));
        let total = skipped
            .checked_add(seed_count)
            .ok_or(OeisTextError::IndexOverflow)?;
        let mut rows = recurrence.generate_rows_rational(initial_rows, first_index, total)?;
        let seeds = rows.split_off(skipped);
        let shift = i64::try_from(first_row)
            .ok()
            .and_then(|row| row.checked_sub(reading.first_row_index))
            .ok_or(OeisTextError::IndexOverflow)?;
        i64::try_from(seed_count)
            .ok()
            .and_then(|count| reading.first_row_index.checked_add(count))
            .ok_or(OeisTextError::IndexOverflow)?;
        Ok(Self {
            displayed: recurrence.shift_index(shift),
            recurrence: recurrence.clone(),
            first_row,
            seeds,
            reading,
        })
    }

    /// [`Self::new`] for a `polytool.recurrence.v1` record.  `first_row`
    /// defaults to the record's first index.
    ///
    /// # Errors
    ///
    /// [`OeisTextError::Recurrence`] for a malformed record, otherwise as
    /// [`Self::new`].
    pub fn from_recurrence_json(
        record: &RecurrenceJson,
        first_row: Option<usize>,
        reading: OeisTriangleReading,
    ) -> Result<Self, OeisTextError> {
        let (recurrence, first_index, initial_rows) = record.to_recurrence_parts()?;
        Self::new(
            &recurrence,
            first_index,
            &initial_rows,
            first_row.unwrap_or(first_index),
            reading,
        )
    }

    /// The reading of the rows.
    pub fn reading(&self) -> &OeisTriangleReading {
        &self.reading
    }

    /// The recurrence with `n` the OEIS row number.
    #[deprecated(
        since = "0.6.0",
        note = "unused; it will be removed (`OeisText` renders the displayed recurrence)"
    )]
    pub fn displayed_recurrence(&self) -> &Recurrence {
        &self.displayed
    }

    /// The seed rows, ascending coefficients; the first is OEIS row
    /// [`OeisTriangleReading::first_row_index`].
    pub fn seed_rows(&self) -> &[Vec<BigRational>] {
        &self.seeds
    }

    /// OEIS row number from which the recurrence applies (the first row after
    /// the seeds).
    pub fn recurrence_start(&self) -> i64 {
        // Checked in the constructor.
        self.reading.first_row_index + self.seeds.len() as i64
    }

    /// Lazy iterator over the rows `P(o,x), P(o+1,x), ...` (ascending exact
    /// coefficients, `o` the first OEIS row number), seeds first; see
    /// [`RecurrenceRowIterator`].
    pub fn rows_iter(&self) -> RecurrenceRowIterator {
        self.recurrence
            .generate_rows_rational_iter(&self.seeds, self.first_row)
    }

    /// The first `count` rows as integer coefficient rows (ascending, not yet
    /// read as OEIS rows).
    ///
    /// # Errors
    ///
    /// [`OeisTextError::Evaluation`] when a step fails,
    /// [`OeisTextError::NonIntegralCoefficient`] for a non-integral row.
    pub fn generate_integer_rows(&self, count: usize) -> Result<Vec<Vec<BigInt>>, OeisTextError> {
        let mut rows = Vec::with_capacity(count);
        for (position, row) in self.rows_iter().take(count).enumerate() {
            let row = row?;
            let row_index = self.reading.first_row_index + position as i64;
            rows.push(
                integer_row(&row).ok_or(OeisTextError::NonIntegralCoefficient { row_index })?,
            );
        }
        Ok(rows)
    }

    /// The first `count` OEIS rows (read with [`Self::reading`]; fewer when a
    /// regular triangle reading stops early).
    ///
    /// # Errors
    ///
    /// As [`Self::generate_integer_rows`].
    pub fn triangle_rows(&self, count: usize) -> Result<Vec<Vec<BigInt>>, OeisTextError> {
        Ok(self.reading.read_rows(&self.generate_integer_rows(count)?))
    }

    /// Every text field, for the first `example_rows` OEIS rows (fewer when a
    /// regular triangle reading stops early; the programs list the same rows).
    ///
    /// # Errors
    ///
    /// As [`Self::generate_integer_rows`].
    ///
    /// ```
    /// use polytool::oeis::{OeisRecurrencePresentation, OeisTriangleReading};
    /// use polytool::recurrence::RecurrenceJson;
    ///
    /// // P(n) = (2 + x) P(n-1) + (1 + x) P(n-2), P(0) = 1, P(1) = 1 + x.
    /// let record: RecurrenceJson = serde_json::from_str(r#"{
    ///   "schema": "polytool.recurrence.v1", "first_index": 0,
    ///   "initial_polynomials": [["1"], ["1", "1"]],
    ///   "recurrence": {"terms": [
    ///     {"offset": 1, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["2", "1"]]}},
    ///     {"offset": 2, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["1", "1"]]}}]}
    /// }"#).unwrap();
    /// let presentation =
    ///     OeisRecurrencePresentation::from_recurrence_json(&record, None, OeisTriangleReading::default())
    ///         .unwrap();
    /// let text = presentation.oeis_text(4).unwrap();
    /// assert_eq!(
    ///     text.recurrence_sentence,
    ///     "The row polynomials satisfy P(n,x) = (2+x)*P(n-1,x) + (1+x)*P(n-2,x) for n >= 2, \
    ///      with P(0,x) = 1, P(1,x) = 1 + x."
    /// );
    /// assert!(text.triangle_example.starts_with("Triangle begins:\n  0: 1;\n  1: 1,  1;\n"));
    /// ```
    pub fn oeis_text(&self, example_rows: usize) -> Result<OeisText, OeisTextError> {
        let rows = self.triangle_rows(example_rows)?;
        let last_row = self.reading.first_row_index + rows.len() as i64 - 1;
        Ok(OeisText {
            schema: OEIS_TEXT_SCHEMA.to_string(),
            first_row_index: self.reading.first_row_index,
            rows_shown: rows.len(),
            triangle_example: format_triangle_example(&rows, self.reading.first_row_index),
            recurrence_sentence: self.format_recurrence_sentence(),
            mathematica: self.format_mathematica_program(last_row),
            pari: self.format_pari_program(last_row),
        })
    }

    /// The `%F` sentence, for example `The row polynomials satisfy P(n,x) =
    /// (2+x)*P(n-1,x) + (1+x)*P(n-2,x) for n >= 2, with P(0,x) = 1, P(1,x) =
    /// 1 + x.`  A reading other than ascending powers without dropped low
    /// zeros is stated first (`Row n lists the coefficients of P(n,x) in
    /// descending powers of x, where ...`).
    pub fn format_recurrence_sentence(&self) -> String {
        let equation = format!(
            "{} = {} for n >= {}",
            self.left_side(Syntax::OeisText),
            self.right_side(Syntax::OeisText),
            self.recurrence_start()
        );
        let seeds = self
            .seeds
            .iter()
            .enumerate()
            .map(|(position, row)| {
                format!(
                    "P({},x) = {}",
                    self.reading.first_row_index + position as i64,
                    univariate_text(row, Syntax::OeisText)
                )
            })
            .collect::<Vec<_>>();
        let with_seeds = if seeds.is_empty() {
            String::new()
        } else {
            format!(", with {}", seeds.join(", "))
        };
        let trims = self.reading.trims_low_zeros();
        if !self.reading.reversed() && !trims {
            return format!("The row polynomials satisfy {equation}{with_seeds}.");
        }
        let order = match (self.reading.reversed(), trims) {
            (false, _) => "ascending powers of x, starting at the lowest nonzero coefficient",
            (true, false) => "descending powers of x",
            (true, true) => "descending powers of x, down to the lowest nonzero coefficient",
        };
        format!("Row n lists the coefficients of P(n,x) in {order}, where {equation}{with_seeds}.")
    }

    /// The `%t` Mathematica program listing the flattened rows from the first
    /// OEIS row up to row `last_row`.
    pub fn format_mathematica_program(&self, last_row: i64) -> String {
        let mut lines: Vec<String> = self
            .seeds
            .iter()
            .enumerate()
            .map(|(position, row)| {
                format!(
                    "P[{}, x_] := {};",
                    self.reading.first_row_index + position as i64,
                    univariate_text(row, Syntax::Mathematica)
                )
            })
            .collect();
        let right = self.right_side(Syntax::Mathematica);
        let body = match &self.displayed.denominator {
            Some(factor) => format!(
                "Expand[Together[({right})/({})]]",
                bivariate_text(factor, Syntax::Mathematica)
            ),
            None => format!("Expand[{right}]"),
        };
        lines.push(format!("P[n_, x_] := P[n, x] = {body};"));
        let mut row = if self.reading.trims_low_zeros() {
            "CoefficientList[Expand[P[n, x]/x^Exponent[P[n, x], x, Min]], x]".to_string()
        } else {
            "CoefficientList[P[n, x], x]".to_string()
        };
        if self.reading.regular_triangle {
            row = format!(
                "PadRight[{row}, {}]",
                plus("n", self.regular_length_shift())
            );
        }
        if self.reading.reversed() {
            row = format!("Reverse[{row}]");
        }
        lines.push(format!(
            "Table[{row}, {{n, {}, {last_row}}}] // Flatten",
            self.reading.first_row_index
        ));
        lines.join("\n")
    }

    /// The `%o` PARI/GP program listing the flattened rows from the first OEIS
    /// row up to row `last_row` (call `T(N)` for other row ranges).
    pub fn format_pari_program(&self, last_row: i64) -> String {
        let o = self.reading.first_row_index;
        let mut fill = format!("my(v = vector({}));", plus("N", 1 - o));
        for (position, row) in self.seeds.iter().enumerate() {
            let _ = write!(
                fill,
                " v[{}] = {};",
                position + 1,
                univariate_text(row, Syntax::Pari)
            );
        }
        let right = self.right_side(Syntax::Pari);
        let value = match &self.displayed.denominator {
            Some(factor) => format!("({right})/({})", bivariate_text(factor, Syntax::Pari)),
            None => right,
        };
        let _ = write!(
            fill,
            " for(n = {}, N, my(i = {}); v[i] = {value}); v",
            self.recurrence_start(),
            plus("n", 1 - o)
        );
        let mut row_steps = vec![if self.reading.trims_low_zeros() {
            "my(c = Vecrev(p/x^valuation(p, x)));".to_string()
        } else {
            "my(c = Vecrev(p));".to_string()
        }];
        if self.reading.regular_triangle {
            row_steps.push(format!(
                " c = concat(c, vector({} - #c));",
                plus("n", self.regular_length_shift())
            ));
        }
        if self.reading.reversed() {
            row_steps.push(" c = Vecrev(c);".to_string());
        }
        row_steps.push(" c".to_string());
        [
            format!("P(N) = {{{fill}}};"),
            format!("row(n, p) = {{{}}};", row_steps.concat()),
            format!(
                "T(N) = {{my(v = P(N)); concat(vector(#v, i, row({}, v[i])))}};",
                plus("i", o - 1)
            ),
            format!("T({last_row})"),
        ]
        .join("\n")
    }

    /// `L0 - o`: row `n` of a regular triangle has `n + (L0 - o)` entries.
    fn regular_length_shift(&self) -> i64 {
        let first_length = self.seeds.first().map_or_else(
            || {
                self.rows_iter()
                    .next()
                    .and_then(Result::ok)
                    .map_or(0, |row| kept_length(&row, self.reading.zero_policy))
            },
            |row| kept_length(row, self.reading.zero_policy),
        );
        first_length as i64 - self.reading.first_row_index
    }

    fn left_side(&self, syntax: Syntax) -> String {
        let reference = reference_text(syntax, 0, 0);
        match &self.displayed.denominator {
            None => reference,
            Some(factor) => match product_term(factor, false, &reference, syntax) {
                Some((negative, body)) => join_signed(&[(negative, body)], true),
                None => format!("0*{reference}"),
            },
        }
    }

    fn right_side(&self, syntax: Syntax) -> String {
        let mut terms: Vec<(bool, String)> = self
            .displayed
            .terms
            .iter()
            .filter_map(|term| recurrence_term_text(term, syntax))
            .collect();
        if let Some(inhomogeneous) = &self.displayed.inhomogeneous {
            terms.extend(bivariate_terms(inhomogeneous, syntax.inner_spaced()));
        }
        join_signed(&terms, true)
    }
}

/// The `%e` block "Triangle begins:" for rows already read as OEIS rows, the
/// first being row `first_row_index`: one line per row, the row number on the
/// left, columns right-aligned, entries separated by commas and each row ended
/// by a semicolon.
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::oeis::format_triangle_example;
///
/// let rows: Vec<Vec<BigInt>> = [vec![1], vec![2, 1], vec![4, 10, 1]]
///     .iter()
///     .map(|row| row.iter().map(|&c| BigInt::from(c)).collect())
///     .collect();
/// assert_eq!(
///     format_triangle_example(&rows, 0),
///     "Triangle begins:\n  0: 1;\n  1: 2,  1;\n  2: 4, 10, 1;\n"
/// );
/// ```
pub fn format_triangle_example(rows: &[Vec<BigInt>], first_row_index: i64) -> String {
    let cells: Vec<Vec<String>> = rows
        .iter()
        .map(|row| row.iter().map(ToString::to_string).collect())
        .collect();
    let mut widths: Vec<usize> = Vec::new();
    for row in &cells {
        for (column, cell) in row.iter().enumerate() {
            if column == widths.len() {
                widths.push(0);
            }
            widths[column] = widths[column].max(cell.len());
        }
    }
    let labels: Vec<String> = (0..rows.len())
        .map(|position| (first_row_index + position as i64).to_string())
        .collect();
    let label_width = labels.iter().map(String::len).max().unwrap_or(0);
    let mut text = String::from("Triangle begins:\n");
    for (label, row) in labels.iter().zip(&cells) {
        let entries = row
            .iter()
            .enumerate()
            .map(|(column, cell)| format!("{cell:>width$}", width = widths[column]))
            .collect::<Vec<_>>()
            .join(", ");
        let _ = writeln!(text, "  {label:>label_width$}: {entries};");
    }
    text
}

// ---------------------------------------------------------------------------
// Expression text
// ---------------------------------------------------------------------------

/// Output syntax of an expression.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Syntax {
    /// OEIS `%F` text: `P(n-1,x)`, `d/dx P(n-1,x)`, compact parentheses.
    OeisText,
    /// Mathematica: `P[n - 1, x]`, `D[P[n - 1, x], x]`.
    Mathematica,
    /// PARI/GP inside the vector fill: `v[i - 1]`, `deriv(v[i - 1], x)`.
    Pari,
}

impl Syntax {
    /// Whether sums inside the parentheses of a coefficient get spaces.
    fn inner_spaced(self) -> bool {
        self != Self::OeisText
    }
}

fn integer_row(row: &[BigRational]) -> Option<Vec<BigInt>> {
    row.iter()
        .map(|c| c.is_integer().then(|| c.to_integer()))
        .collect()
}

/// Length of an ascending rational row after the zero policy.
fn kept_length(row: &[BigRational], zero_policy: BFileZeroPolicy) -> usize {
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
    end - start
}

/// `base + k` as text: `n + 2`, `n - 1`, `n`.
fn plus(base: &str, k: i64) -> String {
    match k.cmp(&0) {
        std::cmp::Ordering::Equal => base.to_string(),
        std::cmp::Ordering::Greater => format!("{base} + {k}"),
        std::cmp::Ordering::Less => format!("{base} - {}", k.unsigned_abs()),
    }
}

fn rational_text(r: &BigRational) -> String {
    if r.is_integer() {
        r.numer().to_string()
    } else {
        format!("{}/{}", r.numer(), r.denom())
    }
}

fn power_text(variable: &str, exponent: usize) -> Option<String> {
    match exponent {
        0 => None,
        1 => Some(variable.to_string()),
        _ => Some(format!("{variable}^{exponent}")),
    }
}

/// `|c| n^i x^j` for a positive `c`: `(3/2)*n^2*x`, `x`, `1`.
fn monomial_text(coefficient: &BigRational, n_power: usize, x_power: usize) -> String {
    let mut factors = Vec::new();
    let has_variables = n_power > 0 || x_power > 0;
    if !(coefficient.is_one() && has_variables) {
        let text = rational_text(coefficient);
        factors.push(if has_variables && text.contains('/') {
            format!("({text})")
        } else {
            text
        });
    }
    factors.extend(power_text("n", n_power));
    factors.extend(power_text("x", x_power));
    factors.join("*")
}

fn join_signed(terms: &[(bool, String)], spaced: bool) -> String {
    if terms.is_empty() {
        return "0".to_string();
    }
    let (plus_sign, minus_sign) = if spaced { (" + ", " - ") } else { ("+", "-") };
    let mut text = String::new();
    for (position, (negative, body)) in terms.iter().enumerate() {
        match (position, negative) {
            (0, true) => text.push('-'),
            (0, false) => {}
            (_, true) => text.push_str(minus_sign),
            (_, false) => text.push_str(plus_sign),
        }
        text.push_str(body);
    }
    text
}

/// Signed terms of a polynomial in `n` and `x`, grouped by ascending powers
/// of `x`; a coefficient polynomial in `n` (descending powers) is
/// parenthesized, `(2*n-2)*x`, except at `x^0`, whose monomials are listed
/// separately.
fn bivariate_terms(poly: &BivarPoly, inner_spaced: bool) -> Vec<(bool, String)> {
    let width = poly.coeffs.iter().map(Vec::len).max().unwrap_or(0);
    let mut terms = Vec::new();
    for x_power in 0..width {
        let n_terms: Vec<(usize, &BigRational)> = (0..poly.coeffs.len())
            .rev()
            .filter_map(|n_power| {
                poly.coeffs[n_power]
                    .get(x_power)
                    .filter(|c| !c.is_zero())
                    .map(|c| (n_power, c))
            })
            .collect();
        if n_terms.len() == 1 || x_power == 0 {
            for (n_power, c) in n_terms {
                terms.push((c.is_negative(), monomial_text(&c.abs(), n_power, x_power)));
            }
            continue;
        }
        if n_terms.is_empty() {
            continue;
        }
        let negate = n_terms[0].1.is_negative();
        let inner: Vec<(bool, String)> = n_terms
            .iter()
            .map(|(n_power, c)| {
                (
                    c.is_negative() != negate,
                    monomial_text(&c.abs(), *n_power, 0),
                )
            })
            .collect();
        let x_text = power_text("x", x_power).expect("positive power");
        terms.push((
            negate,
            format!("({})*{x_text}", join_signed(&inner, inner_spaced)),
        ));
    }
    terms
}

fn bivariate_text(poly: &BivarPoly, syntax: Syntax) -> String {
    join_signed(&bivariate_terms(poly, syntax.inner_spaced()), true)
}

fn univariate_text(row: &[BigRational], syntax: Syntax) -> String {
    bivariate_text(
        &BivarPoly {
            coeffs: vec![row.to_vec()],
        },
        syntax,
    )
}

/// `c(n,x) * reference` as a signed term; `None` for a zero coefficient.
fn product_term(
    coefficient: &BivarPoly,
    alternating: bool,
    reference: &str,
    syntax: Syntax,
) -> Option<(bool, String)> {
    let terms = bivariate_terms(coefficient, syntax.inner_spaced());
    let (negative, mut body) = match terms.as_slice() {
        [] => return None,
        [(negative, body)] if body == "1" => (*negative, reference.to_string()),
        [(negative, body)] => (*negative, format!("{body}*{reference}")),
        [(first_negative, _), ..] => {
            let negate = *first_negative;
            let inner: Vec<(bool, String)> = terms
                .iter()
                .map(|(negative, body)| (*negative != negate, body.clone()))
                .collect();
            (
                negate,
                format!(
                    "({})*{reference}",
                    join_signed(&inner, syntax.inner_spaced())
                ),
            )
        }
    };
    if alternating {
        body = format!("(-1)^n*{body}");
    }
    Some((negative, body))
}

/// `P(n-r,x)` (or its `d`-th derivative) in `syntax`.
fn reference_text(syntax: Syntax, lag: usize, derivative_order: usize) -> String {
    match syntax {
        Syntax::OeisText => {
            let index = if lag == 0 {
                "n".to_string()
            } else {
                format!("n-{lag}")
            };
            let base = format!("P({index},x)");
            match derivative_order {
                0 => base,
                1 => format!("d/dx {base}"),
                d => format!("d^{d}/dx^{d} {base}"),
            }
        }
        Syntax::Mathematica => {
            let index = plus("n", -(lag as i64));
            let base = format!("P[{index}, x]");
            match derivative_order {
                0 => base,
                1 => format!("D[{base}, x]"),
                d => format!("D[{base}, {{x, {d}}}]"),
            }
        }
        Syntax::Pari => {
            let mut text = format!("v[{}]", plus("i", -(lag as i64)));
            for _ in 0..derivative_order {
                text = format!("deriv({text}, x)");
            }
            text
        }
    }
}

fn recurrence_term_text(term: &RecurrenceTerm, syntax: Syntax) -> Option<(bool, String)> {
    product_term(
        &term.coeff,
        term.sign == RecurrenceSign::AlternatingN,
        &reference_text(syntax, term.offset, term.deriv_order),
        syntax,
    )
}

#[cfg(test)]
mod tests;
