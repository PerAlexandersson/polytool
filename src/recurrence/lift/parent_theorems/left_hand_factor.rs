//! Exact normalization of a left-hand factor `x^e f(n)`.
//!
//! A rule `x^e f(n) P_n = sum_k c_k(n, x) D^(d_k) P_(n-r_k) + g(n, x)`
//! (`e >= 1`, `f` free of `x`), as the adaptive search stores it in
//! [`Recurrence::denominator`], is rewritten as
//! `f(n) P_n = sum_k (c_k / x^e) D^(d_k) P_(n-r_k) + g / x^e` when `x^e`
//! divides every right-hand coefficient `c_k` and `g` in `Q[n, x]` (terms
//! with the same lag, derivative order and sign are added first, so the
//! normalized rule has one term for each).
//!
//! The two rules define the same rows.  Write `R_n` for the right-hand side
//! of the original rule evaluated on the rows `P_(n-1), P_(n-2), ...` and
//! `R'_n` for that of the normalized rule on the same rows; by construction
//! `R_n = x^e R'_n` as polynomials in `x`.  At an index `n` with
//! `f(n) != 0`, the original rule sets `P_n = R_n / (x^e f(n))` and the
//! normalized rule sets `P_n = R'_n / f(n)`; these agree by cancelling `x^e`
//! in `Q[x]`, and both are polynomials (the original division is exact
//! because `x^e` divides `R_n`).  At an index with `f(n) = 0` both left-hand
//! factors vanish identically, so both rules keep that row as a seed
//! ([`Recurrence::generation_initial_count`]).  Induction on `n` gives equal
//! rows for every index.
//!
//! When the left-hand factor is `x^e f(n)` but `x^e` does not divide some
//! right-hand coefficient, the division by `x^e` happens only on the rows
//! (it may be exact for the particular rows without being an identity of
//! the rule), and the rule is not rewritten: the callers report the
//! theorem or the lift as not applicable with that reason.
//!
//! [`LeftHandFactorNormalization`] records the original and the normalized
//! rule; [`LeftHandFactorNormalization::verify`] re-derives the
//! normalization from the original rule.

use super::*;

/// The record of an exact normalization of a left-hand factor `x^e f(n)`
/// (module documentation of `parent_theorems`, "Left-hand factors
/// `x^e f(n)`"): the rule as supplied and the exactly equivalent rule
/// `f(n) P_n = sum (c / x^e) D^d P_(n-r)` used by the check.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LeftHandFactorNormalization {
    /// The power `e >= 1` of `x` divided out of both sides.
    pub x_power: usize,
    /// The rule as supplied, with left-hand factor `x^e f(n)`.
    pub original_rule: RecurrenceJsonData,
    /// The normalized rule, with left-hand factor `f(n)` (omitted when `f`
    /// is the constant 1) and every right-hand coefficient divided by `x^e`.
    pub normalized_rule: RecurrenceJsonData,
    /// The normalization in words.
    pub statement: String,
}

impl LeftHandFactorNormalization {
    /// The normalization of `rule` (module documentation).
    ///
    /// Returns `Ok(None)` when the rule has no left-hand factor, when the
    /// factor is free of `x`, or when it is not of the form `x^e f(n)` (it
    /// has more than one power of `x`); such rules are left as they are.
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::Unsupported`] when the left-hand factor is
    /// `x^e f(n)`, `e >= 1`, but `x^e` does not divide a right-hand
    /// coefficient or the inhomogeneous term: the rule is then not exactly
    /// equivalent to one with a left-hand factor free of `x`.
    pub fn from_recurrence(rule: &Recurrence) -> Result<Option<Self>, RecurrenceError> {
        Ok(normalize_left_hand_factor(rule)?.map(|(_, record)| record))
    }

    /// The normalized rule.
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::InvalidInput`] when the recorded rule does not
    /// parse.
    pub fn normalized_recurrence(&self) -> Result<Recurrence, RecurrenceError> {
        rule_of_json_data(&self.normalized_rule)
    }

    /// Re-derive the normalization from [`Self::original_rule`] and compare
    /// it with the record.
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::VerificationFailed`] when the original rule does
    /// not normalize, or normalizes to a different record.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        let original = rule_of_json_data(&self.original_rule)?;
        let recomputed = Self::from_recurrence(&original).map_err(|error| {
            RecurrenceError::VerificationFailed(format!(
                "the original rule does not normalize: {error}"
            ))
        })?;
        match recomputed {
            Some(recomputed) if &recomputed == self => Ok(()),
            Some(_) => Err(RecurrenceError::VerificationFailed(
                "the left-hand factor normalization differs from its recomputation".into(),
            )),
            None => Err(RecurrenceError::VerificationFailed(
                "the original rule has no left-hand factor x^e f(n) with e >= 1".into(),
            )),
        }
    }
}

/// The rule with a left-hand factor `x^e f(n)` normalized, with its record;
/// `Ok(None)` when there is nothing to normalize (see
/// [`LeftHandFactorNormalization::from_recurrence`]).
pub(in crate::recurrence) fn normalize_left_hand_factor(
    rule: &Recurrence,
) -> Result<Option<(Recurrence, LeftHandFactorNormalization)>, RecurrenceError> {
    let Some(denominator) = &rule.denominator else {
        return Ok(None);
    };
    let coefficients = denominator.normalized_coefficients();
    let powers: Vec<usize> = (0..x_width(&coefficients))
        .filter(|&j| !x_column(&coefficients, j).is_empty())
        .collect();
    let (e, f) = match powers.as_slice() {
        [e] if *e >= 1 => (*e, x_column(&coefficients, *e)),
        _ => return Ok(None),
    };
    let divide = |rows: &Bivariate, what: &str| -> Result<BivarPoly, RecurrenceError> {
        if let Some(lowest) = (0..e.min(x_width(rows))).find(|&j| !x_column(rows, j).is_empty()) {
            return Err(RecurrenceError::Unsupported(format!(
                "the left-hand factor is x^{e} f(n), but x^{e} does not divide {what} (it has a \
                 term in x^{lowest}): the division by x^{e} happens only on the rows, so the \
                 rule is not exactly equivalent to a rule with a left-hand factor free of x"
            )));
        }
        Ok(BivarPoly::from_coefficients_normalized(
            rows.iter()
                .map(|row| row.iter().skip(e).cloned().collect())
                .collect(),
        ))
    };
    // Terms with the same lag, derivative order and sign are added first:
    // divisibility is a property of their sum.
    let mut groups: Vec<(RecurrenceTerm, Bivariate)> = Vec::new();
    for term in &rule.terms {
        let coefficients = term.coeff.normalized_coefficients();
        match groups.iter_mut().find(|(first, _)| {
            (first.offset, first.deriv_order, first.sign)
                == (term.offset, term.deriv_order, term.sign)
        }) {
            Some((_, sum)) => *sum = add_bivariate(sum, &coefficients),
            None => groups.push((term.clone(), coefficients)),
        }
    }
    let mut terms = Vec::with_capacity(groups.len());
    for (term, sum) in &groups {
        let what = format!(
            "the coefficient of D^{} P_(n-{}){}",
            term.deriv_order,
            term.offset,
            if term.is_alternating() {
                " (with the factor (-1)^n)"
            } else {
                ""
            }
        );
        terms.push(RecurrenceTerm {
            offset: term.offset,
            deriv_order: term.deriv_order,
            sign: term.sign,
            coeff: divide(sum, &what)?,
        });
    }
    let inhomogeneous = rule
        .inhomogeneous
        .as_ref()
        .map(|g| divide(&g.normalized_coefficients(), "the inhomogeneous term"))
        .transpose()?;
    let f_is_one = f.len() == 1 && f[0].is_one();
    let normalized = Recurrence {
        terms,
        denominator: (!f_is_one).then(|| {
            BivarPoly::from_coefficients_normalized(f.iter().map(|c| vec![c.clone()]).collect())
        }),
        inhomogeneous,
    };
    let record = LeftHandFactorNormalization {
        x_power: e,
        original_rule: json_data_of_rule(rule),
        normalized_rule: json_data_of_rule(&normalized),
        statement: format!(
            "the left-hand factor x^{e} f(n), f(n) = {}, is normalized to f(n): x^{e} divides \
             every right-hand coefficient in Q[n, x], so f(n) P_n = sum (c / x^{e}) D^d P_(n-r) \
             defines the same rows (x^{e} cancels in Q[x] at every index with f(n) != 0, and \
             both rules keep a row with f(n) = 0 as a seed)",
            format_index_polynomial(&f, "n")
        ),
    };
    Ok(Some((normalized, record)))
}

/// `rule`, or its normalization with the record (module documentation).
///
/// # Errors
///
/// As [`LeftHandFactorNormalization::from_recurrence`].
pub(in crate::recurrence) fn rule_with_x_free_left_hand_factor(
    rule: &Recurrence,
) -> Result<(Recurrence, Option<LeftHandFactorNormalization>), RecurrenceError> {
    Ok(match normalize_left_hand_factor(rule)? {
        Some((normalized, record)) => (normalized, Some(record)),
        None => (rule.clone(), None),
    })
}

fn json_data_of_rule(rule: &Recurrence) -> RecurrenceJsonData {
    RecurrenceJson::from_recurrence_rational(rule, 0, &[], None).recurrence
}

fn rule_of_json_data(data: &RecurrenceJsonData) -> Result<Recurrence, RecurrenceError> {
    let record = RecurrenceJson {
        schema: RECURRENCE_JSON_SCHEMA.into(),
        first_index: 0,
        initial_polynomials: Vec::new(),
        recurrence: data.clone(),
        search: None,
    };
    Ok(record.to_recurrence_parts()?.0)
}
