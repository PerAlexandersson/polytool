//! Independent exact read-back of the emitted LaTeX definition subset.
//!
//! Coefficients are expanded by the shared bivariate expression parser, so
//! factored powers are compared symbolically, rather than at sample indices.

use super::*;
use crate::recurrence::closed_form::parse_polynomial_with_exponent_limit;
use crate::recurrence::implicit::algebra::{Biv, BivarConversion};
use num_traits::One;
use std::collections::BTreeMap;

fn invalid(message: &str) -> RecurrenceError {
    RecurrenceError::InvalidInput(format!("LaTeX read-back: {message}"))
}

fn polynomial(text: &str, limit: usize) -> Result<Biv, RecurrenceError> {
    let mut text = text
        .replace("\\bigl", "")
        .replace("\\bigr", "")
        .replace("\\,", "");
    while let Some(start) = text.find("\\frac{") {
        let rest = &text[start + 6..];
        let (numerator, rest) = rest
            .split_once("}{")
            .ok_or_else(|| invalid("invalid fraction"))?;
        let (denominator, _) = rest
            .split_once('}')
            .ok_or_else(|| invalid("invalid fraction"))?;
        if !numerator.bytes().all(|c| c.is_ascii_digit())
            || !denominator.bytes().all(|c| c.is_ascii_digit())
            || numerator.is_empty()
            || denominator.is_empty()
        {
            return Err(invalid("fraction requires integer arguments"));
        }
        let end = start + 6 + numerator.len() + 2 + denominator.len() + 1;
        text.replace_range(start..end, &format!("({numerator}/{denominator})"));
    }
    // Curly groups in coefficients are only literal exponents.
    while let Some(start) = text.find("^{") {
        let rest = &text[start + 2..];
        let (exponent, _) = rest
            .split_once('}')
            .ok_or_else(|| invalid("unclosed exponent"))?;
        if exponent.is_empty() || !exponent.bytes().all(|c| c.is_ascii_digit()) {
            return Err(invalid("invalid exponent"));
        }
        let end = start + 2 + exponent.len() + 1;
        text.replace_range(start..end, &format!("^{exponent}"));
    }
    if text.contains(['{', '}']) {
        return Err(invalid("unexpected coefficient group"));
    }
    parse_polynomial_with_exponent_limit(&text.replace('t', "k"), limit)
}

// Only split top-level sums; subscripts, powers and factors keep their signs.
fn summands(text: &str) -> Result<Vec<&str>, RecurrenceError> {
    let mut stack = Vec::new();
    let mut start = 0;
    let mut parts = Vec::new();
    for (index, character) in text.char_indices() {
        match character {
            '(' | '{' => stack.push(character),
            ')' | '}' => {
                let expected = if character == ')' { '(' } else { '{' };
                if stack.pop() != Some(expected) {
                    return Err(invalid("unbalanced groups"));
                }
            }
            '+' | '-'
                if stack.is_empty() && index > start && !text[start..index].trim().is_empty() =>
            {
                parts.push(text[start..index].trim());
                start = index;
            }
            _ => {}
        }
    }
    if !stack.is_empty() {
        return Err(invalid("unbalanced groups"));
    }
    parts.push(text[start..].trim());
    Ok(parts)
}

type TermKey = (usize, usize, bool);
type Terms = BTreeMap<TermKey, Biv>;

fn side(text: &str, limit: usize) -> Result<(Terms, Biv), RecurrenceError> {
    let mut terms = Terms::new();
    let mut forcing = Biv::zero();
    for summand in summands(text)? {
        let (negative, body) = if let Some(rest) = summand.strip_prefix('-') {
            (true, rest.trim())
        } else {
            (false, summand.strip_prefix('+').unwrap_or(summand).trim())
        };
        let Some((coefficient, reference)) = body.split_once("P_{") else {
            let value = polynomial(body, limit)?;
            forcing = if negative {
                forcing.subtract(&value)
            } else {
                forcing.add(&value)
            };
            continue;
        };
        let (index, derivative) = reference
            .split_once('}')
            .ok_or_else(|| invalid("invalid row reference"))?;
        let offset = if index == "n" {
            0
        } else {
            index
                .strip_prefix("n-")
                .ok_or_else(|| invalid("invalid lag"))?
                .parse::<usize>()
                .map_err(|_| invalid("invalid lag"))?
        };
        let derivative = match derivative.trim() {
            "" => 0,
            "'" => 1,
            "''" => 2,
            other => other
                .strip_prefix("^{(")
                .and_then(|s| s.strip_suffix(")}"))
                .ok_or_else(|| invalid("invalid derivative"))?
                .parse::<usize>()
                .map_err(|_| invalid("invalid derivative"))?,
        };
        let (coefficient, alternating) = match coefficient.trim().strip_suffix("(-1)^n") {
            Some(coefficient) => (coefficient.trim(), true),
            None => (coefficient.trim(), false),
        };
        let mut value = if coefficient.is_empty() {
            Biv::constant(BigRational::one())
        } else {
            polynomial(coefficient, limit)?
        };
        if negative {
            value = value.scale(&-BigRational::one());
        }
        let entry = terms
            .entry((offset, derivative, alternating))
            .or_insert_with(Biv::zero);
        *entry = entry.add(&value);
    }
    terms.retain(|_, value| !value.is_zero());
    Ok((terms, forcing))
}

fn expected_terms(rule: &Recurrence) -> Terms {
    let mut terms = Terms::new();
    for term in &rule.terms {
        let key = (
            term.offset,
            term.deriv_order,
            term.sign == RecurrenceSign::AlternatingN,
        );
        let entry = terms.entry(key).or_insert_with(Biv::zero);
        *entry = entry.add(&Biv::from_bivar(&term.coeff));
    }
    terms.retain(|_, value| !value.is_zero());
    terms
}

/// Parse an emitted LaTeX definition and compare it exactly with `record`.
///
/// Checks every seed's index and polynomial, the normalized recurrence identity
/// (including derivatives, signs, forcing and leader), and the starting index.
/// Both factored and expanded writer output are accepted; comments are ignored.
/// The supported grammar is only the subset emitted by [`RecurrenceJson::to_latex`].
/// An identically zero leader prevents normalization; its original algebraic
/// relation is checked instead, retaining the presentation's refusal to give a
/// defining range. This is a rendering check, not a proof of the recurrence.
///
/// # Errors
/// Returns [`RecurrenceError::InvalidInput`] for malformed or mismatching text,
/// or a record whose normalized recurrence metadata disagrees with its rule.
pub fn verify_latex_recurrence_output(
    text: &str,
    record: &RecurrenceJson,
) -> Result<(), RecurrenceError> {
    let (rule, first, seeds) = record.to_recurrence_parts()?;
    let rule = rule.normalized_presentation().unwrap_or(rule);
    let start = first
        .checked_add(seeds.len())
        .ok_or_else(|| invalid("index overflow"))?;
    // Parsing these JSON coefficients also validates normalized_recurrence.
    if let Some(normalized) = &record.normalized_recurrence {
        let mut expected = record.clone();
        expected.recurrence = normalized.clone();
        let (expected, _, _) = expected.to_recurrence_parts()?;
        if expected_terms(&expected) != expected_terms(&rule)
            || expected.denominator != rule.denominator
            || expected.inhomogeneous != rule.inhomogeneous
        {
            return Err(invalid("normalized_recurrence disagrees with the rule"));
        }
    }
    let mut limit = 1usize;
    for polynomial in rule
        .terms
        .iter()
        .map(|t| &t.coeff)
        .chain(rule.denominator.iter())
        .chain(rule.inhomogeneous.iter())
    {
        limit = limit
            .max(polynomial.coeffs.len())
            .max(polynomial.coefficient_width());
    }
    limit = limit.max(seeds.iter().map(Vec::len).max().unwrap_or(0));
    let formulas = text
        .lines()
        .filter(|line| !line.starts_with('%'))
        .collect::<Vec<_>>();
    if formulas.len() != seeds.len() + 1 {
        return Err(invalid("wrong number of formulas"));
    }
    for (offset, row) in seeds.iter().enumerate() {
        let prefix = format!("$P_{{{}}} = ", first + offset);
        let expression = formulas[offset]
            .strip_prefix(&prefix)
            .and_then(|s| s.strip_suffix('$'))
            .ok_or_else(|| invalid("invalid initial row label"))?;
        if polynomial(expression, limit)?
            != Biv::from_bivar(&BivarPoly {
                coeffs: vec![row.clone()],
            })
        {
            return Err(invalid("initial polynomial mismatch"));
        }
    }
    let formula = formulas[seeds.len()]
        .strip_prefix('$')
        .and_then(|s| s.strip_suffix('$'))
        .ok_or_else(|| invalid("missing formula delimiters"))?;
    let formula = if record
        .presentation
        .as_ref()
        .and_then(|p| p.reason.as_ref())
        .is_some()
    {
        formula
    } else {
        let suffix = format!(" \\quad (n \\geq {start})");
        formula
            .strip_suffix(&suffix)
            .ok_or_else(|| invalid("starting index mismatch"))?
    };
    let (lhs, rhs) = formula
        .split_once(" = ")
        .ok_or_else(|| invalid("missing equality"))?;
    let leader = rule
        .denominator
        .as_ref()
        .map(Biv::from_bivar)
        .unwrap_or_else(|| Biv::constant(BigRational::one()));
    let mut expected_lhs = Terms::new();
    if !leader.is_zero() {
        expected_lhs.insert((0, 0, false), leader);
    }
    if side(lhs, limit)? != (expected_lhs, Biv::zero()) {
        return Err(invalid("left-hand side mismatch"));
    }
    let forcing = rule
        .inhomogeneous
        .as_ref()
        .map(Biv::from_bivar)
        .unwrap_or_else(Biv::zero);
    if side(rhs, limit)? != (expected_terms(&rule), forcing) {
        return Err(invalid("right-hand side mismatch"));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    #[test]
    fn symbolic_factors_seeds_and_range_reject_tampering() {
        let record: RecurrenceJson = serde_json::from_value(json!({
            "schema":"polytool.recurrence.v1", "first_index":3,
            "initial_polynomials":[["1/2","-2"]],
            "recurrence":{"terms":[{"offset":1,"deriv_order":0,"sign":"none",
                "coeff":{"coeffs":[["0","0","1"],["0","2"],["1"]]}}]}
        }))
        .unwrap();
        for latex in [
            record.to_latex().unwrap(),
            record.to_latex_expanded().unwrap(),
        ] {
            verify_latex_recurrence_output(&latex, &record).unwrap();
            for tampered in [
                latex.replace("P_{3}", "P_{2}"),
                latex.replace("geq 4", "geq 5"),
                latex.replace("P_{n-1}", "P_{n-2}"),
                latex.replace("frac{1}{2}", "frac{1}{3}"),
            ] {
                assert!(
                    verify_latex_recurrence_output(&tampered, &record).is_err(),
                    "{tampered}"
                );
            }
        }
        let text = record.to_latex().unwrap();
        assert!(text.contains("(n + t)^{2}"), "{text}");
        assert!(verify_latex_recurrence_output(
            &text.replace("(n + t)^{2}", "(n + t)^{3}"),
            &record
        )
        .is_err());
        let checked = record.clone().with_latex().unwrap();
        assert_eq!(checked.latex_checked, Some(true));
        let mut wrong = record;
        wrong.normalized_recurrence = Some(serde_json::from_value(json!({"terms":[]})).unwrap());
        assert!(wrong.to_latex().is_err());
        let mut overflow = checked;
        overflow.first_index = usize::MAX;
        assert!(verify_latex_recurrence_output(&text, &overflow).is_err());
    }

    #[test]
    fn zero_leader_preserves_the_explicit_nondefining_relation() {
        let record: RecurrenceJson = serde_json::from_value(json!({
            "schema":"polytool.recurrence.v1", "first_index":0,
            "initial_polynomials":[["0"]], "recurrence":{
                "denominator":{"coeffs":[["0"]]},
                "terms":[{"offset":1,"deriv_order":0,"sign":"none",
                    "coeff":{"coeffs":[["1"]]}}]}
        }))
        .unwrap();
        let record = record
            .with_checked_presentation(&[vec![BigRational::from_integer(0.into())]])
            .unwrap();
        let checked = record.with_latex().unwrap();
        assert!(checked.normalized_recurrence.is_none());
        assert_eq!(checked.latex_checked, Some(true));
        assert!(!checked.latex.as_ref().unwrap().contains("geq"));
        verify_latex_recurrence_output(checked.latex.as_ref().unwrap(), &checked).unwrap();
    }

    #[test]
    fn expanded_and_factored_derivatives_alternation_leader_and_forcing() {
        let record: RecurrenceJson = serde_json::from_value(json!({
            "schema":"polytool.recurrence.v1", "first_index":0,
            "initial_polynomials":[["0"]],
            "recurrence":{"denominator":{"coeffs":[["2","2"],["2","2"]]},
                "inhomogeneous":{"coeffs":[["-1/2","1/2"],["0","-2"]]},
                "terms":[{"offset":1,"deriv_order":3,"sign":"alternating_n",
                    "coeff":{"coeffs":[["0","0","-1/2"],["0","-1"],["-1/2"]]}}]}
        }))
        .unwrap();
        for latex in [
            record.to_latex().unwrap(),
            record.to_latex_expanded().unwrap(),
        ] {
            verify_latex_recurrence_output(&latex, &record).unwrap();
            assert!(verify_latex_recurrence_output(&latex.replace("(-1)^n", ""), &record).is_err());
            assert!(
                verify_latex_recurrence_output(&latex.replace("^{(3)}", "'"), &record).is_err()
            );
            assert!(verify_latex_recurrence_output(&(latex + "\n$0 = 0$"), &record).is_err());
        }
        for invalid in [
            "\\frac{1}{0}",
            "(n+t",
            "n/0",
            "n;",
            "\\sin(n)",
            "n^{999999}",
            "n^{2",
            "n^{2}}",
        ] {
            assert!(polynomial(invalid, 4).is_err());
        }
        // The shared parser's usual degree-64 bound is not a writer limit.
        assert_eq!(
            polynomial("t^{65}", 65).unwrap(),
            Biv::monomial(0, 65, BigRational::one())
        );
    }
}
