//! Lean definitions are independently read back before they can be emitted.

use polytool::recurrence::{verify_lean_recurrence_output, BigRational, RecurrenceJson};
use serde_json::json;

fn record(value: serde_json::Value) -> RecurrenceJson {
    serde_json::from_value(value).unwrap()
}

fn singular() -> (RecurrenceJson, Vec<Vec<BigRational>>) {
    let record = record(json!({
        "schema":"polytool.recurrence.v1", "first_index":0,
        "initial_polynomials":[["0"],["1","1"]],
        "recurrence":{"denominator":{"coeffs":[["-1"],["1"]]}, "terms":[{
            "offset":1,"deriv_order":0,"sign":"none",
            "coeff":{"coeffs":[["0","3"],["3","1"],["1"]]}}]}
    }));
    let (rule, first, seeds) = record.to_recurrence_parts().unwrap();
    let rows = rule.generate_rows_rational(&seeds, first, 8).unwrap();
    (record.with_checked_presentation(&rows).unwrap(), rows)
}

#[test]
fn singular_rule_emits_two_bases_and_four_expanded_examples() {
    let (record, rows) = singular();
    let text = record.to_lean("A216118", &rows).unwrap();
    assert!(text.contains("def A216118 : ℕ → ℝ[X]"));
    assert!(text.contains("  | 0 => 0\n  | 1 => 1 + X"));
    assert!(text.contains("| n + 2 =>"));
    assert_eq!(text.matches("example :").count(), 4);
    assert_eq!(text.matches("  rfl").count(), 2);
    assert_eq!(text.matches("  sorry").count(), 2);
    assert!(text.lines().all(|line| line.chars().count() <= 100));
    verify_lean_recurrence_output(&text, "A216118", &rows).unwrap();
    let shifted = record.reindexed(7).unwrap();
    assert_eq!(shifted.to_lean("A216118", &rows).unwrap(), text);
}

#[test]
fn safeguard_rejects_missing_base_wrong_shift_sign_and_precedence() {
    let (record, rows) = singular();
    let good = record.to_lean("A216118", &rows).unwrap();
    let bad_shift = good.replace("(n + 1)", "(n + 0)");
    assert_ne!(good, bad_shift);
    assert!(verify_lean_recurrence_output(&bad_shift, "A216118", &rows).is_err());
    let coefficient_shift = good.replace("(n : ℝ)", "((n : ℝ) + 1)");
    assert!(verify_lean_recurrence_output(&coefficient_shift, "A216118", &rows).is_err());
    let bad_sign = good.replacen("=> C", "=> -C", 1);
    assert_ne!(good, bad_sign);
    assert!(verify_lean_recurrence_output(&bad_sign, "A216118", &rows).is_err());
    let bad_precedence = good.replacen("=> C", "=> 1 + C", 1);
    assert_ne!(good, bad_precedence);
    assert!(verify_lean_recurrence_output(&bad_precedence, "A216118", &rows).is_err());
    let missing_base = good.replace("  | 1 => 1 + X\n", "");
    assert!(verify_lean_recurrence_output(&missing_base, "A216118", &rows).is_err());
    // Lean's total real division yields row 1 = 0 with this one-base definition.
    let one_base = "def A216118 : ℕ → ℝ[X]\n  | 0 => 0\n  | n + 1 =>
        C (((n : ℝ) + 4) / (n : ℝ)) * (C ((n : ℝ) + 1) + X) * A216118 n
        \nexample : A216118 0 = 0 := by\n  rfl
        \nexample : A216118 1 = 1 + X := by\n  sorry
        \nexample : A216118 2 = 5 * (2 + X) * (1 + X) := by\n  sorry";
    let error = verify_lean_recurrence_output(one_base, "A216118", &rows)
        .unwrap_err()
        .to_string();
    assert!(error.contains("supplied row 1"), "{error}");
}

#[test]
fn derivatives_alternation_forcing_rationals_and_index_shifts_replay() {
    let record = record(json!({"schema":"polytool.recurrence.v1","first_index":3,
        "initial_polynomials":[["1","1/2"]], "recurrence":{
            "terms":[
                {"offset":1,"deriv_order":0,"sign":"alternating_n", "coeff":{"coeffs":[["1/2"],["1"]]}},
                {"offset":1,"deriv_order":1,"sign":"none", "coeff":{"coeffs":[["0","1"]]}}],
            "inhomogeneous":{"coeffs":[["-1/2","1"],["1"]]}}}));
    let (rule, first, seeds) = record.to_recurrence_parts().unwrap();
    let rows = rule.generate_rows_rational(&seeds, first, 9).unwrap();
    let record = record.with_checked_presentation(&rows).unwrap();
    let text = record.to_lean("D", &rows).unwrap();
    assert!(text.contains("(-1 : ℝ) ^ (n + 1)"));
    assert!(text.contains("(D n).derivative"));
    let unshifted_sign = text.replace("(-1 : ℝ) ^ (n + 1)", "(-1 : ℝ) ^ n");
    assert!(verify_lean_recurrence_output(&unshifted_sign, "D", &rows).is_err());
    let bad_derivative = text.replace("(D n).derivative", "D n.derivative");
    assert!(verify_lean_recurrence_output(&bad_derivative, "D", &rows).is_err());
    let wrong_fraction = text.replace("C (1 / 2)", "C ((1 / 2 : ℕ))");
    assert_ne!(text, wrong_fraction);
    assert!(verify_lean_recurrence_output(&wrong_fraction, "D", &rows).is_err());
}

#[test]
fn unsupported_leaders_names_rows_and_examples_refuse_output() {
    let (record, rows) = singular();
    let mut t_dependent = record.clone();
    t_dependent.recurrence.denominator.as_mut().unwrap().coeffs[0].push("1".into());
    assert!(t_dependent
        .to_lean("A", &rows)
        .unwrap_err()
        .to_string()
        .contains("depends on t"));
    let mut too_few_bases = record.clone();
    too_few_bases.initial_polynomials.pop();
    assert!(too_few_bases
        .to_lean("A", &rows)
        .unwrap_err()
        .to_string()
        .contains("vanishes"));
    for name in [
        "",
        "n",
        "X",
        "C",
        "def",
        "Bad Name",
        "A; theorem",
        "A.B",
        "a\n",
    ] {
        assert!(record.to_lean(name, &rows).is_err(), "{name:?}");
    }
    let good = record.to_lean("A", &rows).unwrap();
    let unbound_index = good.replace("example : A 0 = 0", "example : A 0 = C (n : ℝ)");
    assert!(verify_lean_recurrence_output(&unbound_index, "A", &rows).is_err());
    let bad_example = good.replace("example : A 0 = 0", "example : A 0 = 1");
    assert!(verify_lean_recurrence_output(&bad_example, "A", &rows).is_err());
    let polynomial_division = good.replacen("| n + 2 =>", "| n + 2 => X / (X + 1) +", 1);
    assert!(verify_lean_recurrence_output(&polynomial_division, "A", &rows).is_err());
    let mut wrong_rows = rows.clone();
    wrong_rows[7][0] += BigRational::from_integer(1.into());
    assert!(record.to_lean("A", &wrong_rows).is_err());
}

#[test]
fn natural_subtraction_truncates_and_natural_division_floors() {
    let output = "def F : ℕ → ℝ[X]\n | 0 => 0\n | n + 1 => C (((n - 1 : ℕ) : ℝ))
        \nexample : F 0 = 0 := by rfl
        \nexample : F 1 = 0 := by sorry
        \nexample : F 2 = 0 := by sorry";
    let zero = vec![BigRational::from_integer(0.into())];
    verify_lean_recurrence_output(output, "F", &[zero.clone(), zero.clone(), zero.clone()])
        .unwrap();
    let real_subtraction = output.replace("((n - 1 : ℕ) : ℝ)", "((n : ℝ) - 1)");
    assert!(
        verify_lean_recurrence_output(&real_subtraction, "F", &[zero.clone(), zero.clone()])
            .is_err()
    );
    let division = output.replace("((n - 1 : ℕ) : ℝ)", "((n / 2 : ℕ) : ℝ)");
    verify_lean_recurrence_output(&division, "F", &[zero.clone(), zero.clone(), zero.clone()])
        .unwrap();
    let real_division = division.replace("((n / 2 : ℕ) : ℝ)", "((n : ℝ) / 2)");
    assert!(verify_lean_recurrence_output(
        &real_division,
        "F",
        &[zero.clone(), zero.clone(), zero]
    )
    .is_err());
}

#[test]
fn readable_coefficients_examples_and_atom_preserving_layout() {
    let source = record(json!({"schema":"polytool.recurrence.v1","first_index":0,
        "initial_polynomials":[["1","-2","1"]],"recurrence":{"terms":[
            {"offset":1,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["1"],["0","2"]]}},
            {"offset":1,"deriv_order":1,"sign":"none","coeff":{"coeffs":[["0","1","-1"]]}}],
            "inhomogeneous":{"coeffs":[["-3"],["-1"]]}}}));
    let (rule, first, seeds) = source.to_recurrence_parts().unwrap();
    let rows = rule.generate_rows_rational(&seeds, first, 7).unwrap();
    let source = source.with_checked_presentation(&rows).unwrap();
    let text = source.to_lean("ReadableRows", &rows).unwrap();
    assert!(text.contains("| 0 => 1 - 2 * X + X ^ 2"), "{text}");
    assert!(text.contains("C (2 + 2 * (n : ℝ))"), "{text}");
    assert!(
        text.contains("(X - X ^ 2) * (ReadableRows n).derivative"),
        "{text}"
    );
    assert!(text.contains("- C (4 + (n : ℝ))"), "{text}");
    for forbidden in ["(1)", "(2)", "1 * X", "+ (-", "/ (1)", " / 1)", "(n + 0)"] {
        assert!(!text.contains(forbidden), "{forbidden}: {text}");
    }
    for line in text.lines() {
        assert!(line.chars().count() <= 100, "{line}");
        assert_eq!(
            line.matches('(').count(),
            line.matches(')').count(),
            "split atom: {line}"
        );
        if line.starts_with("      ") {
            assert!(
                ["+ ", "- ", "* "]
                    .iter()
                    .any(|op| line.trim_start().starts_with(op)),
                "{line}"
            );
        }
    }
    verify_lean_recurrence_output(&text, "ReadableRows", &rows).unwrap();
}

#[test]
fn structurally_invalid_references_are_refused_even_when_values_agree() {
    let source = record(json!({"schema":"polytool.recurrence.v1","first_index":0,
        "initial_polynomials":[["1"],["1"]],"recurrence":{"terms":[
            {"offset":1,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["1"]]}}]}}));
    let rows = vec![vec![BigRational::from_integer(1.into())]; 8];
    let text = source.to_lean("F", &rows).unwrap();
    for argument in [
        "n + 1 - 1",
        "n + 2 - 1",
        "n + 0",
        "n - 0",
        "n / 1",
        "1 + n",
        "0",
        "n + 2",
    ] {
        let wrong = text.replace("F (n + 1)", &format!("F ({argument})"));
        assert_ne!(text, wrong);
        let error = verify_lean_recurrence_output(&wrong, "F", &rows).unwrap_err();
        assert!(
            error.to_string().contains("recursive reference"),
            "{argument}: {error}"
        );
    }
    let zero_reference = text.replace("F (n + 1)", "F n");
    verify_lean_recurrence_output(&zero_reference, "F", &rows).unwrap();
}

#[test]
fn keyword_and_open_polynomial_collisions_have_typed_explanations() {
    use polytool::recurrence::{validate_lean_recurrence_name, RecurrenceError};
    for name in [
        "have",
        "show",
        "at",
        "from",
        "do",
        "calc",
        "fun",
        "let",
        "in",
        "if",
        "then",
        "else",
        "match",
        "with",
        "where",
        "end",
        "Type",
        "Prop",
        "Sort",
        "this",
        "structure",
        "instance",
    ] {
        let error = validate_lean_recurrence_name(name).unwrap_err();
        assert!(matches!(error, RecurrenceError::InvalidInput(_)));
        assert!(
            error.to_string().contains("keyword or reserved word"),
            "{name}: {error}"
        );
    }
    for name in [
        "X",
        "C",
        "n",
        "derivative",
        "degree",
        "natDegree",
        "eval",
        "coeff",
        "Polynomial",
        "Nat",
        "Real",
    ] {
        let error = validate_lean_recurrence_name(name).unwrap_err();
        assert!(matches!(error, RecurrenceError::InvalidInput(_)));
        assert!(error.to_string().contains("clashes"), "{name}: {error}");
    }
    for name in ["Rows", "family_rows", "A216118", "Rows'"] {
        validate_lean_recurrence_name(name).unwrap();
    }
}

#[test]
fn scalar_coefficients_are_reduced_before_choosing_constant_embedding() {
    let source = record(json!({"schema":"polytool.recurrence.v1","first_index":0,
        "initial_polynomials":[["1"],["3"]],"recurrence":{"denominator":{"coeffs":[["-1"],["1"]]},
            "terms":[{"offset":1,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["-2"],["2"]]}}]}}));
    let (rule, first, seeds) = source.to_recurrence_parts().unwrap();
    let rows = rule.generate_rows_rational(&seeds, first, 8).unwrap();
    let source = source.with_checked_presentation(&rows).unwrap();
    let text = source.to_lean("ReducedRows", &rows).unwrap();
    assert!(
        text.contains("| n + 2 => 2 * ReducedRows (n + 1)"),
        "{text}"
    );
    assert!(!text.contains("C ("), "{text}");
}
