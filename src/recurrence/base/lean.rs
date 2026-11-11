//! Lean polynomial definitions and an independent read-back evaluator.
//!
//! Only the emitted expression grammar is accepted. Evaluation distinguishes
//! naturals, reals and polynomials; real division by zero is zero and natural
//! subtraction truncates. Polynomial division is deliberately unsupported.

use super::*;
use crate::vec_poly;
use num_traits::{One, Signed, ToPrimitive, Zero};

mod readback;
pub use readback::verify_lean_recurrence_output;

fn invalid(message: impl Into<String>) -> RecurrenceError {
    RecurrenceError::InvalidInput(message.into())
}

/// Check a simple Lean identifier accepted by the recurrence renderer.
/// Names are ASCII identifiers with optional trailing apostrophes, at most 60
/// characters, excluding syntax keywords and names used by the emitted subset.
/// Namespaces and quoted identifiers are not part of this output subset.
pub fn validate_lean_recurrence_name(name: &str) -> Result<(), RecurrenceError> {
    if name.is_empty()
        || name.len() > 60
        || !name.as_bytes()[0].is_ascii_alphabetic()
        || !name
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b == b'_' || b == b'\'')
    {
        return Err(invalid(
            "invalid Lean recurrence name: use a simple ASCII identifier of at most 60 characters",
        ));
    }
    if [
        "n",
        "X",
        "C",
        "derivative",
        "degree",
        "natDegree",
        "eval",
        "coeff",
        "Polynomial",
        "Nat",
        "Real",
    ]
    .contains(&name)
    {
        return Err(invalid(format!(
            "Lean recurrence name '{name}' clashes with a binder or an open Polynomial name"
        )));
    }
    if [
        "def",
        "example",
        "by",
        "rfl",
        "sorry",
        "open",
        "import",
        "namespace",
        "end",
        "match",
        "with",
        "let",
        "in",
        "where",
        "if",
        "then",
        "else",
        "fun",
        "theorem",
        "abbrev",
        "noncomputable",
        "termination_by",
        "decreasing_by",
        "have",
        "show",
        "at",
        "from",
        "do",
        "calc",
        "Type",
        "Prop",
        "Sort",
        "this",
        "forall",
        "exists",
        "instance",
        "structure",
        "inductive",
        "class",
        "axiom",
        "constant",
        "constants",
        "opaque",
        "private",
        "protected",
        "public",
        "export",
        "attribute",
        "deriving",
        "variable",
        "variables",
        "universe",
        "universes",
        "section",
        "set_option",
        "syntax",
        "macro",
        "macro_rules",
        "elab",
        "elab_rules",
        "notation",
        "infix",
        "infixl",
        "infixr",
        "prefix",
        "postfix",
        "scoped",
        "local",
        "mutual",
        "return",
        "for",
        "while",
        "repeat",
        "unless",
        "try",
        "catch",
        "finally",
        "break",
        "continue",
        "partial",
        "unsafe",
        "initialize",
        "builtin_initialize",
        "as",
        "renaming",
        "hiding",
        "include",
        "omit",
        "extends",
        "using",
        "and",
        "or",
        "of",
        "nomatch",
        "sorryAx",
        "admit",
    ]
    .contains(&name)
    {
        return Err(invalid(format!(
            "Lean recurrence name '{name}' is a keyword or reserved word"
        )));
    }
    Ok(())
}

fn integer(value: &BigInt) -> String {
    let text = value.to_string();
    if text.len() <= 60 {
        return text;
    }
    let base = BigInt::from(10u8).pow(30);
    format!(
        "({} * 10 ^ 30 + {})",
        integer(&(value / &base)),
        integer(&(value % base))
    )
}

fn real_rational(value: &BigRational) -> String {
    if value.denom().is_one() {
        integer(value.numer())
    } else {
        format!("{} / {}", integer(value.numer()), integer(value.denom()))
    }
}

#[derive(Clone)]
struct SignedProduct {
    negative: bool,
    factors: Vec<String>,
}

fn product_text(product: &SignedProduct) -> String {
    product.factors.join(" * ")
}

fn sum_text(terms: &[SignedProduct]) -> String {
    let mut text = String::new();
    for term in terms {
        text.push_str(match (text.is_empty(), term.negative) {
            (true, false) => "",
            (true, true) => "-",
            (false, false) => " + ",
            (false, true) => " - ",
        });
        text.push_str(&product_text(term));
    }
    if text.is_empty() {
        "0".into()
    } else {
        text
    }
}

fn index_polynomial(coefficients: &[BigRational]) -> String {
    let terms: Vec<_> = coefficients
        .iter()
        .enumerate()
        .filter(|(_, c)| !c.is_zero())
        .map(|(i, c)| {
            let magnitude = c.abs();
            let mut factors = Vec::new();
            if i == 0 || !magnitude.is_one() {
                factors.push(real_rational(&magnitude));
            }
            if i > 0 {
                factors.push(if i == 1 {
                    "(n : ℝ)".into()
                } else {
                    format!("(n : ℝ) ^ {i}")
                });
            }
            SignedProduct {
                negative: c.is_negative(),
                factors,
            }
        })
        .collect();
    sum_text(&terms)
}

fn row_products(row: &[BigRational]) -> Vec<SignedProduct> {
    row.iter()
        .enumerate()
        .filter(|(_, c)| !c.is_zero())
        .map(|(j, c)| {
            let magnitude = c.abs();
            let mut factors = Vec::new();
            if j == 0 || !magnitude.is_one() {
                factors.push(if magnitude.denom().is_one() {
                    integer(magnitude.numer())
                } else {
                    format!("C ({})", real_rational(&magnitude))
                });
            }
            if j > 0 {
                factors.push(if j == 1 {
                    "X".into()
                } else {
                    format!("X ^ {j}")
                });
            }
            SignedProduct {
                negative: c.is_negative(),
                factors,
            }
        })
        .collect()
}

fn grouped_sum(text: String) -> String {
    if text.contains(" + ") || text.contains(" - ") {
        format!("({text})")
    } else {
        text
    }
}

fn coefficient_products(coefficient: &BivarPoly, leader: &[BigRational]) -> Vec<SignedProduct> {
    let mut terms = Vec::new();
    for j in 0..coefficient.coefficient_width() {
        let column: Vec<_> = coefficient
            .coeffs
            .iter()
            .map(|r| r.get(j).cloned().unwrap_or_else(BigRational::zero))
            .collect();
        if column.iter().all(Zero::is_zero) {
            continue;
        }
        let numerator = crate::Polynomial::new(column);
        let denominator = crate::Polynomial::new(leader.to_vec());
        let common = numerator.gcd(&denominator);
        let numerator = numerator.div_rem(&common).0;
        let denominator = denominator.div_rem(&common).0;
        let mut column = numerator.coeffs().to_vec();
        let leader = denominator.coeffs();
        let constant_leader = leader.iter().skip(1).all(Zero::is_zero);
        if constant_leader {
            for c in &mut column {
                *c /= &leader[0];
            }
        }
        let negative = column
            .iter()
            .find(|c| !c.is_zero())
            .is_some_and(Signed::is_negative);
        if negative {
            for c in &mut column {
                *c = -c.clone();
            }
        }
        let constant = column.iter().skip(1).all(Zero::is_zero);
        let scalar = if constant_leader && constant {
            if column[0].denom().is_one() {
                integer(column[0].numer())
            } else {
                format!("C ({})", real_rational(&column[0]))
            }
        } else if constant_leader {
            format!("C ({})", index_polynomial(&column))
        } else {
            format!(
                "C ({} / {})",
                grouped_sum(index_polynomial(&column)),
                grouped_sum(index_polynomial(leader))
            )
        };
        let mut factors = Vec::new();
        if j == 0 || scalar != "1" {
            factors.push(scalar);
        }
        if j > 0 {
            factors.push(if j == 1 {
                "X".into()
            } else {
                format!("X ^ {j}")
            });
        }
        terms.push(SignedProduct { negative, factors });
    }
    terms
}

fn append_recurrence_products(
    output: &mut Vec<SignedProduct>,
    mut coefficients: Vec<SignedProduct>,
    sign: Option<String>,
    reference: String,
) {
    if coefficients.is_empty() {
        return;
    }
    // Keep a short polynomial coefficient together; distribute a longer one
    // so every C(...) and every cast remains an unbroken atom.
    let negative = coefficients[0].negative;
    let mut positive = coefficients.clone();
    if negative {
        for term in &mut positive {
            term.negative = !term.negative;
        }
    }
    let grouped = sum_text(&positive);
    if coefficients.len() > 1 && grouped.chars().count() <= 78 {
        coefficients = vec![SignedProduct {
            negative,
            factors: vec![format!("({grouped})")],
        }];
    }
    for mut coefficient in coefficients {
        if coefficient.factors == ["1"] {
            coefficient.factors.clear();
        }
        if let Some(sign) = &sign {
            coefficient.factors.insert(0, sign.clone());
        }
        coefficient.factors.push(reference.clone());
        output.push(coefficient);
    }
}

fn render_line(
    prefix: &str,
    terms: &[SignedProduct],
    suffix: &str,
) -> Result<String, RecurrenceError> {
    let zero = [SignedProduct {
        negative: false,
        factors: vec!["0".into()],
    }];
    let terms = if terms.is_empty() { &zero[..] } else { terms };
    let mut lines = Vec::new();
    let mut current = prefix.to_string();
    for (term_index, term) in terms.iter().enumerate() {
        let sign = match (term_index == 0, term.negative) {
            (true, false) => "",
            (true, true) => "-",
            (false, false) => " + ",
            (false, true) => " - ",
        };
        let final_term = term_index + 1 == terms.len();
        let whole_length = current.chars().count()
            + sign.len()
            + product_text(term).chars().count()
            + if final_term {
                suffix.chars().count()
            } else {
                0
            };
        if term_index > 0 && whole_length > 100 {
            lines.push(current);
            current = format!("      {} ", if term.negative { "-" } else { "+" });
        } else {
            current.push_str(sign);
        }
        for (factor_index, factor) in term.factors.iter().enumerate() {
            let separator = if factor_index == 0 { "" } else { " * " };
            let final_factor = final_term && factor_index + 1 == term.factors.len();
            let reserve = if final_factor {
                suffix.chars().count()
            } else {
                0
            };
            if current.chars().count() + separator.len() + factor.chars().count() + reserve > 100 {
                if factor_index == 0 {
                    return Err(invalid("Lean atom exceeds the 100-column output subset"));
                }
                lines.push(current);
                current = format!("      * {factor}");
            } else {
                current.push_str(separator);
                current.push_str(factor);
            }
            if current.chars().count() + reserve > 100 {
                return Err(invalid("Lean atom exceeds the 100-column output subset"));
            }
        }
    }
    current.push_str(suffix);
    lines.push(current);
    Ok(lines.join("\n"))
}

impl RecurrenceJson {
    /// Emit a Lean `ℕ → ℝ[X]` definition and `B + 2` expanded examples.
    ///
    /// `rows` contains every supplied row, starting at this record's
    /// `first_index`; Lean renumbers these from zero. Seed count B must already
    /// have been selected with [`Self::with_checked_presentation`]. The rule is
    /// shifted in every coefficient and alternating sign. Division occurs only
    /// on reals inside `C`; t-dependent leaders and leaders with zeros at or
    /// beyond the recursive boundary are rejected. The emitted text is parsed
    /// back and all supplied rows and examples are checked exactly, with Lean's
    /// totalized real division and natural arithmetic. This safeguard is finite
    /// evaluation, not Lean elaboration or a proof of an all-index identity.
    pub fn to_lean(
        &self,
        name: &str,
        rows: &[Vec<BigRational>],
    ) -> Result<String, RecurrenceError> {
        validate_lean_recurrence_name(name)?;
        let (rule, first, seeds) = self.to_recurrence_parts()?;
        if let Some(reason) = self.presentation.as_ref().and_then(|p| p.reason.as_ref()) {
            return Err(invalid(reason.clone()));
        }
        let shift = i64::try_from(first).map_err(|_| invalid("Lean index shift exceeds i64"))?;
        let rule = rule.shift_index(shift).normalized_presentation()?;
        let base = seeds.len();
        if base < rule.max_offset() {
            return Err(invalid("not enough Lean base cases"));
        }
        let one = BivarPoly {
            coeffs: vec![vec![BigRational::one()]],
        };
        let leader = rule.denominator.as_ref().unwrap_or(&one);
        if leader
            .coeffs
            .iter()
            .any(|r| r.iter().skip(1).any(|c| !c.is_zero()))
        {
            return Err(invalid(
                "Lean output unavailable: left-hand factor depends on t",
            ));
        }
        if rule
            .left_hand_nonnegative_integer_zeros()?
            .iter()
            .any(|n| n >= &BigInt::from(base))
        {
            return Err(invalid(
                "Lean output unavailable: left-hand factor vanishes at or beyond B",
            ));
        }
        if rule.terms.iter().any(|term| term.offset == 0) {
            return Err(invalid(
                "Lean output requires strictly earlier recursive references",
            ));
        }
        let examples_count = base
            .checked_add(2)
            .ok_or_else(|| invalid("Lean base count overflow"))?;
        let examples = rule
            .generate_rows_rational(&seeds, 0, examples_count)
            .map_err(|e| invalid(e.to_string()))?;
        let base_shift = i64::try_from(base).map_err(|_| invalid("Lean base count exceeds i64"))?;
        // Do not shift alternating signs here; they are emitted explicitly at n+B.
        let shifted_leader = leader.shift_index(base_shift);
        let denominator: Vec<_> = shifted_leader.coeffs.iter().map(|r| r[0].clone()).collect();
        let mut summands = Vec::new();
        for term in &rule.terms {
            let offset = base - term.offset;
            let mut reference = if offset == 0 {
                format!("{name} n")
            } else {
                format!("{name} (n + {offset})")
            };
            for _ in 0..term.deriv_order {
                reference = format!("({reference}).derivative");
            }
            let coefficients =
                coefficient_products(&term.coeff.shift_index(base_shift), &denominator);
            let sign = (term.sign == RecurrenceSign::AlternatingN)
                .then(|| format!("C ((-1 : ℝ) ^ (n + {base}))"));
            append_recurrence_products(&mut summands, coefficients, sign, reference);
        }
        if let Some(forcing) = &rule.inhomogeneous {
            summands.extend(coefficient_products(
                &forcing.shift_index(base_shift),
                &denominator,
            ));
        }
        let mut lines = vec![
            "-- Requires Mathlib.".into(),
            "noncomputable section".into(),
            "open Polynomial".into(),
            "-- Finite fit: validity for all indices is unproved.".into(),
            format!("def {name} : ℕ → ℝ[X]"),
        ];
        for (k, row) in seeds.iter().enumerate() {
            lines.push(render_line(
                &format!("  | {k} => "),
                &row_products(row),
                "",
            )?);
        }
        lines.push(render_line(&format!("  | n + {base} => "), &summands, "")?);
        lines.push(String::new());
        for (k, row) in examples.iter().enumerate() {
            lines.push(render_line(
                &format!("example : {name} {k} = "),
                &row_products(row),
                " := by",
            )?);
            lines.push(if k < base { "  rfl" } else { "  sorry" }.into());
        }
        lines.push(String::new());
        lines.push("end".into());
        let text = lines.join("\n");
        verify_lean_recurrence_output(&text, name, rows)?;
        Ok(text)
    }

    /// Attach checked Lean output, or a reason for refusing it. The original
    /// fitted recurrence and every other presentation remain available.
    pub fn with_lean(mut self, name: &str, rows: &[Vec<BigRational>]) -> Self {
        match self.to_lean(name, rows) {
            Ok(text) => {
                self.lean = Some(text);
                self.lean_reason = None;
            }
            Err(error) => {
                self.lean = None;
                self.lean_reason = Some(error.to_string());
            }
        }
        self
    }
}
