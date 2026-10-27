//! Export of recurrences as Mathematica, Sage and Python definitions and as
//! LaTeX.

use super::*;
use num_traits::One;
use num_traits::Zero;

#[derive(Copy, Clone)]
enum CodeStyle {
    Mathematica,
    Sage,
}

fn recurrence_sign_code(style: CodeStyle, sign: RecurrenceSign) -> Option<&'static str> {
    match (style, sign) {
        (_, RecurrenceSign::None) => None,
        (CodeStyle::Mathematica, RecurrenceSign::AlternatingN) => Some("(-1)^n"),
        (CodeStyle::Sage, RecurrenceSign::AlternatingN) => Some("(-1)**n"),
    }
}

fn fmt_rational_code(style: CodeStyle, r: &BigRational) -> String {
    match style {
        CodeStyle::Mathematica => {
            if r.denom() == &BigInt::one() {
                format!("{}", r.numer())
            } else {
                format!("{}/{}", r.numer(), r.denom())
            }
        }
        CodeStyle::Sage => {
            if r.denom() == &BigInt::one() {
                format!("{}", r.numer())
            } else if r.numer() < &BigInt::zero() {
                format!("-QQ({})/{}", -r.numer().clone(), r.denom())
            } else {
                format!("QQ({})/{}", r.numer(), r.denom())
            }
        }
    }
}

fn fmt_power_code(style: CodeStyle, var: &str, pow: usize) -> Option<String> {
    match pow {
        0 => None,
        1 => Some(var.to_string()),
        _ => Some(match style {
            CodeStyle::Mathematica => format!("{var}^{pow}"),
            CodeStyle::Sage => format!("{var}**{pow}"),
        }),
    }
}

fn fmt_monomial_code_abs(
    style: CodeStyle,
    coeff_abs: &BigRational,
    n_pow: usize,
    t_pow: usize,
) -> String {
    let has_vars = n_pow > 0 || t_pow > 0;
    let mut factors = Vec::new();

    if !(coeff_abs.is_one() && has_vars) {
        let coeff_text = fmt_rational_code(style, coeff_abs);
        let needs_wrap = has_vars && coeff_text.contains('/');
        factors.push(if needs_wrap {
            format!("({coeff_text})")
        } else {
            coeff_text
        });
    }

    if let Some(n_term) = fmt_power_code(style, "n", n_pow) {
        factors.push(n_term);
    }
    if let Some(t_term) = fmt_power_code(style, "t", t_pow) {
        factors.push(t_term);
    }

    if factors.is_empty() {
        "1".to_string()
    } else {
        factors.join("*")
    }
}

fn join_signed_terms(terms: &[(bool, String)]) -> String {
    if terms.is_empty() {
        return "0".to_string();
    }

    let mut result = String::new();
    for (idx, (negative, body)) in terms.iter().enumerate() {
        if idx == 0 {
            if *negative {
                result.push('-');
            }
            result.push_str(body);
        } else if *negative {
            result.push_str(" - ");
            result.push_str(body);
        } else {
            result.push_str(" + ");
            result.push_str(body);
        }
    }
    result
}

fn wrap_if_sum(expr: &str) -> String {
    if expr.contains(" + ") || expr.contains(" - ") {
        format!("({expr})")
    } else {
        expr.to_string()
    }
}

fn fmt_univariate_poly_rational_code(style: CodeStyle, coeffs: &[BigRational]) -> String {
    let mut terms = Vec::new();
    for (pow, coeff) in coeffs.iter().enumerate() {
        if coeff.is_zero() {
            continue;
        }
        let negative = coeff < &BigRational::zero();
        let abs_coeff = if negative {
            -coeff.clone()
        } else {
            coeff.clone()
        };
        let body = fmt_monomial_code_abs(style, &abs_coeff, 0, pow);
        terms.push((negative, body));
    }
    join_signed_terms(&terms)
}

fn fmt_python_rational(r: &BigRational) -> String {
    if r.denom() == &BigInt::one() {
        format!("{}", r.numer())
    } else {
        format!("Fraction({}, {})", r.numer(), r.denom())
    }
}

fn fmt_python_poly_literal(coeffs: &[BigRational]) -> String {
    format!(
        "[{}]",
        coeffs
            .iter()
            .map(fmt_python_rational)
            .collect::<Vec<_>>()
            .join(", ")
    )
}

fn fmt_python_bivar_literal(poly: &BivarPoly) -> String {
    format!(
        "[{}]",
        poly.coeffs
            .iter()
            .map(|row| fmt_python_poly_literal(row))
            .collect::<Vec<_>>()
            .join(", ")
    )
}

fn python_sign_name(sign: RecurrenceSign) -> &'static str {
    match sign {
        RecurrenceSign::None => "none",
        RecurrenceSign::AlternatingN => "alternating_n",
    }
}

fn python_term_literal(term: &RecurrenceTerm) -> String {
    format!(
        "{{\"offset\": {}, \"deriv_order\": {}, \"sign\": \"{}\", \"coeff\": {}}}",
        term.offset,
        term.deriv_order,
        python_sign_name(term.sign),
        fmt_python_bivar_literal(&term.coeff)
    )
}

impl Recurrence {
    fn term_to_mathematica_code(&self, term: &RecurrenceTerm) -> String {
        let coeff = term.coeff.to_mathematica_code();
        let idx = if term.offset == 0 {
            "n".to_string()
        } else {
            format!("n - {}", term.offset)
        };
        let pref = match term.deriv_order {
            0 => format!("P[{idx}, t]"),
            1 => format!("D[P[{idx}, t], t]"),
            d => format!("D[P[{idx}, t], {{t, {d}}}]"),
        };
        let sign = recurrence_sign_code(CodeStyle::Mathematica, term.sign);

        if term.coeff.is_one() {
            if let Some(sign) = sign {
                format!("{sign}*{pref}")
            } else {
                pref
            }
        } else if coeff == "-1" {
            if let Some(sign) = sign {
                format!("-{sign}*{pref}")
            } else {
                format!("-{pref}")
            }
        } else if let Some(sign) = sign {
            format!("{sign}*{}*{pref}", wrap_if_sum(&coeff))
        } else {
            format!("{}*{pref}", wrap_if_sum(&coeff))
        }
    }

    fn term_to_sage_code(&self, term: &RecurrenceTerm) -> String {
        let coeff = term.coeff.to_sage_code();
        let idx = if term.offset == 0 {
            "n".to_string()
        } else {
            format!("n - {}", term.offset)
        };
        let pref = match term.deriv_order {
            0 => format!("P({idx})"),
            1 => format!("P({idx}).derivative(t)"),
            d => format!("P({idx}).derivative(t, {d})"),
        };
        let sign = recurrence_sign_code(CodeStyle::Sage, term.sign);

        if term.coeff.is_one() {
            if let Some(sign) = sign {
                format!("{sign}*{pref}")
            } else {
                pref
            }
        } else if coeff == "-1" {
            if let Some(sign) = sign {
                format!("-{sign}*{pref}")
            } else {
                format!("-{pref}")
            }
        } else if let Some(sign) = sign {
            format!("{sign}*{}*{pref}", wrap_if_sum(&coeff))
        } else {
            format!("{}*{pref}", wrap_if_sum(&coeff))
        }
    }

    fn rhs_to_mathematica_code(&self) -> String {
        let mut pieces: Vec<String> = self
            .terms
            .iter()
            .map(|term| self.term_to_mathematica_code(term))
            .collect();
        if let Some(ref inh) = self.inhomogeneous {
            pieces.push(inh.to_mathematica_code());
        }
        let rhs = if pieces.is_empty() {
            "0".to_string()
        } else {
            pieces.join(" + ")
        };
        if let Some(ref denom) = self.denominator {
            format!(
                "Expand[Together[({})/({})]]",
                rhs,
                denom.to_mathematica_code()
            )
        } else {
            format!("Expand[{rhs}]")
        }
    }

    fn rhs_to_sage_code(&self) -> String {
        let mut pieces: Vec<String> = self
            .terms
            .iter()
            .map(|term| self.term_to_sage_code(term))
            .collect();
        if let Some(ref inh) = self.inhomogeneous {
            pieces.push(inh.to_sage_code());
        }
        let rhs = if pieces.is_empty() {
            "0".to_string()
        } else {
            pieces.join(" + ")
        };
        if let Some(ref denom) = self.denominator {
            format!("R(K({rhs}) / K({}))", denom.to_sage_code())
        } else {
            format!("R({rhs})")
        }
    }

    fn to_mathematica_definition_rational_indexed_prefix(
        &self,
        initial_polys: &[Vec<BigRational>],
        first_index: i64,
        base_count: usize,
    ) -> String {
        let start_n = i128::from(first_index) + base_count as i128;
        let mut lines = vec!["ClearAll[P];".to_string()];
        for (idx, coeffs) in initial_polys.iter().take(base_count).enumerate() {
            lines.push(format!(
                "P[{}, t_] := {};",
                i128::from(first_index) + idx as i128,
                fmt_univariate_poly_rational_code(CodeStyle::Mathematica, coeffs)
            ));
        }
        lines.push(format!(
            "P[n_Integer /; n >= {start_n}, t_] := P[n, t] = {};",
            self.rhs_to_mathematica_code()
        ));
        lines.push(String::new());
        lines.push(format!(
            "(* Example: Table[P[n, t], {{n, {}, {}}}] *)",
            first_index,
            i128::from(first_index) + 9
        ));
        lines.join("\n")
    }

    /// Export a Mathematica definition of `P[n, t]` with first index 1.
    ///
    /// The first `max_offset()` entries of `initial_polys` (all of them, if fewer are
    /// supplied) become base cases and the recurrence applies from the next index. Use
    /// [`Self::to_mathematica_definition_rational_indexed`] to keep every supplied row and
    /// choose the first index.
    pub fn to_mathematica_definition_rational(&self, initial_polys: &[Vec<BigRational>]) -> String {
        let base_count = self.max_offset().min(initial_polys.len());
        self.to_mathematica_definition_rational_indexed_prefix(initial_polys, 1, base_count)
    }

    /// Export an indexed definition which retains every supplied initial row.
    pub fn to_mathematica_definition_rational_indexed(
        &self,
        initial_polys: &[Vec<BigRational>],
        first_index: i64,
    ) -> String {
        self.to_mathematica_definition_rational_indexed_prefix(
            initial_polys,
            first_index,
            initial_polys.len(),
        )
    }

    /// Integer-coefficient wrapper for [`Self::to_mathematica_definition_rational`].
    pub fn to_mathematica_definition(&self, initial_polys: &[Vec<i64>]) -> String {
        let rational_polys = i64_polys_to_rational(initial_polys);
        self.to_mathematica_definition_rational(&rational_polys)
    }

    fn to_sage_definition_rational_indexed_prefix(
        &self,
        initial_polys: &[Vec<BigRational>],
        first_index: i64,
        base_count: usize,
    ) -> String {
        let mut lines = vec![
            "R.<t> = PolynomialRing(QQ)".to_string(),
            "K = R.fraction_field()".to_string(),
            "_P_cache = {".to_string(),
        ];
        for (idx, coeffs) in initial_polys.iter().take(base_count).enumerate() {
            let coeff_list = coeffs
                .iter()
                .map(|coeff| fmt_rational_code(CodeStyle::Sage, coeff))
                .collect::<Vec<_>>()
                .join(", ");
            lines.push(format!(
                "    {}: R([{}]),",
                i128::from(first_index) + idx as i128,
                coeff_list
            ));
        }
        lines.push("}".to_string());
        lines.push(String::new());
        lines.push("def P(n):".to_string());
        lines.push(format!("    if n < {first_index}:"));
        lines.push(format!(
            "        raise ValueError(\"n must be at least {first_index}\")"
        ));
        lines.push("    if n in _P_cache:".to_string());
        lines.push("        return _P_cache[n]".to_string());
        lines.push(format!("    value = {}", self.rhs_to_sage_code()));
        lines.push("    _P_cache[n] = value".to_string());
        lines.push("    return value".to_string());
        lines.push(String::new());
        lines.push(format!(
            "# Example: [P(n) for n in range({}, {})]",
            first_index,
            i128::from(first_index) + 10
        ));
        lines.join("\n")
    }

    /// Export a Sage definition of `P(n)` over `QQ[t]` with first index 1.
    ///
    /// The first `max_offset()` entries of `initial_polys` (all of them, if fewer are
    /// supplied) become base cases and the recurrence applies from the next index. Use
    /// [`Self::to_sage_definition_rational_indexed`] to keep every supplied row and choose
    /// the first index.
    pub fn to_sage_definition_rational(&self, initial_polys: &[Vec<BigRational>]) -> String {
        let base_count = self.max_offset().min(initial_polys.len());
        self.to_sage_definition_rational_indexed_prefix(initial_polys, 1, base_count)
    }

    /// Export an indexed definition which retains every supplied initial row.
    pub fn to_sage_definition_rational_indexed(
        &self,
        initial_polys: &[Vec<BigRational>],
        first_index: i64,
    ) -> String {
        self.to_sage_definition_rational_indexed_prefix(
            initial_polys,
            first_index,
            initial_polys.len(),
        )
    }

    /// Integer-coefficient wrapper for [`Self::to_sage_definition_rational`].
    pub fn to_sage_definition(&self, initial_polys: &[Vec<i64>]) -> String {
        let rational_polys = i64_polys_to_rational(initial_polys);
        self.to_sage_definition_rational(&rational_polys)
    }

    fn to_python_definition_rational_indexed_prefix(
        &self,
        initial_polys: &[Vec<BigRational>],
        first_index: i64,
        base_count: usize,
    ) -> String {
        let mut lines = vec![
            "from fractions import Fraction".to_string(),
            String::new(),
            "def _trim(p):".to_string(),
            "    p = list(p)".to_string(),
            "    while len(p) > 1 and p[-1] == 0:".to_string(),
            "        p.pop()".to_string(),
            "    return p or [0]".to_string(),
            String::new(),
            "def _is_zero(p):".to_string(),
            "    return all(coeff == 0 for coeff in p)".to_string(),
            String::new(),
            "def _add(a, b):".to_string(),
            "    length = max(len(a), len(b))".to_string(),
            "    out = [0] * length".to_string(),
            "    for i in range(length):".to_string(),
            "        out[i] = (a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)"
                .to_string(),
            "    return _trim(out)".to_string(),
            String::new(),
            "def _scale(p, c):".to_string(),
            "    return _trim([c * coeff for coeff in p])".to_string(),
            String::new(),
            "def _mul(a, b):".to_string(),
            "    if _is_zero(a) or _is_zero(b):".to_string(),
            "        return [0]".to_string(),
            "    out = [0] * (len(a) + len(b) - 1)".to_string(),
            "    for i, ai in enumerate(a):".to_string(),
            "        for j, bj in enumerate(b):".to_string(),
            "            out[i + j] += ai * bj".to_string(),
            "    return _trim(out)".to_string(),
            String::new(),
            "def _derivative(p, order):".to_string(),
            "    out = list(p)".to_string(),
            "    for _ in range(order):".to_string(),
            "        if len(out) <= 1:".to_string(),
            "            return [0]".to_string(),
            "        out = [(i + 1) * out[i + 1] for i in range(len(out) - 1)]".to_string(),
            "    return _trim(out)".to_string(),
            String::new(),
            "def _div_exact(numerator, denominator):".to_string(),
            "    denominator = _trim(denominator)".to_string(),
            "    if _is_zero(denominator):".to_string(),
            "        raise ZeroDivisionError(\"zero recurrence denominator\")".to_string(),
            "    remainder = _trim(numerator)".to_string(),
            "    if _is_zero(remainder):".to_string(),
            "        return [0]".to_string(),
            "    if len(remainder) < len(denominator):".to_string(),
            "        raise ValueError(\"non-exact polynomial quotient\")".to_string(),
            "    quotient = [0] * (len(remainder) - len(denominator) + 1)".to_string(),
            "    denominator_degree = len(denominator) - 1".to_string(),
            "    denominator_lc = denominator[-1]".to_string(),
            "    while not _is_zero(remainder) and len(remainder) - 1 >= denominator_degree:"
                .to_string(),
            "        shift = len(remainder) - 1 - denominator_degree".to_string(),
            "        factor = remainder[-1] / denominator_lc".to_string(),
            "        quotient[shift] += factor".to_string(),
            "        for i, coeff in enumerate(denominator):".to_string(),
            "            remainder[i + shift] -= factor * coeff".to_string(),
            "        remainder = _trim(remainder)".to_string(),
            "    if not _is_zero(remainder):".to_string(),
            "        raise ValueError(\"non-exact polynomial quotient\")".to_string(),
            "    return _trim(quotient)".to_string(),
            String::new(),
            "def _bivar(coeffs, n):".to_string(),
            "    width = max((len(row) for row in coeffs), default=0)".to_string(),
            "    out = [0] * width".to_string(),
            "    n_power = 1".to_string(),
            "    for row in coeffs:".to_string(),
            "        for t_power, coeff in enumerate(row):".to_string(),
            "            out[t_power] += coeff * n_power".to_string(),
            "        n_power *= n".to_string(),
            "    return _trim(out)".to_string(),
            String::new(),
            "_P_cache = {".to_string(),
        ];
        for (idx, coeffs) in initial_polys.iter().take(base_count).enumerate() {
            lines.push(format!(
                "    {}: {},",
                i128::from(first_index) + idx as i128,
                fmt_python_poly_literal(coeffs)
            ));
        }
        lines.push("}".to_string());
        lines.push(format!(
            "_terms = [{}]",
            self.terms
                .iter()
                .map(python_term_literal)
                .collect::<Vec<_>>()
                .join(", ")
        ));
        lines.push(format!(
            "_denominator = {}",
            self.denominator
                .as_ref()
                .map(fmt_python_bivar_literal)
                .unwrap_or_else(|| "None".to_string())
        ));
        lines.push(format!(
            "_inhomogeneous = {}",
            self.inhomogeneous
                .as_ref()
                .map(fmt_python_bivar_literal)
                .unwrap_or_else(|| "None".to_string())
        ));
        lines.push(String::new());
        lines.push("def P(n):".to_string());
        lines.push(format!("    if n < {first_index}:"));
        lines.push(format!(
            "        raise ValueError(\"n must be at least {first_index}\")"
        ));
        lines.push("    if n in _P_cache:".to_string());
        lines.push("        return _P_cache[n]".to_string());
        lines.push(format!(
            "    start = max(_P_cache) + 1 if _P_cache else {first_index}"
        ));
        lines.push("    for k in range(start, n + 1):".to_string());
        lines.push("        rhs = [0]".to_string());
        lines.push("        for term in _terms:".to_string());
        lines.push(
            "            ref = _derivative(P(k - term[\"offset\"]), term[\"deriv_order\"])"
                .to_string(),
        );
        lines.push("            coeff = _bivar(term[\"coeff\"], k)".to_string());
        lines.push("            product = _mul(coeff, ref)".to_string());
        lines
            .push("            if term[\"sign\"] == \"alternating_n\" and k % 2 == 1:".to_string());
        lines.push("                product = _scale(product, -1)".to_string());
        lines.push("            rhs = _add(rhs, product)".to_string());
        lines.push("        if _inhomogeneous is not None:".to_string());
        lines.push("            rhs = _add(rhs, _bivar(_inhomogeneous, k))".to_string());
        lines.push("        if _denominator is not None:".to_string());
        lines.push("            rhs = _div_exact(rhs, _bivar(_denominator, k))".to_string());
        lines.push("        _P_cache[k] = rhs".to_string());
        lines.push("    return _P_cache[n]".to_string());
        lines.push(String::new());
        lines.push(format!(
            "# Example: [P(n) for n in range({}, {})]",
            first_index,
            i128::from(first_index) + 10
        ));
        lines.join("\n")
    }

    /// Export standalone Python code defining `P(n)` with exact `Fraction` arithmetic and
    /// first index 1.
    ///
    /// The first `max_offset()` entries of `initial_polys` (all of them, if fewer are
    /// supplied) become base cases and the recurrence applies from the next index. Use
    /// [`Self::to_python_definition_rational_indexed`] to keep every supplied row and
    /// choose the first index.
    pub fn to_python_definition_rational(&self, initial_polys: &[Vec<BigRational>]) -> String {
        let base_count = self.max_offset().min(initial_polys.len());
        self.to_python_definition_rational_indexed_prefix(initial_polys, 1, base_count)
    }

    /// Export an indexed definition which retains every supplied initial row.
    pub fn to_python_definition_rational_indexed(
        &self,
        initial_polys: &[Vec<BigRational>],
        first_index: i64,
    ) -> String {
        self.to_python_definition_rational_indexed_prefix(
            initial_polys,
            first_index,
            initial_polys.len(),
        )
    }

    /// Integer-coefficient wrapper for [`Self::to_python_definition_rational`].
    pub fn to_python_definition(&self, initial_polys: &[Vec<i64>]) -> String {
        let rational_polys = i64_polys_to_rational(initial_polys);
        self.to_python_definition_rational(&rational_polys)
    }

    /// Format the recurrence as a LaTeX expression.
    pub fn to_latex(&self) -> String {
        let mut s = String::new();
        if let Some(ref denom) = self.denominator {
            let dl = denom.to_latex();
            if denom.num_terms() > 1 {
                s.push_str(&format!("\\bigl({dl}\\bigr) "));
            } else {
                s.push_str(&format!("{dl} \\cdot "));
            }
        }
        s.push_str("P(n) = ");

        let mut first = true;
        for term in &self.terms {
            let cl = term.coeff.to_latex();
            let pref = fmt_poly_ref_latex_with_sign(term.sign, term.offset, term.deriv_order);

            if first {
                first = false;
                if cl == "1" {
                    s.push_str(&pref);
                } else if cl == "-1" {
                    s.push_str(&format!("-{pref}"));
                } else if term.coeff.num_terms() > 1 {
                    s.push_str(&format!("\\bigl({cl}\\bigr) {pref}"));
                } else {
                    s.push_str(&format!("{cl} \\cdot {pref}"));
                }
            } else if cl.starts_with('-') && term.coeff.is_monomial() {
                let pos = &cl[1..];
                if pos == "1" {
                    s.push_str(&format!(" - {pref}"));
                } else {
                    s.push_str(&format!(" - {pos} \\cdot {pref}"));
                }
            } else if term.coeff.num_terms() > 1 {
                s.push_str(&format!(" + \\bigl({cl}\\bigr) {pref}"));
            } else if cl == "1" {
                s.push_str(&format!(" + {pref}"));
            } else {
                s.push_str(&format!(" + {cl} \\cdot {pref}"));
            }
        }

        if let Some(ref inh) = self.inhomogeneous {
            if !first {
                s.push_str(" + ");
            }
            s.push_str(&inh.to_latex());
        }

        s
    }
}

impl BivarPoly {
    fn to_code(&self, style: CodeStyle) -> String {
        let mut terms: Vec<(bool, String)> = Vec::new();
        let width = self.coefficient_width();
        if width == 0 {
            return "0".into();
        }
        let max_j = width - 1;
        let max_i = self.coeffs.len().saturating_sub(1);
        for j in 0..=max_j {
            for i in 0..=max_i {
                let Some(c) = self.coefficient(i, j) else {
                    continue;
                };
                if c.is_zero() {
                    continue;
                }
                let negative = c < &BigRational::zero();
                let coeff_abs = if negative { -c.clone() } else { c.clone() };
                terms.push((negative, fmt_monomial_code_abs(style, &coeff_abs, i, j)));
            }
        }
        join_signed_terms(&terms)
    }

    /// Format as a Mathematica expression in `n` and `t` with exact rational coefficients;
    /// the zero polynomial is `0`.
    pub fn to_mathematica_code(&self) -> String {
        self.to_code(CodeStyle::Mathematica)
    }

    /// Format as a Sage expression in `n` and `t`; non-integer coefficients are written
    /// `QQ(a)/b`, and the zero polynomial is `0`.
    pub fn to_sage_code(&self) -> String {
        self.to_code(CodeStyle::Sage)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mathematica_definition_contains_initial_values_and_rule() {
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1],
            vec![2],
            vec![3],
            vec![5],
            vec![8],
            vec![13],
        ];
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 2,
            homogeneous: true,
            ..Default::default()
        };
        let rec = find_polynomial_recurrence(&polys, &opts).expect("should find recurrence");
        let code = rec.to_mathematica_definition(&polys);
        assert!(code.contains("P[1, t_] := 1;"));
        assert!(code.contains("P[2, t_] := 1;"));
        assert!(code.contains(
            "P[n_Integer /; n >= 3, t_] := P[n, t] = Expand[P[n - 1, t] + P[n - 2, t]];"
        ));
    }

    #[test]
    fn sage_definition_contains_derivative_rule() {
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1],
            vec![1, 1],
            vec![1, 4, 1],
            vec![1, 11, 11, 1],
            vec![1, 26, 66, 26, 1],
        ];
        let opts = RecurrenceOptions {
            var_deg: 2,
            idx_deg: 1,
            diff_deg: 1,
            rec_len: 1,
            homogeneous: true,
            ..Default::default()
        };
        let rec = find_polynomial_recurrence(&polys, &opts).expect("should find recurrence");
        let code = rec.to_sage_definition(&polys);
        assert!(code.contains("R.<t> = PolynomialRing(QQ)"));
        assert!(code.contains("_P_cache = {"));
        assert!(code.contains("1: R([1]),"));
        assert!(code.contains("value = R("));
        assert!(code.contains("P(n - 1).derivative(t)"));
    }

    #[test]
    fn python_definition_contains_exact_fraction_rule() {
        let polys: Vec<Vec<i64>> = vec![vec![1], vec![1], vec![2], vec![3], vec![5], vec![8]];
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 2,
            homogeneous: true,
            ..Default::default()
        };
        let rec = find_polynomial_recurrence(&polys, &opts).expect("should find recurrence");
        let code = rec.to_python_definition(&polys);
        assert!(code.contains("from fractions import Fraction"));
        assert!(code.contains("def P(n):"));
        assert!(code.contains("_P_cache = {"));
        assert!(code.contains("1: [1],"));
        assert!(code.contains("\"offset\": 1"));
        assert!(code.contains("_div_exact"));
    }

    #[test]
    fn python_definition_contains_derivative_metadata() {
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1],
            vec![1, 1],
            vec![1, 4, 1],
            vec![1, 11, 11, 1],
            vec![1, 26, 66, 26, 1],
        ];
        let opts = RecurrenceOptions {
            var_deg: 2,
            idx_deg: 1,
            diff_deg: 1,
            rec_len: 1,
            homogeneous: true,
            ..Default::default()
        };
        let rec = find_polynomial_recurrence(&polys, &opts).expect("should find recurrence");
        let code = rec.to_python_definition(&polys);
        assert!(code.contains("\"deriv_order\": 1"));
        assert!(code.contains("ref = _derivative(P(k - term[\"offset\"]),"));
        assert!(code.contains("coeff = _bivar(term[\"coeff\"], k)"));
    }
}
