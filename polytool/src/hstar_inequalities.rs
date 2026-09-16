//! Named inequalities for Ehrhart h*-vectors.

use num_bigint::BigInt;
use num_rational::Ratio;
use num_traits::{One, Signed, Zero};

type Q = Ratio<BigInt>;

/// One named h*-condition or inequality.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct HStarInequalityCheck {
    pub family: String,
    pub name: String,
    pub formula: String,
    pub reference: String,
    pub url: Option<String>,
    pub index: Option<usize>,
    pub applicable: bool,
    pub holds: bool,
    pub lhs: Option<BigInt>,
    pub rhs: Option<BigInt>,
    pub value: Option<String>,
    pub details: String,
}

/// Complete h*-inequality report for one vector.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct HStarInequalityReport {
    pub hstar: Vec<BigInt>,
    pub dimension: usize,
    pub degree: usize,
    pub codegree: Option<usize>,
    pub all_applicable_hold: bool,
    pub checks: Vec<HStarInequalityCheck>,
}

fn check(
    family: &str,
    name: &str,
    formula: &str,
    reference: &str,
    url: Option<&str>,
    applicable: bool,
    holds: bool,
) -> HStarInequalityCheck {
    HStarInequalityCheck {
        family: family.to_string(),
        name: name.to_string(),
        formula: formula.to_string(),
        reference: reference.to_string(),
        url: url.map(ToString::to_string),
        index: None,
        applicable,
        holds,
        lhs: None,
        rhs: None,
        value: None,
        details: String::new(),
    }
}

fn coeff(hstar: &[BigInt], index: usize) -> BigInt {
    hstar.get(index).cloned().unwrap_or_else(BigInt::zero)
}

fn pad_to_dimension(mut hstar: Vec<BigInt>, dimension: usize) -> Vec<BigInt> {
    if hstar.len() <= dimension {
        hstar.resize(dimension + 1, BigInt::zero());
    }
    hstar
}

fn degree(hstar: &[BigInt]) -> usize {
    hstar.iter().rposition(|c| !c.is_zero()).unwrap_or(0)
}

fn sum_range(hstar: &[BigInt], start: usize, end: usize) -> BigInt {
    if start > end {
        return BigInt::zero();
    }
    (start..=end).fold(BigInt::zero(), |acc, i| acc + coeff(hstar, i))
}

fn q(value: i64) -> Q {
    Q::from_integer(BigInt::from(value))
}

fn format_q(value: &Q) -> String {
    if value.denom() == &BigInt::one() {
        value.to_integer().to_string()
    } else {
        format!("{}/{}", value.numer(), value.denom())
    }
}

fn poly_mul(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![BigInt::zero(); a.len() + b.len() - 1];
    for (i, ca) in a.iter().enumerate() {
        if ca.is_zero() {
            continue;
        }
        for (j, cb) in b.iter().enumerate() {
            if !cb.is_zero() {
                out[i + j] += ca * cb;
            }
        }
    }
    while out.last().is_some_and(|c| c.is_zero()) {
        out.pop();
    }
    out
}

fn solve_linear_rational(rows: Vec<Vec<Q>>, variables: usize) -> Option<Vec<Q>> {
    let mut matrix = rows;
    let mut pivot_row = 0usize;
    let mut pivot_cols = Vec::new();

    for col in 0..variables {
        let Some(found) = (pivot_row..matrix.len()).find(|&row| !matrix[row][col].is_zero()) else {
            continue;
        };
        matrix.swap(pivot_row, found);
        let pivot = matrix[pivot_row][col].clone();
        for entry in matrix[pivot_row].iter_mut().skip(col) {
            *entry = entry.clone() / pivot.clone();
        }
        let normalized_pivot = matrix[pivot_row].clone();
        for (row_index, row) in matrix.iter_mut().enumerate() {
            if row_index == pivot_row || row[col].is_zero() {
                continue;
            }
            let factor = row[col].clone();
            for (entry, pivot_entry) in row.iter_mut().zip(&normalized_pivot).skip(col) {
                *entry = entry.clone() - factor.clone() * pivot_entry.clone();
            }
        }
        pivot_cols.push(col);
        pivot_row += 1;
        if pivot_row == matrix.len() {
            break;
        }
    }

    for row in &matrix {
        if row[..variables].iter().all(|entry| entry.is_zero()) && !row[variables].is_zero() {
            return None;
        }
    }

    let mut solution = vec![q(0); variables];
    for (row, &col) in pivot_cols.iter().enumerate() {
        solution[col] = matrix[row][variables].clone();
    }
    Some(solution)
}

fn stapledon_decomposition_shifted(
    hstar: &[BigInt],
    dimension: usize,
    degree: usize,
) -> Option<(Vec<Q>, Vec<Q>)> {
    if degree > dimension {
        return None;
    }
    let codegree = dimension + 1 - degree;
    let multiplier = vec![BigInt::one(); codegree];
    let mut g = poly_mul(hstar, &multiplier);
    g.resize(dimension + 1, BigInt::zero());

    let a_vars = dimension / 2 + 1;
    let b_degree = degree.checked_sub(1);
    let b_vars = b_degree.map_or(0, |bd| bd / 2 + 1);
    let variables = a_vars + b_vars;
    let mut rows = Vec::with_capacity(dimension + 1);

    for (i, coefficient) in g.iter().enumerate().take(dimension + 1) {
        let mut row = vec![q(0); variables + 1];
        let a_col = i.min(dimension - i);
        row[a_col] = q(1);
        if let Some(bd) = b_degree {
            if i >= codegree {
                let j = i - codegree;
                if j <= bd {
                    let b_col = a_vars + j.min(bd - j);
                    row[b_col] = row[b_col].clone() + q(1);
                }
            }
        }
        row[variables] = Q::from_integer(coefficient.clone());
        rows.push(row);
    }

    let solution = solve_linear_rational(rows, variables)?;
    let mut a = vec![q(0); dimension + 1];
    for i in 0..=dimension {
        a[i] = solution[i.min(dimension - i)].clone();
    }

    let b = if let Some(bd) = b_degree {
        let mut b = vec![q(0); bd + 1];
        for j in 0..=bd {
            b[j] = solution[a_vars + j.min(bd - j)].clone();
        }
        b
    } else {
        Vec::new()
    };
    Some((a, b))
}

fn add_basic_checks(checks: &mut Vec<HStarInequalityCheck>, hstar: &[BigInt], dimension: usize) {
    let mut h0 = check(
        "basic",
        "constant term",
        "h*_0 = 1",
        "Balletti-Higashitani, §1.1; standard Ehrhart h*-vector properties",
        Some("https://arxiv.org/pdf/1703.09600"),
        true,
        coeff(hstar, 0) == BigInt::one(),
    );
    h0.lhs = Some(coeff(hstar, 0));
    h0.rhs = Some(BigInt::one());
    checks.push(h0);

    for i in 0..=dimension {
        let c = coeff(hstar, i);
        let mut nonnegative = check(
            "basic",
            "nonnegativity",
            "h*_i >= 0",
            "Balletti-Higashitani, §1.1, citing Stanley nonnegativity",
            Some("https://arxiv.org/pdf/1703.09600"),
            true,
            !c.is_negative(),
        );
        nonnegative.index = Some(i);
        nonnegative.lhs = Some(c);
        nonnegative.rhs = Some(BigInt::zero());
        if !nonnegative.holds {
            nonnegative.details = format!("negative coefficient at index {i}");
        }
        checks.push(nonnegative);
    }
}

fn add_known_inequalities(
    checks: &mut Vec<HStarInequalityCheck>,
    hstar: &[BigInt],
    dimension: usize,
    degree: usize,
) {
    if dimension >= 1 {
        let lhs = coeff(hstar, dimension);
        let rhs = coeff(hstar, 1);
        let mut c = check(
            "Ehrhart interpretation",
            "interior/lattice-point inequality",
            "h*_d <= h*_1",
            "Balletti-Higashitani, Prop. 2.1",
            Some("https://arxiv.org/pdf/1703.09600"),
            true,
            lhs <= rhs,
        );
        c.lhs = Some(lhs);
        c.rhs = Some(rhs);
        checks.push(c);
    }

    if dimension >= 2 {
        for i in 1..=dimension / 2 {
            let lhs = sum_range(hstar, dimension - i, dimension - 1);
            let rhs = sum_range(hstar, 2, i + 1);
            let mut c = check(
                "Hibi",
                "tail-sum inequality",
                "h*_{d-1}+...+h*_{d-i} <= h*_2+...+h*_{i+1}",
                "Balletti-Higashitani, §2.1, citing Hibi",
                Some("https://arxiv.org/pdf/1703.09600"),
                true,
                lhs <= rhs,
            );
            c.index = Some(i);
            c.lhs = Some(lhs);
            c.rhs = Some(rhs);
            checks.push(c);
        }
    }

    for i in 0..=degree / 2 {
        let lhs = sum_range(hstar, 0, i);
        let rhs = sum_range(hstar, degree.saturating_sub(i), degree);
        let mut c = check(
            "Stanley",
            "partial-sum inequality",
            "h*_0+...+h*_i <= h*_s+...+h*_{s-i}",
            "Balletti-Higashitani, §2.1, citing Stanley [Sta93]",
            Some("https://arxiv.org/pdf/1703.09600"),
            true,
            lhs <= rhs,
        );
        c.index = Some(i);
        c.lhs = Some(lhs);
        c.rhs = Some(rhs);
        checks.push(c);
    }

    let has_interior_points = coeff(hstar, dimension).is_positive();
    if has_interior_points && dimension >= 2 {
        for i in 1..=dimension - 1 {
            let lhs = coeff(hstar, 1);
            let rhs = coeff(hstar, i);
            let mut c = check(
                "Hibi",
                "interior-point inequality",
                "if h*_d > 0, then h*_1 <= h*_i",
                "Balletti-Higashitani, §2.1, citing Hibi",
                Some("https://arxiv.org/pdf/1703.09600"),
                true,
                lhs <= rhs,
            );
            c.index = Some(i);
            c.lhs = Some(lhs);
            c.rhs = Some(rhs);
            checks.push(c);
        }
    } else {
        let mut c = check(
            "Hibi",
            "interior-point inequality",
            "if h*_d > 0, then h*_1 <= h*_i",
            "Balletti-Higashitani, §2.1, citing Hibi",
            Some("https://arxiv.org/pdf/1703.09600"),
            false,
            true,
        );
        c.details = "not applicable because h*_d = 0 or dimension < 2".to_string();
        checks.push(c);
    }
}

fn add_balletti_higashitani(checks: &mut Vec<HStarInequalityCheck>, hstar: &[BigInt]) {
    let h1 = coeff(hstar, 1);
    let h2 = coeff(hstar, 2);
    let h3 = coeff(hstar, 3);
    let applicable = h3.is_zero();
    let scott_i = h2.is_zero();
    let scott_ii_rhs = BigInt::from(3) * &h2 + BigInt::from(3);
    let scott_ii = h1 <= scott_ii_rhs;
    let scott_iii = h1 == BigInt::from(7) && h2 == BigInt::one();
    let holds = !applicable || scott_i || scott_ii || scott_iii;
    let mut c = check(
        "Balletti-Higashitani",
        "universal Scott inequality",
        "if h*_3 = 0, then h*_2 = 0, or h*_1 <= 3 h*_2 + 3, or (h*_1,h*_2)=(7,1)",
        "Balletti-Higashitani, Theorem 1.4",
        Some("https://arxiv.org/pdf/1703.09600"),
        applicable,
        holds,
    );
    c.lhs = Some(h1);
    c.rhs = Some(scott_ii_rhs);
    c.details = if applicable {
        format!("h*_2=0: {scott_i}; h*_1 <= 3h*_2+3: {scott_ii}; exceptional (7,1): {scott_iii}")
    } else {
        "not applicable because h*_3 != 0".to_string()
    };
    checks.push(c);
}

fn add_stapledon(
    checks: &mut Vec<HStarInequalityCheck>,
    hstar: &[BigInt],
    dimension: usize,
    degree: usize,
) {
    if degree > dimension {
        let mut c = check(
            "Stapledon",
            "shifted decomposition",
            "(1+...+t^{l-1})h*(t)=a(t)+t^l b(t), with a,b nonnegative symmetric",
            "Stapledon inequalities, summarized in Balletti-Higashitani, §2.1",
            Some("https://arxiv.org/pdf/1703.09600"),
            false,
            true,
        );
        c.details = "not applicable because degree exceeds dimension".to_string();
        checks.push(c);
        return;
    }

    let Some((a, b)) = stapledon_decomposition_shifted(hstar, dimension, degree) else {
        let mut c = check(
            "Stapledon",
            "shifted decomposition",
            "(1+...+t^{l-1})h*(t)=a(t)+t^l b(t)",
            "Stapledon inequalities, summarized in Balletti-Higashitani, §2.1",
            Some("https://arxiv.org/pdf/1703.09600"),
            true,
            false,
        );
        c.details = "could not solve the shifted Stapledon decomposition".to_string();
        checks.push(c);
        return;
    };

    let families = [("a", a), ("b", b)];
    for (label, coefficients) in families {
        if coefficients.is_empty() {
            continue;
        }
        for (i, value) in coefficients.iter().enumerate() {
            let integer = value.denom() == &BigInt::one();
            let nonnegative = integer && !value.numer().is_negative();
            let mut c = check(
                "Stapledon",
                &format!("{label}-coefficient nonnegativity"),
                "(1+...+t^{l-1})h*(t)=a(t)+t^l b(t), with a,b nonnegative symmetric",
                "Stapledon inequalities, summarized in Balletti-Higashitani, §2.1",
                Some("https://arxiv.org/pdf/1703.09600"),
                true,
                nonnegative,
            );
            c.index = Some(i);
            c.value = Some(format_q(value));
            c.details = format!("{label}_{i} = {}", format_q(value));
            checks.push(c);
        }
    }
}

/// Check the standard named h*-inequalities for a supplied dimension.
pub fn hstar_inequality_report_bigint(
    hstar: &[BigInt],
    dimension: Option<usize>,
) -> HStarInequalityReport {
    let inferred_dimension = hstar.len().saturating_sub(1);
    let dimension = dimension.unwrap_or(inferred_dimension);
    let mut h = pad_to_dimension(hstar.to_vec(), dimension);
    let degree = degree(&h);
    let mut checks = Vec::new();

    add_basic_checks(&mut checks, &h, dimension);

    if degree > dimension {
        let mut c = check(
            "basic",
            "dimension bound",
            "deg h*(t) <= d",
            "Standard Ehrhart h*-vector convention",
            None,
            true,
            false,
        );
        c.lhs = Some(BigInt::from(degree));
        c.rhs = Some(BigInt::from(dimension));
        checks.push(c);
    } else {
        add_known_inequalities(&mut checks, &h, dimension, degree);
        add_balletti_higashitani(&mut checks, &h);
        add_stapledon(&mut checks, &h, dimension, degree);
    }

    let codegree = (degree <= dimension).then_some(dimension + 1 - degree);
    let all_applicable_hold = checks
        .iter()
        .filter(|check| check.applicable)
        .all(|check| check.holds);
    while h.len() > 1 && h.last().is_some_and(|c| c.is_zero()) {
        h.pop();
    }
    HStarInequalityReport {
        hstar: h,
        dimension,
        degree,
        codegree,
        all_applicable_hold,
        checks,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn b(values: &[i64]) -> Vec<BigInt> {
        values.iter().map(|&v| BigInt::from(v)).collect()
    }

    #[test]
    fn simplex_hstar_passes_basic_named_checks() {
        let report = hstar_inequality_report_bigint(&b(&[1, 0, 0]), Some(2));
        assert!(report.all_applicable_hold);
    }

    #[test]
    fn eulerian_square_hstar_passes() {
        let report = hstar_inequality_report_bigint(&b(&[1, 4, 1]), Some(3));
        assert!(report.all_applicable_hold);
    }

    #[test]
    fn reports_stanley_failure() {
        let report = hstar_inequality_report_bigint(&b(&[1, 1, 0, 1]), Some(3));
        assert!(report
            .checks
            .iter()
            .any(|check| check.family == "Stanley" && !check.holds));
    }

    #[test]
    fn reports_balletti_higashitani_failure() {
        let report = hstar_inequality_report_bigint(&b(&[1, 20, 1, 0]), Some(3));
        assert!(report.checks.iter().any(|check| {
            check.family == "Balletti-Higashitani" && check.applicable && !check.holds
        }));
    }

    #[test]
    fn explicit_dimension_does_not_discard_higher_coefficients() {
        let report = hstar_inequality_report_bigint(&b(&[1, 0, 1]), Some(1));
        assert_eq!(report.hstar, b(&[1, 0, 1]));
        assert_eq!(report.degree, 2);
        assert_eq!(report.codegree, None);
        assert!(!report.all_applicable_hold);
        assert!(report
            .checks
            .iter()
            .any(|check| check.name == "dimension bound" && !check.holds));
    }
}
