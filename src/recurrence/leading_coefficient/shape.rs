//! Recurrence shapes and sampled degree-law construction.

use super::*;

// ---------------------------------------------------------------------------
// The rule
// ---------------------------------------------------------------------------

pub(super) struct ShapeTerm {
    pub(super) lag: usize,
    pub(super) derivative_order: usize,
    pub(super) alternating: bool,
    pub(super) coefficients: Bivariate,
}

pub(super) struct Shape {
    pub(super) terms: Vec<ShapeTerm>,
    /// `f(n)`.
    pub(super) left_factor: Vec<Q>,
    pub(super) max_lag: usize,
    pub(super) alternating: bool,
}

pub(super) fn rule_shape(rule: &Recurrence) -> Result<Shape, String> {
    if rule.inhomogeneous.is_some() {
        return Err("the rule has an inhomogeneous term".into());
    }
    let left_factor = match &rule.denominator {
        None => vec![Q::one()],
        Some(denominator) => {
            let coefficients = denominator.normalized_coefficients();
            if (1..x_width(&coefficients)).any(|j| !column(&coefficients, j).is_empty()) {
                return Err("the left-hand factor f(n, x) depends on x".into());
            }
            let f = column(&coefficients, 0);
            if f.is_empty() {
                return Err("the left-hand factor is zero".into());
            }
            f
        }
    };
    let mut terms = Vec::new();
    for term in &rule.terms {
        let coefficients = term.coeff.normalized_coefficients();
        if x_degree(&coefficients).is_none() {
            continue;
        }
        if term.offset == 0 {
            return Err("a term refers to P_n itself (lag 0)".into());
        }
        terms.push(ShapeTerm {
            lag: term.offset,
            derivative_order: term.deriv_order,
            alternating: term.sign == RecurrenceSign::AlternatingN,
            coefficients,
        });
    }
    if terms.is_empty() {
        return Err("the rule has no nonzero term".into());
    }
    let max_lag = terms.iter().map(|t| t.lag).max().unwrap_or(1);
    let alternating = terms.iter().any(|t| t.alternating);
    Ok(Shape {
        terms,
        left_factor,
        max_lag,
        alternating,
    })
}

pub(super) fn x_width(b: &Bivariate) -> usize {
    b.iter().map(Vec::len).max().unwrap_or(0)
}

/// `[x^j] b(n, x)` as a polynomial in `n`, without trailing zeros.
pub(super) fn column(b: &Bivariate, j: usize) -> Vec<Q> {
    let mut column: Vec<Q> = b
        .iter()
        .map(|row| row.get(j).cloned().unwrap_or_else(Q::zero))
        .collect();
    vec_poly::strip_trailing_zeros(&mut column);
    column
}

pub(super) fn x_degree(b: &Bivariate) -> Option<usize> {
    (0..x_width(b)).rev().find(|&j| !column(b, j).is_empty())
}

/// `p(c + period m)` as a polynomial in `m`.
pub(super) fn on_class(p: &[Q], c: usize, period: usize) -> Vec<Q> {
    if p.is_empty() {
        return Vec::new();
    }
    vec_poly::compose(p, &[integer(c as i64), integer(period as i64)]).stripped()
}

/// Polynomials without trailing zeros (the zero polynomial is empty).
pub(super) trait StripTrailingZeros {
    fn stripped(self) -> Self;
}

impl StripTrailingZeros for Vec<Q> {
    fn stripped(mut self) -> Self {
        vec_poly::strip_trailing_zeros(&mut self);
        self
    }
}

pub(super) fn integer(value: i64) -> Q {
    Q::from_integer(BigInt::from(value))
}

/// `(base + growth m)_d = prod_(i < d) (base - i + growth m)` in `m`.
pub(super) fn falling_factorial_on_class(base: usize, growth: usize, order: usize) -> Vec<Q> {
    let mut product = vec![Q::one()];
    for i in 0..order {
        product = vec_poly::mul(
            &product,
            &[integer(base as i64 - i as i64), integer(growth as i64)],
        );
    }
    product
}

pub(super) fn leading_coefficient(row: &[Q]) -> Option<&Q> {
    row.iter().rev().find(|c| !c.is_zero())
}

pub(super) fn row_degree(row: &[Q]) -> Option<usize> {
    row.iter().rposition(|c| !c.is_zero())
}

// ---------------------------------------------------------------------------
// Degree-law inference
// ---------------------------------------------------------------------------

/// Candidate laws read off the rows `rows[i] = P_(first + i)`: by increasing
/// period, the earliest start that fits every later sampled row and the next
/// few starts, with at least two periods and the largest lag of evidence.
pub(super) fn infer_degree_laws(
    rows: &[Vec<Q>],
    first: usize,
    maximum_period: usize,
    max_lag: usize,
    alternating: bool,
) -> Vec<RecurrenceDegreeLaw> {
    let degrees: Vec<Option<usize>> = rows.iter().map(|r| row_degree(r)).collect();
    let signs: Vec<Option<LeadingCoefficientSign>> = rows
        .iter()
        .map(|r| leading_coefficient(r).and_then(LeadingCoefficientSign::of))
        .collect();
    let earliest = degrees
        .iter()
        .rposition(Option::is_none)
        .map_or(0, |zero| zero + 1);
    let mut laws = Vec::new();
    for period in 1..=maximum_period {
        if alternating && period % 2 == 1 {
            continue;
        }
        let evidence = 2 * period + max_lag + 2;
        if rows.len() < earliest + evidence {
            continue;
        }
        let last = rows.len() - 1;
        let growth = degrees[last].unwrap() as i64 - degrees[last - period].unwrap() as i64;
        if growth < 0 {
            continue;
        }
        // The earliest start of the suffix on which the law fits.
        let mut start = last - period;
        while start > earliest {
            let i = start - 1;
            let fits = degrees[i + period].unwrap() as i64 - degrees[i].unwrap() as i64 == growth
                && signs[i + period] == signs[i];
            if !fits {
                break;
            }
            start = i;
        }
        if start + evidence > rows.len() {
            continue;
        }
        for start in start..=(start + EXTRA_LAW_STARTS).min(rows.len() - evidence) {
            laws.push(RecurrenceDegreeLaw {
                start: first + start,
                period,
                base_degrees: (0..period).map(|r| degrees[start + r].unwrap()).collect(),
                growth_per_period: growth as usize,
                leading_coefficient_signs: (0..period).map(|r| signs[start + r].unwrap()).collect(),
            });
        }
    }
    laws
}
