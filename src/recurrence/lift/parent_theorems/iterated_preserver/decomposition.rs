//! Step-operator construction and preserver decompositions.

use super::*;

// ---------------------------------------------------------------------------
// The step operator and its decompositions
// ---------------------------------------------------------------------------

/// The step `f(n) T_n = sum_d c_d(n, x) D^d` and the lag `L` of a
/// single-lag rule `P_n = T_n(P_(n-L))`, or why the rule does not have that
/// form.
pub(super) fn step_operator(
    rule: &Recurrence,
) -> Result<(StepCoefficients, Vec<Q>, usize), String> {
    let shape = recurrence_shape(rule)?;
    let Some(&(lag, _)) = shape.terms.keys().next() else {
        return Err("the rule has no nonzero term".into());
    };
    if lag == 0 || shape.terms.keys().any(|(offset, _)| *offset != lag) {
        return Err(
            "the rule is not first order: P_n = T_n(P_(n-1)), or T_n(P_(n-L)) with one lag L"
                .into(),
        );
    }
    let terms = shape
        .terms
        .into_iter()
        .map(|((_, order), coefficients)| (order, coefficients))
        .collect();
    Ok((terms, shape.denominator, lag))
}

/// `delta = max_d (deg_x c_d - d)`.
pub(super) fn degree_step(terms: &StepCoefficients) -> Option<i64> {
    terms
        .iter()
        .filter_map(|(d, c)| x_degree(c).map(|degree| degree as i64 - *d as i64))
        .max()
}

/// `sum_d c_d D^d` as an Ore operator.
pub(super) fn ore_operator(terms: &StepCoefficients) -> DOp {
    let mut op = DOp::default();
    for (d, c) in terms {
        op.add_term((0, *d), Biv::new(c.clone()));
    }
    op
}

/// A recognized step `T_n = g(x) h(n, x) S_n`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct Decomposition {
    /// `g`, nonzero, free of `n`.
    pub(super) fixed_factor: Vec<Q>,
    /// `h_j(n)`, ascending in `x`.
    pub(super) index_factor: Vec<IndexFraction>,
    /// `S_n`; `None` for the identity.
    pub(super) preserver: Option<PreserverOperator>,
    pub(super) factors: Vec<IteratedPreserverFactor>,
}

impl Decomposition {
    pub(super) fn from_certificate(
        certificate: &IteratedPreserverCertificate,
    ) -> Result<Self, RecurrenceError> {
        let fixed_factor = certificate
            .multiplier_fixed_factor
            .iter()
            .map(|value| {
                parse_rational_coeff(value)
                    .map_err(|error| RecurrenceError::from(error).with_context("multiplier"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let index_factor = certificate
            .multiplier_index_factor
            .iter()
            .map(IndexRationalFunction::parse)
            .collect::<Result<Vec<_>, _>>()?;
        if vec_poly::is_zero(&fixed_factor) || index_factor.is_empty() {
            return Err(RecurrenceError::InvalidInput(
                "the recorded multiplier is empty or zero".into(),
            ));
        }
        Ok(Self {
            fixed_factor,
            index_factor,
            preserver: certificate.preserver.clone(),
            factors: certificate.factors.clone(),
        })
    }

    /// `S_n = numerator / denominator(n)`.
    fn preserver_parts(&self) -> Option<(DOp, Vec<Q>)> {
        if !self.factors.is_empty() {
            if self.preserver.is_some() {
                return None;
            }
            return factorization::product_parts(&self.factors);
        }
        match &self.preserver {
            None => Some((DOp::identity(), vec![Q::one()])),
            Some(operator) => operator.weyl_parts(),
        }
    }

    /// `h = H / q` with a common denominator `q`.
    fn index_factor_over_common_denominator(&self) -> (Biv, Vec<Q>) {
        let common = self
            .index_factor
            .iter()
            .fold(vec![Q::one()], |acc, (_, q)| vec_poly::mul(&acc, q));
        let mut numerator = Biv::zero();
        for (j, (p, _)) in self.index_factor.iter().enumerate() {
            let others = self
                .index_factor
                .iter()
                .enumerate()
                .filter(|(i, _)| *i != j)
                .fold(p.clone(), |acc, (_, (_, q))| vec_poly::mul(&acc, q));
            numerator = numerator
                .add(&Biv::from_index_polynomial(&others).multiply(&Biv::monomial(0, j, Q::one())));
        }
        (numerator, common)
    }

    /// `q(n) delta_S(n) sum_d c_d D^d = f(n) (g H) o S_num` in the Ore
    /// algebra (the step `T_n = (1/f) sum c_d D^d` equals `g h S_n`).
    pub(super) fn identity_holds(&self, terms: &StepCoefficients, f: &[Q]) -> bool {
        let Some((numerator, denominator)) = self.preserver_parts() else {
            return false;
        };
        let (h, common) = self.index_factor_over_common_denominator();
        let multiplier = Biv::from_variable_polynomial(&self.fixed_factor).multiply(&h);
        if multiplier.is_zero() || vec_poly::is_zero(&common) || vec_poly::is_zero(&denominator) {
            return false;
        }
        let left = ore_operator(terms).left_multiply(&Biv::from_index_polynomial(&vec_poly::mul(
            &common,
            &denominator,
        )));
        let right = numerator
            .left_multiply(&multiplier)
            .left_multiply(&Biv::from_index_polynomial(f));
        left.subtract(&right).is_zero()
    }

    pub(super) fn index_factor_functions(&self, reflect: bool) -> Vec<IndexRationalFunction> {
        self.index_factor
            .iter()
            .enumerate()
            .map(|(j, (p, q))| {
                let p = if reflect && j % 2 == 1 {
                    vec_poly::neg(p)
                } else {
                    p.clone()
                };
                IndexRationalFunction::from_polynomials(&p, q)
            })
            .collect()
    }

    pub(super) fn multiplier_text(&self) -> String {
        let index_part: Vec<String> = self
            .index_factor_functions(false)
            .iter()
            .enumerate()
            .filter(|(_, c)| c.text() != "0")
            .map(|(j, c)| match j {
                0 => format!("({})", c.text()),
                1 => format!("({}) x", c.text()),
                _ => format!("({}) x^{j}", c.text()),
            })
            .collect();
        format!(
            "({}) ({})",
            format_index_polynomial(&self.fixed_factor, "x"),
            index_part.join(" + ")
        )
    }
}

/// The monic gcd over `Q[x]` of every coefficient of `n^i` in every `c_d`.
pub(super) fn fixed_content(terms: &StepCoefficients) -> Vec<Q> {
    let mut content: Option<Polynomial<Q>> = None;
    for c in terms.values() {
        for row in c {
            let row = vec_poly::trim(row);
            if vec_poly::is_zero(&row) {
                continue;
            }
            let row = Polynomial::new(row);
            content = Some(match content {
                None => row.gcd(&row),
                Some(g) => g.gcd(&row),
            });
        }
    }
    content.map_or_else(|| vec![Q::one()], |g| g.coeffs().to_vec())
}

/// `c / g` row by row (exact: `g` divides every row).
pub(super) fn divide_rows(c: &Bivariate, g: &[Q]) -> Bivariate {
    let divisor = Polynomial::new(g.to_vec());
    c.iter()
        .map(|row| {
            let row = vec_poly::trim(row);
            if vec_poly::is_zero(&row) {
                return Vec::new();
            }
            Polynomial::new(row).div_rem(&divisor).0.coeffs().to_vec()
        })
        .collect()
}

/// `x^k c`.
fn times_x_power(c: &Bivariate, k: usize) -> Bivariate {
    c.iter()
        .map(|row| {
            let row = vec_poly::trim(row);
            if vec_poly::is_zero(&row) {
                Vec::new()
            } else {
                vec_poly::shift_by(&row, k)
            }
        })
        .collect()
}

/// The coefficient of `y^j` of a [`Biv`] as a polynomial in `n`.
fn biv_column(b: &Biv, j: usize) -> Vec<Q> {
    let mut column: Vec<Q> = b
        .rows()
        .iter()
        .map(|row| row.get(j).cloned().unwrap_or_else(Q::zero))
        .collect();
    vec_poly::strip_trailing_zeros(&mut column);
    if column.is_empty() {
        vec![Q::zero()]
    } else {
        column
    }
}

fn over(p: Vec<Q>, f: &[Q]) -> IndexRationalFunction {
    let p = if p.is_empty() { vec![Q::zero()] } else { p };
    IndexRationalFunction::from_polynomials(&p, f)
}

/// The single power `m` of `x` in `c`, with its coefficient in `n`.
fn single_power(c: &Bivariate) -> Option<(usize, Vec<Q>)> {
    let powers: Vec<usize> = (0..x_width(c))
        .filter(|&j| !x_column(c, j).is_empty())
        .collect();
    match powers.as_slice() {
        [m] => Some((*m, x_column(c, *m))),
        _ => None,
    }
}

/// The menu shapes of `s = sum_d s_d D^d` over `f(n)`: pairs of the index
/// factor `h` and the operator (module documentation, "Recognition").
fn menu_shapes(s: &StepCoefficients, f: &[Q]) -> Vec<(Vec<IndexFraction>, PreserverOperator)> {
    let unit = || vec![(vec![Q::one()], vec![Q::one()])];
    let mut out = Vec::new();
    let max_order = s.keys().copied().max().unwrap_or(0);
    if max_order == 0 {
        return out;
    }
    // phi(n, D).
    if s.values().all(|c| x_degree(c) == Some(0)) {
        let coefficients = (0..=max_order)
            .map(|d| over(s.get(&d).map(|c| x_column(c, 0)).unwrap_or_default(), f))
            .collect();
        out.push((
            unit(),
            PreserverOperator::DerivativePolynomial { coefficients },
        ));
    }
    // f(n, theta) and the Euler step, from the Euler form x^J s.
    let (j, euler) = to_theta(&ore_operator(s));
    let levels: Vec<usize> = euler.keys().map(|(_, level)| *level).collect();
    let theta_coefficients = |level: usize, count: usize| -> Vec<IndexRationalFunction> {
        let b = euler.get(&(0, level)).cloned().unwrap_or_else(Biv::zero);
        (0..count).map(|t| over(biv_column(&b, t), f)).collect()
    };
    let theta_degree = |level: usize| {
        euler
            .get(&(0, level))
            .and_then(Biv::variable_degree)
            .unwrap_or(0)
    };
    if euler.keys().all(|(shift, _)| *shift == 0) {
        if levels.iter().all(|&level| level == j) {
            let count = theta_degree(j) + 1;
            out.push((
                unit(),
                PreserverOperator::EulerPolynomial {
                    coefficients: theta_coefficients(j, count),
                },
            ));
        } else if levels.iter().all(|&level| level == j || level == j + 1)
            && levels.contains(&(j + 1))
            && theta_degree(j) == 2
            && theta_degree(j + 1) <= 1
        {
            out.push((
                unit(),
                PreserverOperator::EulerStep {
                    diagonal: theta_coefficients(j, 3),
                    raising: theta_coefficients(j + 1, 2),
                },
            ));
        }
    }
    // kappa(n) x^m D^j and kappa(n) D^j x^m.
    if let Some((power, kappa)) = s.get(&max_order).and_then(single_power) {
        let factor = vec![(kappa, f.to_vec())];
        if s.len() == 1 {
            out.push((
                factor.clone(),
                PreserverOperator::PowerTimesDerivative {
                    power,
                    derivative_order: max_order,
                },
            ));
        }
        if power >= 1 && s.len() > 1 {
            out.push((
                factor,
                PreserverOperator::DerivativeOfPowerMultiple {
                    power,
                    derivative_order: max_order,
                },
            ));
        }
    }
    // Branden's Weyl step A + B D + C D^2 (last: the operators above have
    // their own theorems, and the Euler step is a half-line theorem).
    if max_order <= 2 {
        out.push((
            unit(),
            PreserverOperator::WeylOperator {
                terms: s
                    .iter()
                    .map(|(order, c)| OperatorTermJson {
                        deriv_order: *order,
                        coeff: biv_json(&Biv::new(c.clone())),
                    })
                    .collect(),
                denominator: vec_poly::trim(f)
                    .iter()
                    .map(format_rational_coeff)
                    .collect(),
            },
        ));
    }
    out
}

/// Every decomposition `T_n = g_k(x) h(n, x) S_n` whose operator identity
/// holds exactly, in a fixed order.
pub(super) fn decompositions(terms: &StepCoefficients, f: &[Q]) -> Vec<Decomposition> {
    let g = fixed_content(terms);
    let quotient: StepCoefficients = terms
        .iter()
        .map(|(d, c)| (*d, divide_rows(c, &g)))
        .collect();
    let mut candidates = Vec::new();
    if terms.keys().all(|&d| d == 0) {
        let c0 = quotient.get(&0).cloned().unwrap_or_default();
        let index_factor = (0..=x_degree(&c0).unwrap_or(0))
            .map(|j| {
                let column = x_column(&c0, j);
                (
                    if column.is_empty() {
                        vec![Q::zero()]
                    } else {
                        column
                    },
                    f.to_vec(),
                )
            })
            .collect();
        candidates.push(Decomposition {
            fixed_factor: g.clone(),
            index_factor,
            preserver: None,
            factors: Vec::new(),
        });
    } else {
        let valuation = g.iter().take_while(|c| c.is_zero()).count();
        for k in 0..=valuation {
            let fixed_factor = g[k..].to_vec();
            let shifted: StepCoefficients = quotient
                .iter()
                .map(|(d, c)| (*d, times_x_power(c, k)))
                .collect();
            for (index_factor, preserver) in menu_shapes(&shifted, f) {
                candidates.push(Decomposition {
                    fixed_factor: fixed_factor.clone(),
                    index_factor,
                    preserver: Some(preserver),
                    factors: Vec::new(),
                });
            }
        }
    }
    let mut out: Vec<Decomposition> = Vec::new();
    for candidate in candidates {
        if candidate.identity_holds(terms, f) && !out.contains(&candidate) {
            out.push(candidate);
        }
    }
    out
}
