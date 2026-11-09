//! Row-hypergeometric families: rows whose coefficient ratio
//! `[x^k] P_n / [x^(k-1)] P_n` is a rational function of `(n, k)` that
//! factors into linear factors in `k`, proved real-rooted for every `n` by a
//! product of multiplier sequences applied to `(1 + x)^M(n)`
//! ([`check_hypergeometric_row_theorem`]).
//!
//! # The row relation
//!
//! From the rows `P_n`, `n >= s` (`s` the first base where the rule holds
//! for the generated rows), a relation
//!
//! ```text
//! p(n, k) [x^k] P_n + q(n, k) [x^(k-1)] P_n = 0      for every integer k
//! ```
//!
//! with polynomials `p, q` of total degree at most 4 is fitted exactly
//! (common linear factors in `k` cancelled).  In operator form it is
//! `H(n) P_n = 0` with `H(n) = p(n, theta) + x q(n, theta + 1)`,
//! `theta = x D`.  The fit is only a candidate.  `H(n) P_n = 0` is proved
//! for every `n >= s` by an invariant row-annihilator system of the rule
//! (the certificate of the Hadamard-factor lifts, module documentation of
//! `hadamard_factor`, "Row annihilators"): exact identities
//! `nu_l(n, x) N(A^(l)) = sum_j B_(l,j)(n) o A^(j)` under the rule's
//! companion step, the generators checked on the rows at the start base and
//! after every restart base, and `mu(n, x) H(n) = sum_l Z_l(n) o A^(l)(n)`
//! with `mu` nonzero; at the integer bases where `mu(n, .)` vanishes,
//! `H(n) P_n = 0` is checked exactly.
//!
//! # From the relation to a multiplier sequence
//!
//! Write `p = c_p(n) prod_i (k + u_i(n))` and `q = c_q(n) prod_j (k + w_j(n))`
//! with `u_i, w_j` affine in `n` (other shapes are not applicable).  The
//! coefficients are read in one of two coordinates, `U_j = [x^(j + v(n))] P_n`
//! (the identity, `v(n)` the valuation of the rows, affine with integer
//! coefficients) or `U_j = [x^(d(n) - j)] P_n` (the reversal at the degree
//! `d(n)`, affine with integer coefficients).  In both,
//! `P(n, j) U_j + Q(n, j) U_(j-1) = 0` for every integer `j`, with
//! `R(n, j) = -Q / P = C(n) prod (j + a_l(n)) / prod (j + b_m(n))`
//! (identity: `a = w + v`, `b = u + v`, `C = -c_q / c_p`; reversal:
//! `a = -d - 1 - u`, `b = -d - 1 - w`, `C = -(-1)^(|p| + |q|) c_p / c_q`).
//! The support of the rows is proved as follows:
//!
//! - `[x^k] P_n = 0` for `k < v(n)`: by induction on the rule when every
//!   term `c x^e D^d P_(n-r)` has `e - d - v_1 r >= 0` (rows before the
//!   induction checked exactly), or, for a constant valuation `v`, by the
//!   relation with `p(n, k) != 0` for `0 <= k < v` (then each coefficient
//!   below `x^v` vanishes, from `[x^(-1)] P_n = 0` upwards);
//! - (reversal) `[x^k] P_n = 0` for `k > d(n)`: by induction on the rule
//!   when every term has `e - d - d_1 r <= 0`, or by the relation: `q` has
//!   the factor `k - d(n) - 1`, `c_p(n) != 0`, and every root `-u(n)` of
//!   `p` satisfies `u(n) + d(n) + 1 > 0` or lies above the crude degree
//!   bound `D(n) = D_0 + beta n` of the rule (`beta = max (e - d) / r`), so
//!   `p(n, k) != 0` for `d(n) < k <= D(n)`.
//!
//! The decided conditions start at a structure start at most
//! `MAX_STRUCTURE_DELAY` rows after the first base where the relation holds;
//! every row before the induction start is checked exactly.
//!
//! The factors are matched with multiplier sequences, every condition
//! decided for all `n` by [`crate::decide_index_rational_function_sign`]:
//!
//! - the **base**: a numerator factor `j - M(n) - 1` (`M` affine with integer
//!   coefficients, `M(n) >= 0`; for the reversal `M(n) <= d(n) - v`) and the
//!   denominator factor `j`, giving `binom(M, j)`;
//! - **complements**: further numerator factors `j - M'(n) - 1` with integer
//!   coefficients and `M' >= M`, giving `k -> M'! / (M' - k)!`, a multiplier
//!   sequence by the Pólya--Schur theorem (its exponential generating
//!   function is `(1 + x)^M'`);
//! - **pairs**: a numerator factor `j + a` with a denominator factor `j + b`,
//!   `a - b` affine with integer coefficients, `a - b >= 0` and `b > -1`,
//!   giving the polynomial `k -> (k + b + 1) ... (k + a)` with roots `<= 0`
//!   (Laguerre's theorem);
//! - **reciprocals**: the other denominator factors `j + b`, `b > -1`,
//!   giving `k -> 1 / (b + 1)_k`;
//! - the **scale** `C_eff(n) = C(n) (-1)^c > 0` (`c` complements, the base
//!   included), giving `k -> C_eff^k`.
//!
//! Every denominator factor is positive at `j >= 1`, so `P(n, j) != 0`
//! there and, for `j >= 0`, `U_j = U_0 prod_(t <= j) R(n, t)`, which is
//! `U_0 binom(M, j) mu_j` with `mu` the product of the sequences above (zero
//! for `j > M`, where the base factor vanishes).  So the polynomial
//! `sum_j U_j y^j` is `U_0` times the image of `(1 + y)^M` under the
//! nonnegative multiplier sequence `mu`, hence real-rooted with every root
//! `< 0` when `U_0 != 0`.  `U_0` is `[x^v] P_n` (identity) or `[x^d(n)] P_n`
//! (reversal); the rule gives a scalar recurrence
//! `f(n) U_0(n) = sum_r sigma_r(n) U_0(n - r)` (each coefficient of the rule
//! that reaches `U_0` is `U_0(n - r)` times a fixed number of ratio factors
//! `R(n - r, t)`), and `U_0(n) != 0` with one sign follows from
//! `sigma_r(n) / f(n) > 0` for every lag that occurs and nonzero seeds of one
//! sign (or from a hypergeometric solution of that recurrence); when neither
//! applies, the shared module
//! ([`crate::recurrence::certify_leading_coefficient`]) may prove it for the
//! scalar recurrence read as a family of constant rows (for example a
//! growth-ratio invariant for `U_0(n) = 3 U_0(n - 1) - 2 U_0(n - 2)`).  The
//! rows before the induction start are checked exactly.
//!
//! The conclusion: every `P_n` is nonzero and real-rooted with every root in
//! `(-inf, 0)`, or in `(-inf, 0]` when the rows are divisible by `x` (or the
//! reversal does not show `[x^v] P_n != 0`).  No interlacing is claimed.

use super::*;
use crate::factorization::irreducible_factorization_monic_bivariate_with_budget;
use crate::recurrence::implicit::algebra::{
    content_and_primitive, from_theta, theta_add, Biv, DOp, ThetaOp,
};
use crate::recurrence::lift::hadamard_factor::annihilator::vanishing_bases;
use crate::recurrence::lift::hadamard_factor::annihilator_system::find_invariant_annihilator_system;
use crate::recurrence::lift::hadamard_factor::certify::{
    check_system_rows, MAX_HADAMARD_OPERATOR_TERMS,
};
use crate::recurrence::lift::hadamard_factor::prepare::prepare_lift_rule;
use num_traits::ToPrimitive;
use ore_algebra::exact_kernel_vectors_by_group;

/// Theorem statement and references of [`check_hypergeometric_row_theorem`].
const REFERENCE: &str = "Row-hypergeometric rows: an exact row relation p(n, k) [x^k] P_n + \
     q(n, k) [x^(k-1)] P_n = 0, proved for every n by an invariant row-annihilator system of the \
     rule, writes every row as U_0(n) x^v (or its reversal) times the image of (1 + x)^M(n) under \
     a product of nonnegative multiplier sequences: C^k (C > 0); k -> M'!/(M' - k)!, whose \
     exponential generating function is (1 + x)^M' (Polya--Schur, Theorem 7.1 (3) => (1) of \
     P. Branden, Unimodality, log-concavity, real-rootedness and beyond, arXiv:1410.6601: a \
     sequence is a multiplier sequence when G(x) = sum_k lambda_k x^k / k! or G(-x) can be \
     written C x^m e^(a x) prod_k (1 + alpha_k x) with a >= 0, alpha_k >= 0); k -> phi(k) for a \
     polynomial phi with only real roots, all <= 0 (Laguerre's theorem, RealRooted \
     isPFMultiplierSequence_eval_of_roots_nonpos); k -> 1/(alpha)_k for alpha > 0 (RealRooted \
     isPFMultiplierSequence_inv_ascPochhammer, via the generalized Laguerre polynomials); \
     products of PF multiplier sequences are PF multiplier sequences (RealRooted \
     IsPFMultiplierSequence.mul) and map real-rooted polynomials with nonnegative coefficients \
     to such polynomials (IsPFMultiplierSequence.diagonalOperator_isPF)";

/// Largest total degree in `(n, k)` of a fitted row relation.
const MAX_RELATION_DEGREE: usize = 4;

/// Rows (from the first base where the rule holds) used to fit the
/// relation.
const FIT_ROWS: usize = 14;

/// Largest number of rows by which the start of the decided conditions may
/// move past the first base (the rows before it are checked exactly).
const MAX_STRUCTURE_DELAY: usize = 4;

/// Largest degree of the numerator and denominator of `rho(n)` in
/// `U_0(n) = rho(n) U_0(n - 1)`.
const MAX_FIRST_COEFFICIENT_DEGREE: usize = 3;

/// Recombination budget of the bivariate factorizations.
const FACTORIZATION_BUDGET: usize = 256;

/// `a_0 + a_1 n`.
type Affine = (Q, Q);

fn affine_function(a: &Affine) -> IndexFunction {
    IndexFunction::polynomial(vec![a.0.clone(), a.1.clone()])
}

fn affine_add(a: &Affine, b: &Affine) -> Affine {
    (&a.0 + &b.0, &a.1 + &b.1)
}

fn affine_neg(a: &Affine) -> Affine {
    (-a.0.clone(), -a.1.clone())
}

fn affine_constant(c: i64) -> Affine {
    (Q::from_integer(c.into()), Q::zero())
}

fn affine_is_integral(a: &Affine) -> bool {
    a.0.is_integer() && a.1.is_integer()
}

fn affine_text(a: &Affine) -> String {
    IndexFunction::polynomial(vec![a.0.clone(), a.1.clone()]).text()
}

/// `j + a` as text.
fn factor_text(a: &Affine) -> String {
    let text = affine_text(a);
    if text == "0" {
        "j".into()
    } else if let Some(rest) = text.strip_prefix('-') {
        format!("j - {rest}")
    } else {
        format!("j + {text}")
    }
}

fn times(a: &IndexFunction, b: &IndexFunction) -> IndexFunction {
    IndexFunction::new(
        vec_poly::mul(&a.numerator, &b.numerator),
        vec_poly::mul(&a.denominator, &b.denominator),
    )
}

fn holds(function: &IndexFunction, first: usize, condition: IndexSignCondition) -> bool {
    matches!(
        function.decide(first, condition),
        IndexSignDecision::Holds { .. }
    )
}

/// `p(n, k) = scalar(n) prod_i (k + roots[i](n))`.
#[derive(Debug, Clone, PartialEq, Eq)]
struct LinearFactors {
    scalar: Vec<Q>,
    roots: Vec<Affine>,
}

impl LinearFactors {
    fn biv(&self) -> Biv {
        let mut out = Biv::from_index_polynomial(&self.scalar);
        for root in &self.roots {
            let factor = Biv::new(vec![vec![root.0.clone(), Q::one()], vec![root.1.clone()]]);
            out = out.multiply(&factor);
        }
        out
    }

    /// `p(n, k0)` as a polynomial in `n`.
    fn at(&self, k0: &Affine) -> Vec<Q> {
        let mut out = self.scalar.clone();
        for root in &self.roots {
            let value = affine_add(root, k0);
            out = vec_poly::mul(&out, &[value.0, value.1]);
        }
        vec_poly::trim(&out)
    }

    fn text(&self, variable: &str) -> String {
        let mut parts = vec![format!(
            "({})",
            IndexFunction::polynomial(self.scalar.clone()).text()
        )];
        for root in &self.roots {
            parts.push(format!(
                "({})",
                factor_text(root).replacen('j', variable, 1)
            ));
        }
        parts.join(" ")
    }
}

/// The factorization of `p(n, k)` (`Biv` in `(n, k)`) into a scalar in `n`
/// and linear factors `k + u(n)`, `u` affine; `Err` naming the shape
/// otherwise.
fn linear_factors(p: &Biv) -> Result<LinearFactors, String> {
    let Some(degree) = p.variable_degree() else {
        return Err("the relation has a zero coefficient".into());
    };
    if degree == 0 {
        return Ok(LinearFactors {
            scalar: vec_poly::trim(&p.variable_coefficient(0)),
            roots: Vec::new(),
        });
    }
    let (content, primitive) = content_and_primitive(p);
    let leading = vec_poly::trim(&primitive.variable_coefficient(degree));
    if leading.len() != 1 {
        return Err(
            "the leading coefficient in k of the relation depends on n after removing \
                    the content"
                .into(),
        );
    }
    let lead = leading[0].clone();
    let coefficients: Vec<Vec<Q>> = (0..=degree)
        .map(|i| {
            vec_poly::trim(&vec_poly::scale(
                &primitive.variable_coefficient(i),
                &lead.recip(),
            ))
        })
        .collect();
    let factors =
        irreducible_factorization_monic_bivariate_with_budget(&coefficients, FACTORIZATION_BUDGET)
            .map_err(|error| format!("the relation was not factored: {error}"))?;
    let mut roots = Vec::new();
    for factor in factors {
        if factor.coefficients.len() != 2 || factor.coefficients[0].len() > 2 {
            return Err("a factor of the relation is not linear in k and affine in n".into());
        }
        let root = (
            factor.coefficients[0]
                .first()
                .cloned()
                .unwrap_or_else(Q::zero),
            factor.coefficients[0]
                .get(1)
                .cloned()
                .unwrap_or_else(Q::zero),
        );
        for _ in 0..factor.multiplicity {
            roots.push(root.clone());
        }
    }
    roots.sort();
    Ok(LinearFactors {
        scalar: vec_poly::trim(&vec_poly::scale(&content, &lead)),
        roots,
    })
}

/// The relation `p(n, k) T_k + q(n, k) T_(k-1) = 0` of smallest total degree
/// (at most [`MAX_RELATION_DEGREE`]) satisfied by every coefficient of the
/// rows `rows[i] = P_(first + i)` with index `>= from`; `None` when there is
/// none.  Checked exactly on every row used.
fn fit_row_relation(rows: &[Vec<Q>], first: usize, from: usize) -> Option<(Biv, Biv)> {
    let used: Vec<(usize, &Vec<Q>)> = rows
        .iter()
        .enumerate()
        .map(|(i, row)| (first + i, row))
        .filter(|(n, _)| *n >= from)
        .collect();
    let coefficient = |row: &Vec<Q>, k: i64| -> Q {
        if k < 0 {
            Q::zero()
        } else {
            row.get(k as usize).cloned().unwrap_or_else(Q::zero)
        }
    };
    for degree in 1..=MAX_RELATION_DEGREE {
        let monomials: Vec<(usize, usize)> = (0..=degree)
            .flat_map(|a| (0..=degree - a).map(move |b| (a, b)))
            .collect();
        let width = 2 * monomials.len();
        let mut equations: Vec<Vec<Q>> = Vec::new();
        for (n, row) in &used {
            let n_value = Q::from_integer(BigInt::from(*n));
            for k in 0..=row.len() as i64 {
                let k_value = Q::from_integer(BigInt::from(k));
                let (current, previous) = (coefficient(row, k), coefficient(row, k - 1));
                if current.is_zero() && previous.is_zero() {
                    continue;
                }
                let mut equation = Vec::with_capacity(width);
                for value in [&current, &previous] {
                    for (a, b) in &monomials {
                        equation.push(value * n_value.pow(*a as i32) * k_value.pow(*b as i32));
                    }
                }
                equations.push(equation);
            }
        }
        if equations.len() < width + 8 {
            return None;
        }
        let integral: Vec<Vec<BigInt>> = equations
            .iter()
            .map(|equation| {
                let denominator = equation
                    .iter()
                    .fold(BigInt::one(), |acc, value| acc.lcm(value.denom()));
                equation
                    .iter()
                    .map(|value| value.numer() * (&denominator / value.denom()))
                    .collect()
            })
            .collect();
        let Some(vector) = exact_kernel_vectors_by_group(&integral, &vec![0; width], 1)
            .into_iter()
            .next()
            .flatten()
        else {
            continue;
        };
        let solves = equations.iter().all(|equation| {
            equation
                .iter()
                .zip(&vector)
                .fold(Q::zero(), |acc, (a, b)| acc + a * b)
                .is_zero()
        });
        if !solves {
            continue;
        }
        let build = |offset: usize| {
            let mut rows = vec![vec![Q::zero(); degree + 1]; degree + 1];
            for (index, (a, b)) in monomials.iter().enumerate() {
                rows[*a][*b] = vector[offset + index].clone();
            }
            Biv::new(rows)
        };
        let (p, q) = (build(0), build(monomials.len()));
        if p.is_zero() || q.is_zero() {
            continue;
        }
        return Some((p, q));
    }
    None
}

/// `H(n) = p(n, theta) + x q(n, theta + 1)`: `[x^k] H P = p(n, k) [x^k] P +
/// q(n, k) [x^(k-1)] P`.
fn relation_operator(p: &Biv, q: &Biv) -> DOp {
    let mut op = ThetaOp::new();
    theta_add(&mut op, (0, 0), p.clone());
    theta_add(&mut op, (0, 1), q.shift(&Q::zero(), &Q::one()));
    from_theta(&op)
}

/// The coordinates of [`check_hypergeometric_row_theorem`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Coordinates {
    Identity,
    Reversal,
}

/// The support facts (module documentation) decided from `first`: the
/// valuation bound when the valuation is not zero, and for the reversal
/// the degree bound.
fn support_hypotheses(
    shape: &RecurrenceShape,
    p: &LinearFactors,
    q: &LinearFactors,
    valuation: &Affine,
    degree: Option<&Affine>,
    (rows, first_index, rule_start): (&[Vec<Q>], usize, usize),
    first: usize,
) -> Vec<LiftHypothesis> {
    // Rows checked exactly before the bounds hold by induction on the rule.
    let max_lag = shape.terms.keys().map(|(lag, _)| *lag).max().unwrap_or(1);
    let earlier: Vec<(usize, &Vec<Q>)> = (first_index.max(first.saturating_sub(max_lag))
        ..rule_start.max(first))
        .filter_map(|n| rows.get(n - first_index).map(|row| (n, row)))
        .collect();
    let mut out = Vec::new();
    if *valuation != affine_constant(0) {
        out.push(lower_support_hypothesis(
            shape, p, valuation, &earlier, first,
        ));
    }
    if let Some(degree) = degree {
        out.extend(upper_support_hypotheses(
            shape, p, q, degree, &earlier, first,
        ));
    }
    out
}

/// A matching of the ratio factors with multiplier sequences (module
/// documentation): `numerators[i]` is `j + num[i]`, `denominators[m]` is
/// `j + den[m]`.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Matching {
    base_numerator: usize,
    complements: Vec<usize>,
    pairs: Vec<(usize, usize)>,
    reciprocals: Vec<usize>,
}

/// The ratio `R(n, j) = C(n) prod (j + num) / prod (j + den)` in the chosen
/// coordinates.
#[derive(Clone)]
struct Ratio {
    scale: IndexFunction,
    numerators: Vec<Affine>,
    denominators: Vec<Affine>,
}

impl Ratio {
    /// `R(n, t)` for a fixed integer `t >= 1`.
    fn at(&self, t: usize) -> IndexFunction {
        let t = affine_constant(t as i64);
        let mut out = self.scale.clone();
        for a in &self.numerators {
            out = times(&out, &affine_function(&affine_add(a, &t)));
        }
        for b in &self.denominators {
            out = out.div(&affine_function(&affine_add(b, &t)));
        }
        out
    }

    fn text(&self) -> String {
        let product = |factors: &[Affine]| {
            if factors.is_empty() {
                "1".to_string()
            } else {
                factors
                    .iter()
                    .map(|a| format!("({})", factor_text(a)))
                    .collect::<Vec<_>>()
                    .join(" ")
            }
        };
        format!(
            "({}) {} / ({})",
            self.scale.text(),
            product(&self.numerators),
            product(&self.denominators)
        )
    }
}

/// `M` of a numerator factor `j - M - 1` with integer coefficients.
fn complement_size(a: &Affine) -> Option<Affine> {
    let size = affine_add(&affine_neg(a), &affine_constant(-1));
    affine_is_integral(&size).then_some(size)
}

/// The first matching (in a fixed search order) whose conditions all hold
/// for `n >= first`; `upper` is `d(n) - v` for the reversal.
fn find_matching(ratio: &Ratio, first: usize, upper: Option<&Affine>) -> Option<Matching> {
    let zero = affine_constant(0);
    let one = affine_constant(1);
    let nonnegative =
        |a: &Affine| holds(&affine_function(a), first, IndexSignCondition::Nonnegative);
    let positive = |a: &Affine| holds(&affine_function(a), first, IndexSignCondition::Positive);
    for (base_numerator, a) in ratio.numerators.iter().enumerate() {
        let Some(size) = complement_size(a) else {
            continue;
        };
        if !nonnegative(&size) {
            continue;
        }
        if let Some(upper) = upper {
            if !nonnegative(&affine_add(upper, &affine_neg(&size))) {
                continue;
            }
        }
        let Some(base_denominator) = ratio.denominators.iter().position(|b| *b == zero) else {
            continue;
        };
        let mut used_denominators = vec![false; ratio.denominators.len()];
        used_denominators[base_denominator] = true;
        let mut complements = Vec::new();
        let mut pairs = Vec::new();
        let mut ok = true;
        for (index, a) in ratio.numerators.iter().enumerate() {
            if index == base_numerator {
                continue;
            }
            if let Some(other) = complement_size(a) {
                if nonnegative(&affine_add(&other, &affine_neg(&size))) {
                    complements.push(index);
                    continue;
                }
            }
            let partner = ratio.denominators.iter().enumerate().position(|(m, b)| {
                !used_denominators[m] && {
                    let difference = affine_add(a, &affine_neg(b));
                    affine_is_integral(&difference)
                        && nonnegative(&difference)
                        && positive(&affine_add(b, &one))
                }
            });
            match partner {
                Some(m) => {
                    used_denominators[m] = true;
                    pairs.push((index, m));
                }
                None => {
                    ok = false;
                    break;
                }
            }
        }
        if !ok {
            continue;
        }
        let reciprocals: Vec<usize> = (0..ratio.denominators.len())
            .filter(|m| !used_denominators[*m])
            .collect();
        if reciprocals
            .iter()
            .any(|m| !positive(&affine_add(&ratio.denominators[*m], &one)))
        {
            continue;
        }
        let sign = if (complements.len() + 1) % 2 == 0 {
            Q::one()
        } else {
            -Q::one()
        };
        if !holds(
            &ratio.scale.scale(&sign),
            first,
            IndexSignCondition::Positive,
        ) {
            continue;
        }
        return Some(Matching {
            base_numerator,
            complements,
            pairs,
            reciprocals,
        });
    }
    None
}

/// Check the row-hypergeometric theorem (module documentation) for the
/// family of `rule`, `first_index` and `seeds`.
///
/// Fits an exact row relation `p(n, k) [x^k] P_n + q(n, k) [x^(k-1)] P_n = 0`
/// (total degree at most 4), proves it for every `n` by an invariant
/// row-annihilator system of the rule, and decides, for all `n` with the
/// exact sign checker, that the ratio `-q/p` is the ratio of `U_0(n)` times
/// the image of `(1 + x)^M(n)` under a product of nonnegative multiplier
/// sequences (Pólya--Schur, Laguerre, `1/(alpha)_k`), in the identity
/// coordinates or after the reversal at the degree.  Covers, for example,
/// `sum_k binom(n, k) (2n + 1)_k / (n + 2)_k x^k` (OEIS A330798) and
/// `sum_k binom(n, k) (k + 1)(k + 2)(k + 3) x^k` (A124848).  A failed
/// hypothesis says nothing against real-rootedness.
///
/// # Errors
///
/// Errors parsing or normalizing the recurrence and seeds.
pub fn check_hypergeometric_row_theorem(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<Q>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    check_real_rootedness_theorem(
        RealRootednessTheorem::HypergeometricRow,
        rule,
        first_index,
        seeds,
        options,
    )
}

/// The registry entry of the row-hypergeometric theorem ([`super::registry::ROUTES`]).
pub(super) const ROUTE: TheoremRoute = TheoremRoute {
    theorem: RealRootednessTheorem::HypergeometricRow,
    check_at_supplied_seeds: hypergeometric_row_report,
    start_policy: StartPolicy::MovesStart,
    is_preferred: any_proof,
    certificate_name: None,
    discovery: DiscoveryPolicy::Always,
    in_residue_class_menu: true,
    in_index_section_menu: true,
    parent: ParentPolicy {
        stage: ParentStage::ReplacesEarly,
        region: ParentRegion::Unrestricted,
        record_only_if_applicable: true,
    },
};

/// The theorem at the supplied seeds.
fn hypergeometric_row_report(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<Q>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    let mut report = RealRootednessTheoremReport::new(
        RealRootednessTheorem::HypergeometricRow,
        REFERENCE,
        rule,
        first_index,
        seeds,
        options,
    );
    let shape = match recurrence_shape(rule) {
        Ok(shape) => shape,
        Err(reason) => return Ok(report.not_applicable(reason)),
    };
    if seeds.is_empty() || shape.terms.is_empty() {
        return Ok(report.not_applicable("the rule needs seed rows and a lag term"));
    }
    let total = seeds.len() + FIT_ROWS;
    let rows: Vec<Vec<Q>> = match rule.generate_rows_rational(seeds, first_index, total) {
        Ok(rows) => rows.iter().map(|row| vec_poly::trim(row)).collect(),
        Err(error) => {
            return Ok(report.not_applicable(format!("the rule cannot generate rows: {error}")))
        }
    };
    // The prepared rule; when the supplied seeds do not all follow the rule
    // (more seeds than its order), the family is prepared from its last
    // identifying seeds on, and the earlier rows are checked exactly below.
    let prepared = match prepare_lift_rule(rule, &rows, first_index, "input") {
        Ok(prepared) => prepared,
        Err(error) => {
            let skip = seeds.len().saturating_sub(rule.max_offset());
            match (skip > 0)
                .then(|| prepare_lift_rule(rule, &rows[skip..], first_index + skip, "input"))
            {
                Some(Ok(prepared)) => prepared,
                _ => return Ok(report.not_applicable(error.to_string())),
            }
        }
    };
    let rule_start = first_index + seeds.len();
    let base = prepared.valid_from_base();
    // The first start (up to MAX_STRUCTURE_DELAY rows after the first base)
    // with a fitted relation of hypergeometric shape whose reduced form
    // kills every generated row from there on.
    let mut shape_failure = None;
    let attempt = |from: usize| -> Option<(LinearFactors, LinearFactors, DOp)> {
        let (p, q) = fit_row_relation(&rows, first_index, from)?;
        let (mut p, mut q) = (linear_factors(&p).ok()?, linear_factors(&q).ok()?);
        // Cancel the linear factors common to p and q.
        let mut kept = Vec::new();
        for root in std::mem::take(&mut p.roots) {
            if let Some(position) = q.roots.iter().position(|r| *r == root) {
                q.roots.remove(position);
            } else {
                kept.push(root);
            }
        }
        p.roots = kept;
        let operator = relation_operator(&p.biv(), &q.biv());
        rows.iter()
            .enumerate()
            .map(|(i, row)| (first_index + i, row))
            .filter(|(n, _)| *n >= from)
            .all(|(n, row)| operator.apply_row(row, n).iter().all(Zero::is_zero))
            .then_some((p, q, operator))
    };
    let mut found = None;
    for from in base..=base + MAX_STRUCTURE_DELAY {
        if let Some(result) = attempt(from) {
            found = Some((from, result));
            break;
        }
        if shape_failure.is_none() {
            if let Some((p, q)) = fit_row_relation(&rows, first_index, from) {
                if let Err(reason) = linear_factors(&p).and_then(|_| linear_factors(&q)) {
                    shape_failure = Some(reason);
                }
            }
        }
    }
    let Some((from, (p, q, operator))) = found else {
        return Ok(report.not_applicable(match shape_failure {
            Some(reason) => format!("the row relation is not of hypergeometric shape: {reason}"),
            None => format!(
                "no row relation p(n, k) [x^k] P_n + q(n, k) [x^(k-1)] P_n = 0 of total degree \
                 at most {MAX_RELATION_DEGREE} with linear factors in k holds on the rows from \
                 n = {base} + delay on, delay <= {MAX_STRUCTURE_DELAY}"
            ),
        }));
    };
    // Valuation and degree of the rows, affine in n from some row on.
    let fitted: Vec<(usize, &Vec<Q>)> = rows
        .iter()
        .enumerate()
        .map(|(i, row)| (first_index + i, row))
        .filter(|(n, _)| *n >= from)
        .collect();
    if fitted.iter().any(|(_, row)| vec_poly::is_zero(row)) {
        return Ok(report.not_applicable("a row after the first base is zero"));
    }
    let as_q = |value: usize| Q::from_integer(BigInt::from(value));
    let valuations: Vec<(usize, Q)> = fitted
        .iter()
        .map(|(n, row)| (*n, as_q(row.iter().position(|c| !c.is_zero()).unwrap_or(0))))
        .collect();
    let degrees: Vec<(usize, Q)> = fitted
        .iter()
        .map(|(n, row)| (*n, as_q(row.len() - 1)))
        .collect();
    let Some((shape_from, valuation, degree)) = (0..=MAX_STRUCTURE_DELAY).find_map(|delay| {
        let valuation = fit_affine(&valuations[delay.min(valuations.len())..])?;
        Some((
            from + delay,
            valuation,
            fit_affine(&degrees[delay.min(degrees.len())..]),
        ))
    }) else {
        return Ok(report.not_applicable(
            "the valuation of the rows is not an affine function of n with integer coefficients",
        ));
    };
    report.normal_form = format!(
        "row relation p(n, k) [x^k] P_n + q(n, k) [x^(k-1)] P_n = 0, n >= {from}, with p = {}, \
         q = {}; valuation v(n) = {}{}",
        p.text("k"),
        q.text("k"),
        affine_text(&valuation),
        degree
            .as_ref()
            .map(|d| format!(", degree d(n) = {}", affine_text(d)))
            .unwrap_or_default()
    );
    report.lean_declarations = vec![
        "RealRooted.isPFMultiplierSequence_eval_of_roots_nonpos".into(),
        "RealRooted.isPFMultiplierSequence_inv_ascPochhammer".into(),
        "RealRooted.IsPFMultiplierSequence.mul".into(),
        "RealRooted.IsPFMultiplierSequence.diagonalOperator_isPF".into(),
    ];

    // The relation for every n >= from.
    report.hypotheses.push(relation_hypothesis(
        &operator,
        &prepared,
        from,
        &rows,
        first_index,
    ));

    let mut chosen = None;
    let mut fallback = None;
    'search: for coordinates in [Coordinates::Identity, Coordinates::Reversal] {
        let ratio = match (coordinates, &degree) {
            (Coordinates::Identity, _) => Ratio {
                scale: IndexFunction::new(vec_poly::neg(&q.scalar), p.scalar.clone()),
                numerators: q.roots.iter().map(|w| affine_add(w, &valuation)).collect(),
                denominators: p.roots.iter().map(|u| affine_add(u, &valuation)).collect(),
            },
            (Coordinates::Reversal, Some(degree)) => {
                let shifted = |r: &Affine| {
                    affine_neg(&affine_add(&affine_add(degree, &affine_constant(1)), r))
                };
                let sign = if (p.roots.len() + q.roots.len()) % 2 == 0 {
                    -Q::one()
                } else {
                    Q::one()
                };
                Ratio {
                    scale: IndexFunction::new(vec_poly::scale(&p.scalar, &sign), q.scalar.clone()),
                    numerators: p.roots.iter().map(shifted).collect(),
                    denominators: q.roots.iter().map(shifted).collect(),
                }
            }
            (Coordinates::Reversal, None) => continue,
        };
        let upper = (coordinates == Coordinates::Reversal)
            .then(|| {
                degree
                    .as_ref()
                    .map(|d| affine_add(d, &affine_neg(&valuation)))
            })
            .flatten();
        for structure_start in shape_from.max(from)..=from + MAX_STRUCTURE_DELAY {
            let Some(matching) = find_matching(&ratio, structure_start, upper.as_ref()) else {
                continue;
            };
            let support = support_hypotheses(
                &shape,
                &p,
                &q,
                &valuation,
                (coordinates == Coordinates::Reversal)
                    .then_some(degree.as_ref())
                    .flatten(),
                (&rows, first_index, rule_start),
                structure_start,
            );
            let verified = support
                .iter()
                .all(|h| h.status == LiftHypothesisStatus::Verified);
            if verified || fallback.is_none() {
                fallback = Some((
                    coordinates,
                    ratio.clone(),
                    matching,
                    upper.clone(),
                    structure_start,
                    support,
                ));
            }
            if verified {
                chosen = fallback.take();
                break 'search;
            }
        }
    }
    let Some((coordinates, ratio, matching, upper, s, support)) = chosen.or(fallback) else {
        report.hypotheses.push(unresolved_hypothesis(
            "multiplier_matching",
            "the ratio factors match a product of multiplier sequences applied to (1 + x)^M",
            format!(
                "no matching of the ratio factors (identity: R(n, j) = -q(n, j + v(n)) / p(n, j \
                 + v(n))) with a base binom(M, j), complements, pairs and reciprocal rising \
                 factorials holds for every n >= {from} + delay, delay <= \
                 {MAX_STRUCTURE_DELAY}{}",
                if degree.is_some() {
                    ", also after the reversal at the degree"
                } else {
                    ""
                }
            ),
        ));
        return Ok(report.finish(RealRootednessTheoremConclusion {
            first_index,
            zero_residue_classes: None,
            interlacing_from: None,
            interlacing: None,
            simple_roots: false,
            root_region: None,
            statement: String::new(),
        }));
    };

    report.hypotheses.extend(support);
    report
        .hypotheses
        .extend(matching_hypotheses(&ratio, &matching, s, upper.as_ref()));

    // The first coefficient U_0(n): its scalar recurrence
    // f(n) U_0(n) = sum_r sigma_r(n) U_0(n - r) from the rule.
    let anchor: Affine = match (coordinates, &degree) {
        (Coordinates::Reversal, Some(degree)) => degree.clone(),
        _ => valuation.clone(),
    };
    let slope = match coordinates {
        Coordinates::Identity => valuation.1.clone(),
        Coordinates::Reversal => -anchor.1.clone(),
    };
    let mut lag_coefficients: BTreeMap<usize, IndexFunction> = BTreeMap::new();
    for (&(lag, order), coefficient) in &shape.terms {
        for power in 0..x_width(coefficient) {
            let c = x_column(coefficient, power);
            if c.is_empty() {
                continue;
            }
            // [x^kappa] (x^power D^order P_(n - lag)) is the falling factorial
            // of kappa - power + order of length order times
            // [x^(kappa - power + order)] P_(n - lag), whose index in the
            // coordinates of row n - lag is j below.
            let offset = order as i64 - power as i64;
            let source = affine_add(&anchor, &affine_constant(offset));
            let sign = if coordinates == Coordinates::Identity {
                1
            } else {
                -1
            };
            let Some(j) = (&slope * as_q(lag))
                .to_integer()
                .to_i64()
                .map(|shift| shift + sign * offset)
                .filter(|j| *j >= 0)
            else {
                continue;
            };
            let mut product = IndexFunction::polynomial(c);
            for t in 0..order {
                product = times(
                    &product,
                    &affine_function(&affine_add(&source, &affine_constant(-(t as i64)))),
                );
            }
            for t in 1..=j as usize {
                product = times(&product, &ratio.at(t).shifted(-(lag as i64)));
            }
            let entry = lag_coefficients
                .entry(lag)
                .or_insert_with(|| IndexFunction::constant(Q::zero()));
            *entry = entry.add(&product);
        }
    }
    lag_coefficients.retain(|_, f| !vec_poly::is_zero(&f.numerator));
    // Every lag of the rule: a skipped term (zero by the support) also needs
    // the support facts at its row.
    let max_lag = shape.terms.keys().map(|(lag, _)| *lag).max().unwrap_or(1);
    let start = rule_start.max(s + max_lag);
    let values: Vec<(usize, Q)> = (start - max_lag..first_index + rows.len())
        .map(|n| {
            let row = &rows[n - first_index];
            let index = &anchor.0 + &anchor.1 * as_q(n);
            let value = index
                .to_integer()
                .to_usize()
                .and_then(|k| row.get(k).cloned())
                .unwrap_or_else(Q::zero);
            (n, value)
        })
        .collect();
    let mut first_coefficient = first_coefficient_hypothesis(
        &lag_coefficients,
        &IndexFunction::polynomial(shape.denominator.clone()),
        &values,
        start,
        max_lag,
    );
    // Otherwise the shared leading-coefficient module, applied to U_0 as a
    // family of constants: its scalar recurrence with seeds U_0(n),
    // start - max_lag <= n < start.
    if first_coefficient.status != LiftHypothesisStatus::Verified {
        if let Some(hypothesis) = shared_first_coefficient_hypothesis(
            &lag_coefficients,
            &shape.denominator,
            &values,
            start,
            max_lag,
        ) {
            first_coefficient = hypothesis;
        }
    }
    report.hypotheses.push(first_coefficient);

    // The region and the rows before the induction start.
    let full_support = coordinates == Coordinates::Identity || {
        let size = complement_size(&ratio.numerators[matching.base_numerator]);
        size.as_ref() == upper.as_ref()
    };
    let closed = valuation != affine_constant(0) || !full_support;
    let region = RealInterval {
        lower: None,
        lower_closed: false,
        upper: Some(Q::zero()),
        upper_closed: closed,
    };
    let in_region = |row: &[Q]| {
        let row = vec_poly::trim(row);
        !vec_poly::is_zero(&row)
            && is_real_rooted_rational(&row)
            && if closed {
                row.iter().all(|c| !c.is_negative()) || row.iter().all(|c| !c.is_positive())
            } else {
                has_real_roots_below(&row, &Q::zero())
            }
    };
    let failing = (first_index..start).find(|n| !in_region(&rows[n - first_index]));
    report.hypotheses.push(exact_hypothesis(
        "earlier_rows",
        &format!(
            "the rows P_{first_index}, ..., P_{} before the induction start are nonzero and \
             real-rooted with every root in {}",
            start - 1,
            RootInterval::from_real(&region)
        ),
        match failing {
            None => Ok("checked exactly".into()),
            Some(n) => Err((Some(n), format!("P_{n} fails the exact check"))),
        },
    ));
    if let Some(requested) = &options.root_region {
        let contains = requested.lower.is_none()
            && requested
                .upper
                .as_ref()
                .is_none_or(|upper| upper.is_positive() || (upper.is_zero() && !closed));
        report.hypotheses.push(exact_hypothesis(
            "requested_region",
            "the requested root region contains the certified region",
            if contains {
                Ok(format!("{} is contained", RootInterval::from_real(&region)))
            } else {
                Err((
                    None,
                    format!(
                        "{} is not contained in the requested region",
                        RootInterval::from_real(&region)
                    ),
                ))
            },
        ));
    }
    let size = complement_size(&ratio.numerators[matching.base_numerator])
        .map(|size| affine_text(&size))
        .unwrap_or_default();
    let statement = format!(
        "every P_n, n >= {first_index}, is nonzero and real-rooted with every root in {}: for \
         n >= {s} the coefficients U_j ({}) satisfy U_j = U_0(n) prod_(t <= j) R(n, t), \
         R(n, j) = {}, so sum_j U_j y^j is U_0(n) times the image of (1 + y)^({size}) under a \
         nonnegative multiplier sequence, and U_0(n) != 0 for n >= {} by its scalar recurrence",
        RootInterval::from_real(&region),
        match coordinates {
            Coordinates::Identity => format!(
                "U_j = [x^(j + v(n))] P_n, v(n) = {}",
                affine_text(&valuation)
            ),
            Coordinates::Reversal => format!("U_j = [x^(d(n) - j)] P_n, d(n) = {}", {
                affine_text(&anchor)
            }),
        },
        ratio.text(),
        start - max_lag
    );
    Ok(report.finish(RealRootednessTheoremConclusion {
        first_index,
        zero_residue_classes: None,
        interlacing_from: None,
        interlacing: None,
        simple_roots: false,
        root_region: Some(RootInterval::from_real(&region)),
        statement,
    }))
}

/// `a_0 + a_1 n` with integer coefficients through every point, from the
/// first two; `None` when there are fewer than two points or no such line.
fn fit_affine(points: &[(usize, Q)]) -> Option<Affine> {
    let (n0, y0) = points.first()?;
    let (n1, y1) = points.get(1)?;
    let slope = (y1 - y0) / Q::from_integer(BigInt::from(n1 - n0));
    let intercept = y0 - &slope * Q::from_integer(BigInt::from(*n0));
    let line = (intercept, slope);
    (affine_is_integral(&line)
        && points
            .iter()
            .all(|(n, y)| *y == &line.0 + &line.1 * Q::from_integer(BigInt::from(*n))))
    .then_some(line)
}

/// The valuation bound `val P_n >= v(n)` (so `[x^k] P_n = 0` for
/// `k < v(n)`): by induction on the rule when every term `c x^e D^d P_(n-r)`
/// has `e - d - v_1 r >= 0` and the rows before the rule start are checked
/// exactly; otherwise, for a constant valuation `v`, by the relation with
/// `p(n, k) != 0` for `0 <= k < v`.
fn lower_support_hypothesis(
    shape: &RecurrenceShape,
    p: &LinearFactors,
    valuation: &Affine,
    earlier: &[(usize, &Vec<Q>)],
    first: usize,
) -> LiftHypothesis {
    let name = "lower_support";
    let statement = "[x^k] P_n = 0 for k < v(n), the valuation";
    let terms_hold = shape.terms.iter().all(|(&(lag, order), coefficient)| {
        (0..x_width(coefficient))
            .filter(|power| !x_column(coefficient, *power).is_empty())
            .all(|power| {
                Q::from_integer(BigInt::from(power))
                    - Q::from_integer(BigInt::from(order))
                    - &valuation.1 * Q::from_integer(BigInt::from(lag))
                    >= Q::zero()
            })
    });
    let rows_hold = earlier.iter().all(|(n, row)| {
        let bound = &valuation.0 + &valuation.1 * Q::from_integer(BigInt::from(*n));
        (0..row.len()).all(|k| row[k].is_zero() || Q::from_integer(BigInt::from(k)) >= bound)
    });
    if terms_hold && rows_hold {
        return exact_hypothesis(
            name,
            statement,
            Ok(format!(
                "induction on the rule: every term x^e D^d P_(n-r) has e - d - ({}) r >= 0, and \
                 the rows before the rule start are checked exactly",
                affine_text(&(valuation.1.clone(), Q::zero()))
            )),
        );
    }
    let Some(v) = (valuation.1.is_zero())
        .then(|| valuation.0.to_integer().to_usize())
        .flatten()
    else {
        return unresolved_hypothesis(
            name,
            statement,
            "neither the rule nor the relation bounds the valuation".into(),
        );
    };
    let mut product = IndexFunction::constant(Q::one());
    for k in 0..v {
        product = times(
            &product,
            &IndexFunction::polynomial(p.at(&affine_constant(k as i64))),
        );
    }
    sign_hypothesis(
        name,
        &format!("{statement}: p(n, k) != 0 for 0 <= k < {v}"),
        "prod_(k < v) p(n, k)",
        &product,
        first,
        IndexSignCondition::Nonzero,
    )
}

/// The degree bound `deg P_n <= d(n)` (so `[x^k] P_n = 0` for `k > d(n)`):
/// by induction on the rule when every term `c x^e D^d P_(n-r)` has
/// `e - d - d_1 r <= 0` and the rows before the rule start are checked
/// exactly; otherwise by the relation, with the factor `k - d(n) - 1` of `q`,
/// `c_p(n) != 0` and `u(n) + d(n) + 1 > 0` for every factor `k + u` of `p`.
fn upper_support_hypotheses(
    shape: &RecurrenceShape,
    p: &LinearFactors,
    q: &LinearFactors,
    degree: &Affine,
    earlier: &[(usize, &Vec<Q>)],
    first: usize,
) -> Vec<LiftHypothesis> {
    let name = "upper_support";
    let statement = "[x^k] P_n = 0 for k > d(n), the degree";
    let terms_hold = shape.terms.iter().all(|(&(lag, order), coefficient)| {
        x_degree(coefficient).is_none_or(|power| {
            Q::from_integer(BigInt::from(power))
                - Q::from_integer(BigInt::from(order))
                - &degree.1 * Q::from_integer(BigInt::from(lag))
                <= Q::zero()
        })
    });
    let rows_hold = earlier.iter().all(|(n, row)| {
        Q::from_integer(BigInt::from(row.len().saturating_sub(1)))
            <= &degree.0 + &degree.1 * Q::from_integer(BigInt::from(*n))
    });
    if terms_hold && rows_hold {
        return vec![exact_hypothesis(
            name,
            statement,
            Ok(format!(
                "induction on the rule: every term x^e D^d P_(n-r) has e - d - ({}) r <= 0, and \
                 the rows before the rule start are checked exactly",
                affine_text(&(degree.1.clone(), Q::zero()))
            )),
        )];
    }
    let top = affine_neg(&affine_add(degree, &affine_constant(1)));
    let mut out = vec![exact_hypothesis(
        &format!("{name}_factor"),
        "q(n, k) has the factor k - d(n) - 1",
        if q.roots.contains(&top) {
            Ok(format!("d(n) = {}", affine_text(degree)))
        } else {
            Err((
                None,
                format!("no factor k - d(n) - 1, d(n) = {}", affine_text(degree)),
            ))
        },
    )];
    out.push(sign_hypothesis(
        &format!("{name}_scalar"),
        "c_p(n) != 0, the scalar factor of p",
        "c_p(n)",
        &IndexFunction::polynomial(p.scalar.clone()),
        first,
        IndexSignCondition::Nonzero,
    ));
    // A crude degree bound D(n) = D_0 + beta n from the rule: beta =
    // max (e - d) / r over the terms, D_0 covering the rows checked exactly.
    let beta = shape
        .terms
        .iter()
        .filter_map(|(&(lag, order), coefficient)| {
            x_degree(coefficient).map(|power| {
                (Q::from_integer(BigInt::from(power)) - Q::from_integer(BigInt::from(order)))
                    / Q::from_integer(BigInt::from(lag))
            })
        })
        .fold(
            degree.1.clone(),
            |acc, value| if value > acc { value } else { acc },
        );
    let crude = earlier
        .iter()
        .map(|(n, row)| {
            Q::from_integer(BigInt::from(row.len().saturating_sub(1)))
                - &beta * Q::from_integer(BigInt::from(*n))
        })
        .reduce(|a, b| if a > b { a } else { b })
        .map(|intercept| (intercept, beta.clone()));
    for (index, u) in p.roots.iter().enumerate() {
        let below = sign_hypothesis(
            &format!("{name}_root_{index}"),
            "u(n) + d(n) + 1 > 0, so p(n, k) != 0 for every integer k > d(n) and [x^k] P_n = 0 \
             there",
            "u(n) + d(n) + 1",
            &affine_function(&affine_add(&affine_add(u, degree), &affine_constant(1))),
            first,
            IndexSignCondition::Positive,
        );
        if below.status == LiftHypothesisStatus::Verified {
            out.push(below);
            continue;
        }
        let above = crude.as_ref().map(|crude| {
            sign_hypothesis(
                &format!("{name}_root_{index}"),
                &format!(
                    "u(n) + D(n) < 0 for the degree bound D(n) = {} of the rule (deg P_n <= D(n) \
                     by induction: D(n) >= e - d + D(n - r) for every term x^e D^d P_(n-r), the \
                     rows before checked exactly), so the root -u(n) of p lies above every \
                     degree and p(n, k) != 0 for d(n) < k <= D(n)",
                    IndexFunction::polynomial(vec![crude.0.clone(), crude.1.clone()]).text()
                ),
                "u(n) + D(n)",
                &affine_function(&affine_add(u, crude)),
                first,
                IndexSignCondition::Negative,
            )
        });
        match above {
            Some(above) if above.status == LiftHypothesisStatus::Verified => out.push(above),
            _ => out.push(below),
        }
    }
    out
}

/// The hypotheses of a matching (module documentation), decided for
/// `n >= first`.
fn matching_hypotheses(
    ratio: &Ratio,
    matching: &Matching,
    first: usize,
    upper: Option<&Affine>,
) -> Vec<LiftHypothesis> {
    let mut out = Vec::new();
    let one = affine_constant(1);
    let size = complement_size(&ratio.numerators[matching.base_numerator])
        .unwrap_or_else(|| affine_constant(0));
    out.push(sign_hypothesis(
        "base_size",
        &format!(
            "the base factor {} and the factor j give binom(M, j) with M(n) >= 0",
            factor_text(&ratio.numerators[matching.base_numerator])
        ),
        "M(n)",
        &affine_function(&size),
        first,
        IndexSignCondition::Nonnegative,
    ));
    if let Some(upper) = upper {
        out.push(sign_hypothesis(
            "base_within_degree",
            "M(n) <= d(n) - v (reversal at the degree)",
            "d(n) - v - M(n)",
            &affine_function(&affine_add(upper, &affine_neg(&size))),
            first,
            IndexSignCondition::Nonnegative,
        ));
    }
    for index in &matching.complements {
        let other =
            complement_size(&ratio.numerators[*index]).unwrap_or_else(|| affine_constant(0));
        out.push(sign_hypothesis(
            &format!("complement_{index}"),
            &format!(
                "the factor {} gives k -> M'!/(M' - k)! with M'(n) >= M(n) (Polya--Schur)",
                factor_text(&ratio.numerators[*index])
            ),
            "M'(n) - M(n)",
            &affine_function(&affine_add(&other, &affine_neg(&size))),
            first,
            IndexSignCondition::Nonnegative,
        ));
    }
    for (numerator, denominator) in &matching.pairs {
        let (a, b) = (
            &ratio.numerators[*numerator],
            &ratio.denominators[*denominator],
        );
        out.push(sign_hypothesis(
            &format!("pair_{numerator}_{denominator}_length"),
            &format!(
                "({}) / ({}) pairs to the polynomial (k + b + 1) ... (k + a): a - b is an \
                 integer >= 0 (Laguerre)",
                factor_text(a),
                factor_text(b)
            ),
            "a(n) - b(n)",
            &affine_function(&affine_add(a, &affine_neg(b))),
            first,
            IndexSignCondition::Nonnegative,
        ));
        out.push(sign_hypothesis(
            &format!("pair_{numerator}_{denominator}_roots"),
            "b(n) + 1 > 0: every root of the paired polynomial is negative",
            "b(n) + 1",
            &affine_function(&affine_add(b, &one)),
            first,
            IndexSignCondition::Positive,
        ));
    }
    for index in &matching.reciprocals {
        let b = &ratio.denominators[*index];
        out.push(sign_hypothesis(
            &format!("reciprocal_{index}"),
            &format!(
                "the factor {} gives k -> 1/(b + 1)_k with b(n) + 1 > 0",
                factor_text(b)
            ),
            "b(n) + 1",
            &affine_function(&affine_add(b, &one)),
            first,
            IndexSignCondition::Positive,
        ));
    }
    let sign = if (matching.complements.len() + 1) % 2 == 0 {
        Q::one()
    } else {
        -Q::one()
    };
    out.push(sign_hypothesis(
        "scale_positive",
        "C_eff(n) = C(n) (-1)^(complements) > 0 (the scaling x -> C_eff x)",
        "C_eff(n)",
        &ratio.scale.scale(&sign),
        first,
        IndexSignCondition::Positive,
    ));
    out
}

/// `H(n) P_n = 0` for every `n >= from`: an invariant annihilator system
/// of the rule whose module contains `H` (module documentation).
fn relation_hypothesis(
    operator: &DOp,
    prepared: &crate::recurrence::lift::hadamard_factor::prepare::PreparedLiftRule,
    from: usize,
    rows: &[Vec<Q>],
    first_index: usize,
) -> LiftHypothesis {
    let name = "row_relation";
    let statement = "H(n) P_n = 0 for every n >= the first base, H = p(n, theta) + x q(n, theta + \
                     1): an invariant row-annihilator system of the rule contains H";
    // A cheap rejection on the rows first.
    if let Some((n, _)) = rows
        .iter()
        .enumerate()
        .map(|(i, row)| (first_index + i, row))
        .filter(|(n, _)| *n >= from)
        .find(|(n, row)| !operator.apply_row(row, *n).iter().all(Zero::is_zero))
    {
        return exact_hypothesis(
            name,
            statement,
            Err((Some(n), format!("H does not kill P_{n}"))),
        );
    }
    let decomposition = match find_invariant_annihilator_system(
        operator,
        prepared,
        from,
        MAX_HADAMARD_OPERATOR_TERMS,
    ) {
        Ok(Some(decomposition)) => decomposition,
        Ok(None) => {
            return unresolved_hypothesis(
                name,
                statement,
                "no invariant row-annihilator system of the rule containing H was found".into(),
            )
        }
        Err(error) => {
            return unresolved_hypothesis(
                name,
                statement,
                format!("the annihilator search stopped: {error}"),
            )
        }
    };
    if decomposition.multiplier.is_zero() {
        return unresolved_hypothesis(name, statement, "the multiplier mu is zero".into());
    }
    let restart = decomposition.system.restart_bases(from);
    if let Err(error) =
        check_system_rows(&decomposition.system.generators, prepared, from, &restart)
    {
        return unresolved_hypothesis(name, statement, error.to_string());
    }
    let exceptional = vanishing_bases(&decomposition.multiplier, from);
    let last = exceptional.iter().copied().max().unwrap_or(from);
    if last >= first_index + rows.len() {
        return unresolved_hypothesis(
            name,
            statement,
            format!("the exceptional base {last} lies beyond the generated rows"),
        );
    }
    if let Some(n) = exceptional.iter().find(|n| {
        !operator
            .apply_row(&rows[*n - first_index], **n)
            .iter()
            .all(Zero::is_zero)
    }) {
        return exact_hypothesis(
            name,
            statement,
            Err((Some(*n), format!("H does not kill P_{n} at a zero of mu"))),
        );
    }
    exact_hypothesis(
        name,
        statement,
        Ok(format!(
            "invariant system with {} generator(s) from base {from} (restart bases {:?}), \
             mu H = sum_l Z_l A^(l), H P_m = 0 checked exactly at the zeros {:?} of mu",
            decomposition.system.generators.len(),
            restart,
            exceptional
        )),
    )
}

/// `U_0(n) != 0` with one sign for every `n >= start - max_lag` from the
/// shared leading-coefficient module: the scalar recurrence
/// `f(n) U_0(n) = sum_r sigma_r(n) U_0(n - r)` (valid for `n >= start`),
/// cleared of the denominators `L(n)` of the `sigma_r`, is a recurrence of
/// constant rows `f(n) L(n) U_0(n) = sum_r (L sigma_r)(n) U_0(n - r)` with
/// the seeds `U_0(start - max_lag), ..., U_0(start - 1)`, and a proved
/// certificate with degree zero and one sign from the first seed on is the
/// claim (the module also checks `f(n) L(n) != 0`).  `None` when it is not
/// proved.
fn shared_first_coefficient_hypothesis(
    lag_coefficients: &BTreeMap<usize, IndexFunction>,
    left: &[Q],
    values: &[(usize, Q)],
    start: usize,
    max_lag: usize,
) -> Option<LiftHypothesis> {
    let first = start.checked_sub(max_lag)?;
    let seeds: Vec<Vec<Q>> = values
        .iter()
        .filter(|(n, _)| *n < start)
        .map(|(_, value)| vec![value.clone()])
        .collect();
    if seeds.len() != max_lag || lag_coefficients.is_empty() {
        return None;
    }
    let common = lag_coefficients
        .values()
        .fold(vec![Q::one()], |acc, c| vec_poly::mul(&acc, &c.denominator));
    let constant_in_x = |p: Vec<Q>| {
        BivarPoly::from_coefficients_normalized(p.into_iter().map(|c| vec![c]).collect())
    };
    let terms = lag_coefficients
        .iter()
        .map(|(&lag, c)| {
            let others = lag_coefficients
                .iter()
                .filter(|(other, _)| **other != lag)
                .fold(vec![Q::one()], |acc, (_, o)| {
                    vec_poly::mul(&acc, &o.denominator)
                });
            RecurrenceTerm::new(lag, 0, constant_in_x(vec_poly::mul(&c.numerator, &others)))
        })
        .collect();
    let scalar = Recurrence {
        terms,
        denominator: Some(constant_in_x(vec_poly::mul(left, &common))),
        inhomogeneous: None,
    };
    shared_degree_law_hypothesis(
        "first_coefficient",
        &format!(
            "the first coefficient U_0(n) is nonzero with one sign for every n >= {first} (its \
             scalar recurrence, as a family of constant rows, by the shared leading-coefficient \
             module)"
        ),
        (&scalar, first, &seeds),
        (first, Some(0)),
        |step| step == 0,
        |_| false,
    )
}

/// `U_0(n) != 0` with one sign for every `n >= start - max_lag`: a rational
/// `rho = B / A` (degrees at most [`MAX_FIRST_COEFFICIENT_DEGREE`]) with
/// `A(m) U_0(m) = B(m) U_0(m - 1)` on every computed pair (an exact fit),
/// the identity `sum_r sigma_r(n) prod_(i < r) A(n - i) / B(n - i) = f(n)`
/// of rational functions (so the hypergeometric sequence `H(m) = rho(m)
/// H(m - 1)` started at `U_0(start - max_lag)` solves the scalar recurrence
/// of `U_0`, whose leader `f(n)` does not vanish, and agrees with `U_0` on
/// the seeds, hence everywhere), `rho(m) > 0` for `m > start - max_lag`
/// and `U_0(start - max_lag) != 0`.
fn first_coefficient_hypothesis(
    lag_coefficients: &BTreeMap<usize, IndexFunction>,
    left: &IndexFunction,
    values: &[(usize, Q)],
    start: usize,
    max_lag: usize,
) -> LiftHypothesis {
    let name = "first_coefficient";
    let statement = "the first coefficient U_0(n) is nonzero with one sign for every n from the \
                     induction seeds on (a hypergeometric solution U_0(n) = rho(n) U_0(n - 1) of \
                     its scalar recurrence with rho(n) > 0)";
    let fail = |evidence: String| exact_hypothesis(name, statement, Err((None, evidence)));
    if lag_coefficients.is_empty() {
        return fail("every lag coefficient of the scalar recurrence of U_0 vanishes".into());
    }
    let Some((first_seed, first_value)) = values.first() else {
        return fail("no computed first coefficients".into());
    };
    if first_value.is_zero() {
        return fail(format!("U_0({first_seed}) = 0"));
    }
    let pairs: Vec<(Q, &Q, &Q)> = values
        .windows(2)
        .map(|w| (Q::from_integer(BigInt::from(w[1].0)), &w[1].1, &w[0].1))
        .collect();
    let mut rho = None;
    for degree in 0..=MAX_FIRST_COEFFICIENT_DEGREE {
        let width = 2 * (degree + 1);
        if pairs.len() < width + 3 {
            break;
        }
        // A(m) U_0(m) - B(m) U_0(m - 1) = 0, unknowns (A_0..A_d, B_0..B_d).
        let equations: Vec<Vec<Q>> = pairs
            .iter()
            .map(|(m, current, previous)| {
                let mut equation = Vec::with_capacity(width);
                for i in 0..=degree {
                    equation.push(*current * m.pow(i as i32));
                }
                for i in 0..=degree {
                    equation.push(-(*previous * m.pow(i as i32)));
                }
                equation
            })
            .collect();
        let integral: Vec<Vec<BigInt>> = equations
            .iter()
            .map(|equation| {
                let denominator = equation
                    .iter()
                    .fold(BigInt::one(), |acc, value| acc.lcm(value.denom()));
                equation
                    .iter()
                    .map(|value| value.numer() * (&denominator / value.denom()))
                    .collect()
            })
            .collect();
        let Some(vector) = exact_kernel_vectors_by_group(&integral, &vec![0; width], 1)
            .into_iter()
            .next()
            .flatten()
        else {
            continue;
        };
        let solves = equations.iter().all(|equation| {
            equation
                .iter()
                .zip(&vector)
                .fold(Q::zero(), |acc, (a, b)| acc + a * b)
                .is_zero()
        });
        let a = vec_poly::trim(&vector[..=degree]);
        let b = vec_poly::trim(&vector[degree + 1..]);
        if solves && !vec_poly::is_zero(&a) && !vec_poly::is_zero(&b) {
            rho = Some(IndexFunction::new(b, a));
            break;
        }
    }
    // Without a hypergeometric solution: positive lag coefficients and
    // seeds of one sign keep U_0 of that sign.
    let positivity = || -> Option<String> {
        let seeds: Vec<&Q> = values
            .iter()
            .filter(|(n, _)| *n < start)
            .map(|(_, value)| value)
            .collect();
        let one_sign = seeds.len() == max_lag
            && (seeds.iter().all(|v| v.is_positive()) || seeds.iter().all(|v| v.is_negative()));
        let positive = lag_coefficients.values().all(|coefficient| {
            matches!(
                coefficient
                    .div(left)
                    .decide(start, IndexSignCondition::Positive),
                IndexSignDecision::Holds { .. }
            )
        });
        (one_sign && positive).then(|| {
            format!(
                "every sigma_r(n) / f(n) > 0 for n >= {start} (exact sign checker: {}) and the \
                 seeds U_0({}), ..., U_0({}) are nonzero of one sign",
                lag_coefficients
                    .iter()
                    .map(|(lag, c)| format!("sigma_{lag}/f = {}", c.div(left).text()))
                    .collect::<Vec<_>>()
                    .join(", "),
                start - max_lag,
                start - 1
            )
        })
    };
    let Some(rho) = rho else {
        if let Some(evidence) = positivity() {
            return exact_hypothesis(name, statement, Ok(evidence));
        }
        return unresolved_hypothesis(
            name,
            statement,
            format!(
                "no ratio U_0(n) / U_0(n - 1) of degree at most {MAX_FIRST_COEFFICIENT_DEGREE} \
                 fits the computed first coefficients, and the lag coefficients are not all \
                 positive"
            ),
        );
    };
    // sum_r sigma_r(n) prod_(i < r) 1 / rho(n - i) - f(n) = 0.
    let mut sum = left.negated();
    for (lag, coefficient) in lag_coefficients {
        let mut term = coefficient.clone();
        for i in 0..*lag {
            term = term.div(&rho.shifted(-(i as i64)));
        }
        sum = sum.add(&term);
    }
    if !vec_poly::is_zero(&sum.numerator) {
        if let Some(evidence) = positivity() {
            return exact_hypothesis(name, statement, Ok(evidence));
        }
        return unresolved_hypothesis(
            name,
            statement,
            format!(
                "the fitted ratio rho(n) = {} does not solve the scalar recurrence of U_0",
                rho.text()
            ),
        );
    }
    let from = start - max_lag + 1;
    let mut hypothesis = sign_hypothesis(
        name,
        statement,
        "rho(n)",
        &rho,
        from,
        IndexSignCondition::Positive,
    );
    hypothesis.evidence = format!(
        "U_0(n) = rho(n) U_0(n - 1) from U_0({first_seed}) = {first_value}, the identity \
         sum_r sigma_r(n) prod_(i < r) 1/rho(n - i) = f(n) holds exactly; {}",
        hypothesis.evidence
    );
    hypothesis
}
