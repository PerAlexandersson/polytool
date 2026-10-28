//! The bounded candidate search of [`find_symmetric_companion_lifts`]:
//! schedule inference, the constant-term fit of `a_n`, and the exact
//! linearization of the second coupled equation in the remaining parameters.

use super::certificate::{
    build_certificate, verify_symmetric_companion_certificate, SymmetricCompanionIdentityStatus,
};
use super::groebner::{rational_points, MultivariatePolynomial, SolverStop};
use super::operators::{evaluate_at_index, CompanionParameters};
use super::report::*;
use super::source::PreparedSource;
use crate::linalg::{rational_kernel_with_certificate, RationalKernelError, RationalKernelLimits};
use crate::recurrence::*;
use crate::vec_poly;
use num_traits::{One, Zero};
use std::collections::{BTreeMap, BTreeSet};

type Q = BigRational;

/// Largest supported degree in `n` of the parameter ansatz.
pub const MAX_SYMMETRIC_COMPANION_PARAMETER_DEGREE: usize = 1;

/// Bounds of [`find_symmetric_companion_lifts`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SymmetricCompanionLiftOptions {
    /// The symmetry degree `d_(n0)` at the first index; the schedule is
    /// `d_n = d_(n0) + (n - n0)`.  `None` infers it from the first nonzero row
    /// (`d = valuation + degree`) as a candidate that the certificate then
    /// checks.  Default `None`.
    pub symmetry_degree_at_first_index: Option<usize>,
    /// Largest degree in `n` of `a_n`, `c_n` and `epsilon_n` (0: constants,
    /// 1: affine); at most [`MAX_SYMMETRIC_COMPANION_PARAMETER_DEGREE`].
    /// Degrees are tried in increasing order.  Default 1.
    pub max_parameter_degree_in_n: usize,
    /// Number of rows `A_(n0), ..., A_(n0 + fit_rows - 1)` used by the
    /// search; at least 4.  Default 10.
    pub fit_rows: usize,
    /// Number of indices replayed exactly by the verifier's defensive row
    /// self-check (not part of the proof); at least 2.  Default 8.
    pub self_check_rows: usize,
    /// Largest number of unknown monomials of a linearized system; a larger
    /// ansatz is skipped and sets
    /// [`exhausted_budget`](SymmetricCompanionSearchSummary::exhausted_budget).
    /// Default 64.
    pub max_linearized_monomials: usize,
    /// Step budget (polynomial reductions) of the exact Gröbner solve used
    /// when the linearized kernel has dimension other than one.  Default
    /// 20000.
    pub max_groebner_steps: usize,
    /// Largest number of candidates certified per search.  Default 8.
    pub max_candidates: usize,
}

impl Default for SymmetricCompanionLiftOptions {
    fn default() -> Self {
        Self {
            symmetry_degree_at_first_index: None,
            max_parameter_degree_in_n: 1,
            fit_rows: 10,
            self_check_rows: 8,
            max_linearized_monomials: 64,
            max_groebner_steps: 20_000,
            max_candidates: 8,
        }
    }
}

const KERNEL_LIMITS: RationalKernelLimits = RationalKernelLimits {
    max_entries: 400_000,
    max_operations: 50_000_000,
    max_coefficient_bits: 20_000,
};

// ---------------------------------------------------------------------------
// Polynomials in the unknown parameter coefficients
// ---------------------------------------------------------------------------

/// `sum_mu m_mu(x) theta^mu`: a polynomial in the unknown coefficients
/// `theta` with coefficients in `Q[x]`, keyed by exponent vectors.
#[derive(Debug, Clone, Default)]
struct UnknownPolynomial(BTreeMap<Vec<u8>, Vec<Q>>);

impl UnknownPolynomial {
    fn add_term(&mut self, monomial: Vec<u8>, coefficient: &[Q]) {
        if vec_poly::is_zero(coefficient) {
            return;
        }
        let entry = self
            .0
            .entry(monomial.clone())
            .or_insert_with(|| vec![Q::zero()]);
        *entry = vec_poly::add(entry, coefficient);
        if vec_poly::is_zero(entry) {
            self.0.remove(&monomial);
        }
    }

    fn constant(variables: usize, coefficient: &[Q]) -> Self {
        let mut p = Self::default();
        p.add_term(vec![0; variables], coefficient);
        p
    }

    fn add(&self, other: &Self) -> Self {
        let mut sum = self.clone();
        for (monomial, coefficient) in &other.0 {
            sum.add_term(monomial.clone(), coefficient);
        }
        sum
    }

    fn negate(&self) -> Self {
        let mut out = Self::default();
        for (monomial, coefficient) in &self.0 {
            out.add_term(monomial.clone(), &vec_poly::neg(coefficient));
        }
        out
    }

    fn sub(&self, other: &Self) -> Self {
        self.add(&other.negate())
    }

    fn mul(&self, other: &Self) -> Self {
        let mut product = Self::default();
        for (left, a) in &self.0 {
            for (right, b) in &other.0 {
                let monomial = left.iter().zip(right).map(|(x, y)| x + y).collect();
                product.add_term(monomial, &vec_poly::mul(a, b));
            }
        }
        product
    }

    /// Apply a `Q`-linear map on `Q[x]` to every coefficient.
    fn map_x(&self, map: impl Fn(&[Q]) -> Vec<Q>) -> Self {
        let mut out = Self::default();
        for (monomial, coefficient) in &self.0 {
            out.add_term(monomial.clone(), &map(coefficient));
        }
        out
    }
}

/// One parameter `p(n) = sum_i p_i n^i` of the ansatz.
#[derive(Debug, Clone)]
enum ParameterAnsatz {
    /// Known exactly (fitted from constant terms).
    Known(Vec<Q>),
    /// Unknown: the variable index of `p_i` for `i = 0..=degree`.
    Unknown(Vec<usize>),
}

impl ParameterAnsatz {
    fn at(&self, n: usize, variables: usize) -> UnknownPolynomial {
        match self {
            Self::Known(p) => UnknownPolynomial::constant(variables, &[evaluate_at_index(p, n)]),
            Self::Unknown(indices) => {
                let mut out = UnknownPolynomial::default();
                let mut power = Q::one();
                for &index in indices {
                    let mut monomial = vec![0; variables];
                    monomial[index] = 1;
                    out.add_term(monomial, &[power.clone()]);
                    power *= Q::from_integer(n.into());
                }
                out
            }
        }
    }
}

// ---------------------------------------------------------------------------
// The constant-term fit of a_n
// ---------------------------------------------------------------------------

enum ConstantTermFit {
    Known(Vec<Q>),
    /// Every constant term vanishes: `a_n` stays unknown.
    Unknown,
    Incompatible(String),
    Underdetermined(String),
}

/// Exact Lagrange interpolation through `points`.
fn interpolate(points: &[(Q, Q)]) -> Vec<Q> {
    let mut result = vec![Q::zero()];
    for (i, (xi, yi)) in points.iter().enumerate() {
        let mut basis = vec![Q::one()];
        let mut denominator = Q::one();
        for (j, (xj, _)) in points.iter().enumerate() {
            if i != j {
                basis = vec_poly::mul(&basis, &[-xj.clone(), Q::one()]);
                denominator *= xi - xj;
            }
        }
        result = vec_poly::add(&result, &vec_poly::scale(&basis, &(yi / denominator)));
    }
    vec_poly::trim(&result)
}

/// `R_n(0) = 0` forces `a_n = A_(n+1)(0) / A_n(0)` wherever `A_n(0) != 0`.
fn fit_constant_terms(rows: &[Vec<Q>], first_index: usize, degree: usize) -> ConstantTermFit {
    let mut points = Vec::new();
    for (offset, pair) in rows.windows(2).enumerate() {
        let n = first_index + offset;
        let (current, next) = (&pair[0][0], &pair[1][0]);
        if current.is_zero() {
            if !next.is_zero() {
                return ConstantTermFit::Incompatible(format!(
                    "A_{n}(0) = 0 but A_{}(0) != 0, so R_{n}(0) != 0 for every a_{n}",
                    n + 1
                ));
            }
        } else {
            points.push((Q::from_integer(n.into()), next / current));
        }
    }
    if points.is_empty() {
        return ConstantTermFit::Unknown;
    }
    if points.len() <= degree {
        return ConstantTermFit::Underdetermined(format!(
            "only {} nonzero constant terms for a degree-{degree} fit of a_n",
            points.len()
        ));
    }
    let fitted = interpolate(&points[..=degree]);
    for (n, value) in &points {
        if &vec_poly::evaluate(&fitted, n) != value {
            return ConstantTermFit::Incompatible(format!(
                "the ratios A_(n+1)(0) / A_n(0) are not a polynomial of degree <= {degree} in n \
                 (they disagree at n = {n}), so R_n(0) = 0 fails for every such a_n"
            ));
        }
    }
    ConstantTermFit::Known(fitted)
}

// ---------------------------------------------------------------------------
// The linearized system
// ---------------------------------------------------------------------------

/// `d x f + x (1 - x) f'`.
fn derivative_part(f: &[Q], d: &Q) -> Vec<Q> {
    let x_f = vec_poly::shift(f);
    let derivative = vec_poly::derivative(f);
    vec_poly::add(
        &vec_poly::scale(&x_f, d),
        &vec_poly::mul(&derivative, &[Q::zero(), Q::one(), -Q::one()]),
    )
}

struct Ansatz {
    variables: usize,
    a: ParameterAnsatz,
    c: ParameterAnsatz,
    epsilon: ParameterAnsatz,
    epsilon_indices: Vec<usize>,
}

impl Ansatz {
    fn new(a: Option<Vec<Q>>, degree: usize) -> Self {
        let mut next = 0;
        let mut take = |count: usize| {
            let indices: Vec<usize> = (next..next + count).collect();
            next += count;
            indices
        };
        let a = match a {
            Some(fitted) => ParameterAnsatz::Known(fitted),
            None => ParameterAnsatz::Unknown(take(degree + 1)),
        };
        let c = ParameterAnsatz::Unknown(take(degree + 1));
        let epsilon_indices = take(degree + 1);
        Self {
            variables: next,
            a,
            c,
            epsilon: ParameterAnsatz::Unknown(epsilon_indices.clone()),
            epsilon_indices,
        }
    }

    /// `R_n = A_(n+1) - a_n (1 + x) A_n - c_n (d_n x A_n + x (1 - x) A_n')`.
    fn remainder(&self, rows: &[Vec<Q>], offset: usize, n: usize, d: &Q) -> UnknownPolynomial {
        let v = self.variables;
        let one_plus_x = vec_poly::mul(&rows[offset], &[Q::one(), Q::one()]);
        UnknownPolynomial::constant(v, &rows[offset + 1])
            .sub(
                &self
                    .a
                    .at(n, v)
                    .mul(&UnknownPolynomial::constant(v, &one_plus_x)),
            )
            .sub(&self.c.at(n, v).mul(&UnknownPolynomial::constant(
                v,
                &derivative_part(&rows[offset], d),
            )))
    }

    /// `G_n = eps_n R_(n+1) - eps_(n+1) Vt_n R_n - eps_n eps_(n+1) (c_n - eps_n) x A_n`,
    /// which vanishes for the true parameters (it is `(E A)_n` with `B`
    /// eliminated) and has a factor `eps` in every term.
    fn equation(&self, rows: &[Vec<Q>], offset: usize, n: usize, d: &Q) -> UnknownPolynomial {
        let v = self.variables;
        let r_n = self.remainder(rows, offset, n, d);
        let r_next = self.remainder(rows, offset + 1, n + 1, &(d + Q::one()));
        let (a, c) = (self.a.at(n, v), self.c.at(n, v));
        let (eps, eps_next) = (self.epsilon.at(n, v), self.epsilon.at(n + 1, v));
        // Vt_n f = (a + c - eps)(1 + x) f + c [(d - 1) x f + x (1 - x) f' - (1 - x) f].
        let d_minus_one = d - Q::one();
        let multiplier = a.add(&c).sub(&eps);
        let v_tilde_r = multiplier
            .mul(&r_n.map_x(|f| vec_poly::mul(f, &[Q::one(), Q::one()])))
            .add(&c.mul(&r_n.map_x(|f| {
                vec_poly::sub(
                    &derivative_part(f, &d_minus_one),
                    &vec_poly::mul(f, &[Q::one(), -Q::one()]),
                )
            })));
        let x_a = UnknownPolynomial::constant(v, &vec_poly::shift(&rows[offset]));
        eps.mul(&r_next)
            .sub(&eps_next.mul(&v_tilde_r))
            .sub(&eps.mul(&eps_next).mul(&c.sub(&eps)).mul(&x_a))
    }
}

enum Linearization {
    /// Candidate points `theta` (coefficients of the ansatz unknowns).
    Points(Vec<Vec<Q>>, String),
    NoSolution(String),
    Underdetermined(String),
    OverBudget(String),
}

/// Solve the homogeneous linear system in the monomials of the equations;
/// recover a point from a one-dimensional kernel, and otherwise solve the
/// independent equations exactly with a lex Gröbner basis (with
/// `z epsilon(n0) = 1` excluding vanishing couplings at `n0`).
fn solve_linearized(
    ansatz: &Ansatz,
    equations: &[UnknownPolynomial],
    first_index: usize,
    options: &SymmetricCompanionLiftOptions,
) -> (Linearization, usize, usize, Option<usize>) {
    let monomials: BTreeSet<Vec<u8>> = equations.iter().flat_map(|e| e.0.keys().cloned()).collect();
    let columns: Vec<Vec<u8>> = monomials.into_iter().collect();
    if columns.len() > options.max_linearized_monomials {
        let reason = format!(
            "{} monomials exceed max_linearized_monomials = {}",
            columns.len(),
            options.max_linearized_monomials
        );
        return (Linearization::OverBudget(reason), columns.len(), 0, None);
    }
    let position: BTreeMap<&Vec<u8>, usize> =
        columns.iter().enumerate().map(|(i, m)| (m, i)).collect();
    let mut matrix: Vec<Vec<Q>> = Vec::new();
    for equation in equations {
        let width = equation.0.values().map(Vec::len).max().unwrap_or(0);
        for power in 0..width {
            let mut row = vec![Q::zero(); columns.len()];
            for (monomial, coefficient) in &equation.0 {
                if let Some(value) = coefficient.get(power) {
                    row[position[monomial]] = value.clone();
                }
            }
            if row.iter().any(|c| !c.is_zero()) && !matrix.contains(&row) {
                matrix.push(row);
            }
        }
    }
    let equation_count = matrix.len();
    let kernel = match rational_kernel_with_certificate(&matrix, columns.len(), KERNEL_LIMITS) {
        Ok(kernel) => kernel,
        Err(RationalKernelError::BudgetExceeded(what)) => {
            let reason = format!("kernel computation exceeded its {what} limit");
            return (
                Linearization::OverBudget(reason),
                columns.len(),
                equation_count,
                None,
            );
        }
        Err(error) => {
            let reason = error.to_string();
            return (
                Linearization::OverBudget(reason),
                columns.len(),
                equation_count,
                None,
            );
        }
    };
    let dimension = kernel.basis.len();
    let outcome = match kernel.basis.as_slice() {
        [] => Linearization::NoSolution(
            "the linearized system has only the zero solution: no parameters with epsilon != 0 \
             in this ansatz"
                .into(),
        ),
        [vector] if recover_point(ansatz, &columns, vector).is_some() => {
            let theta = recover_point(ansatz, &columns, vector).unwrap_or_default();
            Linearization::Points(vec![theta], "one-dimensional kernel".into())
        }
        _ => {
            let independent: Vec<&Vec<Q>> = kernel
                .independent_rows
                .iter()
                .map(|&i| &matrix[i])
                .collect();
            solve_exactly(
                ansatz,
                &columns,
                &independent,
                first_index,
                options,
                dimension,
            )
        }
    };
    (outcome, columns.len(), equation_count, Some(dimension))
}

/// The exact fallback: rational points of the independent equations together
/// with `z epsilon(n0) - 1 = 0` (`z` the first, eliminated, variable).
fn solve_exactly(
    ansatz: &Ansatz,
    columns: &[Vec<u8>],
    rows: &[&Vec<Q>],
    first_index: usize,
    options: &SymmetricCompanionLiftOptions,
    dimension: usize,
) -> Linearization {
    let variables = ansatz.variables + 1;
    let lift_monomial = |monomial: &[u8], z: u32| -> Vec<u32> {
        std::iter::once(z)
            .chain(monomial.iter().map(|&e| u32::from(e)))
            .collect()
    };
    let mut system = Vec::new();
    for row in rows {
        let mut p = MultivariatePolynomial::zero(variables);
        for (monomial, coefficient) in columns.iter().zip(row.iter()) {
            p.add_term(lift_monomial(monomial, 0), coefficient);
        }
        system.push(p);
    }
    let mut nonvanishing = MultivariatePolynomial::zero(variables);
    let mut power = Q::one();
    for &index in &ansatz.epsilon_indices {
        let mut monomial = vec![0u8; ansatz.variables];
        monomial[index] = 1;
        nonvanishing.add_term(lift_monomial(&monomial, 1), &power);
        power *= Q::from_integer(first_index.into());
    }
    nonvanishing.add_term(vec![0; variables], &-Q::one());
    system.push(nonvanishing);
    match rational_points(&system, variables, options.max_groebner_steps) {
        Ok(points) if points.is_empty() => Linearization::NoSolution(format!(
            "the linearized kernel has dimension {dimension}; the exact Gröbner solve of the \
             {} independent equations has no rational point with epsilon(n0) != 0",
            rows.len()
        )),
        Ok(points) => {
            let count = points.len();
            let thetas = points.into_iter().map(|p| p[1..].to_vec()).collect();
            Linearization::Points(
                thetas,
                format!(
                    "linearized kernel of dimension {dimension}; exact Gröbner solve of the {} \
                     independent equations: {count} rational point(s)",
                    rows.len()
                ),
            )
        }
        Err(SolverStop::PositiveDimensional) => Linearization::Underdetermined(format!(
            "the linearized kernel has dimension {dimension} and the exact system has infinitely \
             many solutions; more rows (fit_rows) may determine the parameters"
        )),
        Err(SolverStop::Budget(reason)) => Linearization::OverBudget(format!(
            "linearized kernel of dimension {dimension}; exact solve stopped: {reason}"
        )),
    }
}

/// From a kernel vector `k`, find `t` and `theta` with `t k = m(theta)`.
fn recover_point(ansatz: &Ansatz, columns: &[Vec<u8>], vector: &[Q]) -> Option<Vec<Q>> {
    let v = ansatz.variables;
    let lookup = |exponents: &[(usize, u8)]| {
        let mut monomial = vec![0u8; v];
        for (index, power) in exponents {
            monomial[*index] += *power;
        }
        columns
            .iter()
            .position(|m| *m == monomial)
            .map_or_else(Q::zero, |i| vector[i].clone())
    };
    for &pivot in &ansatz.epsilon_indices {
        let linear = lookup(&[(pivot, 1)]);
        let square = lookup(&[(pivot, 2)]);
        if linear.is_zero() || square.is_zero() {
            continue;
        }
        let t = &square / (&linear * &linear);
        let theta: Vec<Q> = (0..v)
            .map(|index| {
                if index == pivot {
                    &t * &linear
                } else {
                    lookup(&[(index, 1), (pivot, 1)]) / &linear
                }
            })
            .collect();
        // Every column must be the monomial of theta, scaled by t.
        let consistent = columns.iter().zip(vector).all(|(monomial, value)| {
            let mut product = Q::one();
            for (index, power) in monomial.iter().enumerate() {
                for _ in 0..*power {
                    product *= &theta[index];
                }
            }
            product == &t * value
        });
        return consistent.then_some(theta);
    }
    None
}

fn parameter_values(ansatz: &ParameterAnsatz, theta: &[Q]) -> Vec<Q> {
    match ansatz {
        ParameterAnsatz::Known(p) => p.clone(),
        ParameterAnsatz::Unknown(indices) => vec_poly::trim(
            &indices
                .iter()
                .map(|&i| theta[i].clone())
                .collect::<Vec<_>>(),
        ),
    }
}

// ---------------------------------------------------------------------------
// Schedule inference
// ---------------------------------------------------------------------------

/// `d_(n0)` from the first nonzero row (`d = valuation + degree`), checked
/// against every nonzero row; `Err` with a reason otherwise.
pub(super) fn infer_schedule(rows: &[Vec<Q>]) -> Result<usize, String> {
    let mut inferred: Option<i64> = None;
    for (offset, row) in rows.iter().enumerate() {
        let Some(degree) = vec_poly::degree(row) else {
            continue;
        };
        let valuation = row.iter().position(|c| !c.is_zero()).unwrap_or(0);
        let d_here = i64::try_from(valuation + degree).unwrap_or(i64::MAX);
        let d0 = d_here - i64::try_from(offset).unwrap_or(i64::MAX);
        match inferred {
            None => inferred = Some(d0),
            Some(previous) if previous != d0 => {
                return Err(format!(
                    "the rows are not symmetric with a unit-step schedule: row offset {offset} \
                     has centre {d_here}, expected {}",
                    previous + i64::try_from(offset).unwrap_or(i64::MAX)
                ))
            }
            Some(_) => {}
        }
    }
    match inferred {
        None => Err("every searched row is zero; no schedule can be inferred".into()),
        Some(d0) if d0 < 0 => Err(format!(
            "the inferred symmetry degree d_(n0) = {d0} is negative"
        )),
        Some(d0) => Ok(usize::try_from(d0).unwrap_or(usize::MAX)),
    }
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// Recognize the palindromic sequence `A_n` defined by a scalar recurrence as
/// the symmetric part of a parent with a first-order differential rule, with
/// an exact all-index certificate.
///
/// The parent is `P_(n+1) = (a_n + b_n x) P_n + c_n x (1 - x) D P_n` with
/// `P_n = A_n + x B_n`, `A_n` symmetric at `d_n = d_(n0) + (n - n0)` and `B_n`
/// at `d_n - 1`; with `delta_n = b_n - a_n - c_n d_n`,
/// `epsilon_n = c_n - delta_n` the companion is forced:
/// `B_n = [A_(n+1) - U_n A_n] / (epsilon_n x)`,
/// `U_n = [a_n + (a_n + c_n d_n) x] + c_n x (1 - x) D`.  The search, the
/// certificate obligations and the preservation theorem behind
/// [`SymmetricCompanionRealRootednessStatus::TheoremApplies`] (with the
/// directed weak-interlacing convention) are described in `doc/algorithms.md`,
/// section "Symmetric-companion lifts", and in the source documentation of
/// this module group.
///
/// `recurrence` (`f(n) A_n = sum c(n,x) D^d A_(n-r)`, `f` free of `x`, no
/// `(-1)^n` factor, homogeneous) with `initial_rows[i] = A_(first_index + i)`
/// defines `A_n` for `n >= first_index`; supplied rows beyond the seeds must
/// agree with the recurrence.  The domain of the lift is `n >= first_index`.
///
/// Every candidate (one per ansatz degree at most) is turned into a
/// [`SymmetricCompanionCertificate`] and replayed by
/// [`verify_symmetric_companion_certificate`]; the reported statuses are the
/// verifier's.  The search stops after the first ansatz degree that yields a
/// certified lift.
///
/// # Errors
///
/// [`RecurrenceError::InvalidInput`] for invalid options or rows that disagree
/// with the recurrence; [`RecurrenceError::Unsupported`] for a recurrence
/// outside the supported class.  A search without candidates is not an error.
pub fn find_symmetric_companion_lifts(
    recurrence: &Recurrence,
    initial_rows: &[Vec<BigRational>],
    first_index: usize,
    options: &SymmetricCompanionLiftOptions,
) -> Result<SymmetricCompanionLiftReport, RecurrenceError> {
    if options.max_parameter_degree_in_n > MAX_SYMMETRIC_COMPANION_PARAMETER_DEGREE {
        return Err(RecurrenceError::InvalidInput(format!(
            "max_parameter_degree_in_n = {} exceeds the supported {}",
            options.max_parameter_degree_in_n, MAX_SYMMETRIC_COMPANION_PARAMETER_DEGREE
        )));
    }
    if options.fit_rows < 4 {
        return Err(RecurrenceError::InvalidInput(
            "fit_rows must be at least 4".into(),
        ));
    }
    if options.self_check_rows < 2 {
        return Err(RecurrenceError::InvalidInput(
            "self_check_rows must be at least 2".into(),
        ));
    }
    let source = PreparedSource::new(recurrence, initial_rows, first_index)?;
    let rows = source.rows(options.fit_rows)?;
    let mut report = SymmetricCompanionLiftReport::new(source.to_json(), options.fit_rows);

    let (d0, origin) = match options.symmetry_degree_at_first_index {
        Some(d0) => (d0, SymmetryScheduleOrigin::Supplied),
        None => match infer_schedule(&rows) {
            Ok(d0) => (d0, SymmetryScheduleOrigin::InferredFromRows),
            Err(reason) => {
                report.obstructions.push(reason);
                report.symmetry_schedule = SymmetryScheduleRecord::new(first_index, None, None);
                return Ok(report);
            }
        },
    };
    report.symmetry_schedule = SymmetryScheduleRecord::new(first_index, Some(d0), Some(origin));

    let mut seen: Vec<CompanionParameters> = Vec::new();
    for degree in 0..=options.max_parameter_degree_in_n {
        let (attempt, candidates) = run_attempt(&source, &rows, d0, degree, options);
        report.search.attempts.push(attempt);
        let mut certified_here = false;
        for parameters in candidates {
            if seen.contains(&parameters) || seen.len() >= options.max_candidates {
                continue;
            }
            seen.push(parameters.clone());
            match build_certificate(&source, &parameters) {
                Ok(certificate) => {
                    let verification = verify_symmetric_companion_certificate(
                        &certificate,
                        options.self_check_rows,
                    )?;
                    certified_here |=
                        verification.identity_status == SymmetricCompanionIdentityStatus::Certified;
                    report.lifts.push(SymmetricCompanionLift::new(
                        certificate,
                        verification,
                        &source,
                        &parameters,
                    ));
                }
                Err(reason) => report.obstructions.push(format!(
                    "degree-{degree} candidate a_n = {}, c_n = {}, epsilon_n = {} has no \
                     certificate: {reason}",
                    text(&parameters.a),
                    text(&parameters.c),
                    text(&parameters.epsilon)
                )),
            }
        }
        if certified_here {
            break;
        }
    }
    report.finish();
    Ok(report)
}

/// One ansatz degree: the attempt record and the candidates.
fn run_attempt(
    source: &PreparedSource,
    rows: &[Vec<Q>],
    d0: usize,
    degree: usize,
    options: &SymmetricCompanionLiftOptions,
) -> (SymmetricCompanionAnsatzAttempt, Vec<CompanionParameters>) {
    let n0 = source.first_index;
    let mut attempt = SymmetricCompanionAnsatzAttempt {
        parameter_degree_in_n: degree,
        a_from_constant_terms: false,
        unknowns: 0,
        monomials: 0,
        equations: 0,
        kernel_dimension: None,
        outcome: SymmetricCompanionAttemptOutcome::NoSolution,
        detail: String::new(),
    };
    let fitted = match fit_constant_terms(rows, n0, degree) {
        ConstantTermFit::Known(a) => Some(a),
        ConstantTermFit::Unknown => None,
        ConstantTermFit::Incompatible(reason) => {
            attempt.outcome = SymmetricCompanionAttemptOutcome::ConstantTermsIncompatible;
            attempt.detail = reason;
            return (attempt, Vec::new());
        }
        ConstantTermFit::Underdetermined(reason) => {
            attempt.outcome = SymmetricCompanionAttemptOutcome::Underdetermined;
            attempt.detail = reason;
            return (attempt, Vec::new());
        }
    };
    attempt.a_from_constant_terms = fitted.is_some();
    let ansatz = Ansatz::new(fitted, degree);
    attempt.unknowns = ansatz.variables;
    let equations: Vec<UnknownPolynomial> = (0..rows.len() - 2)
        .map(|offset| {
            let d = Q::from_integer((d0 + offset).into());
            ansatz.equation(rows, offset, n0 + offset, &d)
        })
        .collect();
    let (outcome, monomials, equation_count, dimension) =
        solve_linearized(&ansatz, &equations, n0, options);
    attempt.monomials = monomials;
    attempt.equations = equation_count;
    attempt.kernel_dimension = dimension;
    let mut candidates = Vec::new();
    match outcome {
        Linearization::Points(points, how) => {
            let mut texts = Vec::new();
            for theta in points {
                let parameters = CompanionParameters::new(
                    parameter_values(&ansatz.a, &theta),
                    parameter_values(&ansatz.c, &theta),
                    parameter_values(&ansatz.epsilon, &theta),
                    n0,
                    d0,
                );
                texts.push(format!(
                    "(a_n = {}, c_n = {}, epsilon_n = {})",
                    text(&parameters.a),
                    text(&parameters.c),
                    text(&parameters.epsilon)
                ));
                candidates.push(parameters);
            }
            attempt.outcome = if candidates.is_empty() {
                SymmetricCompanionAttemptOutcome::NoSolution
            } else {
                SymmetricCompanionAttemptOutcome::CandidateFound
            };
            attempt.detail = format!("{how}: {}", texts.join(", "));
        }
        Linearization::NoSolution(reason) => attempt.detail = reason,
        Linearization::Underdetermined(reason) => {
            attempt.outcome = SymmetricCompanionAttemptOutcome::Underdetermined;
            attempt.detail = reason;
        }
        Linearization::OverBudget(reason) => {
            attempt.outcome = SymmetricCompanionAttemptOutcome::BudgetExceeded;
            attempt.detail = reason;
        }
    }
    (attempt, candidates)
}

fn text(p: &[Q]) -> String {
    crate::recurrence::implicit::algebra::poly_text(p, "n")
}
