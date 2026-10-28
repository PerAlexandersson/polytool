//! The factor menu of the Hadamard-factor lifts: exact values `u(n, k)`, the
//! Gamma-function form behind the conjugation ratios, the support conditions
//! under which the conjugation is valid on the whole index domain, and the
//! transfer theorem for each factor kind.
//!
//! Every factor is a product of one or two [`HadamardFactorTerm`]s.  In
//! Gamma form a term is `C prod Gamma(l_j)^(e_j)` with linear arguments
//! `l_j` of three kinds (`n + c`, `k + c`, `n - k + c`) and the convention
//! `1 / Gamma(m) = 0` at the poles `m = 0, -1, -2, ...`; for example
//! `binom(n + a, k + b) = Gamma(n + a + 1) / (Gamma(k + b + 1) Gamma(n + a -
//! k - b + 1))` for `n + a >= 0`.  The ratios `u(n + i, k - s) / u(n + R, k)`
//! of the conjugation are then products of Pochhammer symbols in `(n, k)`.

use crate::recurrence::implicit::algebra::{binomial, factorial, q, Biv};
use crate::recurrence::{parse_rational_coeff, RecurrenceError};
use crate::BigRational;
use num_traits::{One, Signed, Zero};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

type Q = BigRational;

/// One elementary Hadamard factor `u(n, k)` (the coefficient of `x^k` in row
/// `n` is multiplied by `u(n, k)`).
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum HadamardFactorTerm {
    /// `binom(n + index_offset, k + coefficient_offset)`; requires
    /// `n + index_offset >= 0` on the index domain.
    Binomial {
        /// `a` of `binom(n + a, k + b)`.
        index_offset: i64,
        /// `b` of `binom(n + a, k + b)`.
        coefficient_offset: i64,
    },
    /// `k!`.
    Factorial,
    /// `1 / k!`.
    InverseFactorial,
    /// `(n - k)!`, undefined for `k > n`.
    ComplementFactorial,
    /// `1 / (n - k)!`, zero for `k > n`.
    InverseComplementFactorial,
    /// The rising factorial `(alpha)_k = alpha (alpha + 1) ... (alpha + k - 1)`
    /// with `alpha` an exact rational string, not a nonpositive integer.
    RisingFactorial {
        /// `alpha`.
        parameter: String,
    },
    /// `1 / (alpha)_k`, with `alpha` not a nonpositive integer.
    InverseRisingFactorial {
        /// `alpha`.
        parameter: String,
    },
}

/// A Hadamard factor: the product of one or two [`HadamardFactorTerm`]s.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HadamardFactor {
    /// The terms of the product (one or two).
    pub terms: Vec<HadamardFactorTerm>,
}

/// The argument kind of a Gamma atom.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub(super) enum GammaArgument {
    /// `Gamma(n + c)`.
    Index(i64),
    /// `Gamma(k + c)`, `c` rational.
    Coefficient(Q),
    /// `Gamma(n - k + c)`.
    Complement(i64),
}

fn parse_parameter(text: &str) -> Result<Q, RecurrenceError> {
    parse_rational_coeff(text)
        .map_err(|error| RecurrenceError::from(error).with_context("rising factorial parameter"))
}

fn is_nonpositive_integer(value: &Q) -> bool {
    value.is_integer() && !value.is_positive()
}

/// `(alpha)_k`.
fn rising(alpha: &Q, k: usize) -> Q {
    (0..k).fold(Q::one(), |acc, j| acc * (alpha + q(j as i64)))
}

impl HadamardFactorTerm {
    fn text(&self) -> String {
        let shifted = |variable: &str, offset: i64| match offset {
            0 => variable.to_string(),
            o if o > 0 => format!("{variable}+{o}"),
            o => format!("{variable}{o}"),
        };
        match self {
            Self::Binomial {
                index_offset,
                coefficient_offset,
            } => format!(
                "binom({}, {})",
                shifted("n", *index_offset),
                shifted("k", *coefficient_offset)
            ),
            Self::Factorial => "k!".into(),
            Self::InverseFactorial => "1/k!".into(),
            Self::ComplementFactorial => "(n-k)!".into(),
            Self::InverseComplementFactorial => "1/(n-k)!".into(),
            Self::RisingFactorial { parameter } => format!("({parameter})_k"),
            Self::InverseRisingFactorial { parameter } => format!("1/({parameter})_k"),
        }
    }

    fn validate(&self, first_index: usize) -> Result<(), RecurrenceError> {
        match self {
            Self::Binomial { index_offset, .. } => {
                if (first_index as i64)
                    .checked_add(*index_offset)
                    .is_none_or(|a| a < 0)
                {
                    return Err(RecurrenceError::InvalidInput(format!(
                        "{} needs n + {index_offset} >= 0 for every n >= {first_index}",
                        self.text()
                    )));
                }
                Ok(())
            }
            Self::RisingFactorial { parameter } | Self::InverseRisingFactorial { parameter } => {
                if is_nonpositive_integer(&parse_parameter(parameter)?) {
                    return Err(RecurrenceError::InvalidInput(format!(
                        "the rising factorial parameter {parameter} is a nonpositive integer"
                    )));
                }
                Ok(())
            }
            _ => Ok(()),
        }
    }

    /// `u(n, k)`; `None` where the term is undefined (a pole).
    fn value(&self, n: usize, k: usize) -> Option<Q> {
        match self {
            Self::Binomial {
                index_offset,
                coefficient_offset,
            } => {
                let top = usize::try_from(n as i64 + index_offset).ok()?;
                let bottom = k as i64 + coefficient_offset;
                Some(match usize::try_from(bottom) {
                    Ok(bottom) => binomial(top, bottom),
                    Err(_) => Q::zero(),
                })
            }
            Self::Factorial => Some(factorial(k)),
            Self::InverseFactorial => Some(factorial(k).recip()),
            Self::ComplementFactorial => (k <= n).then(|| factorial(n - k)),
            Self::InverseComplementFactorial => Some(if k <= n {
                factorial(n - k).recip()
            } else {
                Q::zero()
            }),
            Self::RisingFactorial { parameter } => {
                Some(rising(&parse_parameter(parameter).ok()?, k))
            }
            Self::InverseRisingFactorial { parameter } => {
                Some(rising(&parse_parameter(parameter).ok()?, k).recip())
            }
        }
    }

    /// The Gamma atoms with their exponents (the constant normalisations,
    /// such as `1/Gamma(alpha)`, cancel in every ratio).
    fn atoms(&self) -> Result<Vec<(GammaArgument, i32)>, RecurrenceError> {
        Ok(match self {
            Self::Binomial {
                index_offset: a,
                coefficient_offset: b,
            } => vec![
                (GammaArgument::Index(a + 1), 1),
                (GammaArgument::Coefficient(q(b + 1)), -1),
                (GammaArgument::Complement(a - b + 1), -1),
            ],
            Self::Factorial => vec![(GammaArgument::Coefficient(Q::one()), 1)],
            Self::InverseFactorial => vec![(GammaArgument::Coefficient(Q::one()), -1)],
            Self::ComplementFactorial => vec![(GammaArgument::Complement(1), 1)],
            Self::InverseComplementFactorial => vec![(GammaArgument::Complement(1), -1)],
            Self::RisingFactorial { parameter } => {
                vec![(GammaArgument::Coefficient(parse_parameter(parameter)?), 1)]
            }
            Self::InverseRisingFactorial { parameter } => {
                vec![(GammaArgument::Coefficient(parse_parameter(parameter)?), -1)]
            }
        })
    }

    /// Why `k -> u(n, k)` is a nonnegative multiplier sequence for every
    /// `n >= n0`, or `None` when no theorem is recorded for this kind.
    fn multiplier_justification(&self) -> Option<String> {
        match self {
            Self::Binomial {
                coefficient_offset: b,
                ..
            } if *b >= 0 => Some(format!(
                "{}: sum_k u(n,k) x^k / k! = L_(A-{b})^({b})(-x) with A = n + a (zero when A < \
                 {b}), a Laguerre polynomial with parameter {b} > -1, whose zeros are real and \
                 positive (Szego, Orthogonal Polynomials, Theorem 3.3.1); third condition of \
                 Theorem 40",
                self.text()
            )),
            Self::Binomial {
                coefficient_offset: b,
                ..
            } => Some(format!(
                "{}: sum_k u(n,k) x^k / k! = x^{c} A! / (A+{c})! L_A^({c})(-x) with A = n + a, a \
                 power of x times a Laguerre polynomial with parameter {c} > -1, whose zeros are \
                 real and positive (Szego, Orthogonal Polynomials, Theorem 3.3.1); third \
                 condition of Theorem 40",
                self.text(),
                c = -b
            )),
            Self::InverseFactorial => Some(
                "1/k!: sum_k binom(m,k) x^k / k! = L_m(-x) is real-rooted with negative zeros \
                 for every m (Branden, Example 41); second condition of Theorem 40"
                    .into(),
            ),
            Self::InverseRisingFactorial { parameter } => {
                let alpha = parse_parameter(parameter).ok()?;
                alpha.is_positive().then(|| {
                    format!(
                        "1/({parameter})_k: sum_k binom(m,k) x^k / ({parameter})_k = m! / \
                         ({parameter})_m L_m^({parameter}-1)(-x) is real-rooted with negative \
                         zeros for every m, since the Laguerre parameter {parameter} - 1 > -1 \
                         (Szego, Orthogonal Polynomials, Theorem 3.3.1); second condition of \
                         Theorem 40"
                    )
                })
            }
            Self::InverseComplementFactorial => Some(
                "1/(n-k)!: sum_k u(n,k) x^k / k! = (1 + x)^n / n!; third condition of Theorem 40"
                    .into(),
            ),
            Self::Factorial | Self::ComplementFactorial | Self::RisingFactorial { .. } => None,
        }
    }
}

/// The transfer theorem of every factor with one.
pub(super) const POLYA_SCHUR_THEOREM: &str = "Polya--Schur characterization of multiplier \
     sequences (G. Polya and I. Schur, J. Reine Angew. Math. 144 (1914) 89--113; stated as \
     Theorem 40 of P. Branden, Unimodality, log-concavity, real-rootedness and beyond, \
     arXiv:1410.6601): a real sequence (lambda_k) is a multiplier sequence (x^k -> lambda_k x^k \
     maps real-rooted polynomials to real-rooted polynomials or zero) iff sum_k binom(m,k) \
     lambda_k x^k is real-rooted with zeros of one sign for every m, iff sum_k lambda_k x^k / k! \
     (or its reflection) is C x^m e^(a x) prod (1 + alpha_j x) with a, alpha_j >= 0; a product \
     of multiplier sequences is one (composition), and a nonzero factor lambda(n) of the row \
     changes nothing";

impl HadamardFactor {
    /// A single-term factor.
    pub fn single(term: HadamardFactorTerm) -> Self {
        Self { terms: vec![term] }
    }

    /// Plain text, for example `binom(n+1, k+1) * 1/k!`.
    pub fn text(&self) -> String {
        self.terms
            .iter()
            .map(HadamardFactorTerm::text)
            .collect::<Vec<_>>()
            .join(" * ")
    }

    /// Check the shape (one or two terms) and every term's parameters for the
    /// index domain `n >= first_index`.
    pub fn validate(&self, first_index: usize) -> Result<(), RecurrenceError> {
        if self.terms.is_empty() || self.terms.len() > 2 {
            return Err(RecurrenceError::InvalidInput(
                "a Hadamard factor has one or two terms".into(),
            ));
        }
        self.terms.iter().try_for_each(|t| t.validate(first_index))
    }

    /// `u(n, k)` for `k >= 0`; `None` where some term is undefined.
    pub fn value(&self, n: usize, k: usize) -> Option<BigRational> {
        self.terms
            .iter()
            .try_fold(Q::one(), |acc, term| Some(acc * term.value(n, k)?))
    }

    /// The Gamma atoms after cancelling equal arguments with opposite
    /// exponents; zero net exponents are dropped.
    pub(super) fn atoms(&self) -> Result<BTreeMap<GammaArgument, i32>, RecurrenceError> {
        let mut atoms = BTreeMap::new();
        for term in &self.terms {
            for (argument, exponent) in term.atoms()? {
                *atoms.entry(argument).or_insert(0) += exponent;
            }
        }
        atoms.retain(|_, exponent| *exponent != 0);
        Ok(atoms)
    }

    /// The justification that `k -> u(n, k)` is a nonnegative multiplier
    /// sequence for every `n` of the domain, or `None` (no transfer theorem).
    pub(super) fn multiplier_justification(&self) -> Option<String> {
        if let [first, second] = self.terms.as_slice() {
            // k! binom(n + a, k - c) with c >= 0: sum_k u x^k / k! = x^c (1 + x)^(n + a).
            for (a, b) in [(first, second), (second, first)] {
                if let (
                    HadamardFactorTerm::Factorial,
                    HadamardFactorTerm::Binomial {
                        coefficient_offset, ..
                    },
                ) = (a, b)
                {
                    if *coefficient_offset <= 0 {
                        return Some(format!(
                            "{}: sum_k u(n,k) x^k / k! = x^{} (1 + x)^(n+a); third condition of \
                             Theorem 40",
                            self.text(),
                            -coefficient_offset
                        ));
                    }
                }
            }
        }
        let parts = self
            .terms
            .iter()
            .map(HadamardFactorTerm::multiplier_justification)
            .collect::<Option<Vec<_>>>()?;
        Some(parts.join("; "))
    }
}

/// The term set of an operator in Laurent--Euler form: shift `i` and power
/// `x^s` (`s` may be negative), with `R` the leader shift.
pub(super) struct TermShape<'a> {
    pub(super) terms: &'a [(usize, i64)],
    pub(super) leader_shift: usize,
}

/// Check the conditions under which the conjugation identity of
/// `certify` holds at every `n >= first_index`, `k >= 0` (module
/// documentation of `hadamard_factor`).  Returns the reason when one fails.
pub(super) fn conjugation_conditions(
    atoms: &BTreeMap<GammaArgument, i32>,
    first_index: usize,
    shape: &TermShape<'_>,
) -> Result<(), String> {
    let r = shape.leader_shift as i64;
    for (argument, exponent) in atoms {
        match argument {
            GammaArgument::Index(c) => {
                if first_index as i64 + c < 1 {
                    return Err(format!(
                        "Gamma(n + {c}) is not finite and nonzero at n = {first_index}"
                    ));
                }
            }
            GammaArgument::Coefficient(c) if !c.is_integer() => {}
            GammaArgument::Coefficient(c) => {
                let c = c.to_integer();
                if *exponent < 0 {
                    if !c.is_positive() && shape.terms.iter().any(|(_, s)| *s < 0) {
                        return Err(format!(
                            "1/Gamma(k + {c}) vanishes at small k and the operator has a term \
                             with a negative power of x"
                        ));
                    }
                } else {
                    let c = i64::try_from(&c).unwrap_or(i64::MIN);
                    if c < 1 {
                        return Err(format!("Gamma(k + {c}) has a pole at k >= 0"));
                    }
                    if let Some((_, s)) = shape.terms.iter().find(|(_, s)| *s > c - 1) {
                        return Err(format!(
                            "Gamma(k + {c}) grows with k and the operator has a term x^{s} with \
                             s > {}",
                            c - 1
                        ));
                    }
                }
            }
            GammaArgument::Complement(c) => {
                if *exponent > 0 {
                    return Err(format!("Gamma(n - k + {c}) has poles for large k"));
                }
                if let Some((i, s)) = shape.terms.iter().find(|(i, s)| *i as i64 + s > r) {
                    return Err(format!(
                        "1/Gamma(n - k + {c}) vanishes above the row support and the operator \
                         has a term x^{s} S^{i} with i + s > {r}"
                    ));
                }
            }
        }
    }
    Ok(())
}

/// A product of polynomial factors in `(n, k)` (keyed by their monic form)
/// times a rational constant.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(super) struct FactorProduct {
    pub(super) constant: Q,
    pub(super) factors: BTreeMap<Biv, usize>,
}

impl FactorProduct {
    pub(super) fn one() -> Self {
        Self {
            constant: Q::one(),
            factors: BTreeMap::new(),
        }
    }

    pub(super) fn push(&mut self, factor: &Biv) {
        let (lead, monic) = factor.monic();
        self.constant *= lead;
        *self.factors.entry(monic).or_insert(0) += 1;
    }

    pub(super) fn expand(&self) -> Biv {
        let mut product = Biv::constant(self.constant.clone());
        for (factor, multiplicity) in &self.factors {
            product = product.mul(&factor.pow(*multiplicity));
        }
        product
    }
}

/// `alpha n + beta k + gamma` as a [`Biv`] in `(n, k)`.
pub(super) fn linear(alpha: i64, beta: i64, gamma: &Q) -> Biv {
    Biv::new(vec![vec![gamma.clone(), q(beta)], vec![q(alpha)]])
}

/// The ratio `u(n + i, k - s) / u(n + R, k)` of the atoms as numerator and
/// denominator products (module documentation of `hadamard_factor`).
pub(super) fn atom_ratio(
    atoms: &BTreeMap<GammaArgument, i32>,
    (i, s): (usize, i64),
    leader_shift: usize,
) -> (FactorProduct, FactorProduct) {
    let mut numerator = FactorProduct::one();
    let mut denominator = FactorProduct::one();
    let (i, r) = (i as i64, leader_shift as i64);
    for (argument, exponent) in atoms {
        // The ratio for exponent +1 as (upper factors, lower factors).
        let mut upper = Vec::new();
        let mut lower = Vec::new();
        match argument {
            GammaArgument::Index(c) => {
                // Gamma(n + i + c) / Gamma(n + R + c) = 1 / prod_(j = i..R-1) (n + j + c).
                for j in i..r {
                    lower.push(linear(1, 0, &q(j + c)));
                }
            }
            GammaArgument::Coefficient(c) => {
                // Gamma(k - s + c) / Gamma(k + c).
                if s >= 0 {
                    for j in 1..=s {
                        lower.push(linear(0, 1, &(c - q(j))));
                    }
                } else {
                    for j in 0..-s {
                        upper.push(linear(0, 1, &(c + q(j))));
                    }
                }
            }
            GammaArgument::Complement(c) => {
                // Gamma(n + i - k + s + c) / Gamma(n + R - k + c), t = R - i - s.
                let t = r - i - s;
                if t >= 0 {
                    for j in 0..t {
                        lower.push(linear(1, -1, &q(i + s + c + j)));
                    }
                } else {
                    for j in 0..-t {
                        upper.push(linear(1, -1, &q(r + c + j)));
                    }
                }
            }
        }
        for _ in 0..exponent.unsigned_abs() {
            let (num, den) = if *exponent > 0 {
                (&upper, &lower)
            } else {
                (&lower, &upper)
            };
            num.iter().for_each(|f| numerator.push(f));
            den.iter().for_each(|f| denominator.push(f));
        }
    }
    (numerator, denominator)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn binomial_term(a: i64, b: i64) -> HadamardFactorTerm {
        HadamardFactorTerm::Binomial {
            index_offset: a,
            coefficient_offset: b,
        }
    }

    /// The ratio of the atoms agrees with the ratio of exact values wherever
    /// the reference value is nonzero, and is zero where the shifted value
    /// is zero.
    #[test]
    fn atom_ratios_match_values() {
        let factors = [
            HadamardFactor::single(binomial_term(1, 1)),
            HadamardFactor::single(binomial_term(0, -1)),
            HadamardFactor::single(HadamardFactorTerm::InverseFactorial),
            HadamardFactor::single(HadamardFactorTerm::InverseComplementFactorial),
            HadamardFactor::single(HadamardFactorTerm::InverseRisingFactorial {
                parameter: "1/2".into(),
            }),
            HadamardFactor {
                terms: vec![binomial_term(0, 0), HadamardFactorTerm::Factorial],
            },
        ];
        let r = 2usize;
        for factor in &factors {
            let atoms = factor.atoms().unwrap();
            for (i, s) in [
                (0usize, 0i64),
                (0, 1),
                (0, 2),
                (1, 0),
                (1, 1),
                (2, 0),
                (1, -1),
            ] {
                let shape = TermShape {
                    terms: &[(i, s)],
                    leader_shift: r,
                };
                if conjugation_conditions(&atoms, 1, &shape).is_err() {
                    continue;
                }
                let (num, den) = atom_ratio(&atoms, (i, s), r);
                let (num, den) = (num.expand(), den.expand());
                for n in 1..7usize {
                    for k in 0..10usize {
                        if (k as i64) < s {
                            continue;
                        }
                        let reference = factor.value(n + r, k).unwrap();
                        let shifted = factor.value(n + i, (k as i64 - s) as usize).unwrap();
                        let at = |p: &Biv| {
                            let column = p.eval_n(&q(n as i64));
                            crate::vec_poly::evaluate(&column, &q(k as i64))
                        };
                        if reference.is_zero() {
                            assert!(shifted.is_zero(), "{} at {n},{k},{i},{s}", factor.text());
                        } else {
                            assert!(!at(&den).is_zero());
                            assert_eq!(
                                shifted,
                                &reference * at(&num) / at(&den),
                                "{}",
                                factor.text()
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn products_cancel_and_text() {
        let factor = HadamardFactor {
            terms: vec![binomial_term(0, 0), HadamardFactorTerm::Factorial],
        };
        let atoms = factor.atoms().unwrap();
        assert_eq!(atoms.len(), 2);
        assert_eq!(factor.text(), "binom(n, k) * k!");
        assert!(factor.multiplier_justification().is_some());
        assert!(HadamardFactor::single(HadamardFactorTerm::Factorial)
            .multiplier_justification()
            .is_none());
        assert!(HadamardFactor::single(binomial_term(-2, 0))
            .validate(1)
            .is_err());
    }
}
