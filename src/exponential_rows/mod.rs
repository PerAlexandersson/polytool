//! Exact rows `n! [z^n] exp(A(z) + t B(z))` for labelled weighted blocks.
//!
//! [`BlockWeightSpec`] describes factorial-normalized block weights, not
//! ordinary power-series coefficients. [`ExponentialFamily::rows_iter`] uses
//! the insertion recurrence; [`ExponentialFamily::row`] extracts a sparse row
//! without constructing earlier polynomial rows (a closed form for a single
//! marked size). [`check_exponential_rows`] streams into the existing exact
//! root and directed interlacing engines. Finite checks are not all-rank proofs.
//! The separate [`exponential_row_routes`](crate::exponential_rows::exponential_row_routes) registry recognizes symbolic seeds
//! and substitutions; it never infers a series identity from a finite prefix.

mod checks;
mod routes;
#[cfg(test)]
mod tests;
pub use checks::*;
pub use routes::*;

use crate::numeric::{binomial_bigint, factorial_bigint};
use crate::recurrence::parse_rational_coeff;
use crate::vec_poly;
use num_rational::BigRational as Q;
use num_traits::{One, Zero};
use serde::{Deserialize, Serialize};
use std::collections::VecDeque;

/// Invalid exponential-family specification or row request.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ExponentialRowError(pub String);
impl std::fmt::Display for ExponentialRowError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}
impl std::error::Error for ExponentialRowError {}

pub(super) fn rational(s: &str) -> Result<Q, ExponentialRowError> {
    parse_rational_coeff(s).map_err(|e| ExponentialRowError(e.to_string()))
}
pub(super) fn strings(row: &[Q]) -> Vec<String> {
    row.iter().map(ToString::to_string).collect()
}

/// A zero-constant formal series, specified by block weights `w_r = r![z^r]W`.
/// Negative rational weights are allowed for generation, but not proof routes.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum BlockWeightSpec {
    /// Finite weights `[w_1, w_2, ...]`; missing sizes have weight zero.
    Finite {
        /// Exact rational strings, starting at size one.
        weights: Vec<String>,
    },
    /// Unit weights on a set of positive sizes (duplicates are ignored).
    Sizes {
        /// Positive block sizes.
        sizes: Vec<usize>,
    },
    /// Unit weights on every size at least `minimum`.
    Tail {
        /// First allowed size, at least one.
        minimum: usize,
    },
    /// `scale * (exp(rate*z) - sum_{r<minimum}(rate*z)^r/r!)`.
    Exponential {
        /// Exact rational multiplier.
        scale: String,
        /// Exact rational exponential rate.
        rate: String,
        /// First allowed size, at least one.
        minimum: usize,
    },
    /// Exact composition `P(seed)-P(0)` with ordinary rational polynomial coefficients.
    PolynomialSubstitution {
        /// Inner zero-constant formal series.
        seed: Box<BlockWeightSpec>,
        /// Ordinary coefficients of P, starting at degree zero.
        coefficients: Vec<String>,
    },
    /// Exact composition with `(1-d*x)Phi'=a+b*x+lambda*Phi`, `Phi(0)=0`.
    RationalOdeSubstitution {
        /// Inner zero-constant series.
        seed: Box<BlockWeightSpec>,
        /// Constant derivative term.
        a: String,
        /// Linear derivative term.
        b: String,
        /// Constant feedback coefficient.
        lambda: String,
        /// Linear derivative-factor coefficient.
        d: String,
    },
    /// Exact composition `Phi(seed)`, with `Phi(0)=0`, `Phi'=a+b*x+lambda*Phi`.
    OdeSubstitution {
        /// Inner zero-constant series.
        seed: Box<BlockWeightSpec>,
        /// Constant derivative term.
        a: String,
        /// Linear derivative term.
        b: String,
        /// Feedback term.
        lambda: String,
    },
}

impl BlockWeightSpec {
    /// Expand block weights through `max_n` (index zero is always zero).
    pub fn block_weights(&self, max_n: usize) -> Result<Vec<Q>, ExponentialRowError> {
        self.validate(0)?;
        self.expand(max_n)
    }
    pub(super) fn validate(&self, depth: usize) -> Result<(), ExponentialRowError> {
        if depth > 16 {
            return Err(ExponentialRowError("substitution depth exceeds 16".into()));
        }
        match self {
            Self::Finite { weights } => {
                for w in weights {
                    rational(w)?;
                }
            }
            Self::Sizes { sizes } if sizes.contains(&0) => {
                return Err(ExponentialRowError("block sizes must be positive".into()))
            }
            Self::Tail { minimum } | Self::Exponential { minimum, .. } if *minimum == 0 => {
                return Err(ExponentialRowError(
                    "minimum block size must be positive".into(),
                ))
            }
            Self::Exponential { scale, rate, .. } => {
                rational(scale)?;
                rational(rate)?;
            }
            Self::PolynomialSubstitution { seed, coefficients } => {
                seed.validate(depth + 1)?;
                for c in coefficients {
                    rational(c)?;
                }
            }
            Self::OdeSubstitution { seed, a, b, lambda }
            | Self::RationalOdeSubstitution {
                seed, a, b, lambda, ..
            } => {
                seed.validate(depth + 1)?;
                rational(a)?;
                rational(b)?;
                rational(lambda)?;
                if let Self::RationalOdeSubstitution { d, .. } = self {
                    rational(d)?;
                }
            }
            _ => {}
        }
        Ok(())
    }
    fn expand(&self, n: usize) -> Result<Vec<Q>, ExponentialRowError> {
        if n == usize::MAX {
            return Err(ExponentialRowError(
                "row index must be below usize::MAX".into(),
            ));
        }
        let mut out = vec![Q::zero(); n + 1];
        match self {
            Self::Finite { weights } => {
                for (i, w) in weights.iter().take(n).enumerate() {
                    out[i + 1] = rational(w)?;
                }
            }
            Self::Sizes { sizes } => {
                for &r in sizes.iter().filter(|&&r| r <= n) {
                    out[r] = Q::one();
                }
            }
            Self::Tail { minimum } => {
                for w in out.iter_mut().skip(*minimum) {
                    *w = Q::one();
                }
            }
            Self::Exponential {
                scale,
                rate,
                minimum,
            } => {
                let (scale, rate) = (rational(scale)?, rational(rate)?);
                let mut power = Q::one();
                for (r, w) in out.iter_mut().enumerate().skip(1) {
                    power *= &rate;
                    if r >= *minimum {
                        *w = &scale * &power;
                    }
                }
            }
            Self::PolynomialSubstitution { seed, coefficients } => {
                let mut inner = seed.expand(n)?;
                for (r, w) in inner.iter_mut().enumerate().skip(1) {
                    *w /= Q::from_integer(factorial_bigint(r));
                }
                let mut power = vec![Q::one()];
                for coefficient in coefficients.iter().skip(1).take(n) {
                    power = truncated_product(&power, &inner, n);
                    let coefficient = rational(coefficient)?;
                    for (r, c) in power.iter().enumerate() {
                        out[r] += c * &coefficient;
                    }
                }
                for (r, w) in out.iter_mut().enumerate().skip(1) {
                    *w *= Q::from_integer(factorial_bigint(r));
                }
            }
            Self::OdeSubstitution { seed, a, b, lambda }
            | Self::RationalOdeSubstitution {
                seed, a, b, lambda, ..
            } => {
                let mut inner = seed.expand(n)?;
                for (r, w) in inner.iter_mut().enumerate().skip(1) {
                    *w /= Q::from_integer(factorial_bigint(r));
                }
                let (a, b, l) = (rational(a)?, rational(b)?, rational(lambda)?);
                let d = match self {
                    Self::RationalOdeSubstitution { d, .. } => rational(d)?,
                    _ => Q::zero(),
                };
                // phi_i: i! [x^i] Phi. Composition uses ordinary coefficients.
                let mut power = vec![Q::one()];
                let mut phi = a;
                let mut factorial = Q::one();
                for i in 1..=n {
                    power = truncated_product(&power, &inner, n);
                    factorial *= Q::from_integer(i.into());
                    for (r, c) in power.iter().enumerate() {
                        out[r] += c * &phi / &factorial;
                    }
                    phi = (&l + &d * Q::from_integer(i.into())) * phi
                        + if i == 1 { b.clone() } else { Q::zero() };
                }
                for (r, w) in out.iter_mut().enumerate().skip(1) {
                    *w *= Q::from_integer(factorial_bigint(r));
                }
            }
        }
        Ok(out)
    }
    fn finite_support(&self) -> Option<usize> {
        match self {
            Self::Finite { weights } => Some(weights.len()),
            Self::Sizes { sizes } => Some(sizes.iter().copied().max().unwrap_or(0)),
            _ => None,
        }
    }
}

fn truncated_product(a: &[Q], b: &[Q], n: usize) -> Vec<Q> {
    let mut out = vec![Q::zero(); (a.len() + b.len() - 1).min(n + 1)];
    for (i, x) in a.iter().enumerate().filter(|(_, x)| !x.is_zero()) {
        for (j, y) in b
            .iter()
            .take(out.len().saturating_sub(i))
            .enumerate()
            .filter(|(_, y)| !y.is_zero())
        {
            out[i + j] += x * y;
        }
    }
    out = vec_poly::trim(&out);
    out
}

/// Exact labelled-block exponential family; `A(0)=B(0)=0` is enforced by shape.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ExponentialFamily {
    /// Factorial-normalized weights of unmarked blocks.
    pub unmarked: BlockWeightSpec,
    /// Factorial-normalized weights of marked blocks.
    pub marked: BlockWeightSpec,
}
impl ExponentialFamily {
    /// Insertion recurrence through inclusive rank `max_n`, starting at `F_0=1`.
    /// Finite block supports retain only the required history window.
    pub fn rows_iter(&self, max_n: usize) -> Result<ExponentialRowIterator, ExponentialRowError> {
        let a = self.unmarked.block_weights(max_n)?;
        let b = self.marked.block_weights(max_n)?;
        let window = match (self.unmarked.finite_support(), self.marked.finite_support()) {
            (Some(a), Some(b)) => a.max(b).max(1),
            _ => max_n.saturating_add(1),
        };
        Ok(ExponentialRowIterator {
            a,
            b,
            history: VecDeque::new(),
            next: 0,
            max_n,
            window,
        })
    }
    /// Dense rows, indexed by rank from zero through `max_n` inclusive.
    pub fn rows(&self, max_n: usize) -> Result<Vec<Vec<Q>>, ExponentialRowError> {
        Ok(self.rows_iter(max_n)?.collect())
    }
    /// Extract one row without constructing preceding polynomial rows.
    /// For one marked size `k` with weight `beta`, use
    /// `n! beta^m b_(n-km) / ((n-km)! m! (k!)^m)`.
    /// Otherwise extract coefficients of `exp(A) B^m / m!` by truncated products.
    pub fn row(&self, n: usize) -> Result<Vec<Q>, ExponentialRowError> {
        self.selected_rows(&[n]).map(|mut rows| rows.remove(0))
    }
    /// Sparse rows in caller order; share the scalar `exp(A)` coefficient table.
    pub fn selected_rows(&self, indices: &[usize]) -> Result<Vec<Vec<Q>>, ExponentialRowError> {
        self.unmarked.validate(0)?;
        self.marked.validate(0)?;
        let n = indices.iter().copied().max().unwrap_or(0);
        let a = self.unmarked.expand(n)?;
        let b = self.marked.expand(n)?;
        let mut scalar = vec![Q::zero(); n + 1];
        scalar[0] = Q::one();
        let support: Vec<_> = (1..=n).filter(|&r| !a[r].is_zero()).collect();
        for j in 1..=n {
            let mut sum = Q::zero();
            for &r in support.iter().take_while(|&&r| r <= j) {
                sum += Q::from_integer(binomial_bigint(j - 1, r - 1)) * &a[r] * &scalar[j - r];
            }
            scalar[j] = sum;
        }
        let marked: Vec<_> = (1..=n).filter(|&r| !b[r].is_zero()).collect();
        if marked.len() <= 1 {
            return Ok(indices
                .iter()
                .map(|&j| {
                    let Some(&k) = marked.first() else {
                        return vec![scalar[j].clone()];
                    };
                    let mut row = vec![scalar[j].clone()];
                    let mut factor = Q::one();
                    for m in 1..=j / k {
                        factor *= Q::new(binomial_bigint(j - k * (m - 1), k), m.into()) * &b[k];
                        row.push(&factor * &scalar[j - k * m]);
                    }
                    // k! already occurs in each binomial; no ordinary/factorial ambiguity.
                    row = vec_poly::trim(&row);
                    row
                })
                .collect());
        }
        let ordinary_b: Vec<_> = b
            .iter()
            .enumerate()
            .map(|(r, w)| w / Q::from_integer(factorial_bigint(r)))
            .collect();
        let ordinary_scalar: Vec<_> = scalar
            .iter()
            .enumerate()
            .map(|(r, w)| w / Q::from_integer(factorial_bigint(r)))
            .collect();
        let mut power = vec![Q::one()];
        let mut factorial = Q::one();
        let mut rows: Vec<Vec<Q>> = indices.iter().map(|_| Vec::new()).collect();
        let minimum = marked[0];
        for m in 0..=n / minimum {
            for (row, &j) in rows.iter_mut().zip(indices) {
                if m <= j / minimum {
                    let value = (0..=j.min(power.len() - 1))
                        .map(|r| &power[r] * &ordinary_scalar[j - r])
                        .sum::<Q>();
                    row.push(value * Q::from_integer(factorial_bigint(j)) / &factorial);
                }
            }
            power = truncated_product(&power, &ordinary_b, n);
            factorial *= Q::from_integer((m + 1).into());
        }
        for row in &mut rows {
            *row = vec_poly::trim(row);
        }
        Ok(rows)
    }
}

/// Bounded insertion-row iterator, retaining only the block-support window.
pub struct ExponentialRowIterator {
    a: Vec<Q>,
    b: Vec<Q>,
    history: VecDeque<Vec<Q>>,
    next: usize,
    max_n: usize,
    window: usize,
}
impl Iterator for ExponentialRowIterator {
    type Item = Vec<Q>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.next > self.max_n {
            return None;
        }
        let n = self.next;
        let mut row = vec![Q::one()];
        if n > 0 {
            row = vec![Q::zero(); n + 1];
            for r in 1..=n.min(self.history.len()) {
                if self.a[r].is_zero() && self.b[r].is_zero() {
                    continue;
                }
                let c = Q::from_integer(binomial_bigint(n - 1, r - 1));
                let previous = &self.history[self.history.len() - r];
                for (j, w) in previous.iter().enumerate() {
                    row[j] += &c * &self.a[r] * w;
                    row[j + 1] += &c * &self.b[r] * w;
                }
            }
            row = vec_poly::trim(&row);
        }
        self.history.push_back(row.clone());
        if self.history.len() > self.window {
            self.history.pop_front();
        }
        self.next += 1;
        Some(row)
    }
    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.max_n.saturating_add(1).saturating_sub(self.next);
        (remaining, Some(remaining))
    }
}
impl ExactSizeIterator for ExponentialRowIterator {}
impl std::iter::FusedIterator for ExponentialRowIterator {}
