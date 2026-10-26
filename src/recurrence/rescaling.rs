//! Complete recurrence spaces after factorial coefficient rescaling.
//!
//! The opt-in field `coefficient_rescaling` of a
//! [`RecurrenceSpaceRequest`](super::RecurrenceSpaceRequest) replaces the
//! source rows `P[n]` by `Q[n] = T_m(P[n])`, where the linear map
//! `T_m(X^k) = X^k/(k!)^m` divides coefficient `k` by `(k!)^m` for a fixed
//! integer `m >= 1`. The kernel, rank and dual certificates of
//! [`find_recurrence_spaces`](super::find_recurrence_spaces) are then computed
//! on the exact rows of `Q`, and
//! [`verify_recurrence_spaces`](super::verify_recurrence_spaces) replays them
//! from the source recurrence, its seeds and origin, the power `m`, the support
//! and the residual domain recorded in the request.
//!
//! # Finite verification theorem
//!
//! Let `P` satisfy a monic pure-lag recurrence `chi(S) P = 0` of order `r`
//! whose coefficients are polynomials in `X` independent of `n`, valid for base
//! indices `n >= origin`. Consider a candidate
//! `L Q[n] = sum c_(i,j)(n,X) D^j Q[n+i]` with `deg_n c <= U`, `deg_X c <= A`
//! and `j <= J` (maxima over the whole support, leader columns included). If
//! the first `r(U + m(A+J) + 1)` consecutive residuals from a domain start
//! `>= origin` vanish, every later residual vanishes.
//!
//! Proof sketch. Multiply the residual by `X^J` and pull back by `T_m^-1`.
//! With `theta = X D`, `X^J D^j = X^(J-j) theta^(j falling)`, and `T_m^-1`
//! commutes with `theta` while `T_m^-1 X^s T_m = X^s prod_(h=1..s) (theta+h)^m`.
//! Hence a term `n^u X^a D^j S^i` becomes
//! `n^u X^s theta^(j falling) (prod_(h=1..s) (theta+h))^m S^i` with
//! `s = a + J - j >= 0`. Its derivative order is `j + m s <= m(A+J)`. As in the
//! untransformed certificate, `chi^(u+d+1)` annihilates `n^u X^v D^d P[n+i]`, so
//! `chi^H` with `H = U + m(A+J) + 1` annihilates the pulled-back residual.
//! This is a monic recurrence of order `rH`; multiplication by `X^J` and
//! `T_m^-1` are injective, so `rH` zero residuals imply all vanish. The kernel
//! of the residual rows at those `rH` base indices is therefore the complete
//! space of all-index identities for `Q` in the prescribed support. This is a
//! sufficient uniform count; sharper counts would need separate proofs.
//!
//! # Hypotheses enforced
//!
//! Rescaling requires a constant-index monic source (`source` or
//! `source_recurrence`); `indexed_source` and rows-only requests are rejected,
//! as are `automatic_variable_degree` shapes (their degree cap would depend on
//! the residual count, which itself depends on the cap) and `m = 0`. Optional
//! consistency `rows` stay in source coordinates (rows of `P`).
//!
//! Related modules: `space` (complete spaces and replay), `propagation`
//! (polynomial-index sources).

use super::{BigRational as Q, ExactRows, RecurrenceMonomial, RecurrenceSpaceLimits};
use num_bigint::BigInt;
use num_traits::{One, Zero};
use serde::{Deserialize, Serialize};

/// An exact invertible diagonal map on coefficient vectors applied to every
/// source row before the complete-space search.
///
/// The JSON form is `{"kind": "factorial_power", "power": m}`.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Hash)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum RecurrenceCoefficientRescaling {
    /// `T_m(X^k) = X^k/(k!)^m`: coefficient `k` is divided by `(k!)^m`.
    FactorialPower {
        /// The exponent `m`; it must be at least one and at most the request's
        /// `limits.max_degree`.
        power: usize,
    },
}

impl RecurrenceCoefficientRescaling {
    /// The exponent `m` of the factorial divisor `(k!)^m`.
    pub fn power(&self) -> usize {
        match self {
            Self::FactorialPower { power } => *power,
        }
    }

    /// Check the theorem's hypothesis `m >= 1` and the resource bound
    /// `m <= limits.max_degree`.
    pub fn validate(&self, limits: &RecurrenceSpaceLimits) -> Result<(), String> {
        let m = self.power();
        if m == 0 {
            return Err("factorial rescaling needs power m >= 1".into());
        }
        if m > limits.max_degree || u32::try_from(m).is_err() {
            return Err("factorial rescaling power exceeds the degree limit".into());
        }
        Ok(())
    }

    /// The exponent `H = U + m(A+J) + 1` of the finite verification theorem for
    /// a support with maximal index degree `U`, variable degree `A` and
    /// derivative order `J`. `None` on arithmetic overflow.
    pub fn annihilator_power(&self, monomials: &[RecurrenceMonomial]) -> Option<usize> {
        let max = |f: fn(&RecurrenceMonomial) -> usize| monomials.iter().map(f).max().unwrap_or(0);
        let (u, a, j) = (
            max(|c| c.index_degree),
            max(|c| c.variable_degree),
            max(|c| c.derivative_order),
        );
        a.checked_add(j)?
            .checked_mul(self.power())?
            .checked_add(u)?
            .checked_add(1)
    }

    /// Apply the map to an ascending coefficient vector: entry `k` is divided by
    /// `(k!)^m`. Fails (with a message containing `budget`) when a result
    /// would exceed `limits.max_coefficient_bits`.
    pub fn rescale_coefficients(
        &self,
        coefficients: &[Q],
        limits: &RecurrenceSpaceLimits,
    ) -> Result<Vec<Q>, String> {
        self.map(coefficients, limits, false)
    }

    /// Apply the inverse map: entry `k` is multiplied by `(k!)^m`. Fails (with a
    /// message containing `budget`) when a result would exceed
    /// `limits.max_coefficient_bits`.
    pub fn unscale_coefficients(
        &self,
        coefficients: &[Q],
        limits: &RecurrenceSpaceLimits,
    ) -> Result<Vec<Q>, String> {
        self.map(coefficients, limits, true)
    }

    fn map(
        &self,
        coefficients: &[Q],
        limits: &RecurrenceSpaceLimits,
        inverse: bool,
    ) -> Result<Vec<Q>, String> {
        self.validate(limits)?;
        let m = self.power();
        let power = u32::try_from(m).map_err(|_| "rescaling power overflow")?;
        let max_bits = limits.max_coefficient_bits;
        let mut factorial = BigInt::one();
        let mut out = Vec::with_capacity(coefficients.len());
        for (k, c) in coefficients.iter().enumerate() {
            if k > 1 {
                factorial *= BigInt::from(k);
            }
            if c.is_zero() || k < 2 {
                out.push(c.clone());
                continue;
            }
            // The reduced result has a denominator (numerator, for the
            // inverse) of at least (k!)^m / |numerator| (resp. / denominator),
            // so reject before forming a huge power.
            let other = if inverse { c.denom() } else { c.numer() };
            let lower = (factorial.bits() - 1)
                .saturating_mul(m as u64)
                .saturating_sub(other.bits());
            if lower > max_bits {
                return Err("rescaling coefficient bit budget".into());
            }
            let scale = Q::from_integer(factorial.pow(power));
            let value = if inverse { c * scale } else { c / scale };
            if value.numer().bits().max(value.denom().bits()) > max_bits {
                return Err("rescaling coefficient bit budget".into());
            }
            out.push(value);
        }
        Ok(out)
    }
}

/// Replayable data of the finite verification theorem for one rescaled shape.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceRescalingBound {
    /// The rescaling that was applied to every source row.
    pub rescaling: RecurrenceCoefficientRescaling,
    /// Order `r` of the constant-index monic source.
    pub source_order: usize,
    /// Largest power `U` of `n` over all support columns, leader included.
    pub index_degree: usize,
    /// Largest power `A` of `X` over all support columns, leader included.
    pub variable_degree: usize,
    /// Largest derivative order `J` over all support columns.
    pub derivative_order: usize,
    /// `H = U + m(A+J) + 1`; the result's `h` and `residual_count = r H`.
    pub annihilator_power: usize,
    /// Rescaled rows `Q[domain_start + j]` for `j < target_lag`, as ascending
    /// rational strings: initial data that, with a representative whose leader
    /// is nonzero at every rank of the domain, determine the later rows.
    pub initial_rows: ExactRows,
    /// Human-readable completeness statement for this shape.
    pub statement: String,
}

impl RecurrenceRescalingBound {
    pub(super) fn new(
        rescaling: RecurrenceCoefficientRescaling,
        source_order: usize,
        monomials: &[RecurrenceMonomial],
        initial_rows: ExactRows,
        domain_start: usize,
    ) -> Result<Self, String> {
        let max = |f: fn(&RecurrenceMonomial) -> usize| monomials.iter().map(f).max().unwrap_or(0);
        let annihilator_power = rescaling
            .annihilator_power(monomials)
            .ok_or("rescaled annihilator power overflow")?;
        let count = annihilator_power
            .checked_mul(source_order)
            .ok_or("rescaled residual bound overflow")?;
        let m = rescaling.power();
        let statement = format!(
            "Q[n] = T(P[n]) with T(X^k) = X^k/(k!)^{m}, where P satisfies the supplied monic \
             order-{source_order} source from its origin. Every identity in this support that \
             vanishes at the {count} base indices n = {domain_start}..{} vanishes for all \
             n >= {domain_start}, since X^J T^-1 maps each residual into a sequence annihilated \
             by chi^{annihilator_power}. The kernel is the complete space of all-index \
             identities for Q in this support.",
            domain_start.saturating_add(count).saturating_sub(1)
        );
        Ok(Self {
            rescaling,
            source_order,
            index_degree: max(|c| c.index_degree),
            variable_degree: max(|c| c.variable_degree),
            derivative_order: max(|c| c.derivative_order),
            annihilator_power,
            initial_rows,
            statement,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(n: i64, d: i64) -> Q {
        Q::new(n.into(), d.into())
    }

    #[test]
    fn factorial_maps_are_mutually_inverse() {
        let t = RecurrenceCoefficientRescaling::FactorialPower { power: 2 };
        let limits = RecurrenceSpaceLimits::default();
        let row = vec![q(1, 1), q(3, 1), q(-8, 1), q(5, 7)];
        let scaled = t.rescale_coefficients(&row, &limits).unwrap();
        assert_eq!(scaled, vec![q(1, 1), q(3, 1), q(-2, 1), q(5, 252)]);
        assert_eq!(t.unscale_coefficients(&scaled, &limits).unwrap(), row);
    }

    #[test]
    fn invalid_powers_and_bit_budget() {
        let limits = RecurrenceSpaceLimits {
            max_coefficient_bits: 16,
            ..Default::default()
        };
        let zero = RecurrenceCoefficientRescaling::FactorialPower { power: 0 };
        assert!(zero.validate(&limits).is_err());
        let t = RecurrenceCoefficientRescaling::FactorialPower { power: 3 };
        let row = vec![Q::one(); 8];
        let error = t.rescale_coefficients(&row, &limits).unwrap_err();
        assert!(error.contains("budget"));
    }

    #[test]
    fn annihilator_power_matches_the_theorem() {
        let t = RecurrenceCoefficientRescaling::FactorialPower { power: 2 };
        let cols = [
            RecurrenceMonomial {
                lag: 1,
                derivative_order: 0,
                index_degree: 1,
                variable_degree: 0,
            },
            RecurrenceMonomial {
                lag: 0,
                derivative_order: 1,
                index_degree: 0,
                variable_degree: 1,
            },
        ];
        // U = 1, A = 1, J = 1, m = 2.
        assert_eq!(t.annihilator_power(&cols), Some(6));
    }
}
