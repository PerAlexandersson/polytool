//! Componentwise maps of vector (multi-component) states and exact transport
//! of first-order vector systems through them.
//!
//! # Componentwise maps
//!
//! A scalar, bijective registry map `T` acts on a state `F_n = (F_(n,1), ...,
//! F_(n,c))` component by component ([`RecurrenceTransformation::map_state`]).
//! A map with a degree schedule acts on component `j` at the schedule shifted
//! by the component's degree offset ([`RecurrenceTransformation::component_degree_offsets`]
//! of the map that produced the vector): the symmetric decomposition
//! `P_n = a_n + x b_n` has offsets `(0, -1)`, so `gamma` or `reversal` with
//! the decomposition's schedule `d_n` acts on `a_n` at `d_n` and on `b_n` at
//! `d_n - 1`, their symmetry indices. Maps that produce vectors and the lossy
//! `derivative` do not act on vector states. The inverse is applied
//! componentwise ([`RecurrenceTransformation::inverse_state`]).
//!
//! # Transport of vector systems
//!
//! For a Möbius-class map (`multiply_power`, `divide_power`, `affine`,
//! `reversal`, `mobius`, `f_polynomial`) write the gauge of component `j` as
//! `G_(n,j) = lambda^(w_j + h n) F_(n,j)(phi)` with `phi = num/den`; the
//! components share `phi`, `lambda` and the step `h`, and `w_j` carries the
//! degree offset. The chain rule gives
//! `F^(k)_(n,j)(phi) = kappa^k lambda^(alpha_k) T_k[G_(n,j)]` with
//! `alpha_0 = -(w_j + h n)`, `alpha_(k+1) = alpha_k + sigma` and
//! `T_(k+1) = lambda T_k' + alpha_k lambda' T_k` (the scalar recursion of
//! module `transform::registry`). Substituting `x -> phi` in row `i` of
//! `q(n,x) F_(n+1,i) = sum_(j,k) M_(i,j,k)(n,x) D^k F_(n,j) + g_i(n,x)`
//! and multiplying by `lambda^(w_i + h (n+1))` gives
//!
//! ```text
//! q(n,phi) G_(n+1,i) = sum_(j,k) M_(i,j,k)(n,phi) kappa^k
//!     lambda^(w_i - w_j + h + sigma k) T_k[G_(n,j)] + g_i(n,phi) lambda^(w_i + h (n+1)).
//! ```
//!
//! Each `c(n,phi)` is `hom(c) den^(-deg c)` with the homogenization `hom(c)`.
//! Every row is multiplied by one common power of `lambda` that clears all
//! negative exponents, so the left factor `hom(q) lambda^(...)` stays shared;
//! a common power of `lambda` is then divided out exactly, and the first
//! nonzero coefficient of the left factor (order `n^0 x^0, n^0 x^1, ...,
//! n^1 x^0, ...`) is normalized to one. The forcing is
//! transported only when its factor is index free (`h = 0`, or a constant
//! `lambda` with `lambda^h = 1`). The result is an identity of rational
//! functions, an exact consequence of the source system wherever it holds;
//! nothing is evaluated at the pole of `phi`.

use super::*;
use crate::vec_poly;
use num_traits::{One, Zero};
use std::collections::BTreeMap;

impl RecurrenceTransformation {
    /// Degree offsets of the image components relative to the map's degree
    /// schedule: `[0, -1]` for the symmetric decomposition (`a_n` symmetric
    /// at `d_n`, `b_n` at `d_n - 1`), `[0]` for every scalar map.
    pub fn component_degree_offsets(&self) -> Vec<i64> {
        match self {
            Self::SymmetricDecomposition { .. } => vec![0, -1],
            _ => vec![0],
        }
    }

    /// Whether the map is of the Möbius class (`multiply_power`,
    /// `divide_power`, `affine`, `reversal`, `mobius`, `f_polynomial`), whose
    /// relations and vector systems are transported by a gauge.
    pub fn is_mobius_class(&self) -> bool {
        matches!(
            self,
            Self::MultiplyPower { .. }
                | Self::DividePower { .. }
                | Self::Affine { .. }
                | Self::Reversal { .. }
                | Self::Mobius { .. }
                | Self::FPolynomial { .. }
        )
    }

    /// This map with its degree schedule shifted to `d_n + offset` (maps
    /// without a schedule are returned unchanged). An error when the shifted
    /// first degree is negative or the shifted schedule is invalid.
    pub fn with_degree_offset(&self, offset: i64) -> Result<Self, TransformError> {
        let mut shifted = self.clone();
        if offset == 0 {
            return Ok(shifted);
        }
        let schedule = match &mut shifted {
            Self::Reversal { schedule }
            | Self::Mobius { schedule, .. }
            | Self::FPolynomial { schedule }
            | Self::Gamma { schedule }
            | Self::SymmetricDecomposition { schedule } => schedule,
            _ => return Ok(shifted),
        };
        let degree = i128::try_from(schedule.first_degree)
            .map_err(|_| TransformError::Limit("degree bound overflow".into()))?
            + i128::from(offset);
        schedule.first_degree = usize::try_from(degree).map_err(|_| {
            TransformError::InvalidSchedule(format!(
                "the degree bound {degree} at index {} shifted by {offset} is negative",
                schedule.first_index
            ))
        })?;
        schedule.validate()?;
        Ok(shifted)
    }

    /// Apply the map to the state `F_index`: a scalar state is mapped by
    /// [`RecurrenceTransformation::map_row`] (possibly to a vector); a vector
    /// state is mapped componentwise, component `j` by the map shifted by
    /// `offsets[j]` ([`RecurrenceTransformation::with_degree_offset`]).
    /// Vector-producing and lossy maps are rejected on vector states.
    pub fn map_state(
        &self,
        index: usize,
        state: &[Vec<Q>],
        offsets: &[i64],
    ) -> Result<Vec<Vec<Q>>, TransformError> {
        if state.len() == 1 {
            return self.map_row(index, &state[0]);
        }
        self.componentwise_check(state.len(), offsets)?;
        state
            .iter()
            .zip(offsets)
            .map(|(component, &offset)| {
                Ok(self
                    .with_degree_offset(offset)?
                    .map_row(index, component)?
                    .remove(0))
            })
            .collect()
    }

    /// Undo [`RecurrenceTransformation::map_state`]: `offsets` are those of
    /// the state the map was applied to. A vector-producing map returns the
    /// scalar row as a one-component state.
    pub fn inverse_state(
        &self,
        index: usize,
        image: &[Vec<Q>],
        offsets: &[i64],
    ) -> Result<Vec<Vec<Q>>, TransformError> {
        if self.component_count() != 1 || image.len() == 1 {
            return Ok(vec![self.inverse_row(index, image)?]);
        }
        self.componentwise_check(image.len(), offsets)?;
        image
            .iter()
            .zip(offsets)
            .map(|(component, &offset)| {
                self.with_degree_offset(offset)?
                    .inverse_row(index, std::slice::from_ref(component))
            })
            .collect()
    }

    fn componentwise_check(
        &self,
        components: usize,
        offsets: &[i64],
    ) -> Result<(), TransformError> {
        if self.component_count() != 1 || self.is_lossy() {
            return Err(TransformError::TransportUnsupported(format!(
                "{} does not act componentwise on a {components}-component state",
                self.kind_name()
            )));
        }
        if offsets.len() != components {
            return Err(TransformError::ComponentCount {
                expected: offsets.len(),
                found: components,
            });
        }
        Ok(())
    }

    /// The gauge of a Möbius-class map.
    fn gauge(&self) -> Result<Gauge, TransformError> {
        self.validate()?;
        match self {
            Self::MultiplyPower { power } | Self::DividePower { power } => {
                let power =
                    i64::try_from(*power).map_err(|_| TransformError::Limit("power".into()))?;
                Ok(Gauge::monomial(
                    if matches!(self, Self::MultiplyPower { .. }) {
                        power
                    } else {
                        -power
                    },
                ))
            }
            Self::Affine { .. } => {
                let (a, b) = self.affine_parameters()?.expect("affine");
                Ok(Gauge::mobius(a, b, q(0), q(1), 0, 0))
            }
            Self::Reversal { .. } | Self::FPolynomial { .. } | Self::Mobius { .. } => {
                let m = self.mobius_parameters()?.expect("mobius");
                Ok(Gauge::mobius(
                    m.a,
                    m.b,
                    m.c,
                    m.e,
                    m.schedule.weight_constant()?,
                    m.schedule.step,
                ))
            }
            _ => Err(TransformError::TransportUnsupported(format!(
                "{} is not a Möbius-class map",
                self.kind_name()
            ))),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
enum VectorPiece {
    Lhs,
    Term { component: usize, derivative: usize },
    Forcing,
}

/// Transport a first-order vector system
/// `q(n,x) F_(n+1) = M(n,x,D) F_n + G(n,x)` (source index `n`) through the
/// Möbius-class map `map` applied componentwise, component `j` at the
/// schedule shifted by `offsets[j]` (see the module documentation).
///
/// The result holds wherever the source system holds and the mapped rows
/// are defined; its left factor is shared by every row.
///
/// ```
/// use polytool::recurrence::*;
///
/// // (a_n, b_n) = (1, n) for P_n = 1 + n x: a_(n+1) = a_n, b_(n+1) = a_n + b_n.
/// let int = |c: i64| BigRational::from_integer(c.into());
/// let one = BivarPoly::from_coefficients_normalized(vec![vec![int(1)]]);
/// let zero = BivarPoly::from_coefficients_normalized(vec![vec![int(0)]]);
/// let system = SharedDenominatorVectorRecurrence {
///     denominator: one.clone(),
///     matrix: vec![
///         vec![WeylOperator::new(vec![one.clone()]), WeylOperator::new(vec![zero])],
///         vec![WeylOperator::new(vec![one.clone()]), WeylOperator::new(vec![one])],
///     ],
///     forcing: None,
/// };
/// // Multiply both components by x^2: the system is unchanged.
/// let map = RecurrenceTransformation::MultiplyPower { power: 2 };
/// let moved = transport_vector_componentwise(&system, &map, &[0, 0]).unwrap();
/// assert_eq!(moved, system);
/// ```
pub fn transport_vector_componentwise(
    system: &SharedDenominatorVectorRecurrence,
    map: &RecurrenceTransformation,
    offsets: &[i64],
) -> Result<SharedDenominatorVectorRecurrence, TransformError> {
    let dimension = system.dimension();
    if dimension == 0
        || offsets.len() != dimension
        || system.matrix.iter().any(|row| row.len() != dimension)
        || system
            .forcing
            .as_ref()
            .is_some_and(|f| f.len() != dimension)
    {
        return Err(TransformError::ComponentCount {
            expected: offsets.len(),
            found: dimension,
        });
    }
    let gauges = offsets
        .iter()
        .map(|&offset| map.with_degree_offset(offset)?.gauge())
        .collect::<Result<Vec<_>, _>>()?;
    let gauge = &gauges[0];
    let tau = gauge.tau;
    let degree_exponent = |degree: usize| -> Result<i64, TransformError> {
        checked(i64::try_from(degree).ok().and_then(|d| d.checked_mul(tau)))
    };
    let lambda_poly = vec_poly::trim(&gauge.lambda);
    let lambda_derivative = gauge.lambda[1].clone();
    let step = gauge.step;

    let mut rows: Vec<Vec<(VectorPiece, Bp, i64)>> = Vec::with_capacity(dimension);
    let (lhs_hom, lhs_degree) = bp_homogenize(
        &bp_from(&system.denominator),
        &gauge.numerator,
        &gauge.denominator,
    );
    if bp_is_zero(&lhs_hom) {
        return Err(TransformError::TransportUnsupported(
            "the left factor of the system is zero".into(),
        ));
    }
    let lhs_exponent = -degree_exponent(lhs_degree)?;
    for (i, row) in system.matrix.iter().enumerate() {
        let w_i = gauges[i].weight_constant;
        let mut pieces = vec![(VectorPiece::Lhs, lhs_hom.clone(), lhs_exponent)];
        for (j, operator) in row.iter().enumerate() {
            let w_j = gauges[j].weight_constant;
            for (k, coefficient) in operator.coefficients.iter().enumerate() {
                let coefficient = bp_from(coefficient);
                if bp_is_zero(&coefficient) {
                    continue;
                }
                let (hom, degree) =
                    bp_homogenize(&coefficient, &gauge.numerator, &gauge.denominator);
                let k_i64 = i64::try_from(k).map_err(|_| TransformError::Limit("order".into()))?;
                // alpha_0(n) = -(w_j + h n).
                let mut alpha = bp_affine_index(q(-w_j), q(-step));
                let mut t_operator: Vec<Bp> = vec![bp_constant(Q::one())];
                for _ in 0..k {
                    let mut next: Vec<Bp> = vec![Vec::new(); t_operator.len() + 1];
                    for (order, t) in t_operator.iter().enumerate() {
                        let from_derivative = bp_mul_x(&bp_dx(t), &lambda_poly);
                        let from_order = bp_mul_x(t, &lambda_poly);
                        let from_weight = bp_scale(&bp_mul(&alpha, t), &lambda_derivative);
                        next[order] = bp_add(&bp_add(&next[order], &from_derivative), &from_weight);
                        next[order + 1] = bp_add(&next[order + 1], &from_order);
                    }
                    t_operator = next;
                    alpha = bp_add(&alpha, &bp_constant(q(gauge.sigma)));
                }
                let scale = rational_power(&gauge.kappa, k_i64);
                let exponent = checked(
                    w_i.checked_sub(w_j)
                        .and_then(|v| v.checked_add(step))
                        .and_then(|v| v.checked_add(gauge.sigma.checked_mul(k_i64)?))
                        .and_then(|v| v.checked_sub(degree_exponent(degree).ok()?)),
                )?;
                for (order, t) in t_operator.into_iter().enumerate() {
                    let poly = bp_scale(&bp_mul(&hom, &t), &scale);
                    if !bp_is_zero(&poly) {
                        pieces.push((
                            VectorPiece::Term {
                                component: j,
                                derivative: order,
                            },
                            poly,
                            exponent,
                        ));
                    }
                }
            }
        }
        if let Some(g) = system.forcing.as_ref().map(|f| bp_from(&f[i])) {
            if !bp_is_zero(&g) {
                let index_independent = step == 0
                    || (gauge.lambda_is_constant()
                        && rational_power(&gauge.lambda[0], step).is_one());
                if !index_independent {
                    return Err(TransformError::InhomogeneousNotTransportable);
                }
                let (hom, degree) = bp_homogenize(&g, &gauge.numerator, &gauge.denominator);
                let exponent = checked(
                    w_i.checked_add(step)
                        .and_then(|v| v.checked_sub(degree_exponent(degree).ok()?)),
                )?;
                pieces.push((VectorPiece::Forcing, hom, exponent));
            }
        }
        rows.push(pieces);
    }

    // Combine: constant lambda gives scalar factors; otherwise one common
    // power of lambda clears every row, keeping the left factor shared.
    let mut combined: Vec<BTreeMap<VectorPiece, Bp>> = vec![BTreeMap::new(); dimension];
    if gauge.lambda_is_constant() {
        for (target, pieces) in combined.iter_mut().zip(rows) {
            for (key, poly, exponent) in pieces {
                let scaled = bp_scale(&poly, &rational_power(&gauge.lambda[0], exponent));
                let entry = target.entry(key).or_default();
                *entry = bp_add(entry, &scaled);
            }
        }
    } else {
        let shift = rows.iter().flatten().map(|(_, _, e)| -e).max().unwrap_or(0);
        for (target, pieces) in combined.iter_mut().zip(rows) {
            for (key, poly, exponent) in pieces {
                let power =
                    usize::try_from(shift + exponent).expect("nonnegative by choice of shift");
                let scaled = bp_mul_x(&poly, &vec_poly::pow(&lambda_poly, power));
                let entry = target.entry(key).or_default();
                *entry = bp_add(entry, &scaled);
            }
        }
        loop {
            let divided = combined
                .iter()
                .map(|row| {
                    row.iter()
                        .map(|(key, poly)| {
                            bp_div_linear(poly, &gauge.lambda).map(|p| (key.clone(), p))
                        })
                        .collect::<Option<BTreeMap<_, _>>>()
                })
                .collect::<Option<Vec<_>>>();
            match divided {
                Some(next)
                    if next[0]
                        .get(&VectorPiece::Lhs)
                        .is_some_and(|p| !bp_is_zero(p)) =>
                {
                    combined = next;
                }
                _ => break,
            }
        }
    }

    // Normalize the first nonzero coefficient of the left factor (in the
    // order n^0 x^0, n^0 x^1, ..., n^1 x^0, ...) to one, as the vector fits do.
    let lhs = combined[0]
        .get(&VectorPiece::Lhs)
        .cloned()
        .unwrap_or_default();
    if let Some(lead) = lhs.iter().flatten().find(|c| !c.is_zero()).cloned() {
        let inverse = lead.recip();
        for row in &mut combined {
            for poly in row.values_mut() {
                *poly = bp_scale(poly, &inverse);
            }
        }
    }
    let lhs = combined[0]
        .get(&VectorPiece::Lhs)
        .cloned()
        .unwrap_or_default();
    let order = combined
        .iter()
        .flat_map(|row| row.keys())
        .filter_map(|key| match key {
            VectorPiece::Term { derivative, .. } => Some(*derivative),
            _ => None,
        })
        .max()
        .unwrap_or(0);
    let zero = || bp_to_bivar(Vec::new());
    let matrix = combined
        .iter()
        .map(|row| {
            (0..dimension)
                .map(|j| {
                    WeylOperator::new(
                        (0..=order)
                            .map(|d| {
                                row.get(&VectorPiece::Term {
                                    component: j,
                                    derivative: d,
                                })
                                .cloned()
                                .map_or_else(zero, bp_to_bivar)
                            })
                            .collect(),
                    )
                })
                .collect()
        })
        .collect();
    let forcing = combined
        .iter()
        .any(|row| {
            row.get(&VectorPiece::Forcing)
                .is_some_and(|g| !bp_is_zero(g))
        })
        .then(|| {
            combined
                .iter()
                .map(|row| {
                    row.get(&VectorPiece::Forcing)
                        .cloned()
                        .map_or_else(zero, bp_to_bivar)
                })
                .collect()
        });
    Ok(SharedDenominatorVectorRecurrence {
        denominator: bp_to_bivar(lhs),
        matrix,
        forcing,
    })
}

#[cfg(test)]
mod tests;
