//! The source specification: a scalar recurrence with seed rows defining the
//! target sequence `A_n`, validated for the supported class.

use crate::recurrence::implicit::algebra::{roots_in_domain, Biv, DOp};
use crate::recurrence::*;
use crate::vec_poly;
use num_traits::{One, Zero};

type Q = BigRational;

/// A validated source: `f(m) A_m = sum c(m,x) D^d A_(m-r)` (target index
/// `m`), with `f` free of `x`, seed rows `A_(n0), ..., A_(n0 + K - 1)` and the
/// forward relation `M` (base `n`, leading shift `R`).
#[derive(Debug, Clone)]
pub(super) struct PreparedSource {
    pub(super) recurrence: Recurrence,
    pub(super) first_index: usize,
    /// The rows taken as given (`K` of them, `K >= R`).
    pub(super) seeds: Vec<Vec<Q>>,
    /// `M`, the source rule in forward form (`ImplicitRelation::from_recurrence`).
    pub(super) relation: ImplicitRelation,
    /// The recurrence order `R`.
    pub(super) order: usize,
}

impl PreparedSource {
    /// Validate the supported class and the supplied rows.
    ///
    /// Supported: homogeneous, no `(-1)^n` factor, every offset positive,
    /// derivatives allowed, and a left factor `f(n)` free of `x` (absent means
    /// one).  The first `K` supplied rows (`K` from
    /// [`Recurrence::generation_initial_count`]) are the seeds; every later
    /// supplied row must agree with the generated one, and `f(m)` must not
    /// vanish at any `m >= n0 + K`, so every later row is generated.
    pub(super) fn new(
        recurrence: &Recurrence,
        initial_rows: &[Vec<Q>],
        first_index: usize,
    ) -> Result<Self, RecurrenceError> {
        if recurrence.terms.is_empty() {
            return Err(RecurrenceError::Unsupported(
                "symmetric-companion lifts need a recurrence with at least one term".into(),
            ));
        }
        for (position, term) in recurrence.terms.iter().enumerate() {
            if term.offset == 0 {
                return Err(RecurrenceError::InvalidInput(format!(
                    "recurrence term {position} has offset 0; offsets must be positive"
                )));
            }
            if term.sign != RecurrenceSign::None {
                return Err(RecurrenceError::Unsupported(format!(
                    "recurrence term {position} has an alternating (-1)^n factor; \
                     symmetric-companion lifts support only sign-free terms"
                )));
            }
        }
        if recurrence
            .inhomogeneous
            .as_ref()
            .is_some_and(|inhomogeneous| !inhomogeneous.is_zero())
        {
            return Err(RecurrenceError::Unsupported(
                "the recurrence has an inhomogeneous term; symmetric-companion lifts support \
                 only homogeneous recurrences"
                    .into(),
            ));
        }
        let leader_in_target_index = match &recurrence.denominator {
            None => vec![Q::one()],
            Some(denominator) => {
                let leader = Biv::from_bivar(denominator);
                if !leader.is_y_free() {
                    return Err(RecurrenceError::Unsupported(
                        "the recurrence's left factor depends on x; symmetric-companion lifts \
                         support only a left factor f(n) free of x"
                            .into(),
                    ));
                }
                let leader = leader.eval_y(&Q::zero());
                if vec_poly::is_zero(&leader) {
                    return Err(RecurrenceError::InvalidInput(
                        "the recurrence's left factor is zero".into(),
                    ));
                }
                leader
            }
        };
        let order = recurrence.max_offset();
        if initial_rows.len() < order {
            return Err(RecurrenceError::InvalidInput(format!(
                "the recurrence has order {order} but only {} initial rows were supplied",
                initial_rows.len()
            )));
        }
        let seed_count = recurrence.generation_initial_count(first_index, initial_rows.len());
        let rows: Vec<Vec<Q>> = initial_rows.iter().map(|row| vec_poly::trim(row)).collect();
        let seeds = rows[..seed_count].to_vec();
        let first_generated = first_index.checked_add(seed_count).ok_or_else(|| {
            RecurrenceError::LimitExceeded("the first generated index overflows".into())
        })?;
        let singular = roots_in_domain(&leader_in_target_index, first_generated, None);
        if let Some(index) = singular.first() {
            return Err(RecurrenceError::Unsupported(format!(
                "the recurrence's left factor vanishes at index {index}, beyond the supplied \
                 rows, so the sequence is not determined there; supply rows through index \
                 {index}"
            )));
        }
        let (relation, _) = ImplicitRelation::from_recurrence(recurrence)?;
        let prepared = Self {
            recurrence: recurrence.clone(),
            first_index,
            seeds,
            relation,
            order,
        };
        let generated = prepared.rows(rows.len())?;
        for (offset, (supplied, generated)) in rows.iter().zip(&generated).enumerate() {
            if supplied != generated {
                return Err(RecurrenceError::InvalidInput(format!(
                    "supplied row {} disagrees with the row generated by the recurrence",
                    first_index + offset
                )));
            }
        }
        Ok(prepared)
    }

    /// Rows `A_(n0), ..., A_(n0 + count - 1)`, canonical (`[0]` for zero).
    pub(super) fn rows(&self, count: usize) -> Result<Vec<Vec<Q>>, RecurrenceError> {
        self.recurrence
            .generate_rows_rational(&self.seeds, self.first_index, count)
            .map(|rows| rows.iter().map(|row| vec_poly::trim(row)).collect())
            .map_err(|error| {
                RecurrenceError::VerificationFailed(format!(
                    "the source recurrence cannot be replayed: {error}"
                ))
            })
    }

    /// First base `n` from which `M` holds at every base: `A_(n+R)` is a
    /// generated row exactly when `n + R >= n0 + K`.
    pub(super) fn first_regular_base(&self) -> usize {
        self.first_index + self.seeds.len() - self.order
    }

    /// The forward leader `f(n + R)` of `M`, as a polynomial in the base `n`.
    pub(super) fn forward_leader(&self) -> Vec<Q> {
        let column = self.relation.op().column(self.order);
        column
            .0
            .get(&(0, 0))
            .map(|c| c.eval_y(&Q::zero()))
            .unwrap_or_else(|| vec![Q::zero()])
    }

    /// The source as a `polytool.recurrence.v1` record with its seeds.
    pub(super) fn to_json(&self) -> RecurrenceJson {
        RecurrenceJson::from_recurrence_rational(
            &self.recurrence,
            self.first_index,
            &self.seeds,
            None,
        )
    }

    /// `(E A)_n` for a shift-2 operator, from rows starting at `n0`.
    pub(super) fn apply(
        operator: &DOp,
        rows: &[Vec<Q>],
        origin: usize,
        base: usize,
    ) -> Option<Vec<Q>> {
        operator
            .apply(rows, origin, base)
            .map(|r| vec_poly::trim(&r))
    }
}
