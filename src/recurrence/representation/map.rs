//! Certified normal forms of composed coordinate changes, used by the
//! representation search to deduplicate representations.
//!
//! A chain of registry maps `T_1, ..., T_s` (applied in this order) is
//! reduced to a list of segments:
//!
//! - **Möbius segments** collect consecutive maps of the Möbius class
//!   (`multiply_power`, `divide_power`, `affine`, `reversal`, `mobius`,
//!   `f_polynomial`). Each such map, and every composite of them, has the form
//!   `Q_n(x) = K_0 K_1^n prod_s (x + s)^(u_s + v_s n) P_n(phi(x))` with a
//!   Möbius transformation `phi(x) = (a x + b)/(c x + e)` (a projective
//!   matrix), finitely many monic linear factors with affine integer
//!   exponents and a rational constant `K_0 K_1^n`. Composition substitutes
//!   `phi` of the outer map into the factors of the inner map, multiplies the
//!   matrices and merges equal factors; the matrix is normalized by its first
//!   nonzero entry. This is an identity of rational functions, so equal
//!   normal forms give equal rows at every index where both chains are
//!   defined.
//! - **Diagonal segments** collect consecutive `diagonal` maps; their weights
//!   multiply. Two products are equal when the powers of `k!` agree and the
//!   rational parts agree as polynomials after cross multiplication.
//! - **Opaque segments** (`gamma`, `symmetric_decomposition`,
//!   `derivative_antecedent`, the lossy `derivative`) are equal only when their
//!   parameters are equal.
//!
//! Maps after a vector map (`symmetric_decomposition`) act componentwise,
//! with schedules shifted by the component's degree offset (`b` at
//! `d_n - 1`). The normal form then keeps the segments up to the vector map
//! and, per component, the reduced normal form of the componentwise suffix;
//! a suffix that reduces to the identity on every component is dropped, so
//! a reversal applied twice at the same bounds deduplicates on vectors too.
//!
//! Identity segments are removed. Equality of normal forms is a sufficient,
//! exact certificate of equal maps; unequal normal forms never imply unequal
//! maps, so such representations are simply kept apart. Rows are never
//! compared for deduplication. A chain with lossy maps is inverted with the
//! lost data recorded for each of them (`LostDataRecord`).

use super::*;
use crate::recurrence::implicit::diagonal_products_equal;
use num_traits::{One, Signed, Zero};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

type Q = BigRational;

/// Largest absolute exponent formed while composing Möbius-class maps.
const MAX_EXPONENT: i128 = 1 << 20;

/// Human-readable record of the normal form of a composed map.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RepresentationMapRecord {
    /// Whether the composite is the identity map (no segments left).
    pub identity: bool,
    /// One formula per segment, in application order (`P` is the source
    /// sequence, `Q` the image).
    pub segments: Vec<String>,
    /// Original index domain on which every map of the chain is defined.
    pub index_domain: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct MobiusSegment {
    /// `(a, b, c, e)`, normalized so that the first nonzero entry is one.
    matrix: [Q; 4],
    /// Monic factor `x + s` (key `s`) with exponent `u + v n` (value `(u, v)`).
    factors: BTreeMap<Q, (i128, i128)>,
    /// `(K_0, K_1)` with constant `K_0 K_1^n`.
    constant: (Q, Q),
}

#[derive(Debug, Clone)]
enum Segment {
    Mobius(Box<MobiusSegment>),
    Diagonal(Vec<DiagonalWeight>),
    Opaque(RecurrenceTransformation),
}

/// Normal form of a chain of registry maps applied to a sequence whose rows
/// start at a fixed original index.
#[derive(Debug, Clone)]
pub(super) struct MapNormalForm {
    segments: Vec<Segment>,
    /// For a chain with a vector map followed by componentwise maps: the
    /// normal form of the componentwise suffix of every component (schedules
    /// shifted by the component's degree offset); empty otherwise, and empty
    /// when every component suffix reduces to the identity.
    components: Vec<Vec<Segment>>,
    /// Names of the components of the vector map, for the record.
    component_names: Vec<String>,
    first_index: usize,
    last_index: Option<usize>,
}

/// Reduce a chain of scalar maps to segments (merging consecutive Möbius and
/// diagonal maps and dropping identity segments).
fn reduce_segments(chain: &[&RecurrenceTransformation]) -> Result<Vec<Segment>, RecurrenceError> {
    let mut segments: Vec<Segment> = Vec::new();
    for t in chain {
        let segment = if let Some(m) = MobiusSegment::from_transformation(t)? {
            Segment::Mobius(Box::new(m))
        } else if let RecurrenceTransformation::Diagonal { weight } = t {
            Segment::Diagonal(vec![weight.clone()])
        } else {
            Segment::Opaque((*t).clone())
        };
        let merged = match (segments.last(), &segment) {
            (Some(Segment::Mobius(inner)), Segment::Mobius(outer)) => {
                Some(Segment::Mobius(Box::new(inner.then(outer)?)))
            }
            (Some(Segment::Diagonal(inner)), Segment::Diagonal(outer)) => {
                let mut weights = inner.clone();
                weights.extend(outer.iter().cloned());
                Some(Segment::Diagonal(weights))
            }
            _ => None,
        };
        let segment = match merged {
            Some(merged) => {
                segments.pop();
                merged
            }
            None => segment,
        };
        let identity = match &segment {
            Segment::Mobius(m) => m.is_identity(),
            Segment::Diagonal(weights) => diagonal_identity(weights)?,
            Segment::Opaque(_) => false,
        };
        if !identity {
            segments.push(segment);
        }
    }
    Ok(segments)
}

fn segments_equal(a: &[Segment], b: &[Segment]) -> bool {
    a.len() == b.len()
        && a.iter().zip(b).all(|pair| match pair {
            (Segment::Mobius(a), Segment::Mobius(b)) => a == b,
            (Segment::Diagonal(a), Segment::Diagonal(b)) => {
                diagonal_products_equal(a, b).unwrap_or(false)
            }
            (Segment::Opaque(a), Segment::Opaque(b)) => a == b,
            _ => false,
        })
}

fn segment_text(segment: &Segment) -> String {
    match segment {
        Segment::Mobius(m) => m.text(),
        Segment::Diagonal(weights) => format!(
            "Q_n(x) = sum_k lambda(n,k) [x^k]P_n x^k with lambda = {}",
            weights
                .iter()
                .map(|w| format!("({})", w.formula()))
                .collect::<Vec<_>>()
                .join(" ")
        ),
        Segment::Opaque(t) => t.record().row_map,
    }
}

/// The degree offsets of the state entering each map of `chain`, and of the
/// final state (`chain.len() + 1` entries): `[0]` for scalar states, the
/// offsets of the vector map for vector states.
pub(super) fn chain_offsets(chain: &[&RecurrenceTransformation]) -> Vec<Vec<i64>> {
    let mut offsets = vec![vec![0]];
    for t in chain {
        let last = offsets.last().expect("nonempty").clone();
        offsets.push(if t.component_count() > 1 {
            t.component_degree_offsets()
        } else {
            last
        });
    }
    offsets
}

fn parse(text: &str, name: &str) -> Result<Q, RecurrenceError> {
    parse_rational_coeff(text).map_err(|e| RecurrenceError::from(e).with_context(name))
}

fn checked_exponent(value: Option<i128>) -> Result<i128, RecurrenceError> {
    value.filter(|v| v.abs() <= MAX_EXPONENT).ok_or_else(|| {
        RecurrenceError::LimitExceeded("map normal form exponent overflow".to_string())
    })
}

fn power(base: &Q, exponent: i128) -> Result<Q, RecurrenceError> {
    if exponent == 0 {
        return Ok(Q::one());
    }
    if base.is_zero() {
        return Err(RecurrenceError::InvalidInput(
            "zero base with a nonzero exponent in a map normal form".into(),
        ));
    }
    let e = u32::try_from(exponent.unsigned_abs())
        .map_err(|_| RecurrenceError::LimitExceeded("exponent overflow".into()))?;
    let value = Q::new(base.numer().pow(e), base.denom().pow(e));
    Ok(if exponent > 0 { value } else { value.recip() })
}

impl MobiusSegment {
    fn identity() -> Self {
        Self {
            matrix: [Q::one(), Q::zero(), Q::zero(), Q::one()],
            factors: BTreeMap::new(),
            constant: (Q::one(), Q::one()),
        }
    }

    fn is_identity(&self) -> bool {
        *self == Self::identity()
    }

    fn normalize_matrix(&mut self) {
        if let Some(lead) = self.matrix.iter().find(|c| !c.is_zero()).cloned() {
            for entry in &mut self.matrix {
                *entry /= &lead;
            }
        }
    }

    /// Multiply by `(x + shift)^(u + v n)`.
    fn add_factor(&mut self, shift: Q, exponent: (i128, i128)) -> Result<(), RecurrenceError> {
        let entry = self.factors.entry(shift.clone()).or_insert((0, 0));
        entry.0 = checked_exponent(entry.0.checked_add(exponent.0))?;
        entry.1 = checked_exponent(entry.1.checked_add(exponent.1))?;
        if *entry == (0, 0) {
            self.factors.remove(&shift);
        }
        Ok(())
    }

    /// Multiply by `base^(u + v n)`.
    fn add_constant(&mut self, base: &Q, exponent: (i128, i128)) -> Result<(), RecurrenceError> {
        self.constant.0 *= power(base, exponent.0)?;
        self.constant.1 *= power(base, exponent.1)?;
        Ok(())
    }

    /// Multiply by `(alpha x + beta)^(u + v n)` with `(alpha, beta) != 0`.
    fn add_linear(
        &mut self,
        alpha: &Q,
        beta: &Q,
        exponent: (i128, i128),
    ) -> Result<(), RecurrenceError> {
        if alpha.is_zero() {
            self.add_constant(beta, exponent)
        } else {
            self.add_constant(alpha, exponent)?;
            self.add_factor(beta / alpha, exponent)
        }
    }

    fn from_transformation(t: &RecurrenceTransformation) -> Result<Option<Self>, RecurrenceError> {
        let mut segment = Self::identity();
        let schedule_exponent = |s: &DegreeSchedule| -> Result<(i128, i128), RecurrenceError> {
            let w = (s.first_degree as i128)
                .checked_sub(
                    (s.step as i128)
                        .checked_mul(s.first_index as i128)
                        .unwrap_or(0),
                )
                .ok_or(RecurrenceError::LimitExceeded(
                    "degree schedule overflow".into(),
                ))?;
            Ok((
                checked_exponent(Some(w))?,
                checked_exponent(Some(s.step as i128))?,
            ))
        };
        let mobius = |a: Q, b: Q, c: Q, e: Q, s: &DegreeSchedule, seg: &mut Self| {
            let exponent = schedule_exponent(s)?;
            seg.add_linear(&c, &e, exponent)?;
            seg.matrix = [a, b, c, e];
            Ok::<(), RecurrenceError>(())
        };
        match t {
            RecurrenceTransformation::MultiplyPower { power } => {
                segment.add_factor(Q::zero(), (*power as i128, 0))?;
            }
            RecurrenceTransformation::DividePower { power } => {
                segment.add_factor(Q::zero(), (-(*power as i128), 0))?;
            }
            RecurrenceTransformation::Affine { a, b } => {
                segment.matrix = [parse(a, "a")?, parse(b, "b")?, Q::zero(), Q::one()];
            }
            RecurrenceTransformation::Reversal { schedule } => {
                mobius(
                    Q::zero(),
                    Q::one(),
                    Q::one(),
                    Q::zero(),
                    schedule,
                    &mut segment,
                )?;
            }
            RecurrenceTransformation::FPolynomial { schedule } => {
                mobius(
                    Q::one(),
                    Q::zero(),
                    Q::one(),
                    Q::one(),
                    schedule,
                    &mut segment,
                )?;
            }
            RecurrenceTransformation::Mobius {
                a,
                b,
                c,
                e,
                schedule,
            } => {
                mobius(
                    parse(a, "a")?,
                    parse(b, "b")?,
                    parse(c, "c")?,
                    parse(e, "e")?,
                    schedule,
                    &mut segment,
                )?;
            }
            _ => return Ok(None),
        }
        segment.normalize_matrix();
        Ok(Some(segment))
    }

    /// The composite `outer o inner` (`inner` applied first):
    /// `R = B(x) A(phi_B(x)) P(phi_A(phi_B(x)))`.
    fn then(&self, outer: &Self) -> Result<Self, RecurrenceError> {
        let [a, b, c, e] = &outer.matrix;
        let mut result = Self {
            matrix: multiply(&self.matrix, &outer.matrix),
            factors: outer.factors.clone(),
            constant: (
                &outer.constant.0 * &self.constant.0,
                &outer.constant.1 * &self.constant.1,
            ),
        };
        for (shift, &(u, v)) in &self.factors {
            // phi_B(x) + s = ((a + s c) x + (b + s e)) / (c x + e).
            let alpha = a + shift * c;
            let beta = b + shift * e;
            result.add_linear(&alpha, &beta, (u, v))?;
            result.add_linear(c, e, (-u, -v))?;
        }
        result.normalize_matrix();
        Ok(result)
    }

    fn text(&self) -> String {
        let mut parts = Vec::new();
        let (k0, k1) = &self.constant;
        if !k0.is_one() {
            parts.push(format!("({k0})"));
        }
        if !k1.is_one() {
            parts.push(format!("({k1})^n"));
        }
        for (shift, (u, v)) in &self.factors {
            let base = if shift.is_zero() {
                "x".to_string()
            } else {
                format!("({})", affine_text(&Q::one(), shift, "x"))
            };
            let exponent = affine_text(
                &Q::from_integer((*v).into()),
                &Q::from_integer((*u).into()),
                "n",
            );
            parts.push(match exponent.as_str() {
                "1" => base,
                e if e.chars().all(|c| c.is_ascii_digit()) => format!("{base}^{e}"),
                e => format!("{base}^({e})"),
            });
        }
        let [a, b, c, e] = &self.matrix;
        let argument = if c.is_zero() {
            affine_text(&(a / e), &(b / e), "x")
        } else {
            format!("({})/({})", affine_text(a, b, "x"), affine_text(c, e, "x"))
        };
        parts.push(format!("P_n({argument})"));
        format!("Q_n(x) = {}", parts.join(" "))
    }
}

/// `coefficient var + constant` with unit coefficients and signs simplified.
fn affine_text(coefficient: &Q, constant: &Q, var: &str) -> String {
    let mut text = if coefficient.is_zero() {
        String::new()
    } else if coefficient.is_one() {
        var.to_string()
    } else if *coefficient == -Q::one() {
        format!("-{var}")
    } else {
        format!("{coefficient} {var}")
    };
    if !constant.is_zero() {
        if text.is_empty() {
            text = constant.to_string();
        } else if constant.is_negative() {
            text = format!("{text} - {}", -constant);
        } else {
            text = format!("{text} + {constant}");
        }
    }
    if text.is_empty() {
        "0".to_string()
    } else {
        text
    }
}

fn multiply(left: &[Q; 4], right: &[Q; 4]) -> [Q; 4] {
    let [a, b, c, d] = left;
    let [e, f, g, h] = right;
    [a * e + b * g, a * f + b * h, c * e + d * g, c * f + d * h]
}

fn diagonal_identity(weights: &[DiagonalWeight]) -> Result<bool, RecurrenceError> {
    diagonal_products_equal(weights, &[]).map_err(RecurrenceError::from)
}

impl MapNormalForm {
    /// The normal form of `chain` (applied in order) on rows starting at
    /// `origin`.
    pub(super) fn from_chain(
        chain: &[&RecurrenceTransformation],
        origin: usize,
    ) -> Result<Self, RecurrenceError> {
        let mut first_index = origin;
        let mut last_index: Option<usize> = None;
        for t in chain {
            if let Some(schedule) = t.schedule() {
                first_index = first_index.max(schedule.first_index);
                last_index = match (last_index, schedule.last_index) {
                    (Some(a), Some(b)) => Some(a.min(b)),
                    (a, b) => a.or(b),
                };
            }
        }
        let split = chain
            .iter()
            .position(|t| t.component_count() > 1)
            .map_or(chain.len(), |p| p + 1);
        let segments = reduce_segments(&chain[..split])?;
        let mut components = Vec::new();
        let mut component_names = Vec::new();
        if split < chain.len() {
            let vector = chain[split - 1];
            component_names = vector
                .record()
                .components
                .into_iter()
                .map(|c| c.name)
                .collect();
            for offset in vector.component_degree_offsets() {
                let shifted = chain[split..]
                    .iter()
                    .map(|t| t.with_degree_offset(offset).map_err(RecurrenceError::from))
                    .collect::<Result<Vec<_>, _>>()?;
                components.push(reduce_segments(&shifted.iter().collect::<Vec<_>>())?);
            }
            if components.iter().all(Vec::is_empty) {
                components.clear();
            }
        }
        Ok(Self {
            segments,
            components,
            component_names,
            first_index,
            last_index,
        })
    }

    /// Number of segments (maps to invert when reconstructing the source).
    pub(super) fn segment_count(&self) -> usize {
        self.segments.len() + self.components.iter().map(Vec::len).max().unwrap_or(0)
    }

    /// Whether the two normal forms are certified equal maps on the same
    /// index domain. `false` means "not certified equal".
    pub(super) fn certified_equal(&self, other: &Self) -> bool {
        self.first_index == other.first_index
            && self.last_index == other.last_index
            && segments_equal(&self.segments, &other.segments)
            && self.components.len() == other.components.len()
            && self
                .components
                .iter()
                .zip(&other.components)
                .all(|(a, b)| segments_equal(a, b))
    }

    /// The readable record.
    pub(super) fn record(&self) -> RepresentationMapRecord {
        let mut segments = self.segments.iter().map(segment_text).collect::<Vec<_>>();
        for (c, component) in self.components.iter().enumerate() {
            let name = self
                .component_names
                .get(c)
                .cloned()
                .unwrap_or_else(|| c.to_string());
            segments.push(if component.is_empty() {
                format!("componentwise on {name}: the identity")
            } else {
                format!(
                    "componentwise on {name}: {}",
                    component
                        .iter()
                        .map(segment_text)
                        .collect::<Vec<_>>()
                        .join("; then ")
                )
            });
        }
        RepresentationMapRecord {
            identity: segments.is_empty(),
            segments,
            index_domain: match self.last_index {
                Some(last) => format!("{} <= n <= {last}", self.first_index),
                None => format!("n >= {}", self.first_index),
            },
        }
    }
}

/// Apply `chain` (in order) to scalar source rows `rows[i] = P_(origin + i)`,
/// returning the image states `[row][component]`. Maps after a vector map
/// act componentwise ([`RecurrenceTransformation::map_state`]).
pub(super) fn apply_chain(
    chain: &[&RecurrenceTransformation],
    origin: usize,
    rows: &[Vec<Q>],
) -> Result<Vec<Vec<Vec<Q>>>, RecurrenceError> {
    let offsets = chain_offsets(chain);
    let mut states = rows.iter().map(|r| vec![r.clone()]).collect::<Vec<_>>();
    for (step, t) in chain.iter().enumerate() {
        let mut next = Vec::with_capacity(states.len());
        for (i, state) in states.iter().enumerate() {
            next.push(
                t.map_state(origin + i, state, &offsets[step])
                    .map_err(|e| {
                        RecurrenceError::from(e)
                            .with_context(format!("map {step} ({})", t.kind_name()))
                    })?,
            );
        }
        states = next;
    }
    Ok(states)
}

/// Undo `chain` on image states, returning the reconstructed source rows.
#[cfg(test)]
pub(super) fn invert_chain(
    chain: &[&RecurrenceTransformation],
    origin: usize,
    states: &[Vec<Vec<Q>>],
) -> Result<Vec<Vec<Q>>, RecurrenceError> {
    invert_chain_with_lost_data(chain, origin, states, &[])
}

/// Undo `chain` on image states with the recorded lost data of its lossy
/// maps, returning the reconstructed source rows. A lossy map without a
/// record (or a record at a bijective map) is an error.
pub(super) fn invert_chain_with_lost_data(
    chain: &[&RecurrenceTransformation],
    origin: usize,
    states: &[Vec<Vec<Q>>],
    lost: &[LostDataRecord],
) -> Result<Vec<Vec<Q>>, RecurrenceError> {
    if lost.iter().any(|r| r.chain_position >= chain.len()) {
        return Err(RecurrenceError::InvalidInput(
            "a lost-data record lies outside the chain".into(),
        ));
    }
    let offsets = chain_offsets(chain);
    let mut current = states.to_vec();
    for (step, t) in chain.iter().enumerate().rev() {
        let record = lost.iter().find(|r| r.chain_position == step);
        let mut previous = Vec::with_capacity(current.len());
        for (i, state) in current.iter().enumerate() {
            let lost_row = record
                .map(|r| {
                    r.rows
                        .get(i)
                        .ok_or_else(|| {
                            RecurrenceError::InvalidInput(format!(
                                "lost data of map {step} misses row {i}"
                            ))
                        })
                        .and_then(|row| crate::recurrence::implicit::decode_poly(row))
                })
                .transpose()?;
            let error = |e: TransformError| {
                RecurrenceError::from(e)
                    .with_context(format!("inverse of map {step} ({})", t.kind_name()))
            };
            if offsets[step].len() > 1 {
                if lost_row.is_some() {
                    return Err(RecurrenceError::InvalidInput(format!(
                        "lost data of map {step} on a vector state"
                    )));
                }
                previous.push(
                    t.inverse_state(origin + i, state, &offsets[step])
                        .map_err(error)?,
                );
            } else {
                previous.push(vec![t
                    .inverse_row_with_lost_data(origin + i, state, lost_row.as_deref())
                    .map_err(error)?]);
            }
        }
        current = previous;
    }
    current
        .into_iter()
        .map(|mut s| {
            if s.len() == 1 {
                Ok(s.remove(0))
            } else {
                Err(RecurrenceError::InvalidInput(
                    "the inverted chain does not end in scalar rows".to_string(),
                ))
            }
        })
        .collect()
}

/// The lost data of every lossy map of `chain` applied to the scalar source
/// rows `rows[i] = P_(origin + i)`.
pub(super) fn chain_lost_data(
    chain: &[&RecurrenceTransformation],
    origin: usize,
    rows: &[Vec<Q>],
) -> Result<Vec<LostDataRecord>, RecurrenceError> {
    let mut records = Vec::new();
    for (position, t) in chain.iter().enumerate() {
        let RecurrenceTransformation::Derivative { order } = t else {
            continue;
        };
        let entering = apply_chain(&chain[..position], origin, rows)?;
        let mut lost = Vec::with_capacity(entering.len());
        for (i, state) in entering.iter().enumerate() {
            let row = state.first().filter(|_| state.len() == 1).ok_or_else(|| {
                RecurrenceError::InvalidInput(format!("map {position} needs scalar rows"))
            })?;
            let data = t
                .lost_data(origin + i, row)
                .map_err(|e| RecurrenceError::from(e).with_context(format!("map {position}")))?
                .unwrap_or_default();
            lost.push(crate::recurrence::implicit::encode_poly(&data));
        }
        records.push(LostDataRecord {
            chain_position: position,
            order: *order,
            rows: lost,
        });
    }
    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn qi(value: i128) -> Q {
        Q::from_integer(BigInt::from(value))
    }

    fn reversal(first_degree: usize) -> RecurrenceTransformation {
        RecurrenceTransformation::Reversal {
            schedule: DegreeSchedule::new(0, first_degree, 1),
        }
    }

    fn form(chain: &[RecurrenceTransformation]) -> MapNormalForm {
        MapNormalForm::from_chain(&chain.iter().collect::<Vec<_>>(), 0).unwrap()
    }

    #[test]
    fn reversal_twice_is_the_identity_only_at_the_same_bound() {
        let same = form(&[reversal(0), reversal(0)]);
        assert!(same.record().identity);
        assert!(same.certified_equal(&form(&[])));
        let shifted = form(&[reversal(0), reversal(1)]);
        assert!(!shifted.record().identity);
        assert_eq!(shifted.record().segments, vec!["Q_n(x) = x P_n(x)"]);
        assert!(
            shifted.certified_equal(&form(&[RecurrenceTransformation::MultiplyPower {
                power: 1
            }]))
        );
    }

    #[test]
    fn mobius_compositions_and_diagonal_products() {
        // f-polynomial twice: x -> x/(1+x) twice is x -> x/(1+2x).
        let f = RecurrenceTransformation::FPolynomial {
            schedule: DegreeSchedule::new(0, 0, 1),
        };
        let twice = form(&[f.clone(), f]);
        let direct = form(&[RecurrenceTransformation::Mobius {
            a: "1".into(),
            b: "0".into(),
            c: "2".into(),
            e: "1".into(),
            schedule: DegreeSchedule::new(0, 0, 1),
        }]);
        assert!(twice.certified_equal(&direct));
        let affine = |a: &str, b: &str| RecurrenceTransformation::Affine {
            a: a.into(),
            b: b.into(),
        };
        assert!(
            form(&[affine("2", "1"), affine("1/2", "-1/2")])
                .record()
                .identity
        );
        let up = RecurrenceTransformation::Diagonal {
            weight: DiagonalWeight::Binomial { m: 1 },
        };
        let down = RecurrenceTransformation::Diagonal {
            weight: DiagonalWeight::InverseBinomial { m: 1 },
        };
        assert!(form(&[up.clone(), down.clone()]).record().identity);
        assert!(!form(&[up.clone(), up]).record().identity);
        // A diagonal map between two Möbius maps blocks their merging.
        assert!(!form(&[reversal(0), down, reversal(0)]).record().identity);
    }

    #[test]
    fn chains_map_and_reconstruct_rows() {
        let rows = vec![vec![qi(1)], vec![qi(1), qi(2)], vec![qi(1), qi(4), qi(4)]];
        let chain = [reversal(0), reversal(1)];
        let refs = chain.iter().collect::<Vec<_>>();
        let states = apply_chain(&refs, 0, &rows).unwrap();
        assert_eq!(states[1][0], vec![qi(0), qi(1), qi(2)]);
        assert_eq!(invert_chain(&refs, 0, &states).unwrap(), rows);
    }
}
