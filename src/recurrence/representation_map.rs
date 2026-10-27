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
//!   `derivative_antecedent`) are equal only when their parameters are equal.
//!
//! Identity segments are removed. Equality of normal forms is a sufficient,
//! exact certificate of equal maps; unequal normal forms never imply unequal
//! maps, so such representations are simply kept apart. Rows are never
//! compared for deduplication.

use super::diagonal::diagonal_products_equal;
use super::*;
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
    first_index: usize,
    last_index: Option<usize>,
}

fn parse(text: &str, name: &str) -> Result<Q, String> {
    parse_rational_coeff(text).map_err(|e| format!("{name}: {e}"))
}

fn checked_exponent(value: Option<i128>) -> Result<i128, String> {
    value
        .filter(|v| v.abs() <= MAX_EXPONENT)
        .ok_or_else(|| "map normal form exponent overflow".to_string())
}

fn power(base: &Q, exponent: i128) -> Result<Q, String> {
    if exponent == 0 {
        return Ok(Q::one());
    }
    if base.is_zero() {
        return Err("zero base with a nonzero exponent in a map normal form".into());
    }
    let e = u32::try_from(exponent.unsigned_abs()).map_err(|_| "exponent overflow")?;
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
    fn add_factor(&mut self, shift: Q, exponent: (i128, i128)) -> Result<(), String> {
        let entry = self.factors.entry(shift.clone()).or_insert((0, 0));
        entry.0 = checked_exponent(entry.0.checked_add(exponent.0))?;
        entry.1 = checked_exponent(entry.1.checked_add(exponent.1))?;
        if *entry == (0, 0) {
            self.factors.remove(&shift);
        }
        Ok(())
    }

    /// Multiply by `base^(u + v n)`.
    fn add_constant(&mut self, base: &Q, exponent: (i128, i128)) -> Result<(), String> {
        self.constant.0 *= power(base, exponent.0)?;
        self.constant.1 *= power(base, exponent.1)?;
        Ok(())
    }

    /// Multiply by `(alpha x + beta)^(u + v n)` with `(alpha, beta) != 0`.
    fn add_linear(&mut self, alpha: &Q, beta: &Q, exponent: (i128, i128)) -> Result<(), String> {
        if alpha.is_zero() {
            self.add_constant(beta, exponent)
        } else {
            self.add_constant(alpha, exponent)?;
            self.add_factor(beta / alpha, exponent)
        }
    }

    fn from_transformation(t: &RecurrenceTransformation) -> Result<Option<Self>, String> {
        let mut segment = Self::identity();
        let schedule_exponent = |s: &DegreeSchedule| -> Result<(i128, i128), String> {
            let w = (s.first_degree as i128)
                .checked_sub(
                    (s.step as i128)
                        .checked_mul(s.first_index as i128)
                        .unwrap_or(0),
                )
                .ok_or("degree schedule overflow")?;
            Ok((
                checked_exponent(Some(w))?,
                checked_exponent(Some(s.step as i128))?,
            ))
        };
        let mobius = |a: Q, b: Q, c: Q, e: Q, s: &DegreeSchedule, seg: &mut Self| {
            let exponent = schedule_exponent(s)?;
            seg.add_linear(&c, &e, exponent)?;
            seg.matrix = [a, b, c, e];
            Ok::<(), String>(())
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
    fn then(&self, outer: &Self) -> Result<Self, String> {
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

fn diagonal_identity(weights: &[DiagonalWeight]) -> Result<bool, String> {
    diagonal_products_equal(weights, &[]).map_err(|e| e.to_string())
}

impl MapNormalForm {
    /// The normal form of `chain` (applied in order) on rows starting at
    /// `origin`.
    pub(super) fn from_chain(
        chain: &[&RecurrenceTransformation],
        origin: usize,
    ) -> Result<Self, String> {
        let mut segments: Vec<Segment> = Vec::new();
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
        Ok(Self {
            segments,
            first_index,
            last_index,
        })
    }

    /// Number of segments (maps to invert when reconstructing the source).
    pub(super) fn segment_count(&self) -> usize {
        self.segments.len()
    }

    /// Whether the two normal forms are certified equal maps on the same
    /// index domain. `false` means "not certified equal".
    pub(super) fn certified_equal(&self, other: &Self) -> bool {
        if self.first_index != other.first_index
            || self.last_index != other.last_index
            || self.segments.len() != other.segments.len()
        {
            return false;
        }
        self.segments
            .iter()
            .zip(&other.segments)
            .all(|pair| match pair {
                (Segment::Mobius(a), Segment::Mobius(b)) => a == b,
                (Segment::Diagonal(a), Segment::Diagonal(b)) => {
                    diagonal_products_equal(a, b).unwrap_or(false)
                }
                (Segment::Opaque(a), Segment::Opaque(b)) => a == b,
                _ => false,
            })
    }

    /// The readable record.
    pub(super) fn record(&self) -> RepresentationMapRecord {
        let segments = self
            .segments
            .iter()
            .map(|segment| match segment {
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
            })
            .collect::<Vec<_>>();
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
/// returning the image states `[row][component]`.
pub(super) fn apply_chain(
    chain: &[&RecurrenceTransformation],
    origin: usize,
    rows: &[Vec<Q>],
) -> Result<Vec<Vec<Vec<Q>>>, String> {
    let mut states = rows.iter().map(|r| vec![r.clone()]).collect::<Vec<_>>();
    for (step, t) in chain.iter().enumerate() {
        let mut next = Vec::with_capacity(states.len());
        for (i, state) in states.iter().enumerate() {
            if state.len() != 1 {
                return Err(format!(
                    "map {step} ({}) needs a scalar row; the state has {} components",
                    t.kind_name(),
                    state.len()
                ));
            }
            next.push(
                t.map_row(origin + i, &state[0])
                    .map_err(|e| format!("map {step} ({}): {e}", t.kind_name()))?,
            );
        }
        states = next;
    }
    Ok(states)
}

/// Undo `chain` on image states, returning the reconstructed source rows.
pub(super) fn invert_chain(
    chain: &[&RecurrenceTransformation],
    origin: usize,
    states: &[Vec<Vec<Q>>],
) -> Result<Vec<Vec<Q>>, String> {
    let mut current = states.to_vec();
    for (step, t) in chain.iter().enumerate().rev() {
        let mut previous = Vec::with_capacity(current.len());
        for (i, state) in current.iter().enumerate() {
            let row = t
                .inverse_row(origin + i, state)
                .map_err(|e| format!("inverse of map {step} ({}): {e}", t.kind_name()))?;
            previous.push(vec![row]);
        }
        current = previous;
    }
    Ok(current.into_iter().map(|mut s| s.remove(0)).collect())
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
