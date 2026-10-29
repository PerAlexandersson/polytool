//! Transport of derivative recurrences into gamma coordinates.
//!
//! For rows symmetric at index `d_n` write `u = 1 + x`, `z = x/u^2` and
//! `P_n(x) = u^(d_n) G_n(z)` with `deg G_n <= d_n/2`.  The substitution
//! `z = x/(1 + x)^2` is not invertible as a polynomial map, so the rule is not
//! transported by composition.  Instead:
//!
//! 1. **Derivatives.**  With `z' = (1 - x)/u^3`, induction on `k` gives
//!    `D^k [u^d G(z)] = sum_i u^(d - k - 2i) T_(k,i)(x) G^(i)(z)` with
//!    `T_(0,0) = 1` and
//!    `T_(k+1,i) = (d - k - 2i) T_(k,i) + u T_(k,i)' + (1 - x) T_(k,i-1)`,
//!    a polynomial in `x` and (through `d = d_(m-r)`) in the index.
//! 2. **Clearing.**  Substituting into `f P_m = sum s c D^k P_(m-r) + g`,
//!    dividing by `u^(d_m)` (with `d_(m-r) = d_m - h r`) and multiplying by the
//!    smallest power `u^M` that clears every negative power gives a polynomial
//!    identity `E_0 G_m = sum E_(r,i) G^(i)_(m-r) + E_g` whose coefficients
//!    `E` are polynomials in `(m, x)`.  An additive term `g` is divided by the
//!    index-dependent `u^(d_m)` and is accepted only for `h = 0`.
//! 3. **Splitting.**  `Q(x)` is a quadratic extension of `Q(z)`: with
//!    `w = (1 - x)/(1 + x)`, `w^2 = 1 - 4z`.  At a split index `N` at least the
//!    `x`-degree of every coefficient, each `E` decomposes uniquely into its
//!    symmetric and antisymmetric parts at index `N`:
//!    `E = u^N A(z) + (1 - x) u^(N - 1) B(z)`, so `E/u^N = A(z) + w B(z)`
//!    (`A` are the gamma coefficients of the symmetric part, `B` those of the
//!    antisymmetric part divided by `1 - x` at index `N - 1`).  The involution
//!    `x -> 1/x` fixes `z` and negates `w`, so the identity is equivalent to
//!    the pair of relations (`even`: all `A`; `odd`: all `B`) among the
//!    `G^(i)_(m-r)`.
//!
//! The split index is chosen where the left factor `E_0` is symmetric or
//! antisymmetric (then one part often vanishes identically), and otherwise at
//! the largest `x`-degree.  Different split indices give pairs related by an
//! invertible `Q(z)`-linear change, so the rank of the pair does not depend on
//! the choice.  If the two parts are proportional (rank one), the part with a
//! nonzero leader is a scalar rule `A_0(m,z) G_m = sum ...` that determines
//! `G_m` from earlier rows; it is the transported rule.  If they are
//! independent (rank two), the source rule is not invariant under reversal
//! at index `d_n`: the odd part is a genuine companion constraint, and no
//! single scalar rule carries the source rule; [`transport_gamma`] then
//! returns both parts and the obstruction, and
//! [`RecurrenceTransformation::transport`] reports it as unsupported. The pair
//! itself is an exact description: [`GammaTransport::relation_system`] states
//! both parts as homogeneous implicit relations of `G` (the even part and its
//! odd companion, which together are equivalent to the source rule on
//! symmetric rows), and transform reports carry them as `transported.system`
//! with evidence `transported_identity`, conditional on the symmetry
//! condition. A rule with an additive term or `(-1)^n` factors has no such
//! homogeneous form and stays unsupported.
//!
//! Every transported relation is an exact consequence of the source rule at
//! targets where the rows `P_(m-r), ..., P_m` are symmetric at their
//! indices; that symmetry is the map's image condition, certified separately
//! (see the transform report's conditions).

use super::*;
use crate::vec_poly;
use num_traits::{One, Zero};
use std::collections::BTreeMap;

/// One part (even or odd) of a gamma-transported identity.
#[derive(Debug, Clone)]
pub struct GammaRelationPart {
    /// Whether every coefficient of the part vanishes.
    pub vanishes: bool,
    /// Whether the coefficient of `G_m` (the leader) vanishes identically.
    pub leader_vanishes: bool,
    /// The part as a relation `A_0(m,z) G_m = sum ...` in the [`Recurrence`]
    /// convention (variable `t` stands for `z`).  When `leader_vanishes`, its
    /// left factor is the zero polynomial: a constraint on earlier rows, not a
    /// successor rule.
    pub relation: Recurrence,
    /// Human-readable form; a vanishing leader is written `0 = ...`.
    pub text: String,
}

/// Result of [`transport_gamma`].
#[derive(Debug, Clone)]
pub struct GammaTransport {
    /// The power `M` of `1 + x` multiplied in after dividing by `(1 + x)^(d_m)`.
    pub cleared_power: usize,
    /// The split index `N`.
    pub split_index: usize,
    /// The relation formed by the symmetric parts (`A` coefficients).
    pub even: GammaRelationPart,
    /// The relation formed by the antisymmetric parts (`B` coefficients).
    pub odd: GammaRelationPart,
    /// Whether the two parts are proportional (rank one).
    pub proportional: bool,
    /// The transported scalar rule in `G`, when the parts have rank one.
    pub scalar: Option<Recurrence>,
    /// Why no scalar rule exists, when the parts are independent.
    pub obstruction: Option<String>,
}

fn one_plus_x() -> Vec<Q> {
    vec![Q::one(), Q::one()]
}

/// `sum_i gamma_i x^i (1 + x)^(n - 2i)` coefficients of a polynomial symmetric
/// at index `n` (`None` if it is not symmetric there).
fn symmetric_gamma(p: &[Q], n: usize) -> Option<Vec<Q>> {
    let p = vec_poly::trim(p);
    if vec_poly::is_zero(&p) {
        return Some(Vec::new());
    }
    if vec_poly::reverse(&p, n)? != p {
        return None;
    }
    let mut remainder = p;
    remainder.resize(n + 1, Q::zero());
    let mut gamma = vec![Q::zero(); n / 2 + 1];
    for (i, slot) in gamma.iter_mut().enumerate() {
        let c = remainder[i].clone();
        if c.is_zero() {
            continue;
        }
        let basis = vec_poly::shift_by(&vec_poly::pow(&one_plus_x(), n - 2 * i), i);
        remainder = vec_poly::sub(&remainder, &vec_poly::scale(&basis, &c));
        remainder.resize(n + 1, Q::zero());
        *slot = c;
    }
    debug_assert!(vec_poly::is_zero(&remainder));
    vec_poly::strip_trailing_zeros(&mut gamma);
    Some(gamma)
}

/// `(A, B)` with `p = u^n A(z) + (1 - x) u^(n-1) B(z)` for `deg p <= n`.
fn split_row(p: &[Q], n: usize) -> (Vec<Q>, Vec<Q>) {
    let p = vec_poly::trim(p);
    if vec_poly::is_zero(&p) {
        return (Vec::new(), Vec::new());
    }
    let reversed = vec_poly::reverse(&p, n).expect("degree at most the split index");
    let half = Q::new(BigInt::one(), BigInt::from(2));
    let symmetric = vec_poly::scale(&vec_poly::add(&p, &reversed), &half);
    let antisymmetric = vec_poly::scale(&vec_poly::sub(&p, &reversed), &half);
    let even = symmetric_gamma(&symmetric, n).expect("symmetric by construction");
    let odd = if vec_poly::is_zero(&antisymmetric) {
        Vec::new()
    } else {
        // antisymmetric = (1 - x) t with t symmetric at n - 1.
        let quotient = vec_poly::exact_div_by_t_minus(&antisymmetric, &Q::one())
            .expect("an antisymmetric polynomial vanishes at 1");
        symmetric_gamma(&vec_poly::neg(&quotient), n - 1).expect("symmetric by construction")
    };
    (even, odd)
}

fn split_bivariate(p: &Bp, n: usize) -> (Bp, Bp) {
    let (even, odd): (Vec<_>, Vec<_>) = p.iter().map(|row| split_row(row, n)).unzip();
    (bp_trim(even), bp_trim(odd))
}

/// The common index `valuation + degree` of every nonzero row, if any.
fn symmetry_center(p: &Bp) -> Option<usize> {
    let mut center = None;
    for row in p {
        let Some(degree) = vec_poly::degree(row) else {
            continue;
        };
        let valuation = row.iter().position(|c| !c.is_zero()).expect("nonzero row");
        match center {
            None => center = Some(valuation + degree),
            Some(c) if c != valuation + degree => return None,
            _ => {}
        }
    }
    center
}

fn part_from(map: &BTreeMap<PieceKey, Bp>) -> GammaRelationPart {
    let vanishes = map.values().all(bp_is_zero);
    let leader_vanishes = map.get(&PieceKey::Lhs).is_none_or(bp_is_zero);
    let relation = assemble_rule(map.clone());
    let text = if leader_vanishes {
        let mut rest = map.clone();
        rest.insert(PieceKey::Lhs, bp_constant(Q::one()));
        let text = assemble_rule(rest).to_string();
        match text.strip_prefix("P(n) = ") {
            Some(rhs) if !rhs.is_empty() => format!("0 = {rhs}"),
            _ => "0 = 0".to_string(),
        }
    } else {
        relation.to_string()
    };
    GammaRelationPart {
        vanishes,
        leader_vanishes,
        relation,
        text,
    }
}

fn proportional(even: &BTreeMap<PieceKey, Bp>, odd: &BTreeMap<PieceKey, Bp>) -> bool {
    let keys: Vec<&PieceKey> = even.keys().chain(odd.keys()).collect();
    let get =
        |map: &BTreeMap<PieceKey, Bp>, key: &PieceKey| map.get(key).cloned().unwrap_or_default();
    let Some(pivot) = keys
        .iter()
        .find(|key| !bp_is_zero(&get(even, key)) || !bp_is_zero(&get(odd, key)))
    else {
        return true;
    };
    let (pivot_even, pivot_odd) = (get(even, pivot), get(odd, pivot));
    keys.iter()
        .all(|key| bp_mul(&get(even, key), &pivot_odd) == bp_mul(&pivot_even, &get(odd, key)))
}

/// Transport a rule `f P_m = sum s c D^k P_(m-r) + g` (target index `m`) into
/// gamma coordinates `P_n = (1 + x)^(d_n) G_n(x/(1 + x)^2)` with symmetry
/// indices `d_n` from `schedule` (see the module documentation for the
/// derivation).
///
/// Both parts of the transported identity are returned; `scalar` is the
/// transported rule when they are proportional.  The result holds wherever
/// the source rule holds and the rows involved are symmetric at their
/// indices.
///
/// ```
/// use polytool::recurrence::*;
///
/// // Eulerian A008292: E_1 = 1, symmetric at n - 1.
/// let (rule, _seeds, _first) = polytool::oeis::by_id("A008292")
///     .unwrap()
///     .recurrence_export_parts()
///     .unwrap();
/// let gamma = transport_gamma(&rule, &DegreeSchedule::new(1, 0, 1)).unwrap();
/// assert!(gamma.odd.vanishes);
/// assert_eq!(
///     gamma.scalar.unwrap().to_string(),
///     "P(n) = (1 - 4t + 2nt) P(n-1) + (t - 4t^2) P'(n-1)"
/// );
/// ```
pub fn transport_gamma(
    rule: &Recurrence,
    schedule: &DegreeSchedule,
) -> Result<GammaTransport, TransformError> {
    schedule.validate()?;
    let w = schedule.weight_constant()?;
    let h = schedule.step;
    let u = one_plus_x();
    let one_minus_x = vec![Q::one(), -Q::one()];

    // Pieces E * u^e of the identity divided by u^(d_m).
    let mut pieces: Vec<(PieceKey, Bp, i64)> = Vec::new();
    let lhs = rule
        .denominator
        .as_ref()
        .map(bp_from)
        .unwrap_or_else(|| bp_constant(Q::one()));
    if bp_is_zero(&lhs) {
        return Err(TransformError::TransportUnsupported(
            "the left factor of the rule is zero".into(),
        ));
    }
    pieces.push((PieceKey::Lhs, lhs, 0));
    for term in &rule.terms {
        let coefficient = bp_from(&term.coeff);
        if bp_is_zero(&coefficient) {
            continue;
        }
        let r = i64::try_from(term.offset).map_err(|_| TransformError::Limit("offset".into()))?;
        let k = i64::try_from(term.deriv_order)
            .map_err(|_| TransformError::Limit("derivative order".into()))?;
        let hr = checked(h.checked_mul(r))?;
        // d_(m-r) = (w - h r) + h m.
        let degree = bp_affine_index(q(checked(w.checked_sub(hr))?), q(h));
        let mut t: Vec<Bp> = vec![bp_constant(Q::one())];
        for order in 0..term.deriv_order {
            let mut next: Vec<Bp> = vec![Vec::new(); t.len() + 1];
            for (i, ti) in t.iter().enumerate() {
                let shift = checked((order as i64).checked_add(2 * i as i64))?;
                let factor = bp_add(&degree, &bp_constant(q(-shift)));
                next[i] = bp_add(&next[i], &bp_mul(&factor, ti));
                next[i] = bp_add(&next[i], &bp_mul_x(&bp_dx(ti), &u));
                next[i + 1] = bp_add(&next[i + 1], &bp_mul_x(ti, &one_minus_x));
            }
            t = next;
        }
        for (i, ti) in t.into_iter().enumerate() {
            let poly = bp_mul(&coefficient, &ti);
            if bp_is_zero(&poly) {
                continue;
            }
            let exponent = checked(
                hr.checked_add(k)
                    .and_then(|v| v.checked_add(2 * i as i64))
                    .map(|v| -v),
            )?;
            pieces.push((
                PieceKey::Term {
                    offset: term.offset,
                    derivative: i,
                    alternating: term.sign == RecurrenceSign::AlternatingN,
                },
                poly,
                exponent,
            ));
        }
    }
    if let Some(g) = &rule.inhomogeneous {
        let g = bp_from(g);
        if !bp_is_zero(&g) {
            if h != 0 {
                return Err(TransformError::InhomogeneousNotTransportable);
            }
            pieces.push((PieceKey::Inhomogeneous, g, -w));
        }
    }

    let shift = pieces.iter().map(|(_, _, e)| -e).max().unwrap_or(0).max(0);
    let mut combined: BTreeMap<PieceKey, Bp> = BTreeMap::new();
    for (key, poly, exponent) in pieces {
        let power = usize::try_from(shift + exponent).expect("nonnegative by choice of shift");
        let scaled = bp_mul_x(&poly, &vec_poly::pow(&u, power));
        let entry = combined.entry(key).or_default();
        *entry = bp_add(entry, &scaled);
    }
    combined.retain(|key, poly| *key == PieceKey::Lhs || !bp_is_zero(poly));
    if combined.get(&PieceKey::Lhs).is_none_or(bp_is_zero) {
        return Err(TransformError::TransportUnsupported(
            "the left factor of the rule is zero".into(),
        ));
    }

    let max_degree = combined.values().filter_map(bp_x_degree).max().unwrap_or(0);
    let mut candidates = Vec::new();
    if let Some(center) = symmetry_center(&combined[&PieceKey::Lhs]) {
        if center >= max_degree {
            candidates.push(center);
        }
    }
    if !candidates.contains(&max_degree) {
        candidates.push(max_degree);
    }
    let split = |n: usize| {
        let mut even = BTreeMap::new();
        let mut odd = BTreeMap::new();
        for (key, poly) in &combined {
            let (a, b) = split_bivariate(poly, n);
            even.insert(key.clone(), a);
            odd.insert(key.clone(), b);
        }
        (even, odd)
    };
    let mut chosen = None;
    for &n in &candidates {
        let (even, odd) = split(n);
        let pure = even.values().all(bp_is_zero) || odd.values().all(bp_is_zero);
        if chosen.is_none() || pure {
            chosen = Some((n, even, odd));
        }
        if pure {
            break;
        }
    }
    let (split_index, even, odd) = chosen.expect("at least one split index");
    let is_proportional = proportional(&even, &odd);
    let even_part = part_from(&even);
    let odd_part = part_from(&odd);
    let (scalar, obstruction) = if is_proportional {
        let rule = if !even_part.leader_vanishes {
            even_part.relation.clone()
        } else if !odd_part.leader_vanishes {
            odd_part.relation.clone()
        } else {
            return Err(TransformError::TransportUnsupported(
                "both parts of the gamma-transported identity have a zero leader".into(),
            ));
        };
        (Some(rule), None)
    } else {
        (
            None,
            Some(format!(
                "the even and odd parts of the gamma-transported identity are independent \
                 relations (even: {}; odd: {}), so the rule is not invariant under reversal \
                 at the symmetry index and a single scalar rule in G does not carry it; the odd \
                 part is a companion constraint that would need a vector or implicit system",
                even_part.text, odd_part.text
            )),
        )
    };
    Ok(GammaTransport {
        cleared_power: usize::try_from(shift).expect("nonnegative shift"),
        split_index,
        even: even_part,
        odd: odd_part,
        proportional: is_proportional,
        scalar,
        obstruction,
    })
}

/// The rank-two gamma case as a coupled system of two relations of `G`:
/// the even part and its odd companion (see [`GammaTransport::relation_system`]).
#[derive(Debug, Clone)]
pub struct GammaRelationSystem {
    /// The even part (`A` coefficients) in implicit form at base `n = m - even_shift`.
    pub even: ImplicitRelation,
    /// The largest offset `R_e` of the even part: its base is `m - R_e`.
    pub even_shift: usize,
    /// The odd companion (`B` coefficients) in implicit form at base `n = m - odd_shift`.
    pub odd: ImplicitRelation,
    /// The largest offset `R_o` of the odd part.
    pub odd_shift: usize,
}

impl GammaTransport {
    /// The even part and the odd companion as two homogeneous implicit
    /// relations of `G`, when the parts are independent (rank two) and can
    /// be stated exactly: the source rule has no additive term and no
    /// `(-1)^n` factors.
    ///
    /// Both relations hold exactly at every target `m` where the source rule
    /// holds and the rows `P_(m-r), ..., P_m` are symmetric at their indices
    /// (the map's image condition); together they are equivalent to the
    /// source rule on symmetric rows, since `1` and `w = (1 - x)/(1 + x)` are
    /// linearly independent over `Q(z)`. An odd part whose leader vanishes
    /// is a constraint on earlier rows (for a first-order rule a row
    /// differential equation of `G_(m-1)`).
    ///
    /// ```
    /// use polytool::recurrence::*;
    ///
    /// // (1 + x)^n from P_m = (1 + 2x - m x) P_(m-1) + (x + x^2) P'_(m-1), a
    /// // rule that is not invariant under reversal at d_n = n.
    /// let json: RecurrenceJson = serde_json::from_str(
    ///     r#"{"schema": "polytool.recurrence.v1", "first_index": 0,
    ///         "initial_polynomials": [["1"]],
    ///         "recurrence": {"terms": [
    ///           {"offset": 1, "deriv_order": 0, "sign": "none",
    ///            "coeff": {"coeffs": [["1", "2"], ["0", "-1"]]}},
    ///           {"offset": 1, "deriv_order": 1, "sign": "none",
    ///            "coeff": {"coeffs": [["0", "1", "1"]]}}]}}"#,
    /// )
    /// .unwrap();
    /// let (rule, _, _) = json.to_recurrence_parts().unwrap();
    /// let gamma = transport_gamma(&rule, &DegreeSchedule::new(0, 0, 1)).unwrap();
    /// assert!(gamma.scalar.is_none());
    /// // The gamma rows are G_n = 1: a successor rule and a row equation.
    /// let system = gamma.relation_system().unwrap();
    /// assert_eq!(system.even.text("Q"), "Q[n+1] - Q[n] = 0");
    /// assert_eq!(system.odd.text("Q"), "-x D Q[n] = 0");
    /// ```
    pub fn relation_system(&self) -> Result<GammaRelationSystem, RecurrenceError> {
        if self.proportional {
            return Err(RecurrenceError::InvalidInput(
                "the parts are proportional: the scalar rule carries the source rule".into(),
            ));
        }
        let (even, even_shift) = ImplicitRelation::from_recurrence(&self.even.relation)?;
        let (odd, odd_shift) = ImplicitRelation::from_recurrence(&self.odd.relation)?;
        if even.is_zero() || odd.is_zero() {
            return Err(RecurrenceError::InvalidInput(
                "a part of the gamma-transported identity vanishes".into(),
            ));
        }
        Ok(GammaRelationSystem {
            even,
            even_shift,
            odd,
            odd_shift,
        })
    }
}

/// Whether two rules define the same row map: their left factors and every
/// term (and additive part) agree after cross multiplication by the other
/// rule's left factor.
pub(in crate::recurrence) fn rules_proportional(left: &Recurrence, right: &Recurrence) -> bool {
    let a = canonical_pieces(left);
    let b = canonical_pieces(right);
    let (Some(lhs_a), Some(lhs_b)) = (a.get(&PieceKey::Lhs), b.get(&PieceKey::Lhs)) else {
        return false;
    };
    if bp_is_zero(lhs_a) || bp_is_zero(lhs_b) {
        return false;
    }
    let keys: std::collections::BTreeSet<PieceKey> = a.keys().chain(b.keys()).cloned().collect();
    keys.iter().all(|key| {
        let pa = a.get(key).cloned().unwrap_or_default();
        let pb = b.get(key).cloned().unwrap_or_default();
        bp_mul(&pa, lhs_b) == bp_mul(&pb, lhs_a)
    })
}

#[cfg(test)]
mod tests;
