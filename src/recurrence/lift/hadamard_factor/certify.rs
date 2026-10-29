//! The certificate of a Hadamard-factor lift: conjugation of the input
//! relation by the factor, pseudo-division by the parent relation, the row
//! annihilator of a product parent, an invariant row annihilator of a
//! first-order parent or an invariant row-annihilator system of a parent of
//! any order, and [`HadamardFactorCertificate`] with its independent replay.

use super::annihilator::{invariance_holds, step_operator, vanishing_bases, InvariantAnnihilator};
use super::annihilator_system::{system_invariance_holds, InvariantAnnihilatorSystem};
use super::factor::{atom_ratio, conjugation_conditions, FactorProduct, HadamardFactor, TermShape};
use super::prepare::{left_multiply, shift_free_from_json, shift_free_json, PreparedRule};
use crate::recurrence::implicit::algebra::{
    content_and_primitive, from_theta, roots_in_domain, theta_add, to_theta, Biv, DOp, ThetaOp,
};
use crate::recurrence::implicit::{biv_from_json, biv_json};
use crate::recurrence::lift::section_certify::{
    n_coefficients, n_polynomial, pseudo_divide_row, x_free_leader, ForwardDivision, RatioGauge,
};
use crate::recurrence::lift::section_operators::{check_term_budget, OperatorMatrix};
use crate::recurrence::{
    decode_poly, encode_poly, BivarPolyJson, ImplicitRelation, ImplicitRelationJson,
    OperatorTermJson, RecurrenceError, RecurrenceJson, VeroneseSectionGauge,
};
use crate::vec_poly;
use crate::BigRational;
use num_traits::{One, Zero};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

type Q = BigRational;

/// Largest number of operator terms in a conjugated relation, a division or
/// a certificate.
pub const MAX_HADAMARD_OPERATOR_TERMS: usize = 20_000;

/// An operator in Laurent--Euler form `sum x^s m_(i,s)(n, theta) S^i`,
/// keyed by `(i, s)` with `s` possibly negative.
pub(super) type LaurentOp = BTreeMap<(usize, i64), Biv>;

/// The Laurent--Euler form of a normal-ordered operator, using
/// `D^d = x^(-d) theta (theta - 1) ... (theta - d + 1)`.
pub(super) fn laurent_euler(op: &DOp) -> LaurentOp {
    let (j, theta) = to_theta(op);
    theta
        .into_iter()
        .map(|((i, s), c)| ((i, s as i64 - j as i64), c))
        .collect()
}

/// `M''`: the input relation conjugated by the factor (module documentation
/// of `hadamard_factor`), multiplied by `x^J''` and normal ordered.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct ConjugatedRelation {
    pub(super) op: DOp,
    /// `J''`, the power of `x` that makes every term polynomial.
    pub(super) x_shift: usize,
    /// `D(n, k)`, the least common multiple of the ratio denominators.
    pub(super) denominator: Biv,
}

/// Conjugate `input` (leader shift `R = input.max_shift()`) by the factor
/// `lambda(n) u(n, k)`: the term `x^s m(n, theta) S^i` becomes
/// `x^(s + J'') e(n, theta + s) m(n, theta) S^i` with
/// `e = D(n, k) u(n + i, k - s) / u(n + R, k)` and `D` the least common
/// multiple of the denominators of these ratios.
pub(super) fn conjugate_input(
    input: &DOp,
    factor: &HadamardFactor,
    gauge: &RatioGauge,
    first_index: usize,
) -> Result<ConjugatedRelation, RecurrenceError> {
    factor.validate(first_index)?;
    let r = input.max_shift();
    let laurent = laurent_euler(input);
    let terms: Vec<(usize, i64)> = laurent.keys().copied().collect();
    let atoms = factor.atoms()?;
    conjugation_conditions(
        &atoms,
        first_index,
        &TermShape {
            terms: &terms,
            leader_shift: r,
        },
    )
    .map_err(RecurrenceError::Unsupported)?;
    let numerator_gauge = n_polynomial(&gauge.numerator);
    let denominator_gauge = n_polynomial(&gauge.denominator);
    let mut ratios = Vec::with_capacity(terms.len());
    let mut lcm: BTreeMap<Biv, usize> = BTreeMap::new();
    for &(i, s) in &terms {
        let (mut numerator, mut denominator) = atom_ratio(&atoms, (i, s), r);
        // lambda(n + i) / lambda(n + R) = prod_(j = i+1..R) b(n + j) / a(n + j).
        for j in i + 1..=r {
            numerator.push(&denominator_gauge.shift_n(j as i64));
            denominator.push(&numerator_gauge.shift_n(j as i64));
        }
        for (factor, multiplicity) in &denominator.factors {
            if factor.as_constant().is_some() {
                continue;
            }
            let entry = lcm.entry(factor.clone()).or_insert(0);
            *entry = (*entry).max(*multiplicity);
        }
        ratios.push((numerator, denominator));
    }
    let shift = terms.iter().map(|(_, s)| -s).max().unwrap_or(0).max(0);
    let mut theta = ThetaOp::new();
    for ((&(i, s), (numerator, denominator)), m) in terms.iter().zip(&ratios).zip(laurent.values())
    {
        let mut e = FactorProduct {
            constant: &numerator.constant / &denominator.constant,
            factors: numerator.factors.clone(),
        };
        for (factor, multiplicity) in &lcm {
            let used = denominator.factors.get(factor).copied().unwrap_or(0);
            *e.factors.entry(factor.clone()).or_insert(0) += multiplicity - used;
        }
        let coefficient = e
            .expand()
            .shift(&Q::zero(), &Q::from_integer(s.into()))
            .mul(m);
        theta_add(&mut theta, (i, (s + shift) as usize), coefficient);
    }
    let denominator = FactorProduct {
        constant: Q::one(),
        factors: lcm,
    }
    .expand();
    Ok(ConjugatedRelation {
        op: from_theta(&theta),
        x_shift: shift as usize,
        denominator,
    })
}

/// A parent `g(n) P_(n+1) = a(n) c(x) P_n` with one seed `f = P_(p0)`: then
/// `P_n = Lambda(n) c^(n - p0) f`, and every row is annihilated by
/// `A(n) = c f D - ((n - p0) c' f + c f')`.
pub(super) struct ProductParent {
    c: Vec<Q>,
    seed: Vec<Q>,
    first_index: usize,
}

impl ProductParent {
    /// Recognise the shape; `None` otherwise.
    pub(super) fn recognise(parent: &PreparedRule) -> Option<Self> {
        if parent.leader_shift != 1 || parent.seeds.len() != 1 {
            return None;
        }
        x_free_leader(&parent.relation).ok()?;
        let column = parent.relation.column(0);
        if column.0.len() != 1 {
            return None;
        }
        let coefficient = column.0.get(&(0, 0))?;
        let (_, primitive) = content_and_primitive(coefficient);
        if !primitive.is_n_free() {
            return None;
        }
        let c = vec_poly::trim(&primitive.rows().first().cloned().unwrap_or_default());
        let seed = vec_poly::trim(&parent.seeds[0]);
        if c.iter().all(Zero::is_zero) || seed.iter().all(Zero::is_zero) {
            return None;
        }
        Some(Self {
            c,
            seed,
            first_index: parent.first_index,
        })
    }

    /// `A(n)`.
    pub(super) fn annihilator(&self) -> DOp {
        let c = &self.c;
        let f = &self.seed;
        let leading = vec_poly::mul(c, f);
        let c_prime_f = vec_poly::mul(&vec_poly::derivative(c), f);
        let c_f_prime = vec_poly::mul(c, &vec_poly::derivative(f));
        let index = n_polynomial(&[-Q::from_integer(self.first_index.into()), Q::one()]);
        let constant = index
            .mul(&Biv::from_y(&c_prime_f))
            .add(&Biv::from_y(&c_f_prime))
            .scale(&-Q::one());
        let mut op = DOp::single(0, 1, Biv::from_y(&leading));
        op.add_term((0, 0), constant);
        op
    }

    /// Check `A(n + 1) o c = c o A(n)` and `A(p0) P_(p0) = 0` exactly.
    pub(super) fn annihilator_is_invariant(&self, annihilator: &DOp) -> bool {
        let mut next = DOp::default();
        for (key, coefficient) in &annihilator.0 {
            next.add_term(*key, coefficient.shift_n(1));
        }
        let c = DOp::single(0, 0, Biv::from_y(&self.c));
        next.compose(&c) == c.compose(annihilator)
            && annihilator
                .apply_row(&self.seed, self.first_index)
                .iter()
                .all(Zero::is_zero)
    }
}

/// `mu f = Z o g + rem` in the Weyl algebra (shift-free operators), with
/// `mu` a power of the leading coefficient of `g` and `ord rem < ord g`.
pub(super) fn right_pseudo_divide(
    f: &DOp,
    g: &DOp,
    max_terms: usize,
) -> Result<(Biv, DOp, DOp), RecurrenceError> {
    let order = g.max_derivative();
    let lead =
        g.0.get(&(0, order))
            .cloned()
            .ok_or_else(|| RecurrenceError::VerificationFailed("the divisor is zero".into()))?;
    let mut remainder = f.clone();
    let mut quotient = DOp::default();
    let mut multiplier = Biv::one();
    while !remainder.is_zero() && remainder.max_derivative() >= order {
        let top = remainder.max_derivative();
        let coefficient = remainder.0.get(&(0, top)).cloned().ok_or_else(|| {
            RecurrenceError::VerificationFailed("a shift-free operator was expected".into())
        })?;
        let term = DOp::single(0, top - order, coefficient);
        remainder = left_multiply(&lead, &remainder).sub(&term.compose(g));
        quotient = left_multiply(&lead, &quotient).add(&term);
        multiplier = multiplier.mul(&lead);
        if remainder.0.contains_key(&(0, top)) {
            return Err(RecurrenceError::VerificationFailed(
                "right division did not cancel the leading derivative".into(),
            ));
        }
        check_term_budget(remainder.0.len() + quotient.0.len(), max_terms)?;
    }
    Ok((multiplier, quotient, remainder))
}

/// `w(n) M'' = sum_nu H_nu S^nu L + R` for a scalar parent relation `L`.
pub(super) fn divide_by_parent(
    op: &DOp,
    parent: &DOp,
    max_terms: usize,
) -> Result<ForwardDivision, RecurrenceError> {
    let (k, leader) = x_free_leader(parent)?;
    let system = OperatorMatrix {
        entries: vec![vec![parent.clone()]],
    };
    pseudo_divide_row(std::slice::from_ref(op), &system, k, &leader, max_terms)
}

/// One quotient `H_nu` of a [`HadamardFactorCertificate`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HadamardDivisionQuotient {
    /// The shift `nu` of `H_nu o S^nu o L`.
    pub shift: usize,
    /// `H_nu = sum_d h_d(n, x) D^d`.
    pub differentiations: Vec<OperatorTermJson>,
}

/// The invariance part of a [`HadamardFactorCertificate`] (module
/// documentation of `hadamard_factor`, section "Row annihilators"), kept
/// in the certificate's `row_annihilator_*` fields.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct RowAnnihilatorInvariance {
    /// `nu(n, x)`.
    pub(super) step_multiplier: Biv,
    /// `B(n)`.
    pub(super) step_quotient: DOp,
    /// `s`.
    pub(super) start_base: usize,
    /// The bases `m >= s` where `nu(m, .)` is the zero polynomial.
    pub(super) restart_bases: Vec<usize>,
}

/// One generator `A^(l)` of a `HadamardAnnihilatorSystem`, with its
/// invariance identity
///
/// ```text
/// nu_l(n, x) N(A^(l)) = sum_j B_(l,j)(n) o A^(j),
/// N(A) = g(n) S o A - A_(k-1)(n + 1) o L,
/// ```
///
/// and its quotient `Z_l` in the certificate identity.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HadamardAnnihilatorGenerator {
    /// `A^(l) = sum_(i < k) A^(l)_i(n) S^i` (`k` the parent order), as terms
    /// `c(n, x) D^d S^i`.
    pub annihilator: ImplicitRelationJson,
    /// `nu_l(n, x)` (`coeffs[i][j]` multiplies `n^i x^j`), nonzero.
    pub step_multiplier: BivarPolyJson,
    /// `B_(l,j)` for every generator `j` in the order of
    /// `HadamardAnnihilatorSystem::generators` (shift-free, possibly
    /// empty for zero).
    pub step_quotients: Vec<Vec<OperatorTermJson>>,
    /// `Z_l` (shift-free) of the certificate identity
    /// `mu w M'' = sum_nu H_nu S^nu L + sum_l Z_l A^(l)`.
    pub quotient: Vec<OperatorTermJson>,
}

/// An invariant row-annihilator system of a parent of order `k`
/// (`g(n) P_(n+k) = sum_(i < k) T_i(n) P_(n+i)`; module documentation of
/// `hadamard_factor`, section "Row annihilators"): every generator kills
/// the parent row windows `(P_n, ..., P_(n+k-1))` from `start_base` on.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HadamardAnnihilatorSystem {
    /// The generators `A^(l)`.
    pub generators: Vec<HadamardAnnihilatorGenerator>,
    /// `mu_2(n, x)`, a further multiplier of the certificate identity
    /// `mu(x) mu_2(n, x) w(n) M'' = sum_nu H_nu S^nu L + sum_l Z_l A^(l)`
    /// (`coeffs[i][j]` multiplies `n^i x^j`); every integer base `m >= s`
    /// where `mu_2(m, .)` is the zero polynomial is an exceptional base.
    pub row_multiplier: BivarPolyJson,
    /// The first base `s` of the induction (the first base where both the
    /// input and the parent relation are used); `A^(l) P = 0` at `s` was
    /// checked exactly for every `l`.
    pub start_base: usize,
    /// The bases `m >= s` where some `nu_l(m, .)` is the zero polynomial;
    /// there `(A^(l) P)_(m+1) = 0` was checked exactly for every `l`.
    pub restart_bases: Vec<usize>,
}

/// The exact operator identity behind a certified Hadamard-factor lift
/// `Q_n = lambda(n) sum_k u(n, k) p(n, k) x^k`:
///
/// ```text
/// mu(x) w(n) M'' = sum_nu H_nu o S^nu o L + Z o A
/// ```
///
/// with `M''` the input relation `M` conjugated by `lambda(n) u(n, k)`
/// (recomputed by [`verify`](Self::verify)), `L` the parent relation, `A`
/// a row annihilator: that of a product parent, or an invariant annihilator
/// of a first-order parent `g(n) P_(n+1) = T(n) P_n` with the invariance
/// identity
///
/// ```text
/// nu(n, x) A(n + 1) o T(n) = B(n) o A(n)
/// ```
///
/// (absent, with `Z = 0` and `mu = 1`, when the division leaves no
/// remainder).  For a parent of any order the term `Z o A` may instead be
/// a sum `sum_l Z_l o A^(l)` over an invariant row-annihilator system
/// (`HadamardAnnihilatorSystem`, field `row_annihilator_system`), and the
/// left side then carries a further factor `mu_2(n, x)`.  What it proves,
/// with the index domains recorded here, is stated in the
/// `hadamard_factor` module documentation.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HadamardFactorCertificate {
    /// The factor `u`.
    pub factor: HadamardFactor,
    /// The row normalisation `lambda(n)`.
    pub gauge: VeroneseSectionGauge,
    /// `M`: the input rule at base index `n`.
    pub input_relation: ImplicitRelationJson,
    /// The largest offset `R` of the input rule.
    pub input_leader_shift: usize,
    /// The parent `P`: rule, first index (equal to the input's) and seeds.
    pub parent: RecurrenceJson,
    /// `L`: the parent rule at base index `n`; its leader `g(n) S^k` is free
    /// of `x` and `D`.
    pub parent_relation: ImplicitRelationJson,
    /// `k`, the order of the parent rule.
    pub parent_leader_shift: usize,
    /// `M''` (multiplied by `x^J''`).
    pub conjugated_input_relation: ImplicitRelationJson,
    /// `J''`.
    pub conjugation_x_shift: usize,
    /// `w(n)`, ascending exact strings in `n`.
    pub multiplier: Vec<String>,
    /// `mu(x)`, ascending exact strings in `x`.
    pub row_multiplier: Vec<String>,
    /// The nonzero quotients `H_nu`.
    pub quotients: Vec<HadamardDivisionQuotient>,
    /// `A(n)`, when used: recomputed from the parent by the replay for a
    /// product parent, certified by `row_annihilator_invariance` otherwise.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub row_annihilator: Option<Vec<OperatorTermJson>>,
    /// `nu(n, x)` of the invariance identity (`coeffs[i][j]` multiplies
    /// `n^i x^j`), present exactly when `A` is an invariant annihilator
    /// rather than the product-parent annihilator.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub row_annihilator_step_multiplier: Option<BivarPolyJson>,
    /// `B(n)` of the invariance identity (empty unless
    /// `row_annihilator_step_multiplier` is present).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub row_annihilator_step_quotient: Vec<OperatorTermJson>,
    /// The first base `s` of the annihilator induction, where `A(s) P_s = 0`
    /// was checked exactly (with `row_annihilator_step_multiplier`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub row_annihilator_start_base: Option<usize>,
    /// The bases `m >= s` where `nu(m, .)` is the zero polynomial; there
    /// `A(m + 1) P_(m+1) = 0` was checked exactly.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub row_annihilator_restart_bases: Vec<usize>,
    /// `Z`.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub row_annihilator_quotient: Vec<OperatorTermJson>,
    /// An invariant row-annihilator system of the parent with its quotients
    /// `Z_l`, used in place of `row_annihilator` (which is then absent).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub row_annihilator_system: Option<HadamardAnnihilatorSystem>,
    /// First base index at which `M` holds for the input rows (`n0 + K - R`).
    pub input_valid_from_base: usize,
    /// First base index at which `L` holds for the parent rows.
    pub parent_valid_from_base: usize,
    /// Bases where the identity gives nothing (integer roots of `w`, bases
    /// below `parent_valid_from_base`); `M` was evaluated exactly on the
    /// candidate rows there and vanished.
    pub exceptional_bases: Vec<usize>,
    /// The identifying input indices `n0, ..., n0 + K - 1`, where
    /// `Q_n = lambda(n) u ⊙ P_n` was checked exactly.
    pub identifying_indices: Vec<usize>,
    /// What the certificate proves.
    pub statement: String,
}

fn relation_op(json: &ImplicitRelationJson) -> Result<DOp, RecurrenceError> {
    Ok(ImplicitRelation::from_json(json)?.op())
}

/// The left-hand and right-hand sides of the certificate identity agree;
/// `annihilators` lists the pairs `(A, Z)` of the terms `Z o A`.
pub(super) fn identity_holds(
    conjugated: &DOp,
    parent: &DOp,
    multiplier: &Biv,
    quotients: &BTreeMap<usize, DOp>,
    annihilators: &[(DOp, DOp)],
) -> bool {
    let mut difference = left_multiply(multiplier, conjugated);
    for (nu, h) in quotients {
        difference = difference.sub(&h.compose(&parent.shifted(*nu)));
    }
    for (a, z) in annihilators {
        difference = difference.sub(&z.compose(a));
    }
    difference.is_zero()
}

impl HadamardFactorCertificate {
    /// Replay the certificate: recompute `M''` from the input relation, the
    /// factor, the gauge and the first index (including the conjugation
    /// conditions), recompute the parent relation and, when used, the row
    /// annihilator with its invariance, check the operator identity exactly
    /// and that every integer root of `w` from the valid bases on is listed
    /// as exceptional.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        let failed = |what: &str| Err(RecurrenceError::VerificationFailed(what.to_string()));
        let Some(&first_index) = self.identifying_indices.first() else {
            return failed("the certificate lists no identifying indices");
        };
        let input = relation_op(&self.input_relation)?;
        if input.max_shift() != self.input_leader_shift
            || input.column(self.input_leader_shift).max_derivative() != 0
        {
            return failed("the input relation does not have the recorded leader");
        }
        let gauge = self.gauge.ratio_form(first_index)?;
        let conjugated = conjugate_input(&input, &self.factor, &gauge, first_index)?;
        if conjugated.op != relation_op(&self.conjugated_input_relation)?
            || conjugated.x_shift != self.conjugation_x_shift
        {
            return failed("the conjugated input relation does not match its recomputation");
        }
        let (rule, parent_first, rows) = self.parent.to_recurrence_parts()?;
        if parent_first != first_index {
            return failed("the parent does not start at the input's first index");
        }
        let parent = super::prepare::prepare_rule(&rule, &rows, parent_first, "parent")?;
        if parent.relation != relation_op(&self.parent_relation)? {
            return failed("the parent relation does not match the parent record");
        }
        let (k, _) = x_free_leader(&parent.relation)?;
        if k != self.parent_leader_shift {
            return failed("the parent leader shift does not match");
        }
        let multiplier = n_polynomial(&decode_poly(&self.multiplier)?);
        let row_multiplier = Biv::from_y(&decode_poly(&self.row_multiplier)?);
        if multiplier.is_zero() || row_multiplier.is_zero() {
            return failed("a multiplier is zero");
        }
        let mut quotients: BTreeMap<usize, DOp> = BTreeMap::new();
        for entry in &self.quotients {
            let h = shift_free_from_json(&entry.differentiations)?;
            let slot = quotients.entry(entry.shift).or_default();
            *slot = slot.add(&h);
        }
        let from = self.input_valid_from_base.max(self.parent_valid_from_base);
        let invariance = match (
            &self.row_annihilator_step_multiplier,
            self.row_annihilator_start_base,
        ) {
            (Some(multiplier), Some(start_base)) => Some(RowAnnihilatorInvariance {
                step_multiplier: biv_from_json(multiplier)?,
                step_quotient: shift_free_from_json(&self.row_annihilator_step_quotient)?,
                start_base,
                restart_bases: self.row_annihilator_restart_bases.clone(),
            }),
            (None, None)
                if self.row_annihilator_step_quotient.is_empty()
                    && self.row_annihilator_restart_bases.is_empty() =>
            {
                None
            }
            _ => return failed("incomplete row-annihilator invariance fields"),
        };
        let annihilator = match (&self.row_annihilator, &invariance) {
            (None, None) => {
                if !self.row_annihilator_quotient.is_empty() {
                    return failed("a row-annihilator quotient without an annihilator");
                }
                None
            }
            _ if self.row_annihilator_system.is_some() => {
                return failed("a row annihilator and a row-annihilator system together");
            }
            (None, Some(_)) => return failed("an invariance certificate without an annihilator"),
            (Some(recorded), None) => {
                let Some(product) = ProductParent::recognise(&parent) else {
                    return failed("the parent is not a product parent");
                };
                let a = product.annihilator();
                if a != shift_free_from_json(recorded)? || !product.annihilator_is_invariant(&a) {
                    return failed("the row annihilator does not match the parent");
                }
                Some((a, shift_free_from_json(&self.row_annihilator_quotient)?))
            }
            (Some(recorded), Some(invariance)) => {
                if self.parent_valid_from_base != parent.valid_from_base().max(first_index) {
                    return failed("the parent base does not match the parent record");
                }
                let a = shift_free_from_json(recorded)?;
                verify_invariance(&a, invariance, &parent, from)?;
                Some((a, shift_free_from_json(&self.row_annihilator_quotient)?))
            }
        };
        let mut annihilators: Vec<(DOp, DOp)> = annihilator.into_iter().collect();
        let mut identity_multiplier = row_multiplier.mul(&multiplier);
        if let Some(system) = &self.row_annihilator_system {
            if self.parent_valid_from_base != parent.valid_from_base().max(first_index) {
                return failed("the parent base does not match the parent record");
            }
            let (decoded, terms) = system_from_json(system)?;
            verify_system(
                &decoded,
                system.start_base,
                &system.restart_bases,
                &parent,
                from,
            )?;
            annihilators.extend(terms);
            let second = biv_from_json(&system.row_multiplier)?;
            if second.is_zero() {
                return failed("the row-annihilator system multiplier is zero");
            }
            if vanishing_bases(&second, from)
                .iter()
                .any(|base| !self.exceptional_bases.contains(base))
            {
                return failed(
                    "a zero of the row-annihilator system multiplier is not an exceptional base",
                );
            }
            identity_multiplier = identity_multiplier.mul(&second);
        }
        if !identity_holds(
            &conjugated.op,
            &parent.relation,
            &identity_multiplier,
            &quotients,
            &annihilators,
        ) {
            return failed("the identity mu w M'' = sum H_nu S^nu L + Z A fails");
        }
        let roots = roots_in_domain(&n_coefficients(&multiplier), from, None);
        if roots
            .iter()
            .any(|root| !self.exceptional_bases.contains(root))
        {
            return failed("a root of w is not listed as an exceptional base");
        }
        Ok(())
    }
}

/// Replay the invariance of the annihilator `a` of the prepared `parent`:
/// the parent has order one, the identity `nu A(n + 1) o T(n) = B o A(n)`
/// holds exactly, the start base is `from`, the restart bases are the
/// integer bases `m >= from` with `nu(m, .) = 0`, and `A` kills the parent
/// rows at the start base and after every restart base (exact).
pub(super) fn verify_invariance(
    a: &DOp,
    invariance: &RowAnnihilatorInvariance,
    parent: &super::prepare::PreparedRule,
    from: usize,
) -> Result<(), RecurrenceError> {
    let failed = |what: &str| Err(RecurrenceError::VerificationFailed(what.to_string()));
    let Some(step) = step_operator(&parent.relation) else {
        return failed("an invariant row annihilator needs a first-order parent");
    };
    let invariant = InvariantAnnihilator {
        annihilator: a.clone(),
        step_multiplier: invariance.step_multiplier.clone(),
        step_quotient: invariance.step_quotient.clone(),
    };
    if invariant.step_multiplier.is_zero() || !invariance_holds(&invariant, &step) {
        return failed("the invariance identity nu A(n + 1) T(n) = B A(n) fails");
    }
    if invariance.start_base != from {
        return failed("the annihilator induction does not start at the first used base");
    }
    if vanishing_bases(&invariant.step_multiplier, from) != invariance.restart_bases {
        return failed("the restart bases are not the zeros of nu");
    }
    check_annihilated_rows(a, parent, from, &invariance.restart_bases)
}

/// `A(m) P_m = 0` exactly at `m = from` and at `m = r + 1` for every restart
/// base `r`.
pub(super) fn check_annihilated_rows(
    a: &DOp,
    parent: &super::prepare::PreparedRule,
    from: usize,
    restart_bases: &[usize],
) -> Result<(), RecurrenceError> {
    let indices: Vec<usize> = std::iter::once(from)
        .chain(restart_bases.iter().map(|base| base + 1))
        .collect();
    let last = indices.iter().copied().max().unwrap_or(from);
    let count = last + 1 - parent.first_index.min(last);
    if count > super::search::MAX_HADAMARD_FACTOR_ROWS + parent.seeds.len() {
        return Err(RecurrenceError::LimitExceeded(format!(
            "the annihilator check needs {count} parent rows"
        )));
    }
    let rows = parent.rows(count)?;
    for index in indices {
        let row = &rows[index - parent.first_index];
        if !a.apply_row(row, index).iter().all(Zero::is_zero) {
            return Err(RecurrenceError::VerificationFailed(format!(
                "the row annihilator does not kill parent row {index}"
            )));
        }
    }
    Ok(())
}

/// The system of its JSON form, with the terms `(A^(l), Z_l)`.
fn system_from_json(
    system: &HadamardAnnihilatorSystem,
) -> Result<(InvariantAnnihilatorSystem, Vec<(DOp, DOp)>), RecurrenceError> {
    let mut decoded = InvariantAnnihilatorSystem {
        generators: Vec::new(),
        step_multipliers: Vec::new(),
        step_quotients: Vec::new(),
    };
    let mut terms = Vec::new();
    for generator in &system.generators {
        let a = relation_op(&generator.annihilator)?;
        decoded.generators.push(a.clone());
        decoded
            .step_multipliers
            .push(biv_from_json(&generator.step_multiplier)?);
        decoded.step_quotients.push(
            generator
                .step_quotients
                .iter()
                .map(|b| shift_free_from_json(b))
                .collect::<Result<_, _>>()?,
        );
        terms.push((a, shift_free_from_json(&generator.quotient)?));
    }
    Ok((decoded, terms))
}

/// Replay an invariant row-annihilator system of the prepared `parent`:
/// every invariance identity holds exactly, the start base is `from`, the
/// restart bases are the integer bases `m >= from` where some `nu_l(m, .)`
/// vanishes, and every generator kills the parent rows at the start base
/// and after every restart base (exact).
pub(super) fn verify_system(
    system: &InvariantAnnihilatorSystem,
    start_base: usize,
    restart_bases: &[usize],
    parent: &super::prepare::PreparedRule,
    from: usize,
) -> Result<(), RecurrenceError> {
    let failed = |what: &str| Err(RecurrenceError::VerificationFailed(what.to_string()));
    if !system_invariance_holds(system, &parent.relation)? {
        return failed("the invariance identities nu_l N(A^(l)) = sum_j B_(l,j) A^(j) fail");
    }
    if start_base != from {
        return failed("the annihilator induction does not start at the first used base");
    }
    if system.restart_bases(from) != restart_bases {
        return failed("the restart bases are not the zeros of the step multipliers");
    }
    check_system_rows(&system.generators, parent, from, restart_bases)
}

/// `(A^(l) P)_m = 0` exactly for every generator at `m = from` and at
/// `m = r + 1` for every restart base `r`.
pub(super) fn check_system_rows(
    generators: &[DOp],
    parent: &super::prepare::PreparedRule,
    from: usize,
    restart_bases: &[usize],
) -> Result<(), RecurrenceError> {
    let indices: Vec<usize> = std::iter::once(from)
        .chain(restart_bases.iter().map(|base| base + 1))
        .collect();
    let width = generators.iter().map(DOp::max_shift).max().unwrap_or(0);
    let last = indices.iter().copied().max().unwrap_or(from) + width;
    let count = last + 1 - parent.first_index.min(last);
    if count > super::search::MAX_HADAMARD_FACTOR_ROWS + parent.seeds.len() {
        return Err(RecurrenceError::LimitExceeded(format!(
            "the annihilator check needs {count} parent rows"
        )));
    }
    let rows = parent.rows(count)?;
    for index in indices {
        for a in generators {
            let value = a.apply(&rows, parent.first_index, index).ok_or_else(|| {
                RecurrenceError::VerificationFailed(format!("parent rows missing at base {index}"))
            })?;
            if !value.iter().all(Zero::is_zero) {
                return Err(RecurrenceError::VerificationFailed(format!(
                    "the row-annihilator system does not kill the parent rows at base {index}"
                )));
            }
        }
    }
    Ok(())
}

/// How the division remainder is removed: `mu r = Z o A` (or
/// `sum_l Z_l o A^(l)`) with `mu(n, x) = c(n) p(x)`; `c` is folded into the
/// multiplier `w(n)`.
pub(super) struct RowAnnihilation {
    /// `p(x)`.
    pub(super) row_multiplier: Vec<Q>,
    /// `c(n)`.
    pub(super) index_multiplier: Vec<Q>,
    /// The annihilators.
    pub(super) kind: RowAnnihilationKind,
}

/// One annihilator, or an invariant system.
pub(super) enum RowAnnihilationKind {
    /// `Z o A`.
    Single {
        /// `A(n)`.
        annihilator: DOp,
        /// `Z(n)`.
        quotient: DOp,
        /// The invariance certificate, absent for a product parent.
        invariance: Option<RowAnnihilatorInvariance>,
    },
    /// `sum_l Z_l o A^(l)`.
    System {
        /// The generators and their invariance identities.
        system: InvariantAnnihilatorSystem,
        /// `mu_2(n, x)`.
        row_multiplier: Biv,
        /// `Z_l`.
        quotients: Vec<DOp>,
        /// `s`.
        start_base: usize,
        /// The bases `m >= s` where some `nu_l(m, .)` vanishes.
        restart_bases: Vec<usize>,
    },
}

impl RowAnnihilation {
    /// The invariance certificate of `invariant` starting at `start_base`.
    pub(super) fn invariance_record(
        invariant: &InvariantAnnihilator,
        start_base: usize,
    ) -> RowAnnihilatorInvariance {
        RowAnnihilatorInvariance {
            step_multiplier: invariant.step_multiplier.clone(),
            step_quotient: invariant.step_quotient.clone(),
            start_base,
            restart_bases: vanishing_bases(&invariant.step_multiplier, start_base),
        }
    }
}

/// The JSON form of an invariant system with its quotients `Z_l`.
fn system_json(
    system: &InvariantAnnihilatorSystem,
    quotients: &[DOp],
    start_base: usize,
    restart_bases: Vec<usize>,
) -> HadamardAnnihilatorSystem {
    HadamardAnnihilatorSystem {
        generators: system
            .generators
            .iter()
            .zip(&system.step_multipliers)
            .zip(&system.step_quotients)
            .zip(quotients)
            .map(|(((a, nu), b), z)| HadamardAnnihilatorGenerator {
                annihilator: ImplicitRelation::from_op(a).to_json(),
                step_multiplier: biv_json(nu),
                step_quotients: b.iter().map(shift_free_json).collect(),
                quotient: shift_free_json(z),
            })
            .collect(),
        row_multiplier: biv_json(&Biv::one()),
        start_base,
        restart_bases,
    }
}

/// Assemble the certificate fields from a successful division; `multiplier`
/// is `w(n)` with the index part of the row multiplier folded in.
#[allow(clippy::too_many_arguments)]
pub(super) fn build_certificate(
    factor: &HadamardFactor,
    gauge: &VeroneseSectionGauge,
    input: &PreparedRule,
    parent: &PreparedRule,
    conjugated: &ConjugatedRelation,
    division: &ForwardDivision,
    multiplier: &Biv,
    annihilation: Option<RowAnnihilation>,
    bases: (usize, usize, Vec<usize>),
    statement: String,
) -> HadamardFactorCertificate {
    let (row_multiplier, index_multiplier, kind) = match annihilation {
        Some(annihilation) => (
            annihilation.row_multiplier,
            annihilation.index_multiplier,
            Some(annihilation.kind),
        ),
        None => (vec![Q::one()], vec![Q::one()], None),
    };
    let mut full_row_multiplier =
        Biv::from_y(&row_multiplier).mul(&n_polynomial(&index_multiplier));
    let (row_annihilator, quotient, invariance, system) = match kind {
        None => (None, Vec::new(), None, None),
        Some(RowAnnihilationKind::Single {
            annihilator,
            quotient,
            invariance,
        }) => (
            Some(shift_free_json(&annihilator)),
            shift_free_json(&quotient),
            invariance,
            None,
        ),
        Some(RowAnnihilationKind::System {
            system,
            row_multiplier: second,
            quotients,
            start_base,
            restart_bases,
        }) => {
            full_row_multiplier = full_row_multiplier.mul(&second);
            let mut json = system_json(&system, &quotients, start_base, restart_bases);
            json.row_multiplier = biv_json(&second);
            (None, Vec::new(), None, Some(json))
        }
    };
    let quotients = division
        .quotients
        .iter()
        .filter_map(|(nu, row)| {
            let h = left_multiply(&full_row_multiplier, &row[0]);
            (!h.is_zero()).then(|| HadamardDivisionQuotient {
                shift: *nu,
                differentiations: shift_free_json(&h),
            })
        })
        .collect();
    let (k, _) = x_free_leader(&parent.relation).unwrap_or((1, Biv::one()));
    HadamardFactorCertificate {
        factor: factor.clone(),
        gauge: gauge.clone(),
        input_relation: ImplicitRelation::from_op(&input.relation).to_json(),
        input_leader_shift: input.leader_shift,
        parent: RecurrenceJson::from_recurrence_rational(
            &parent.rule,
            parent.first_index,
            &parent.seeds,
            None,
        ),
        parent_relation: ImplicitRelation::from_op(&parent.relation).to_json(),
        parent_leader_shift: k,
        conjugated_input_relation: ImplicitRelation::from_op(&conjugated.op).to_json(),
        conjugation_x_shift: conjugated.x_shift,
        multiplier: encode_poly(&n_coefficients(multiplier)),
        row_multiplier: encode_poly(&row_multiplier),
        quotients,
        row_annihilator,
        row_annihilator_step_multiplier: invariance.as_ref().map(|i| biv_json(&i.step_multiplier)),
        row_annihilator_step_quotient: invariance
            .as_ref()
            .map_or_else(Vec::new, |i| shift_free_json(&i.step_quotient)),
        row_annihilator_start_base: invariance.as_ref().map(|i| i.start_base),
        row_annihilator_restart_bases: invariance.map_or_else(Vec::new, |i| i.restart_bases),
        row_annihilator_quotient: quotient,
        row_annihilator_system: system,
        input_valid_from_base: bases.0,
        parent_valid_from_base: bases.1,
        exceptional_bases: bases.2,
        identifying_indices: (input.first_index..input.first_index + input.seeds.len()).collect(),
        statement,
    }
}

#[cfg(test)]
mod tests {
    use super::super::factor::{HadamardFactorTerm, TermShape};
    use super::*;
    use crate::recurrence::implicit::algebra::q;

    fn y(coefficients: &[i64]) -> Biv {
        Biv::from_y(&coefficients.iter().map(|c| q(*c)).collect::<Vec<_>>())
    }

    fn pseudo_random_rows(count: usize, seed: u64) -> Vec<Vec<Q>> {
        let mut state = seed;
        (0..count)
            .map(|n| {
                (0..n + 3)
                    .map(|_| {
                        state = state
                            .wrapping_mul(6364136223846793005)
                            .wrapping_add(1442695040888963407);
                        q(((state >> 33) % 7) as i64 - 3)
                    })
                    .collect()
            })
            .collect()
    }

    /// The conjugation lemma, checked on arbitrary rows: wherever
    /// `U(n + R, k) != 0`, `D(n, k) E(n, k) = U(n + R, k) [x^(k + J'')](M'' P)_n`,
    /// and elsewhere `E(n, k) = 0`.
    #[test]
    fn conjugation_lemma_holds_on_arbitrary_rows() {
        // M = (n + 3) S^2 - (2n + 1)(1 + x) S - x^2 D S + (n + 1)(1 - x)^2 + (1 + x) D.
        let mut m = DOp::default();
        m.add_term((2, 0), Biv::new(vec![vec![q(3)], vec![q(1)]]));
        m.add_term(
            (1, 0),
            Biv::new(vec![vec![q(-1), q(-1)], vec![q(-2), q(-2)]]),
        );
        m.add_term((1, 1), y(&[0, 0, -1]));
        m.add_term(
            (0, 0),
            Biv::new(vec![vec![q(1), q(-2), q(1)], vec![q(1), q(-2), q(1)]]),
        );
        m.add_term((0, 1), y(&[1, 1]));
        let gauges = [
            RatioGauge {
                numerator: vec![Q::one()],
                denominator: vec![Q::one()],
                initial: Q::one(),
            },
            RatioGauge {
                numerator: vec![q(1), q(1)],
                denominator: vec![q(3), q(1)],
                initial: q(2),
            },
        ];
        let factors = [
            HadamardFactorTerm::Binomial {
                index_offset: 1,
                coefficient_offset: 1,
            },
            HadamardFactorTerm::Binomial {
                index_offset: 0,
                coefficient_offset: 0,
            },
            HadamardFactorTerm::InverseFactorial,
            HadamardFactorTerm::InverseComplementFactorial,
            HadamardFactorTerm::RisingFactorial {
                parameter: "1/2".into(),
            },
        ];
        let n0 = 1;
        let r = 2;
        let rows = pseudo_random_rows(12, 7);
        let mut checked = 0;
        for term in factors {
            let factor = HadamardFactor::single(term);
            let terms: Vec<(usize, i64)> = laurent_euler(&m).keys().copied().collect();
            let atoms = factor.atoms().unwrap();
            let shape = TermShape {
                terms: &terms,
                leader_shift: r,
            };
            assert!(super::super::factor::conjugation_conditions(&atoms, n0, &shape).is_ok());
            for gauge in &gauges {
                let conjugated = conjugate_input(&m, &factor, gauge, n0).unwrap();
                let lambdas = gauge.values(n0, rows.len());
                let scaled: Vec<Vec<Q>> = rows
                    .iter()
                    .enumerate()
                    .map(|(offset, row)| {
                        let n = n0 + offset;
                        vec_poly::trim(
                            &row.iter()
                                .enumerate()
                                .map(|(k, p)| factor.value(n, k).unwrap() * p * &lambdas[offset])
                                .collect::<Vec<_>>(),
                        )
                    })
                    .collect();
                for base in n0..n0 + rows.len() - r {
                    let e = m.apply(&scaled, n0, base).unwrap();
                    let conj = conjugated.op.apply(&rows, n0, base).unwrap();
                    let at = |p: &[Q], i: usize| p.get(i).cloned().unwrap_or_else(Q::zero);
                    for k in 0..e.len() + 4 {
                        let reference =
                            factor.value(base + r, k).unwrap() * &lambdas[base + r - n0];
                        if reference.is_zero() {
                            assert!(at(&e, k).is_zero());
                            continue;
                        }
                        let d = vec_poly::evaluate(
                            &conjugated.denominator.eval_n(&q(base as i64)),
                            &q(k as i64),
                        );
                        assert!(!d.is_zero());
                        assert_eq!(
                            &d * at(&e, k),
                            reference * at(&conj, k + conjugated.x_shift)
                        );
                        checked += 1;
                    }
                }
            }
        }
        assert!(checked > 100);
    }

    #[test]
    fn right_division_by_a_first_order_operator() {
        // A = x D - 2; F = (x D - 2)(D + x) = x D^2 + x^2 D - D + 2x... computed by composition.
        let mut a = DOp::single(0, 1, y(&[0, 1]));
        a.add_term((0, 0), y(&[-2]));
        let mut z = DOp::single(0, 1, Biv::one());
        z.add_term((0, 0), y(&[0, 1]));
        let f = z.compose(&a);
        let (mu, quotient, remainder) = right_pseudo_divide(&f, &a, 1000).unwrap();
        assert!(remainder.is_zero());
        assert_eq!(left_multiply(&mu, &f), quotient.compose(&a));
    }
}
