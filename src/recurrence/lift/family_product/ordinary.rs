//! Ordinary product closure through tensor companion systems in the Ore algebra.

use super::*;
use crate::recurrence::implicit::biv_json;
use crate::recurrence::lift::hadamard_factor::annihilator::vanishing_bases;
use crate::recurrence::lift::section_certify::x_free_leader;

/// Exact closure identity for an ordinary row product. The square tensor
/// system is reconstructed from both parent rules during replay.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PolynomialFamilyOrdinaryClosure {
    /// `w(n)` in `w(n) M e_0 = sum_b H_b Phi_b`, nonzero and free of x.
    pub multiplier: BivarPolyJson,
    /// `H_b`, one per tensor-system row, in canonical operator form.
    pub quotients: Vec<ImplicitRelationJson>,
    /// Size of the tensor system (product of the parent orders).
    pub system_size: usize,
    /// First base at which both companion rules apply.
    pub valid_from_base: usize,
    /// Bases below the companion start or at zeros of w; the target relation
    /// was checked on the actual product rows at each one.
    pub exceptional_bases: Vec<usize>,
    /// Indices at which the product equals the target's identifying seeds.
    pub identifying_indices: Vec<usize>,
}

fn companion(parent: &PreparedLiftRule) -> Result<(Biv, Vec<Vec<Biv>>), RecurrenceError> {
    if parent.relation.max_derivative() != 0 {
        return Err(unsupported(
            "tensor companion closure currently requires shift-only parent rules",
        ));
    }
    let (order, leader) = x_free_leader(&parent.relation)?;
    if order == 0 || order > 4 {
        return Err(unsupported(
            "parent order exceeds the tensor companion bound",
        ));
    }
    let mut matrix = vec![vec![Biv::zero(); order]; order];
    for (i, row) in matrix.iter_mut().enumerate().take(order - 1) {
        row[i + 1] = leader.clone();
    }
    for (j, value) in matrix[order - 1].iter_mut().enumerate() {
        *value = parent
            .relation
            .0
            .get(&(j, 0))
            .cloned()
            .unwrap_or_default()
            .scale(&(-Q::one()));
    }
    Ok((leader, matrix))
}

pub(super) fn system(
    left: &PreparedLiftRule,
    right: &PreparedLiftRule,
) -> Result<Vec<Vec<DOp>>, RecurrenceError> {
    let (ga, a) = companion(left)?;
    let (gb, b) = companion(right)?;
    let size = a.len() * b.len();
    if size > 9 {
        return Err(unsupported("tensor companion size exceeds nine"));
    }
    let denominator = ga.multiply(&gb);
    let mut system = vec![vec![DOp::default(); size]; size];
    for i in 0..a.len() {
        for j in 0..b.len() {
            let row = i * b.len() + j;
            system[row][row].add_term((1, 0), denominator.clone());
            for h in 0..a.len() {
                for k in 0..b.len() {
                    system[row][h * b.len() + k]
                        .add_term((0, 0), a[i][h].multiply(&b[j][k]).scale(&(-Q::one())));
                }
            }
        }
    }
    Ok(system)
}

pub(super) fn certify(
    image: &PreparedLiftRule,
    left: &PreparedLiftRule,
    right: &PreparedLiftRule,
) -> Result<PolynomialFamilyOrdinaryClosure, RecurrenceError> {
    let system = system(left, right)?;
    let certificate = ore_algebra::certify_component_annihilator(&image.relation, &system, 0)
        .map_err(|error| unsupported(error.to_string()))?
        .ok_or_else(|| {
            unsupported("the product rule is not a certified tensor-component annihilator")
        })?;
    certificate
        .verify()
        .map_err(|error| unsupported(error.to_string()))?;
    let from = left.valid_from_base().max(right.valid_from_base());
    let image_from = image.valid_from_base();
    let mut exceptional_bases: Vec<_> = (image_from..from).collect();
    exceptional_bases.extend(vanishing_bases(
        &certificate.multiplier,
        image_from.max(from),
    ));
    exceptional_bases.sort_unstable();
    exceptional_bases.dedup();
    let last = exceptional_bases.last().copied().unwrap_or(image_from) + image.leader_shift;
    let count = (last - image.first_index + 1).max(image.seeds.len());
    if count > 64 {
        return Err(unsupported(
            "a tensor exceptional row exceeds the check bound",
        ));
    }
    let rows = image_rows(
        left,
        right,
        PolynomialFamilyProductOperation::Ordinary,
        count,
    )?;
    if rows[..image.seeds.len()] != image.seeds {
        return Err(unsupported("product identifying rows differ"));
    }
    for &base in &exceptional_bases {
        let values = image
            .relation
            .apply(&rows, image.first_index, base)
            .ok_or_else(|| unsupported("product rows missing at an exceptional base"))?;
        if values.iter().any(|v| !v.is_zero()) {
            return Err(unsupported("the product rule fails at an exceptional base"));
        }
    }
    Ok(PolynomialFamilyOrdinaryClosure {
        multiplier: biv_json(&certificate.multiplier),
        quotients: certificate
            .quotients
            .iter()
            .map(|q| ImplicitRelation::from_op(q).to_json())
            .collect(),
        system_size: system.len(),
        valid_from_base: from,
        exceptional_bases,
        identifying_indices: (image.first_index..image.first_index + image.seeds.len()).collect(),
    })
}
