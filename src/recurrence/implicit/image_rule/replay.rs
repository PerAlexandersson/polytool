//! Certificate replay and exact supporting checks.

use super::*;

// ---------------------------------------------------------------------------
// JSON of operators
// ---------------------------------------------------------------------------

pub(super) fn operator_json(op: &DOp) -> ImplicitRelationJson {
    ImplicitRelation::from_op(op).to_json()
}

pub(super) fn operator_from_json(json: &ImplicitRelationJson) -> Result<DOp, RecurrenceError> {
    Ok(ImplicitRelation::from_json(json)?.op())
}

/// The runnable rule of `M` (target form) and its order.
pub(super) fn successor_rule(m: &DOp) -> Result<(Recurrence, usize), RecurrenceError> {
    ImplicitRelation::from_op(m)
        .to_successor_rule()
        .ok_or_else(|| RecurrenceError::VerificationFailed("the image rule is not runnable".into()))
}

// ---------------------------------------------------------------------------
// The construction and its checks
// ---------------------------------------------------------------------------

/// Bases `n >= from` where `p(n, .)` is the zero polynomial.
fn zero_bases(p: &Biv, from: usize) -> Vec<usize> {
    let (content, _) = content_and_primitive(p);
    roots_in_domain(&content, from, None)
}

/// The bases `m` where each member must be checked (module documentation,
/// step 2), image indexing; `first_base` is the first base of `L'`.
pub(super) fn system_checks(
    certificate: &AnnihilatorSystemCertificate,
    n0: usize,
    first_base: usize,
) -> Vec<Vec<usize>> {
    let j = certificate.step;
    let start = n0.max(first_base);
    certificate
        .entries
        .iter()
        .map(|entry| {
            let mut bases: Vec<usize> = (n0..start + j).collect();
            for base in zero_bases(&entry.multiplier, n0) {
                bases.push(base + j);
            }
            bases.sort_unstable();
            bases.dedup();
            bases
        })
        .collect()
}

/// Check `(C P')_m = 0` for every member at its bases, exactly.
pub(super) fn system_holds(
    certificate: &AnnihilatorSystemCertificate,
    checks: &[Vec<usize>],
    prepared: &PreparedImage,
    parent: &Parent,
    n0: usize,
) -> Result<bool, RecurrenceError> {
    let needed = certificate
        .members
        .iter()
        .zip(checks)
        .filter_map(|(c, bases)| bases.last().map(|m| m + c.max_shift() + 1 - n0))
        .max()
        .unwrap_or(0);
    require_row_budget(needed)?;
    let rows = prepared.parent_rows(parent, n0, needed)?;
    Ok(certificate.members.iter().zip(checks).all(|(c, bases)| {
        bases.iter().all(|m| {
            c.apply(&rows, n0, *m)
                .is_some_and(|value| vec_poly::is_zero(&value))
        })
    }))
}

pub(super) fn system_json(
    certificate: &AnnihilatorSystemCertificate,
    checks: Vec<Vec<usize>>,
) -> ImageAnnihilatorSystem {
    ImageAnnihilatorSystem {
        step: certificate.step,
        members: certificate
            .members
            .iter()
            .zip(&certificate.entries)
            .zip(checks)
            .map(|((c, entry), bases)| ImageAnnihilatorSystemMember {
                annihilator: operator_json(c),
                multiplier: biv_json(&entry.multiplier),
                recurrence_quotient: operator_json(&entry.recurrence_quotient),
                member_quotients: entry.member_quotients.iter().map(operator_json).collect(),
                checked_bases: bases,
            })
            .collect(),
    }
}

/// Replay a recorded annihilator system relative to `L'`.
pub(super) fn replay_system(
    json: &ImageAnnihilatorSystem,
    recurrence: &DOp,
    first_base: usize,
    prepared: &PreparedImage,
    parent: &Parent,
    n0: usize,
) -> Result<AnnihilatorSystemCertificate, RecurrenceError> {
    let members = json
        .members
        .iter()
        .map(|m| operator_from_json(&m.annihilator))
        .collect::<Result<Vec<_>, _>>()?;
    let entries = json
        .members
        .iter()
        .map(|m| {
            Ok(AnnihilatorSystemEntry {
                multiplier: biv_from_json(&m.multiplier)?,
                recurrence_quotient: operator_from_json(&m.recurrence_quotient)?,
                member_quotients: m
                    .member_quotients
                    .iter()
                    .map(operator_from_json)
                    .collect::<Result<Vec<_>, RecurrenceError>>()?,
            })
        })
        .collect::<Result<Vec<_>, RecurrenceError>>()?;
    let certificate = AnnihilatorSystemCertificate {
        recurrence: recurrence.clone(),
        members,
        step: json.step,
        entries,
    };
    certificate.verify().map_err(|error| {
        RecurrenceError::VerificationFailed(format!("annihilator system: {error}"))
    })?;
    if certificate
        .members
        .iter()
        .filter(|c| c.max_shift() > 0)
        .count()
        > 1
        || certificate
            .members
            .iter()
            .filter(|c| c.max_shift() == 0)
            .count()
            > 1
    {
        return Err(RecurrenceError::VerificationFailed(
            "the annihilator system has at most one member of each kind".into(),
        ));
    }
    let checks = system_checks(&certificate, n0, first_base);
    let recorded: Vec<Vec<usize>> = json
        .members
        .iter()
        .map(|m| m.checked_bases.clone())
        .collect();
    if checks != recorded {
        return Err(RecurrenceError::VerificationFailed(format!(
            "annihilator system: recorded checks {recorded:?}, recomputed {checks:?}"
        )));
    }
    if !system_holds(&certificate, &checks, prepared, parent, n0)? {
        return Err(RecurrenceError::VerificationFailed(
            "annihilator system: a base row check fails".into(),
        ));
    }
    Ok(certificate)
}

/// The relation in the identity and its first base, and the row
/// annihilator, from an annihilator system (module documentation).
pub(super) fn working_relation(
    system: Option<&AnnihilatorSystemCertificate>,
    relation: &DOp,
    first_base: usize,
    n0: usize,
) -> (DOp, usize, Option<DOp>) {
    match system {
        None => (relation.clone(), first_base, None),
        Some(system) => {
            let row = system.row_annihilator().cloned();
            match system.lower_order_member() {
                Some(b) => (b.clone(), n0, row),
                None => (relation.clone(), first_base, row),
            }
        }
    }
}

/// The seed count of a runnable `M` with first index `n0`: at least `R`,
/// and past every target where the leader vanishes.
pub(super) fn seed_count(m: &DOp, n0: usize) -> usize {
    let order = m.max_shift();
    let leader = m.0.get(&(order, 0)).cloned().unwrap_or_default();
    zero_bases(&leader, n0)
        .into_iter()
        .map(|base| base + order + 1 - n0)
        .max()
        .unwrap_or(0)
        .max(order)
}

/// The content in `n` of an operator: the gcd of the contents of all its
/// coefficients (zero for the zero operator).
fn operator_index_content(op: &DOp) -> Vec<Q> {
    let mut gcd = crate::polynomial::Polynomial::<Q>::zero();
    for c in op.0.values() {
        let (content, _) = content_and_primitive(c);
        gcd = gcd.gcd(&crate::polynomial::Polynomial::new(content));
    }
    vec_poly::trim(gcd.coeffs())
}

/// A runnable rule with seeds for rows `Q_n`, `n >= first_index`, from an
/// operator `M` (base form, `sum m_(i,d)(n, x) D^d S^i`, top shift `R`) that
/// is **proved** to annihilate the rows, `(M Q)_n = 0`, at every base
/// `n >= proved_from_base` except the bases listed in `unproved_bases`.
/// `rows(count)` returns `Q_(first_index), ..., Q_(first_index + count - 1)`
/// computed independently (from the definition of `Q`).
///
/// The normalisation:
///
/// 1. `M = c(n) x^e M_0` with `M_0` = `ore_algebra::normalize_annihilator(M)`
///    (content `c(n)` and power `x^e` removed); `(M_0 Q)_n = 0` wherever
///    `(M Q)_n = 0` and `c(n) != 0` (`Q[x]` has no zero divisors), so the
///    integer roots of `c` join the unproved bases.
/// 2. At every unproved base `b >= first_index` (and every base below
///    `proved_from_base`) whose leader `m_R(b, .)` is not the zero
///    polynomial, `(M_0 Q)_b` is evaluated exactly on the rows; where it is
///    not zero, `M_0` is multiplied on the left by `(n - b)`, so the leader
///    vanishes at `b` and the target `b + R` becomes a seed.  Call the result
///    `M_1`; it vanishes on `Q` at every base `n >= first_index` at which
///    its leader is not the zero polynomial.
/// 3. The seeds are `Q_(first_index), ..., Q_(first_index + K - 1)` with `K`
///    at least `R` and past every target where the leader of `M_1`
///    vanishes identically (the integer roots of its content in `n`, found
///    exactly).  At every later target `t` the leader `m_R(t - R, .)` is a
///    nonzero polynomial and `(M_1 Q)_(t-R) = 0`, so `Q_t` is the unique
///    row the rule generates from the previous ones; by induction the rule
///    with these seeds generates exactly `Q`.
///
/// A replay of the rule against `rows` on eight further rows is a
/// self-check (evidence only; the argument above is the proof).
///
/// # Errors
///
/// [`RecurrenceError::Unsupported`] when `M` is not runnable (its top shift
/// column is not a single term without derivatives),
/// [`RecurrenceError::LimitExceeded`] when more than
/// [`MAX_IMAGE_RECURRENCE_ROWS`] rows are needed, errors of `rows`, and
/// [`RecurrenceError::VerificationFailed`] when the self-check fails.
pub(in crate::recurrence) fn runnable_record_from_annihilator(
    annihilator: &DOp,
    first_index: usize,
    proved_from_base: usize,
    unproved_bases: &[usize],
    rows: impl Fn(usize) -> Result<Vec<Vec<Q>>, RecurrenceError>,
) -> Result<RecurrenceJson, RecurrenceError> {
    let n0 = first_index;
    if !ore_algebra::is_runnable(annihilator) {
        return Err(RecurrenceError::Unsupported(
            "the annihilator is not runnable (its top shift column is not a single term \
             without derivatives)"
                .into(),
        ));
    }
    let content = operator_index_content(annihilator);
    let normalized = ore_algebra::normalize_annihilator(annihilator);
    let order = normalized.max_shift();
    let leader = |op: &DOp| op.0.get(&(op.max_shift(), 0)).cloned().unwrap_or_default();
    let leader_zero_bases = zero_bases(&leader(&normalized), n0);
    let mut checks: Vec<usize> = (n0..proved_from_base.max(n0)).collect();
    checks.extend(unproved_bases.iter().copied().filter(|b| *b >= n0));
    checks.extend(roots_in_domain(&content, n0, None));
    checks.sort_unstable();
    checks.dedup();
    checks.retain(|b| !leader_zero_bases.contains(b));
    let last_base = checks
        .iter()
        .chain(&leader_zero_bases)
        .max()
        .copied()
        .unwrap_or(n0);
    let needed = last_base + order + 1 - n0 + 8;
    require_row_budget(needed)?;
    let computed = rows(needed)?;
    if computed.len() < needed {
        return Err(RecurrenceError::VerificationFailed(format!(
            "{} rows were computed, {needed} are needed",
            computed.len()
        )));
    }
    let mut factor = vec![Q::one()];
    for base in &checks {
        let value = normalized.apply(&computed, n0, *base).ok_or_else(|| {
            RecurrenceError::VerificationFailed(format!("rows missing at base {base}"))
        })?;
        if !vec_poly::is_zero(&value) {
            factor = vec_poly::mul(&factor, &[-q_usize(*base), Q::one()]);
        }
    }
    let adjusted = if factor.len() > 1 {
        normalized.left_multiply(&Biv::new(factor.iter().map(|c| vec![c.clone()]).collect()))
    } else {
        normalized
    };
    let seeds = zero_bases(&leader(&adjusted), n0)
        .into_iter()
        .map(|base| base + order + 1 - n0)
        .max()
        .unwrap_or(0)
        .max(order);
    let (rule, _) = successor_rule(&adjusted)?;
    let record = RecurrenceJson::from_recurrence_rational(&rule, n0, &computed[..seeds], None);
    let replay = rule
        .generate_rows_rational(&computed[..seeds], n0, needed)
        .map_err(|error| {
            RecurrenceError::VerificationFailed(format!("the rule cannot be replayed: {error}"))
        })?;
    if replay
        .iter()
        .zip(&computed)
        .any(|(generated, expected)| vec_poly::trim(generated) != vec_poly::trim(expected))
    {
        return Err(RecurrenceError::VerificationFailed(
            "the normalized rule does not replay the rows".into(),
        ));
    }
    Ok(record)
}

/// The exceptional bases of the identity (module documentation, step 3).
pub(super) fn exceptional_bases(
    multiplier: &Biv,
    base_factor: &[Q],
    n0: usize,
    first_base: usize,
) -> Vec<usize> {
    let mut bases: Vec<usize> = (n0..first_base.max(n0)).collect();
    bases.extend(zero_bases(multiplier, n0));
    bases.extend(roots_in_domain(base_factor, n0, None));
    bases.sort_unstable();
    bases.dedup();
    bases
}

fn require_row_budget(count: usize) -> Result<(), RecurrenceError> {
    if count > MAX_IMAGE_RECURRENCE_ROWS {
        return Err(RecurrenceError::LimitExceeded(format!(
            "the image certificate needs {count} rows, more than {MAX_IMAGE_RECURRENCE_ROWS}"
        )));
    }
    Ok(())
}

/// The multiplier side condition on `[-S, -1]` (module documentation).
pub(super) fn check_multiplier_bound(
    prepared: &PreparedImage,
    bounds: PullbackBounds,
) -> Result<(), RecurrenceError> {
    if let Some((a, _)) = &prepared.multiplier {
        let s = bounds.multiplier_euler_bound() as i64;
        require_no_integer_roots_from(a, -s, "multiplier numerator A(k)")?;
    }
    Ok(())
}

/// The shifted parent relation `L'` and its first base `b_L'` (saturated
/// at zero).
pub(super) fn shifted_relation(parent: &Parent, shift: i64) -> (DOp, usize) {
    let first_base = (parent.first_base as i64 - shift).max(0) as usize;
    (parent.relation.shift_index(shift), first_base)
}

/// The row checks of a candidate at its exceptional bases: returns
/// `(image rows, exceptional bases)` or the failure.
pub(super) fn check_rows(
    certificate: &ImageAnnihilatorCertificate,
    base_factor: &[Q],
    prepared: &PreparedImage,
    parent: &Parent,
    n0: usize,
    first_base: usize,
    replay_rows: usize,
) -> Result<(Vec<Vec<Q>>, Vec<usize>), RecurrenceError> {
    let m = &certificate.annihilator;
    let order = m.max_shift();
    let seeds = seed_count(m, n0);
    let exceptional = exceptional_bases(&certificate.multiplier, base_factor, n0, first_base);
    let needed = (seeds + replay_rows).max(exceptional.last().map_or(0, |b| b + order + 1 - n0));
    require_row_budget(needed)?;
    let images = prepared.image_rows(parent, n0, needed)?;
    for base in &exceptional {
        let value = m.apply(&images, n0, *base).ok_or_else(|| {
            RecurrenceError::VerificationFailed(format!("rows missing at base {base}"))
        })?;
        if !vec_poly::is_zero(&value) {
            return Err(RecurrenceError::VerificationFailed(format!(
                "the image rule fails at the exceptional base {base}"
            )));
        }
    }
    Ok((images, exceptional))
}
