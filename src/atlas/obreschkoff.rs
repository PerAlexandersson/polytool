//! Obreschkoff combinations `Q_n = P_n + P_(n+1)` of a node whose
//! consecutive rows interlace.
//!
//! # Theorem
//!
//! Hermite--Kakeya--Obreschkoff (N. Obreschkoff, *Verteilung und Berechnung
//! der Nullstellen reeller Polynome*, 1963; see also S. Fisk, *Polynomials,
//! roots, and interlacing*, Theorem 1.43): for real-rooted `f`, `g` whose
//! zeros interlace (their degrees differ by at most one), every real linear
//! combination `a f + b g` is zero or real-rooted.  Common zeros are
//! allowed: write `f = h f1`, `g = h g1` with `h = gcd(f, g)`; then `f1`,
//! `g1` interlace strictly, `a f1 + b g1` is real-rooted, and so is
//! `h (a f1 + b g1)`.  A zero row gives `Q_n = P_(n+1)` or `Q_n = P_n`.
//! Applied with `f = P_n`, `g = P_(n+1)`: if every `P_n`, `n >= n0`, is zero
//! or real-rooted and `P_(m-1)`, `P_m` interlace for every `m >= n0 + 1`,
//! then every `Q_n`, `n >= n0`, is zero or real-rooted.  No root region is
//! claimed beyond the whole line: the leading coefficients of `P_n` and
//! `P_(n+1)` may have opposite signs.
//!
//! # Identity
//!
//! For a first-order parent `f(m) P_m = sum_d c_d(m, x) D^d P_(m-1)`, valid
//! at every target `m >= p0 + K` (`K` identifying seeds, the leader decided
//! nonzero there), the next row is a Weyl image of the current one:
//! `P_(n+1) = (1/f(n+1)) sum_d c_d(n+1, x) D^d P_n` for every
//! `n >= p0 + K - 1`, so
//!
//! ```text
//! Q_n = T_n(P_n),   T_n = (1/f(n+1)) (f(n+1) + sum_d c_d(n+1, x) D^d).
//! ```
//!
//! The certificate is a preserver-image certification of the candidate
//! `Q_n = T_n(P_n)` (operator [`PreserverOperator::WeylOperator`], index
//! shift 0, no gauge) whose identity is certified for every `n >= n0` by the
//! exact Ore-algebra division of the preserver machinery, for the rule of
//! `Q` computed by [`crate::recurrence::image_recurrence_of_operator`].
//! [`ObreschkoffCertificate::analyze`] replays it, recomputes `T_n` from
//! the recorded parent and requires `n0 >= p0 + K - 1`.  The premises are
//! the parent's real-rootedness from `n0` on and its interlacing from
//! `n0 + 1` on; both must be proved properties of the parent node.

use super::analysis::PremiseRequirement;
use super::analysis::{canonical_definition, preserver_input, CertificateAnalysis};
use super::error::AtlasError;
use super::rows::whole_line;
use super::types::{AtlasClaim, AtlasPremiseRole};
use crate::recurrence::{certify_preserver_image_lift, ImageRecurrenceOptions};
use crate::recurrence::{
    format_rational_coeff, prepare_rule, BivarPoly, BivarPolyJson, InterlacingOrientation,
    OperatorTermJson, PreserverImageCandidate, PreserverImageCertification, PreserverImageOptions,
    PreserverOperator, RecurrenceJson, RecurrenceSign, RowInterlacing, VeroneseSectionGauge,
};
use crate::BigRational;
use num_traits::Zero;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// The certificate of an Obreschkoff combination (module documentation).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObreschkoffCertificate {
    /// The parent `P` (canonical definition).
    pub parent: RecurrenceJson,
    /// The first index `n0` of `Q`.
    pub first_index: usize,
    /// The certified identity `Q_n = T_n(P_n)`, `n >= n0`.
    pub identity: PreserverImageCertification,
}

fn biv_json(polynomial: &BivarPoly) -> BivarPolyJson {
    BivarPolyJson {
        coeffs: polynomial
            .coeffs
            .iter()
            .map(|row| row.iter().map(format_rational_coeff).collect())
            .collect(),
    }
}

/// The operator `T_n` of the module documentation for a first-order parent,
/// with `p0 + K - 1` (the first `n` for which `P_(n+1) = T_n'(P_n)` holds).
///
/// # Errors
///
/// [`AtlasError::Unsupported`] when the parent is not a runnable first-order
/// rule.
pub fn obreschkoff_operator(
    parent: &RecurrenceJson,
) -> Result<(PreserverOperator, usize), AtlasError> {
    let (rule, first_index, rows) = parent.to_recurrence_parts()?;
    let prepared = prepare_rule(&rule, &rows, first_index, "parent")?;
    if prepared.leader_shift != 1
        || prepared
            .rule
            .terms
            .iter()
            .any(|t| t.offset != 1 || t.sign != RecurrenceSign::None)
    {
        return Err(AtlasError::Unsupported(
            "an Obreschkoff combination needs a first-order parent rule".into(),
        ));
    }
    let shifted = prepared.rule.shift_index(1);
    let denominator = match &shifted.denominator {
        None => {
            BivarPoly::from_coefficients_normalized(vec![vec![BigRational::from_integer(1.into())]])
        }
        Some(d) => d.clone(),
    };
    if denominator
        .coeffs
        .iter()
        .any(|row| row.iter().skip(1).any(|c| !c.is_zero()))
    {
        return Err(AtlasError::Unsupported(
            "the parent leader depends on x".into(),
        ));
    }
    let leader: Vec<BigRational> = denominator
        .coeffs
        .iter()
        .map(|row| row.first().cloned().unwrap_or_else(BigRational::zero))
        .collect();
    let mut by_order: BTreeMap<usize, Vec<Vec<BigRational>>> = BTreeMap::new();
    let mut add = |order: usize, polynomial: &BivarPoly| {
        let entry = by_order.entry(order).or_default();
        for (i, row) in polynomial.coeffs.iter().enumerate() {
            if entry.len() <= i {
                entry.resize(i + 1, Vec::new());
            }
            for (j, c) in row.iter().enumerate() {
                if entry[i].len() <= j {
                    entry[i].resize(j + 1, BigRational::zero());
                }
                entry[i][j] += c;
            }
        }
    };
    add(0, &denominator);
    for term in &shifted.terms {
        add(term.deriv_order, &term.coeff);
    }
    let terms = by_order
        .into_iter()
        .map(|(deriv_order, coeffs)| OperatorTermJson {
            deriv_order,
            coeff: biv_json(&BivarPoly::from_coefficients_normalized(coeffs)),
        })
        .collect();
    let operator = PreserverOperator::WeylOperator {
        terms,
        denominator: leader.iter().map(format_rational_coeff).collect(),
    };
    Ok((operator, prepared.first_index + prepared.seeds.len() - 1))
}

/// The candidate `Q_n = T_n(P_n)` for a parent.
///
/// # Errors
///
/// Those of [`obreschkoff_operator`].
pub fn obreschkoff_candidate(
    parent: &RecurrenceJson,
) -> Result<(PreserverImageCandidate, usize), AtlasError> {
    let parent = canonical_definition(parent)?;
    let (operator, start) = obreschkoff_operator(&parent)?;
    Ok((
        PreserverImageCandidate {
            parent_name: "obreschkoff".into(),
            parent,
            index_shift: 0,
            gauge: VeroneseSectionGauge::None,
            operator,
        },
        start,
    ))
}

/// Compute the rule of `Q_n = P_n + P_(n+1)`, `n >= first_index`, and
/// certify the identity (module documentation).  `first_index` must be at
/// least `p0 + K - 1`.
///
/// # Errors
///
/// [`AtlasError`] when the parent is not first order, the image rule is not
/// found within the budget, or the identity is not certified.
pub fn certify_obreschkoff_combination(
    parent: &RecurrenceJson,
    first_index: usize,
    image_options: &ImageRecurrenceOptions,
) -> Result<ObreschkoffCertificate, AtlasError> {
    let (candidate, start) = obreschkoff_candidate(parent)?;
    if first_index < start {
        return Err(AtlasError::Unsupported(format!(
            "the combination starts at {first_index}, before the parent rule holds at the next \
             index ({start})"
        )));
    }
    let image = candidate.image_recurrence(first_index, image_options)?;
    let (rule, n0, seeds) = image.rule.to_recurrence_parts()?;
    let rows = rule
        .generate_rows_rational(&seeds, n0, 16.max(seeds.len()))
        .map_err(|e| AtlasError::Unsupported(format!("the image rule does not replay: {e}")))?;
    let identity = certify_preserver_image_lift(
        &rule,
        &rows,
        n0,
        &candidate,
        &PreserverImageOptions {
            try_root_domain: false,
            ..PreserverImageOptions::default()
        },
    )?;
    let certificate = ObreschkoffCertificate {
        parent: candidate.parent.clone(),
        first_index: n0,
        identity,
    };
    certificate.analyze()?;
    Ok(certificate)
}

impl ObreschkoffCertificate {
    /// Replay the certificate and read off its claim and premises (module
    /// documentation).
    ///
    /// # Errors
    ///
    /// [`AtlasError::VerificationFailed`] when the identity does not replay,
    /// is not certified, or does not belong to this parent and operator.
    pub fn analyze(&self) -> Result<CertificateAnalysis, AtlasError> {
        self.identity.verify()?;
        let input = preserver_input(&self.identity)?.clone();
        let (operator, start) = obreschkoff_operator(&self.parent)?;
        let candidate = &self.identity.candidate;
        if candidate.parent != self.parent
            || candidate.operator != operator
            || candidate.index_shift != 0
            || candidate.gauge != VeroneseSectionGauge::None
            || self.identity.input_first_index != self.first_index
            || input.first_index != self.first_index
            || self.first_index < start
        {
            return Err(AtlasError::VerificationFailed(
                "the identity does not belong to the Obreschkoff combination of the parent".into(),
            ));
        }
        Ok(CertificateAnalysis {
            input,
            claims: vec![AtlasClaim::RealRooted {
                from_index: self.first_index,
                region: whole_line(),
            }],
            premises: vec![
                PremiseRequirement {
                    definition: self.parent.clone(),
                    role: AtlasPremiseRole::Sequence,
                    uses: AtlasClaim::RealRooted {
                        from_index: self.first_index,
                        region: whole_line(),
                    },
                },
                PremiseRequirement {
                    definition: self.parent.clone(),
                    role: AtlasPremiseRole::Sequence,
                    uses: AtlasClaim::Interlacing {
                        from_index: self.first_index + 1,
                        kind: RowInterlacing::Weak,
                        orientation: InterlacingOrientation::Unordered,
                    },
                },
            ],
            equal_to: None,
            theorem: "Hermite--Kakeya--Obreschkoff: real-rooted f, g with interlacing zeros give \
                      a zero or real-rooted a f + b g for all real a, b (Obreschkoff 1963; Fisk, \
                      Polynomials, roots, and interlacing, Thm 1.43)"
                .into(),
            lean_route: None,
            lean_declarations: Vec::new(),
        })
    }
}
