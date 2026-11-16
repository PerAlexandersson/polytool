//! Exact index subsequences of stationary recurrences, sharing the Newton-sum kernel.
use super::*;

/// An all-index subsequence identity, without a real-rootedness assertion.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ConstantCoefficientIndexSubsequenceCertificate {
    /// Original recurrence record, including every prescribed row.
    pub source: RecurrenceJson,
    /// Section step, between one and sixteen.
    pub stride: usize,
    /// Residue between zero and `stride-1`.
    pub residue: usize,
    /// Source index of section row zero; after every exceptional seed row.
    pub first_source_index: usize,
    /// Recurrence of `Q_m=P_(first_source_index+stride*m)`, starting at zero.
    pub subsequence: RecurrenceJson,
    /// Ascending powers of `z`, then `x`, of `prod_i (z-alpha_i^stride)`;
    /// rational coefficients are exact strings. `alpha_i` are characteristic roots.
    pub characteristic_polynomial: Vec<Vec<String>>,
}

impl ConstantCoefficientIndexSubsequenceCertificate {
    /// Recompute the Newton sums, source anchor and every section seed.
    ///
    /// # Errors
    /// Returns a typed verification error when any recorded field changes.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        let expected = certify_constant_coefficient_index_subsequence(
            &self.source,
            self.stride,
            self.residue,
        )?;
        if &expected != self {
            return Err(RecurrenceError::VerificationFailed(
                "stationary index-subsequence certificate changed".into(),
            ));
        }
        Ok(())
    }
}

/// Certify an index section of a stationary polynomial recurrence exactly.
///
/// For a source with `K` seeds and order `r`, start at its last `r` seeds:
/// `N=first_index+K-r`. The returned family is
/// `Q_m=P_(N+residue+stride*m)`. It has characteristic polynomial
/// `prod_i(z-alpha_i^stride)`, calculated by the same exact Newton sums as
/// index-section lifts. Cayley--Hamilton applied to the stride-th power of
/// the companion matrix proves the identity at every index; seed agreement
/// selects the source's solution. No finite fit is used. This construction
/// does not assert real-rootedness, minimal recurrence order or interlacing.
///
/// `stride` is limited to 16 and source order to 24 to bound construction
/// work. Coefficients must be independent of the index, with no derivatives,
/// alternating signs or forcing. A constant nonzero row scaling is allowed.
///
/// # Errors
/// Returns a typed error for unsupported shapes, invalid residues, too few
/// seeds, arithmetic overflow, construction bounds or row-generation failure.
pub fn certify_constant_coefficient_index_subsequence(
    source: &RecurrenceJson,
    stride: usize,
    residue: usize,
) -> Result<ConstantCoefficientIndexSubsequenceCertificate, RecurrenceError> {
    if stride == 0 || stride > 16 || residue >= stride {
        return Err(RecurrenceError::InvalidInput(
            "index subsequence requires 1 <= stride <= 16 and residue < stride".into(),
        ));
    }
    let (rule, first, seeds) = source.to_recurrence_parts()?;
    let (normalized, _) = rule_with_x_free_left_hand_factor(&rule)?;
    let multipliers =
        constant_coefficient_multipliers(&normalized).map_err(RecurrenceError::Unsupported)?;
    let order = multipliers.len();
    if order > 24 {
        return Err(RecurrenceError::LimitExceeded(
            "index-subsequence order exceeds 24".into(),
        ));
    }
    if seeds.len() < order {
        return Err(RecurrenceError::InvalidInput(
            "too few section source seeds".into(),
        ));
    }
    let prefix = seeds.len() - order;
    let anchor = prefix
        .checked_add(residue)
        .ok_or_else(|| RecurrenceError::InvalidInput("section index overflow".into()))?;
    let last = order
        .checked_sub(1)
        .and_then(|r| r.checked_mul(stride))
        .and_then(|n| n.checked_add(anchor))
        .and_then(|n| n.checked_add(1))
        .ok_or_else(|| RecurrenceError::InvalidInput("section index overflow".into()))?;
    let generated = normalized
        .generate_rows_rational(&seeds, first, last.max(seeds.len()))
        .map_err(|e| RecurrenceError::Unsupported(e.to_string()))?;
    let section_seeds: Vec<Vec<Q>> = (0..order)
        .map(|i| generated[anchor + stride * i].clone())
        .collect();
    let characteristic =
        power_characteristic_polynomial(&characteristic_polynomial(&multipliers), stride);
    let section_rule = super::super::super::veronese_search::parent_recurrence_from_characteristic(
        &characteristic,
    );
    let subsequence =
        RecurrenceJson::from_recurrence_rational(&section_rule, 0, &section_seeds, None);
    let first_source_index = first
        .checked_add(anchor)
        .ok_or_else(|| RecurrenceError::InvalidInput("section index overflow".into()))?;
    Ok(ConstantCoefficientIndexSubsequenceCertificate {
        source: source.clone(),
        stride,
        residue,
        first_source_index,
        subsequence,
        characteristic_polynomial: bivariate_to_coefficients(&characteristic)
            .iter()
            .map(|r| r.iter().map(ToString::to_string).collect())
            .collect(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn even_and_odd_replay_and_exact_generated_rows() {
        let entry = crate::oeis::by_id("A158909").unwrap();
        let (rule, seeds) = entry.recurrence_parts().unwrap();
        let first = entry.recurrence_first_index;
        let source = RecurrenceJson::from_recurrence_rational(&rule, first, &seeds, None);
        let original = rule.generate_rows_rational(&seeds, first, 24).unwrap();
        for residue in 0..2 {
            let c = certify_constant_coefficient_index_subsequence(&source, 2, residue).unwrap();
            c.verify().unwrap();
            let (r, f, s) = c.subsequence.to_recurrence_parts().unwrap();
            let generated = r.generate_rows_rational(&s, f, 10).unwrap();
            for i in 0..10 {
                assert_eq!(generated[i], original[c.first_source_index - first + 2 * i]);
            }
            let mut changed = c;
            changed.first_source_index += 1;
            assert!(changed.verify().is_err());
        }
        assert!(certify_constant_coefficient_index_subsequence(&source, 0, 0).is_err());
        assert!(certify_constant_coefficient_index_subsequence(&source, 2, 2).is_err());
        let mut dependent = source;
        dependent.recurrence.terms[0]
            .coeff
            .coeffs
            .push(vec!["1".into()]);
        assert!(certify_constant_coefficient_index_subsequence(&dependent, 2, 0).is_err());
    }
}
