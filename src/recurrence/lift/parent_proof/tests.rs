//! Tests for the enclosing module.

#[cfg(test)]
mod catalog_flag_tests {
    use super::super::*;

    fn q(value: i64) -> Q {
        Q::from_integer(value.into())
    }

    /// The catalog rule, identifying seeds and first index of a bundled entry.
    fn catalog_parent(id: &str) -> (Recurrence, Vec<Vec<Q>>, usize) {
        let entry = crate::oeis::by_id(id).unwrap();
        let (rule, rows) = entry.recurrence_parts().unwrap();
        let count = rule.generation_initial_count(entry.recurrence_first_index, rows.len());
        (rule, rows[..count].to_vec(), entry.recurrence_first_index)
    }

    fn prove(rule: &Recurrence, seeds: &[Vec<Q>], first_index: usize) -> VeroneseParentProof {
        let record = ParentRecurrenceRecord {
            rule,
            seeds,
            first_index,
        };
        prove_parent_real_rootedness(&record, first_index, &[], false).unwrap()
    }

    fn catalog_attempts(proof: &VeroneseParentProof) -> Vec<&VeroneseParentProofAttempt> {
        proof
            .attempts
            .iter()
            .filter(|attempt| attempt.variant.starts_with("oeis_catalog:"))
            .collect()
    }

    /// The Jabotinsky-type rows A039692 (`P_n = (x + n) P_(n-1) + n (n - 1 +
    /// 2x) P_(n-2)`, `P_0 = 1`, `P_1 = 3 + x`) are flagged `lean_proved`; no
    /// polytool route proves them, so the verdict is conditional on the flag,
    /// with one-signed rows proved.
    #[test]
    fn a_flagged_catalog_parent_is_conditional_on_the_flag() {
        let (rule, seeds, first_index) = catalog_parent("A039692");
        assert_eq!(seeds, vec![vec![q(1)], vec![q(3), q(1)]]);
        let proof = prove(&rule, &seeds, first_index);
        assert_eq!(
            proof.status,
            VeroneseParentProofStatus::ConditionalOnAttestation,
            "{proof:#?}"
        );
        assert_eq!(proof.route, Some(VeroneseParentProofRoute::Attestation));
        assert_eq!(proof.conditional_on, ["oeis_catalog:A039692:lean_proved"]);
        assert!(proof.statement.contains("not proved here"));
        let attempt = proof.attempts.last().unwrap();
        assert!(attempt.success);
        assert!(attempt.detail.contains("lean_proved, as of 2026-10-03"));
        assert!(attempt.detail.contains("coefficients of one sign"));
        // A caller attestation of the same record is applied first.
        let attestation = VeroneseParentAttestation {
            label: "caller".into(),
            recurrence: RecurrenceJson::from_recurrence_rational(&rule, first_index, &seeds, None),
            first_index,
            statement: "attested".into(),
            source: "test".into(),
        };
        let record = ParentRecurrenceRecord {
            rule: &rule,
            seeds: &seeds,
            first_index,
        };
        let proof =
            prove_parent_real_rootedness(&record, first_index, &[attestation], false).unwrap();
        assert_eq!(proof.conditional_on, ["caller"]);
    }

    /// A157153 is not flagged (its Lean proof is open), so no catalog
    /// attestation applies.
    #[test]
    fn an_unflagged_catalog_parent_gets_no_catalog_attestation() {
        assert!(crate::oeis::by_id("A157153")
            .unwrap()
            .real_rootedness()
            .is_none());
        let (rule, seeds, first_index) = catalog_parent("A157153");
        let proof = prove(&rule, &seeds, first_index);
        assert_ne!(
            proof.status,
            VeroneseParentProofStatus::ConditionalOnAttestation,
            "{proof:#?}"
        );
        // Flagged entries with the same seeds differ in the rule.
        assert!(catalog_attempts(&proof).is_empty(), "{proof:#?}");
        assert!(proof
            .attempts
            .iter()
            .filter(|attempt| attempt.route == VeroneseParentProofRoute::Attestation)
            .all(|attempt| !attempt.success && !attempt.detail.contains("A157153")));
    }

    /// Matching is exact: one changed rule coefficient or seed and the flag no
    /// longer applies.
    #[test]
    fn a_parent_differing_in_one_coefficient_or_seed_does_not_match() {
        let (rule, seeds, first_index) = catalog_parent("A039692");
        let mut changed_rule = rule.clone();
        changed_rule.terms[1].coeff.coeffs[1][1] = q(3);
        let proof = prove(&changed_rule, &seeds, first_index);
        assert_ne!(
            proof.status,
            VeroneseParentProofStatus::ConditionalOnAttestation
        );
        let attempts = catalog_attempts(&proof);
        assert!(attempts.is_empty(), "{proof:#?}");
        let summary = proof.attempts.last().unwrap();
        assert_eq!(summary.variant, "oeis_catalog");
        assert!(!summary.success);
        assert!(summary.detail.contains("oeis_catalog:A039692:lean_proved"));

        let mut changed_seeds = seeds.clone();
        changed_seeds[1] = vec![q(2)];
        let proof = prove(&rule, &changed_seeds, first_index);
        assert_ne!(
            proof.status,
            VeroneseParentProofStatus::ConditionalOnAttestation
        );
        assert!(catalog_attempts(&proof).is_empty());

        // The same rule and seeds at another first index do not match either.
        let proof = prove(&rule, &seeds, first_index + 1);
        assert!(catalog_attempts(&proof).is_empty());
    }

    /// The flag gives real-rootedness only: a flagged parent that no route
    /// proves and whose rows are not shown one-signed by the weight test
    /// (A008303) stays open.
    #[test]
    fn a_flagged_parent_without_a_one_sign_proof_stays_open() {
        let (rule, seeds, first_index) = catalog_parent("A008303");
        assert!(crate::oeis::by_id("A008303")
            .unwrap()
            .real_rootedness()
            .is_some());
        let proof = prove(&rule, &seeds, first_index);
        assert_eq!(proof.status, VeroneseParentProofStatus::NotProved);
        let attempts = catalog_attempts(&proof);
        assert_eq!(attempts.len(), 1);
        assert!(attempts[0]
            .detail
            .contains("one-signed coefficients are not shown: the weight"));
    }

    fn term(offset: usize, derivative: usize, coefficients: Vec<Vec<Q>>) -> RecurrenceTerm {
        RecurrenceTerm::new(
            offset,
            derivative,
            BivarPoly::from_coefficients_normalized(coefficients),
        )
    }

    #[test]
    fn one_signed_rows_need_nonnegative_combined_weights() {
        // Eulerian: P_n = (1 - x + n x) P_(n-1) + (x - x^2) P'_(n-1), P_1 = 1.
        // The weight n - 1 - s of x^s into x^(s+1) needs the degree bound
        // deg P_m <= m - 1.
        let eulerian = Recurrence {
            terms: vec![
                term(1, 0, vec![vec![q(1), q(-1)], vec![q(0), q(1)]]),
                term(1, 1, vec![vec![q(0), q(1), q(-1)]]),
            ],
            denominator: None,
            inhomogeneous: None,
        };
        let shown = one_signed_rows_for_all_indices(&eulerian, 1, &[vec![q(1)]]).unwrap();
        assert!(shown.ends_with("deg P_m <= -1 + 1 m"), "{shown}");
        // Nonpositive seeds of a homogeneous rule: the rows are nonpositive.
        assert!(one_signed_rows_for_all_indices(&eulerian, 1, &[vec![q(-1)]]).is_ok());
        assert!(
            one_signed_rows_for_all_indices(&eulerian, 1, &[vec![q(1), q(-1)]])
                .unwrap_err()
                .contains("both signs")
        );
        // A coefficient x - 2 x^2 makes the weight n - 2 - 2 s, negative at
        // s = n - 2 (the top coefficient of P_(n-1)).
        let mut too_negative = eulerian.clone();
        too_negative.terms[1] = term(1, 1, vec![vec![q(0), q(1), q(-2)]]);
        assert!(
            one_signed_rows_for_all_indices(&too_negative, 1, &[vec![q(1)]])
                .unwrap_err()
                .contains("not shown to be >= 0")
        );
        // A leader 3 - n vanishes at n = 3.
        let mut leader = eulerian.clone();
        leader.denominator = Some(BivarPoly::from_coefficients_normalized(vec![
            vec![q(3)],
            vec![q(-1)],
        ]));
        assert!(one_signed_rows_for_all_indices(&leader, 1, &[vec![q(1)]])
            .unwrap_err()
            .contains("leader is not positive"));
        // (-1)^n factors and negative inhomogeneous parts are refused.
        let mut alternating = eulerian.clone();
        alternating.terms[0].sign = RecurrenceSign::AlternatingN;
        assert!(one_signed_rows_for_all_indices(&alternating, 1, &[vec![q(1)]]).is_err());
        let mut inhomogeneous = eulerian.clone();
        inhomogeneous.inhomogeneous = Some(BivarPoly::from_coefficients_normalized(vec![
            vec![q(5)],
            vec![q(-1)],
        ]));
        assert!(
            one_signed_rows_for_all_indices(&inhomogeneous, 1, &[vec![q(1)]])
                .unwrap_err()
                .contains("inhomogeneous")
        );
        // Zero seeds of a homogeneous rule give zero rows.
        assert!(one_signed_rows_for_all_indices(&eulerian, 1, &[vec![]]).is_ok());
    }
}

#[cfg(test)]
mod interlacing_tests {
    use super::super::*;

    fn q(value: i64) -> Q {
        Q::from_integer(value.into())
    }

    fn prove(rule: &Recurrence, seeds: &[Vec<Q>], first_index: usize) -> VeroneseParentProof {
        let record = ParentRecurrenceRecord {
            rule,
            seeds,
            first_index,
        };
        prove_parent_real_rootedness(&record, first_index, &[], false).unwrap()
    }

    /// `P_n = P_(n-1) + x P_(n-2)` with the canonical rows `P_0 = 1`,
    /// `P_1 = 1 + x` of the fixed-delay window theorem: route 2 proves the
    /// parent and records `P_(m-1) << P_m` from `m = 1`.
    #[test]
    fn the_window_theorem_route_records_directed_interlacing() {
        let rule = Recurrence {
            terms: vec![
                RecurrenceTerm::new(
                    1,
                    0,
                    BivarPoly::from_coefficients_normalized(vec![vec![q(1)]]),
                ),
                RecurrenceTerm::new(
                    2,
                    0,
                    BivarPoly::from_coefficients_normalized(vec![vec![q(0), q(1)]]),
                ),
            ],
            denominator: None,
            inhomogeneous: None,
        };
        let proof = prove(&rule, &[vec![q(1)], vec![q(1), q(1)]], 0);
        assert_eq!(proof.route, Some(VeroneseParentProofRoute::PaperTheorem));
        let interlacing = proof.interlacing.as_ref().expect("window interlacing");
        assert_eq!(interlacing.kind, RowInterlacing::Weak);
        assert_eq!(
            interlacing.orientation,
            InterlacingOrientation::PreviousInterlacesNext
        );
        assert_eq!(interlacing.from_index, 1);
        assert!(interlacing
            .justification
            .theorem
            .starts_with("fixed_delay_predecessor_window"));
        // Exact evidence on the first rows.
        let rows = rule
            .generate_rows_rational(&[vec![q(1)], vec![q(1), q(1)]], 0, 10)
            .unwrap();
        for pair in rows.windows(2) {
            let integer = |row: &[Q]| row.iter().map(|c| c.to_integer()).collect::<Vec<_>>();
            assert_eq!(
                crate::check_weak_interlacing_bigint_coeffs(&integer(&pair[0]), &integer(&pair[1])),
                Some(true)
            );
        }
    }

    /// The Eulerian rule `P_n = (1 + (n - 1) x) P_(n-1) + x (1 - x) P'_(n-1)`
    /// is proved by Liu--Wang (route 5) without the root-domain route; the
    /// theorem's weak interlacing is recorded.  A periodic product records
    /// none.
    #[test]
    fn theorem_routes_record_their_conclusion_and_products_record_none() {
        let entry = crate::oeis::by_id("A008292").unwrap();
        let (rule, rows) = entry.recurrence_parts().unwrap();
        let count = rule.generation_initial_count(entry.recurrence_first_index, rows.len());
        let proof = prove(&rule, &rows[..count], entry.recurrence_first_index);
        assert_eq!(proof.route, Some(VeroneseParentProofRoute::LiuWang));
        let interlacing = proof.interlacing.as_ref().expect("Liu--Wang interlacing");
        assert_eq!(interlacing.kind, RowInterlacing::Weak);
        assert_eq!(
            interlacing.orientation,
            InterlacingOrientation::PreviousInterlacesNext
        );
        assert!(interlacing.justification.theorem.starts_with("liu_wang: "));
        let json = serde_json::to_value(&proof).unwrap();
        assert_eq!(
            json["interlacing"]["orientation"],
            "previous_interlaces_next"
        );
        let back: VeroneseParentProof = serde_json::from_value(json).unwrap();
        assert_eq!(back, proof);

        // (1 + x)^n: a periodic product, no interlacing recorded.
        let product = Recurrence {
            terms: vec![RecurrenceTerm::new(
                1,
                0,
                BivarPoly::from_coefficients_normalized(vec![vec![q(1), q(1)]]),
            )],
            denominator: None,
            inhomogeneous: None,
        };
        let proof = prove(&product, &[vec![q(1)]], 0);
        assert_eq!(proof.route, Some(VeroneseParentProofRoute::PeriodicProduct));
        assert!(proof.interlacing.is_none());
        assert!(serde_json::to_value(&proof)
            .unwrap()
            .get("interlacing")
            .is_none());
    }
}

#[cfg(test)]
mod left_hand_factor_tests {
    use super::super::*;

    fn q(value: i64) -> Q {
        Q::from_integer(value.into())
    }

    /// Both sides of `rule` multiplied by `x`.
    fn times_x(rule: &Recurrence) -> Recurrence {
        let shift = |p: &BivarPoly| {
            BivarPoly::from_coefficients_normalized(
                p.coeffs
                    .iter()
                    .map(|row| std::iter::once(q(0)).chain(row.iter().cloned()).collect())
                    .collect(),
            )
        };
        let one = BivarPoly::from_coefficients_normalized(vec![vec![q(1)]]);
        Recurrence {
            terms: rule
                .terms
                .iter()
                .map(|term| RecurrenceTerm {
                    coeff: shift(&term.coeff),
                    ..term.clone()
                })
                .collect(),
            denominator: Some(shift(rule.denominator.as_ref().unwrap_or(&one))),
            inhomogeneous: rule.inhomogeneous.as_ref().map(shift),
        }
    }

    fn prove(rule: &Recurrence, seeds: &[Vec<Q>], first_index: usize) -> VeroneseParentProof {
        let record = ParentRecurrenceRecord {
            rule,
            seeds,
            first_index,
        };
        prove_parent_real_rootedness(&record, first_index, &[], false).unwrap()
    }

    /// The type-B Eulerian rule times `x` is proved by the same route, with
    /// the normalization recorded.
    #[test]
    fn an_x_power_left_factor_is_normalized_for_the_theorem_routes() {
        let rule = Recurrence {
            terms: vec![
                RecurrenceTerm::new(
                    1,
                    0,
                    BivarPoly::from_coefficients_normalized(vec![
                        vec![q(1), q(-1)],
                        vec![q(0), q(2)],
                    ]),
                ),
                RecurrenceTerm::new(
                    1,
                    1,
                    BivarPoly::from_coefficients_normalized(vec![vec![q(0), q(2), q(-2)]]),
                ),
            ],
            denominator: None,
            inhomogeneous: None,
        };
        let seeds = vec![vec![q(1)]];
        let plain = prove(&rule, &seeds, 0);
        assert_eq!(plain.status, VeroneseParentProofStatus::Proved);
        assert!(plain.left_hand_factor_normalization.is_none());
        let lifted = prove(&times_x(&rule), &seeds, 0);
        assert_eq!(
            lifted.status,
            VeroneseParentProofStatus::Proved,
            "{lifted:#?}"
        );
        assert_eq!(lifted.route, plain.route);
        assert_eq!(lifted.interlacing.is_some(), plain.interlacing.is_some());
        let normalization = lifted.left_hand_factor_normalization.as_ref().unwrap();
        assert_eq!(normalization.x_power, 1);
        normalization.verify().unwrap();
    }

    /// The rows of A039692 times `x` still match their `lean_proved`
    /// catalog flag through the normalized rule.
    #[test]
    fn a_catalog_flag_matches_the_normalized_rule() {
        let entry = crate::oeis::by_id("A039692").unwrap();
        let (rule, rows) = entry.recurrence_parts().unwrap();
        let first_index = entry.recurrence_first_index;
        let count = rule.generation_initial_count(first_index, rows.len());
        let proof = prove(&times_x(&rule), &rows[..count], first_index);
        assert_eq!(
            proof.status,
            VeroneseParentProofStatus::ConditionalOnAttestation,
            "{proof:#?}"
        );
        assert_eq!(proof.conditional_on, ["oeis_catalog:A039692:lean_proved"]);
        assert!(proof.left_hand_factor_normalization.is_some());
    }

    /// A rule whose division by `x` happens only on the rows is kept as it
    /// is; the theorem routes name the reason.
    #[test]
    fn a_row_only_division_is_reported_by_the_theorem_routes() {
        // x P_n = (x^2 + n - 1) P_(n-1) - x P'_(n-1), rows x^n.
        let rule = Recurrence {
            terms: vec![
                RecurrenceTerm::new(
                    1,
                    0,
                    BivarPoly::from_coefficients_normalized(vec![
                        vec![q(-1), q(0), q(1)],
                        vec![q(1)],
                    ]),
                ),
                RecurrenceTerm::new(
                    1,
                    1,
                    BivarPoly::from_coefficients_normalized(vec![vec![q(0), q(-1)]]),
                ),
            ],
            denominator: Some(BivarPoly::from_coefficients_normalized(vec![vec![
                q(0),
                q(1),
            ]])),
            inhomogeneous: None,
        };
        let proof = prove(&rule, &[vec![q(1)]], 0);
        assert!(proof.left_hand_factor_normalization.is_none());
        let liu_wang = proof
            .attempts
            .iter()
            .find(|attempt| attempt.route == VeroneseParentProofRoute::LiuWang)
            .unwrap();
        assert!(
            liu_wang.detail.contains("happens only on the rows"),
            "{}",
            liu_wang.detail
        );
    }
}

#[cfg(test)]
mod registry_policy_tests {
    use super::super::cascade::run_theorem_routes;
    use super::super::*;

    #[test]
    fn only_product_proofs_can_be_replaced_by_later_theorems() {
        let entry = crate::oeis::by_id("A026386").unwrap();
        let (rule, rows) = entry.recurrence_parts().unwrap();
        let count = rule.generation_initial_count(entry.recurrence_first_index, rows.len());
        let record = ParentRecurrenceRecord {
            rule: &rule,
            seeds: &rows[..count],
            first_index: entry.recurrence_first_index,
        };
        for (initial, expected) in [
            (
                VeroneseParentProofRoute::ProductRecurrence,
                VeroneseParentProofRoute::IteratedPreserver,
            ),
            (
                VeroneseParentProofRoute::PeriodicProduct,
                VeroneseParentProofRoute::PeriodicProduct,
            ),
            (
                VeroneseParentProofRoute::PaperTheorem,
                VeroneseParentProofRoute::PaperTheorem,
            ),
            (
                VeroneseParentProofRoute::RootDomain,
                VeroneseParentProofRoute::RootDomain,
            ),
        ] {
            let mut proved = Some(initial);
            let mut interlacing = None;
            let mut attempts = Vec::new();
            run_theorem_routes(
                &record,
                ParentStage::ReplacesEarly,
                &|_| true,
                &mut proved,
                &mut interlacing,
                &mut attempts,
            )
            .unwrap();
            assert_eq!(proved, Some(expected));
            assert_eq!(
                attempts.is_empty(),
                initial != VeroneseParentProofRoute::ProductRecurrence
            );
        }
    }

    #[test]
    fn parent_region_fallback_matches_the_requested_then_unrestricted_checks() {
        for id in ["A008292", "A099605", "A201013"] {
            let entry = crate::oeis::by_id(id).unwrap();
            let (rule, rows) = entry.recurrence_parts().unwrap();
            let count = rule.generation_initial_count(entry.recurrence_first_index, rows.len());
            let record = ParentRecurrenceRecord {
                rule: &rule,
                seeds: &rows[..count],
                first_index: entry.recurrence_first_index,
            };
            for theorem in [
                RealRootednessTheorem::ResidueClassSubsequences,
                RealRootednessTheorem::IndexSectionLift,
                RealRootednessTheorem::FactorGauge,
            ] {
                let requested = RealRootednessTheoremOptions {
                    root_region: Some(RealInterval {
                        lower: None,
                        lower_closed: false,
                        upper: Some(Q::zero()),
                        upper_closed: false,
                    }),
                    ..Default::default()
                };
                let ray = check_real_rootedness_theorem(
                    theorem,
                    &rule,
                    record.first_index,
                    record.seeds,
                    &requested,
                )
                .unwrap();
                let expected = if ray.is_proved()
                    || ray.status == RealRootednessTheoremStatus::NotApplicable
                {
                    ray
                } else {
                    check_real_rootedness_theorem(
                        theorem,
                        &rule,
                        record.first_index,
                        record.seeds,
                        &Default::default(),
                    )
                    .unwrap()
                };
                let route = ROUTES.iter().find(|r| r.theorem == theorem).unwrap();
                let (actual, _, applicable) = theorem_attempt(&record, route).unwrap();
                assert_eq!(
                    applicable,
                    expected.status != RealRootednessTheoremStatus::NotApplicable
                );
                let upper = expected
                    .conclusion
                    .as_ref()
                    .and_then(|c| c.root_region.as_ref())
                    .and_then(|region| region.at(0).ok())
                    .and_then(|region| region.upper);
                assert_eq!(
                    actual.success,
                    expected.is_proved() && upper.is_some_and(|u| !u.is_positive())
                );
                assert!(actual.detail.contains(&expected.normal_form));
            }
        }
    }
}
