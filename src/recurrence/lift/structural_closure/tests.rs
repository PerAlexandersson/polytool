//! Structural discovery fixtures and optional test-corpus breadth audit.
use super::*;
use crate::real_rootedness::RealInterval;
use std::io::Write;

fn record(id: &str) -> RecurrenceJson {
    let entry = crate::oeis::by_id(id).unwrap();
    let (rule, seeds) = entry.recurrence_parts().unwrap();
    RecurrenceJson::from_recurrence_rational(&rule, entry.recurrence_first_index, &seeds, None)
}

#[test]
fn inverse_operators_are_exact() {
    let r: Vec<Q> = [3, 7, 4, 1]
        .into_iter()
        .map(|n| Q::from_integer(n.into()))
        .collect();
    for c in [-1i64, 1] {
        let q = inverse_derivative(&r, c);
        let mut reconstructed = q.clone();
        for k in 1..q.len() {
            reconstructed[k - 1] += &q[k] * Q::from_integer((c * k as i64).into());
        }
        assert_eq!(reconstructed, r);
    }
}

#[test]
fn forward_gamma_preimage_replays_and_rejects_tampering() {
    let input = record("A001263");
    let search =
        find_discovery_closure_lifts(&input, &Default::default(), &Default::default()).unwrap();
    assert!(search.is_proved(), "{search:?}");
    search.verify(&Default::default()).unwrap();
    let mut changed = search.clone();
    changed.attempts[0].steps[0].statement.push('!');
    assert!(changed.verify(&Default::default()).is_err());
    let mut changed = search.clone();
    changed.conclusion.as_mut().unwrap().from_index += 1;
    assert!(changed.verify(&Default::default()).is_err());
    let mut changed = search;
    if let DiscoveryClosureCertificate::ForwardGamma { certification } =
        &mut changed.attempts[0].certificate
    {
        certification
            .gamma_transport
            .request
            .source_recurrence
            .as_mut()
            .unwrap()
            .initial_polynomials[0] = vec!["999".into()];
    } else {
        panic!("expected forward gamma expansion");
    }
    assert!(changed.verify(&Default::default()).is_err());
}

#[test]
fn disabled_and_region_requested_searches_do_not_prove() {
    let input = record("A001263");
    let disabled =
        find_discovery_closure_lifts(&input, &Default::default(), &CatalogLiftOptions::disabled())
            .unwrap();
    assert!(!disabled.is_proved());
    disabled.verify(&Default::default()).unwrap();
    let options = RealRootednessTheoremOptions {
        root_region: Some(RealInterval {
            lower: None,
            lower_closed: false,
            upper_closed: false,
            upper: Some(Q::zero()),
        }),
        ..Default::default()
    };
    assert!(
        !find_discovery_closure_lifts(&input, &options, &Default::default())
            .unwrap()
            .is_proved()
    );
}

#[test]
#[ignore = "explicit corpus breadth audit; output path supplied by environment"]
fn structural_closure_corpus_breadth() {
    let path = std::env::var("POLYTOOL_CLOSURE_AUDIT").expect("output JSONL path");
    let mut out = std::io::BufWriter::new(std::fs::File::create(path).unwrap());
    let only = std::env::var("POLYTOOL_CLOSURE_AUDIT_IDS").ok();
    for line in include_str!("../../../../data/oeis-catalog.jsonl").lines() {
        let e: serde_json::Value = serde_json::from_str(line).unwrap();
        let id = e["id"].as_str().unwrap();
        if only
            .as_ref()
            .is_some_and(|ids| !ids.split(',').any(|s| s == id))
        {
            continue;
        }
        let input = record(id);
        for stage in ["hadamard_factor", "row_factor", "gamma", "preserver_image"] {
            let start = std::time::Instant::now();
            let found = search(&input, &Default::default(), 6, MAX_DEPTH, Some(stage));
            let item = match found {
                Ok(found) => {
                    serde_json::json!({"id":id,"stage":stage,"proved":found.is_proved(),"candidates":found.candidates_tried,"seconds":start.elapsed().as_secs_f64(),"search":found})
                }
                Err(error) => {
                    serde_json::json!({"id":id,"stage":stage,"proved":false,"error":error.to_string(),"seconds":start.elapsed().as_secs_f64()})
                }
            };
            writeln!(out, "{item}").unwrap();
            out.flush().unwrap();
            eprintln!(
                "{id} {stage} {} {:.3}s",
                item["proved"],
                start.elapsed().as_secs_f64()
            );
        }
    }
}

#[test]
fn inverse_gamma_and_row_factors_are_structural() {
    for (id, stage) in [("A008303", "gamma"), ("A008518", "row_factor")] {
        let found = search(&record(id), &Default::default(), 6, MAX_DEPTH, Some(stage)).unwrap();
        assert!(found.is_proved(), "{id}: {found:?}");
        assert_eq!(found.attempts[0].stage, stage);
    }
}

#[test]
fn differential_preimage_is_computed_without_a_parent_library() {
    let source: RecurrenceJson = serde_json::from_value(serde_json::json!({
        "schema":"polytool.recurrence.v1", "first_index":1,
        "initial_polynomials":[["1"]],
        "recurrence":{"terms":[
            {"offset":1,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["1","-1"],["0","1"]]}},
            {"offset":1,"deriv_order":1,"sign":"none","coeff":{"coeffs":[["0","1","-1"]]}}
        ]}
    }))
    .unwrap();
    let image = PreserverImageCandidate {
        parent_name: "test input".into(),
        parent: source,
        index_shift: 0,
        gauge: VeroneseSectionGauge::None,
        operator: PreserverOperator::derivative_polynomial(&[Q::one(), Q::one()]),
    }
    .image_recurrence(1, &Default::default())
    .unwrap();
    let found = search(
        &image.rule,
        &Default::default(),
        6,
        MAX_DEPTH,
        Some("preserver_image"),
    )
    .unwrap();
    assert!(found.is_proved(), "{found:?}");
    let full = find_discovery_closure_lifts(&image.rule, &Default::default(), &Default::default())
        .unwrap();
    assert!(full.is_proved());
    full.verify(&Default::default()).unwrap();
    let mut altered = full.clone();
    altered.depth = MAX_DEPTH + 1;
    assert!(altered.verify(&Default::default()).is_err());
    if let StructuralClosureParentProof::Registry { proof } = &full.attempts[0].parent_proof {
        assert!(proof.conditional_on.is_empty());
        assert!(proof
            .attempts
            .iter()
            .all(|a| a.route != VeroneseParentProofRoute::Attestation));
    } else {
        panic!("expected a direct registry proof");
    }
}

#[test]
fn earlier_proofs_keep_priority() {
    let proof = discover_recurrence_proofs(&record("A001263"), &Default::default()).unwrap();
    assert_eq!(proof.status, "proved");
    assert!(proof.closure_lifts.is_none());
    proof.verify().unwrap();
}

#[test]
fn matching_twenty_five_rows_is_not_an_all_index_identity() {
    let one = Q::one();
    let b = |coeffs: Vec<Vec<Q>>| BivarPoly::from_coefficients_normalized(coeffs);
    let source = Recurrence {
        terms: vec![RecurrenceTerm::new(
            1,
            0,
            b(vec![vec![one.clone(), one.clone()]]),
        )],
        denominator: None,
        inhomogeneous: None,
    };
    let parent = RecurrenceJson::from_recurrence_rational(&source, 0, &[vec![one.clone()]], None);
    let candidate = PreserverImageCandidate {
        parent_name: "computed test preimage".into(),
        parent,
        index_shift: 0,
        gauge: VeroneseSectionGauge::None,
        operator: PreserverOperator::derivative_polynomial(&[one.clone(), one.clone()]),
    };
    // Q_n=(1+D)(1+x)^n obeys the repeated-root order-two rule below.
    // A perturbation vanishing at n=0,...,24 makes every finite check of
    // sixteen rows pass, while leaving a nonzero all-index Ore remainder.
    let perturbation = (0i64..=24).fold(vec![one.clone()], |p, j| {
        vec_poly::mul_linear(&p, &Q::from_integer((-j).into()), &one)
    });
    let mut a: Vec<Vec<Q>> = perturbation
        .into_iter()
        .map(|c| vec![c.clone(), Q::zero(), c])
        .collect();
    a[0][0] += Q::from_integer(2.into());
    a[0][1] += Q::from_integer(2.into());
    let target = Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, b(a)),
            RecurrenceTerm::new(
                2,
                0,
                b(vec![vec![
                    Q::from_integer((-1).into()),
                    Q::from_integer((-2).into()),
                    Q::from_integer((-1).into()),
                ]]),
            ),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let certification = certify_preserver_image_lift(
        &target,
        &[vec![one.clone()], vec![Q::from_integer(2.into()), one]],
        0,
        &candidate,
        &PreserverImageOptions {
            allow_attestations: false,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(certification.rows_compared, 16);
    assert!(!matches!(
        certification.identity,
        PreserverImageIdentityStatus::Certified { .. }
    ));
    assert!(!certification.is_proved());
}

#[test]
fn row_factor_then_gamma_needs_two_closure_levels() {
    let image = PreserverImageCandidate {
        parent_name: "test source".into(),
        parent: record("A008303"),
        index_shift: 0,
        gauge: VeroneseSectionGauge::None,
        operator: PreserverOperator::DerivativeOfPowerMultiple {
            power: 1,
            derivative_order: 0,
        },
    }
    .image_recurrence(0, &Default::default())
    .unwrap();
    let shallow = search(&image.rule, &Default::default(), 6, 1, Some("row_factor")).unwrap();
    assert!(!shallow.is_proved());
    let deep = search(&image.rule, &Default::default(), 6, 2, Some("row_factor")).unwrap();
    assert!(deep.is_proved(), "{deep:?}");
    let StructuralClosureParentProof::Closure { search } = &deep.attempts[0].parent_proof else {
        panic!("expected a recursive parent proof");
    };
    assert_eq!(search.depth, 1);
    assert_eq!(search.attempts[0].stage, "gamma");
    search.verify(&Default::default()).unwrap();
    let full = find_discovery_closure_lifts(&image.rule, &Default::default(), &Default::default())
        .unwrap();
    full.verify(&Default::default()).unwrap();
}
