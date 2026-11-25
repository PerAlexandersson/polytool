//! Independent identities, acceptance counterexamples and route hypothesis guards.
use super::*;
use num_traits::Signed;
fn sizes(s: &[usize]) -> BlockWeightSpec {
    BlockWeightSpec::Sizes { sizes: s.to_vec() }
}
fn family(a: &[usize], b: &[usize]) -> ExponentialFamily {
    ExponentialFamily {
        unmarked: sizes(a),
        marked: sizes(b),
    }
}
fn request(
    family: ExponentialFamily,
    rows: Vec<usize>,
    check: ExponentialRowCheck,
) -> ExponentialRowsRequest {
    ExponentialRowsRequest {
        family,
        rows,
        check,
        summary: false,
        finite_rank: vec![],
    }
}
fn q(v: i64) -> Q {
    Q::from_integer(v.into())
}
#[test]
fn touchard_bell_and_big_block_catalogue() {
    let touchard = ExponentialFamily {
        unmarked: sizes(&[]),
        marked: BlockWeightSpec::Tail { minimum: 1 },
    };
    let rows = touchard.rows(14).unwrap();
    assert_eq!(rows[5], vec![q(0), q(1), q(15), q(25), q(10), q(1)]);
    assert_eq!(
        rows.iter()
            .take(11)
            .map(|r| r.iter().sum::<Q>())
            .collect::<Vec<_>>(),
        [1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975].map(q)
    );
    assert_eq!(
        rows,
        touchard
            .selected_rows(&(0..=14).collect::<Vec<_>>())
            .unwrap()
    );
    for j in [2, 3] {
        let f = ExponentialFamily {
            unmarked: sizes(&(1..j).collect::<Vec<_>>()),
            marked: BlockWeightSpec::Tail { minimum: j },
        };
        let r = f.rows(14).unwrap();
        assert_eq!(r, f.selected_rows(&(0..=14).collect::<Vec<_>>()).unwrap());
        assert!(r
            .iter()
            .enumerate()
            .all(|(n, p)| vec_poly::degree(p) == Some(n / j)));
        if j == 3 {
            let catalog = crate::oeis::A355144().generate_rows(15).unwrap();
            for (a, b) in r.iter().zip(catalog) {
                assert_eq!(
                    strings(a),
                    b.iter().map(ToString::to_string).collect::<Vec<_>>()
                );
            }
        } else {
            // Threshold two: G_(n+1)=G_n+tG'_n+n*tG_(n-1).
            for n in 1..14 {
                let derivative = vec_poly::derivative(&r[n]);
                let mut rhs = r[n].clone();
                rhs.resize(r[n + 1].len(), Q::zero());
                for (i, c) in derivative.iter().enumerate() {
                    rhs[i + 1] += c;
                }
                for (i, c) in r[n - 1].iter().enumerate() {
                    rhs[i + 1] += q(n as i64) * c;
                }
                rhs = vec_poly::trim(&rhs);
                assert_eq!(rhs, r[n + 1]);
            }
        }
    }
    let report = check_exponential_rows(&request(
        touchard,
        (0..=14).collect(),
        ExponentialRowCheck::Interlacing,
    ))
    .unwrap();
    assert!(report.first_failure.is_none());
    assert!(report.interlacing.iter().all(|p| p.forward == Some(true)));
}
#[test]
fn single_marked_size_matches_dense_including_rational_weights() {
    let f = ExponentialFamily {
        unmarked: BlockWeightSpec::Finite {
            weights: vec!["1/2".into(), "2/3".into()],
        },
        marked: BlockWeightSpec::Finite {
            weights: vec!["0".into(), "0".into(), "3/5".into()],
        },
    };
    assert_eq!(
        f.rows(30).unwrap(),
        f.selected_rows(&(0..=30).collect::<Vec<_>>()).unwrap()
    );
    let mut iter = f.rows_iter(3).unwrap();
    assert_eq!(iter.len(), 4);
    assert_eq!(iter.by_ref().count(), 4);
    assert!(iter.next().is_none());
    assert!(iter.next().is_none());
    assert_eq!(
        family(&[], &[]).rows(3).unwrap(),
        vec![vec![q(1)], vec![q(0)], vec![q(0)], vec![q(0)]]
    );
}
#[test]
#[cfg_attr(
    not(feature = "full-tests"),
    ignore = "expensive exact late failure scan"
)]
fn acceptance_late_failure_305_and_sparse_1000() {
    let f = family(&[1, 2], &[7]);
    let report = check_exponential_rows(&request(
        f.clone(),
        (0..=305).collect(),
        ExponentialRowCheck::FirstFailure,
    ))
    .unwrap();
    let failure = report.first_failure.unwrap();
    assert_eq!(failure.index, 305);
    assert_eq!(failure.degree, Some(43));
    assert!(failure.distinct_real_roots < failure.squarefree_degree);
    failure.verify().unwrap();
    let sparse = f.selected_rows(&[304, 305, 1000]).unwrap();
    assert_eq!(
        strings(&sparse[0]),
        report.rows[304].coefficients.clone().unwrap()
    );
    assert_eq!(
        strings(&sparse[1]),
        report.rows[305].coefficients.clone().unwrap()
    );
    assert_eq!(sparse[2].len(), 143);
    let route = exponential_row_routes(&f).unwrap();
    assert_eq!(route[0].status, "obstruction");
    assert_eq!(route[0].evidence, ExponentialEvidenceLevel::PaperProof);
}
#[test]
fn invalid_cubic_substitutions_fail_at_six() {
    // B=z^2+z^6 = z^2(1+z^4), with ordinary coefficients converted to weights.
    let f = ExponentialFamily {
        unmarked: sizes(&[1]),
        marked: BlockWeightSpec::Finite {
            weights: vec!["0", "2", "0", "0", "0", "720"]
                .into_iter()
                .map(str::to_string)
                .collect(),
        },
    };
    let report = check_exponential_rows(&request(
        f.clone(),
        (0..=8).collect(),
        ExponentialRowCheck::FirstFailure,
    ))
    .unwrap();
    assert_eq!(report.first_failure.as_ref().unwrap().index, 6);
    assert_eq!(f.row(6).unwrap(), [1, 750, 180, 120].map(q));
    report.first_failure.unwrap().verify().unwrap();
    assert!(exponential_row_routes(&f).unwrap().is_empty());
}
#[test]
fn ode_hypotheses_are_structural_and_exact() {
    let seed = family(&[1], &[2]);
    for (a, b, l, accepted) in [
        ("0", "0", "0", true),
        ("1/2", "2/3", "1/7", true),
        ("-1", "1", "1", false),
        ("1", "-1", "1", false),
        ("1", "1", "-1", false),
    ] {
        let f = ExponentialFamily {
            unmarked: seed.unmarked.clone(),
            marked: BlockWeightSpec::OdeSubstitution {
                seed: Box::new(seed.marked.clone()),
                a: a.into(),
                b: b.into(),
                lambda: l.into(),
            },
        };
        let routes = exponential_row_routes(&f).unwrap();
        assert_eq!(!routes.is_empty(), accepted);
        if accepted {
            assert_eq!(routes[0].route, "ode_substitution");
            assert!(f
                .rows(8)
                .unwrap()
                .iter()
                .all(|p| crate::is_real_rooted_bigint_coeffs(
                    &exact_algebra::denominators::clear_row_denominators(p).0
                )));
        }
    }
    let exp = BlockWeightSpec::Exponential {
        scale: "2/3".into(),
        rate: "3/2".into(),
        minimum: 2,
    };
    let ode = BlockWeightSpec::OdeSubstitution {
        seed: Box::new(sizes(&[1])),
        a: "0".into(),
        b: "3/2".into(),
        lambda: "3/2".into(),
    };
    assert_eq!(
        exp.block_weights(12).unwrap(),
        ode.block_weights(12).unwrap()
    );
    let unknown = ExponentialFamily {
        unmarked: sizes(&[1, 2]),
        marked: BlockWeightSpec::Tail { minimum: 3 },
    };
    assert!(exponential_row_routes(&unknown).unwrap().is_empty());
    let bad_seed = ExponentialFamily {
        unmarked: sizes(&[1]),
        marked: BlockWeightSpec::OdeSubstitution {
            seed: Box::new(sizes(&[2, 6])),
            a: "1".into(),
            b: "0".into(),
            lambda: "1".into(),
        },
    };
    assert!(exponential_row_routes(&bad_seed).unwrap().is_empty());
    let negative = ExponentialFamily {
        unmarked: BlockWeightSpec::Finite {
            weights: vec!["-1".into()],
        },
        marked: ode,
    };
    assert!(exponential_row_routes(&negative).unwrap().is_empty());
    let large_sparse = family(&[2], &[1_000_003]);
    assert_eq!(
        exponential_row_routes(&large_sparse).unwrap()[0].status,
        "obstruction"
    );
    let eight = family(&[3], &[8]);
    let r = exponential_row_routes(&eight).unwrap();
    assert_eq!(r[0].evidence, ExponentialEvidenceLevel::UnderReview);
    assert_eq!(r[0].status, "under_review");
    assert!(BlockWeightSpec::Tail { minimum: 0 }
        .block_weights(3)
        .is_err());
    assert!(sizes(&[0]).block_weights(3).is_err());
    assert!(family(&[1], &[2]).row(0).unwrap()[0].is_positive());
}
#[test]
fn finite_rank_normalization_and_pencil_replay() {
    let reference = FiniteRankRequest {
        index: 6,
        reference: vec!["1".into(), "2".into(), "1".into()],
        pencil: None,
    };
    let c = certify_finite_rank(&[q(1), q(2), q(1)], &reference).unwrap();
    assert!(c.certified);
    assert_eq!(c.witnesses, vec![vec!["1", "2", "1"]]);
    c.verify().unwrap();
    let failed = certify_finite_rank(&[q(1), q(1), q(1)], &reference).unwrap();
    assert!(!failed.certified);
    let mut p = reference.clone();
    p.pencil = Some(PencilReference {
        a: "1".into(),
        b: "1".into(),
        factor: vec!["1".into(), "1".into()],
    });
    let c = certify_finite_rank(&[q(1), q(2), q(1)], &p).unwrap();
    assert!(c.certified);
    c.verify().unwrap();
    p.pencil.as_mut().unwrap().a = "2".into();
    assert!(certify_finite_rank(&[q(1), q(2), q(1)], &p).is_err());
    let mut tampered = c.clone();
    tampered.certified = false;
    assert!(tampered.verify().is_err());
}

#[test]
fn polynomial_paper_route_checks_zeros_and_never_upgrades_evidence() {
    let transform = |coefficients: Vec<&str>| BlockWeightSpec::PolynomialSubstitution {
        seed: Box::new(sizes(&[2])),
        coefficients: coefficients.into_iter().map(str::to_string).collect(),
    };
    let good = ExponentialFamily {
        unmarked: sizes(&[1]),
        marked: transform(vec!["2", "3", "1"]),
    };
    let routes = exponential_row_routes(&good).unwrap();
    assert_eq!(routes[0].evidence, ExponentialEvidenceLevel::PaperProof);
    let composed = ExponentialFamily {
        unmarked: good.unmarked.clone(),
        marked: BlockWeightSpec::OdeSubstitution {
            seed: Box::new(good.marked.clone()),
            a: "1".into(),
            b: "0".into(),
            lambda: "1".into(),
        },
    };
    assert_eq!(
        exponential_row_routes(&composed).unwrap()[0].evidence,
        ExponentialEvidenceLevel::PaperProof
    );
    let invalid = ExponentialFamily {
        unmarked: sizes(&[1]),
        marked: transform(vec!["0", "1", "0", "1"]),
    };
    assert!(exponential_row_routes(&invalid).unwrap().is_empty());
    let report = check_exponential_rows(&request(
        invalid,
        (0..=6).collect(),
        ExponentialRowCheck::FirstFailure,
    ))
    .unwrap();
    assert_eq!(report.first_failure.unwrap().index, 6);
    let altered_a = ExponentialFamily {
        unmarked: sizes(&[1, 2]),
        ..good
    };
    assert!(exponential_row_routes(&altered_a).unwrap().is_empty());
}

#[test]
fn interlacing_failure_retains_both_rows_and_replays() {
    // F_1=1+t; F_2=6+6t+t^2 has both roots strictly left of -1.
    let f = ExponentialFamily {
        unmarked: BlockWeightSpec::Finite {
            weights: vec!["1".into(), "5".into()],
        },
        marked: BlockWeightSpec::Finite {
            weights: vec!["1".into(), "4".into()],
        },
    };
    let report =
        check_exponential_rows(&request(f, vec![1, 2], ExponentialRowCheck::Interlacing)).unwrap();
    assert!(report.rows.iter().all(|r| r.real_rooted == Some(true)));
    let certificate = report.first_failure.unwrap();
    assert_eq!(certificate.property, "interlacing_either_orientation");
    assert_eq!(certificate.index, 2);
    certificate.verify().unwrap();
    let mut tampered = certificate;
    tampered.previous_integer_coefficients = None;
    assert!(tampered.verify().is_err());
}

#[test]
fn finite_coefficient_specs_recognize_full_quadratic_identities() {
    let unsupported = family(&[1, 2], &[1]);
    assert!(exponential_row_routes(&unsupported).unwrap().is_empty());
    let failure = check_exponential_rows(&request(
        unsupported,
        vec![0, 1, 2],
        ExponentialRowCheck::FirstFailure,
    ))
    .unwrap();
    assert_eq!(failure.first_failure.unwrap().index, 2);
    // (z^2+z^3)^2 and (z^2+z^4)+3*(z^2+z^4)^2, in block weights.
    for weights in [
        vec!["0", "0", "0", "24", "240", "720"],
        vec!["0", "2", "0", "96", "0", "4320", "0", "120960"],
        vec!["0", "0", "0", "0", "0", "0", "0", "40320"],
    ] {
        let f = ExponentialFamily {
            unmarked: sizes(&[1]),
            marked: BlockWeightSpec::Finite {
                weights: weights.into_iter().map(str::to_string).collect(),
            },
        };
        let routes = exponential_row_routes(&f).unwrap();
        assert_eq!(routes[0].route, "ode_substitution");
        assert_eq!(routes[0].evidence, ExponentialEvidenceLevel::TheoremBacked);
        assert!(!routes[0].consecutive_weak_interlacing);
        assert!(f
            .rows(12)
            .unwrap()
            .iter()
            .all(|p| crate::is_real_rooted_bigint_coeffs(
                &exact_algebra::denominators::clear_row_denominators(p).0
            )));
    }
    let near_miss = ExponentialFamily {
        unmarked: sizes(&[1]),
        marked: BlockWeightSpec::Finite {
            weights: vec!["0", "0", "0", "24", "240", "1440"]
                .into_iter()
                .map(str::to_string)
                .collect(),
        },
    };
    assert!(exponential_row_routes(&near_miss).unwrap().is_empty());
    let quadratic = ExponentialFamily {
        unmarked: sizes(&[1]),
        marked: BlockWeightSpec::Finite {
            weights: vec!["1".into(), "2".into()],
        },
    };
    let report = check_exponential_rows(&request(
        quadratic,
        (0..=10).collect(),
        ExponentialRowCheck::Interlacing,
    ))
    .unwrap();
    assert!(report.first_failure.is_none());
    assert!(report.routes[0].consecutive_weak_interlacing);
}

#[test]
fn rational_ode_closure_is_exact_and_theorem_backed() {
    for (lambda, geometric) in [("0", false), ("1", true)] {
        let spec = BlockWeightSpec::RationalOdeSubstitution {
            seed: Box::new(sizes(&[1])),
            a: "1".into(),
            b: "0".into(),
            lambda: lambda.into(),
            d: "1".into(),
        };
        let weights = spec.block_weights(10).unwrap();
        for (r, w) in weights.iter().enumerate().skip(1) {
            assert_eq!(
                *w,
                Q::from_integer(factorial_bigint(if geometric { r } else { r - 1 }))
            );
        }
        let f = ExponentialFamily {
            unmarked: sizes(&[1]),
            marked: spec,
        };
        let routes = exponential_row_routes(&f).unwrap();
        assert_eq!(routes[0].route, "rational_ode_substitution");
        assert_eq!(routes[0].evidence, ExponentialEvidenceLevel::TheoremBacked);
        assert_eq!(routes[0].status, "proved");
        assert!(routes[0].consecutive_weak_interlacing);
        let report = check_exponential_rows(&request(
            f,
            (0..=8).collect(),
            ExponentialRowCheck::Interlacing,
        ))
        .unwrap();
        assert!(report.first_failure.is_none());
    }
    let spec = |d: &str| BlockWeightSpec::RationalOdeSubstitution {
        seed: Box::new(sizes(&[2])),
        a: "1".into(),
        b: "2".into(),
        lambda: "0".into(),
        d: d.into(),
    };
    let formal = ExponentialFamily {
        unmarked: sizes(&[1]),
        marked: spec("0"),
    };
    assert_eq!(
        exponential_row_routes(&formal).unwrap()[0].evidence,
        ExponentialEvidenceLevel::TheoremBacked
    );
    let negative = ExponentialFamily {
        marked: spec("-1"),
        ..formal.clone()
    };
    assert!(exponential_row_routes(&negative).unwrap().is_empty());
    let wrong_unmarked = ExponentialFamily {
        unmarked: sizes(&[1, 2]),
        marked: spec("1"),
    };
    assert!(exponential_row_routes(&wrong_unmarked).unwrap().is_empty());
}

#[test]
fn divisibility_filter_handles_infinite_unmarked_supports_exactly() {
    for unmarked in [
        BlockWeightSpec::Tail { minimum: 2 },
        BlockWeightSpec::Exponential {
            scale: "1/2".into(),
            rate: "2".into(),
            minimum: 10,
        },
    ] {
        let f = ExponentialFamily {
            unmarked,
            marked: sizes(&[7]),
        };
        let routes = exponential_row_routes(&f).unwrap();
        assert_eq!(routes[0].status, "obstruction");
        assert_eq!(routes[0].evidence, ExponentialEvidenceLevel::PaperProof);
    }
    let zero = ExponentialFamily {
        unmarked: BlockWeightSpec::Exponential {
            scale: "0".into(),
            rate: "-1".into(),
            minimum: 1,
        },
        marked: sizes(&[7]),
    };
    assert_eq!(
        exponential_row_routes(&zero).unwrap()[0].status,
        "necessary_condition_passed"
    );
}

fn expression(s: &str) -> BlockWeightSpec {
    BlockWeightSpec::SeriesExpression {
        expression: s.into(),
    }
}
#[test]
fn series_expression_acceptance_and_sparse_equivalence() {
    for (a, b, structured) in [
        (
            "0",
            "exp(z)-1",
            ExponentialFamily {
                unmarked: sizes(&[]),
                marked: BlockWeightSpec::Tail { minimum: 1 },
            },
        ),
        (
            "z+z^2/2",
            "exp(z)-1-z-z^2/2",
            ExponentialFamily {
                unmarked: sizes(&[1, 2]),
                marked: BlockWeightSpec::Tail { minimum: 3 },
            },
        ),
        ("z", "z^2/2", family(&[1], &[2])),
        (
            "z",
            "z^2+2*z^3+z^4/3",
            ExponentialFamily {
                unmarked: sizes(&[1]),
                marked: BlockWeightSpec::Finite {
                    weights: vec!["0".into(), "2".into(), "12".into(), "8".into()],
                },
            },
        ),
        (
            "z",
            "log(1/(1-z))",
            ExponentialFamily {
                unmarked: sizes(&[1]),
                marked: BlockWeightSpec::RationalOdeSubstitution {
                    seed: Box::new(sizes(&[1])),
                    a: "1".into(),
                    b: "0".into(),
                    lambda: "0".into(),
                    d: "1".into(),
                },
            },
        ),
        (
            "z",
            "1/(1-z)-1",
            ExponentialFamily {
                unmarked: sizes(&[1]),
                marked: BlockWeightSpec::RationalOdeSubstitution {
                    seed: Box::new(sizes(&[1])),
                    a: "1".into(),
                    b: "0".into(),
                    lambda: "1".into(),
                    d: "1".into(),
                },
            },
        ),
    ] {
        let f = ExponentialFamily {
            unmarked: expression(a),
            marked: expression(b),
        };
        assert_eq!(f.rows(14).unwrap(), structured.rows(14).unwrap(), "{b}");
        assert_eq!(
            f.selected_rows(&[12, 0, 3, 14, 12]).unwrap(),
            [12, 0, 3, 14, 12].map(|n| structured.row(n).unwrap()),
            "{b}"
        );
        if b.contains("log") || b.contains("1/(1-z)") {
            let routes = exponential_row_routes(&f).unwrap();
            assert_eq!(
                routes[0].evidence,
                ExponentialEvidenceLevel::TheoremBacked,
                "{b}"
            );
        } else if !b.contains("z^2/2") || a == "z" {
            assert!(!exponential_row_routes(&f).unwrap().is_empty(), "{b}");
        }
    }
    for b in ["z^2*(1+z^4)", "z^2+z^6"] {
        let f = ExponentialFamily {
            unmarked: expression("z"),
            marked: expression(b),
        };
        let report = check_exponential_rows(&request(
            f,
            (0..=8).collect(),
            ExponentialRowCheck::FirstFailure,
        ))
        .unwrap();
        assert_eq!(report.first_failure.as_ref().unwrap().index, 6);
        report.first_failure.unwrap().verify().unwrap();
        assert!(report.routes.is_empty());
    }
    assert_eq!(
        ExponentialFamily {
            unmarked: expression("z"),
            marked: expression("z^2*(1+z)^3")
        }
        .row(8)
        .unwrap(),
        ExponentialFamily {
            unmarked: sizes(&[1]),
            marked: BlockWeightSpec::Finite {
                weights: [0, 2, 18, 72, 120].map(|x| x.to_string()).to_vec()
            }
        }
        .row(8)
        .unwrap()
    );
}

#[test]
fn proved_substitution_faces_cite_pinned_splitness_and_interlacing_theorems() {
    for (outer, route, splits, interlacing) in [
        (
            "x/(1-x)",
            "rational_ode_substitution",
            "insertionRow_splits_subst_of_rational_ode",
            "insertionRow_interl_subst_of_rational_ode",
        ),
        (
            "-log(1-x)",
            "rational_ode_substitution",
            "insertionRow_splits_subst_of_rational_ode",
            "insertionRow_interl_subst_of_rational_ode",
        ),
        (
            "x^2",
            "power_two_substitution",
            "insertionRow_splits_pow_two_pow",
            "insertionRow_interl_pow_two_pow",
        ),
        (
            "x^4",
            "power_two_substitution",
            "insertionRow_splits_pow_two_pow",
            "insertionRow_interl_pow_two_pow",
        ),
        (
            "x^8",
            "power_two_substitution",
            "insertionRow_splits_pow_two_pow",
            "insertionRow_interl_pow_two_pow",
        ),
        (
            "x^2+x^4",
            "even_quartic_substitution",
            "insertionRow_splits_even_quartic_subst",
            "insertionRow_interl_even_quartic_subst",
        ),
        (
            "2*x^2+3*x^4",
            "even_quartic_substitution",
            "insertionRow_splits_even_quartic_subst",
            "insertionRow_interl_even_quartic_subst",
        ),
    ] {
        // The exponential seed has a theorem-backed consecutive-interlacing route.
        let f = ExponentialFamily {
            unmarked: sizes(&[1]),
            marked: expression(&format!("compose({outer},exp(z)-1)")),
        };
        let report = check_exponential_rows(&request(
            f,
            (0..=10).collect(),
            ExponentialRowCheck::Interlacing,
        ))
        .unwrap();
        let r = &report.routes[0];
        assert_eq!(r.route, route, "{outer}");
        assert_eq!(
            r.evidence,
            ExponentialEvidenceLevel::TheoremBacked,
            "{outer}"
        );
        assert_eq!(r.status, "proved");
        assert!(r.consecutive_weak_interlacing, "{outer}");
        assert!(r.references[0].contains("/blob/b09ae7ff6562bcaada3f998622737f7091feb9f3/"));
        for declaration in [splits, interlacing] {
            assert!(r
                .references
                .contains(&["RealRooted", "ExponentialRows", declaration].join(".")));
        }
        assert!(report.first_failure.is_none(), "{outer}");
    }
    for power in [2, 4, 8] {
        let routes = exponential_row_routes(&ExponentialFamily {
            unmarked: sizes(&[1]),
            marked: expression(&format!("(exp(z)-1)^{power}")),
        })
        .unwrap();
        assert_eq!(routes[0].route, "power_two_substitution");
        assert!(routes[0].consecutive_weak_interlacing);
    }
    // Finite monomials retain their support for the separate necessity route.
    let routes = exponential_row_routes(&ExponentialFamily {
        unmarked: sizes(&[2]),
        marked: expression("z^7"),
    })
    .unwrap();
    assert_eq!(routes[0].route, "single_marked_size_divisibility");
    assert_eq!(routes[0].status, "obstruction");
}

#[test]
fn substitution_closures_preserve_seed_evidence_and_interlacing_scope() {
    let polynomial = |seed, coefficients: &[&str]| BlockWeightSpec::PolynomialSubstitution {
        seed: Box::new(seed),
        coefficients: coefficients.iter().map(|c| c.to_string()).collect(),
    };
    let paper_seed = polynomial(sizes(&[1]), &["0", "0", "0", "1"]);
    for coefficients in [
        vec![],
        vec!["0"],
        vec!["-7"],
        vec!["-7", "0", "0", "0", "2"],
        vec!["0", "0", "2", "0", "0"],
    ] {
        let routes = exponential_row_routes(&ExponentialFamily {
            unmarked: sizes(&[1]),
            marked: polynomial(sizes(&[1]), &coefficients),
        })
        .unwrap();
        assert_eq!(routes[0].route, "even_quartic_substitution");
        assert_eq!(routes[0].evidence, ExponentialEvidenceLevel::TheoremBacked);
        assert!(routes[0].consecutive_weak_interlacing);
    }
    for marked in [
        paper_seed.clone(),
        polynomial(sizes(&[1]), &["0", "0", "1", "3", "1"]),
        polynomial(paper_seed.clone(), &["0", "0", "0", "0", "1"]),
        polynomial(paper_seed.clone(), &["0", "0", "1", "0", "1"]),
        BlockWeightSpec::RationalOdeSubstitution {
            seed: Box::new(paper_seed),
            a: "1".into(),
            b: "0".into(),
            lambda: "1".into(),
            d: "1".into(),
        },
    ] {
        let routes = exponential_row_routes(&ExponentialFamily {
            unmarked: sizes(&[1]),
            marked,
        })
        .unwrap();
        assert_eq!(routes[0].evidence, ExponentialEvidenceLevel::PaperProof);
        assert_eq!(routes[0].status, "paper_proof");
    }
    for outer in ["x/(1-x)", "x^4", "x^2+x^4"] {
        let routes = exponential_row_routes(&ExponentialFamily {
            unmarked: sizes(&[1]),
            marked: expression(&format!("compose({outer},z^2+z^3)")),
        })
        .unwrap();
        assert_eq!(routes[0].evidence, ExponentialEvidenceLevel::TheoremBacked);
        assert!(!routes[0].consecutive_weak_interlacing);
        assert!(!routes[0]
            .references
            .iter()
            .any(|s| s.contains(".insertionRow_interl_")));
    }
    // The new theorems also permit non-unit singleton A, including zero rows.
    for outer in ["x/(1-x)", "x^4", "x^2+x^4"] {
        for a in ["0", "2*z"] {
            let routes = exponential_row_routes(&ExponentialFamily {
                unmarked: expression(a),
                marked: expression(&format!("compose({outer},z)")),
            })
            .unwrap();
            assert_eq!(routes[0].evidence, ExponentialEvidenceLevel::TheoremBacked);
        }
    }
    for outer in ["x^3", "x^6", "x^2+3*x^3+x^4", "x^2+x^3+x^4"] {
        let routes = exponential_row_routes(&ExponentialFamily {
            unmarked: sizes(&[1]),
            marked: expression(&format!("compose({outer},exp(z)-1)")),
        })
        .unwrap();
        assert_eq!(routes[0].route, "polynomial_substitution_paper");
        assert_eq!(routes[0].evidence, ExponentialEvidenceLevel::PaperProof);
        assert!(!routes[0].consecutive_weak_interlacing);
    }
    for outer in ["-x^2+x^4", "x^2-x^4", "x^2+x^4+x^6"] {
        assert!(exponential_row_routes(&ExponentialFamily {
            unmarked: sizes(&[1]),
            marked: expression(&format!("compose({outer},exp(z)-1)")),
        })
        .unwrap()
        .is_empty());
    }
}
#[test]
fn expressions_check_full_identities_and_sign_hypotheses() {
    for b in [
        "exp(z)-1",
        "exp(z)-1-z",
        "compose(exp(x)-1-x,z^2)",
        "compose(x*(1+x),z^2)",
        "exp(z^2)-1-z^2",
        "exp(z^2+z^3)-1",
    ] {
        let f = ExponentialFamily {
            unmarked: expression("z"),
            marked: expression(b),
        };
        assert!(!exponential_row_routes(&f).unwrap().is_empty(), "{b}");
    }
    // A late perturbation cannot pass a parameter fit or a short prefix check.
    for b in [
        "exp(z)-1+z^20",
        "compose(exp(x)-1+x^20,z^2)",
        "log(1/(1-z))+z^20",
        "exp(z^2)-1+z",
        "exp(z^2)-1-z^3",
    ] {
        assert!(
            exponential_row_routes(&ExponentialFamily {
                unmarked: expression("z"),
                marked: expression(b)
            })
            .unwrap()
            .is_empty(),
            "{b}"
        );
    }
    let f = ExponentialFamily {
        unmarked: expression("z-z^2/3"),
        marked: expression("exp(z)-1"),
    };
    let report =
        check_exponential_rows(&request(f, vec![0, 2, 4], ExponentialRowCheck::RealRooted))
            .unwrap();
    assert!(report.routes.is_empty());
    assert_eq!(report.coefficient_signs[0].status, "negative");
    assert_eq!(
        report.coefficient_signs[0]
            .first_negative_coefficient
            .as_deref(),
        Some("-1/3")
    );
    assert_eq!(report.coefficient_signs[0].first_negative_degree, Some(2));
    let f = ExponentialFamily {
        unmarked: expression("z"),
        marked: expression("exp(z)-1-z^20"),
    };
    let report = check_exponential_rows(&request(f, vec![3], ExponentialRowCheck::None)).unwrap();
    assert_eq!(report.coefficient_signs[1].status, "unknown");
    assert!(report.routes.is_empty());
    for b in ["1+z", "exp(1)", "compose(z,1+z)", "log(z)", "1/z"] {
        assert!(expression(b).block_weights(0).is_err(), "{b}");
    }
    assert_eq!(
        expression("coeffs:exponential:[0,1,2]")
            .block_weights(3)
            .unwrap(),
        vec![q(0), q(1), q(2), q(0)]
    );
    assert_eq!(
        expression("coeffs:[0,1,2]").block_weights(3).unwrap(),
        vec![q(0), q(1), q(4), q(0)]
    );
}

#[test]
fn expression_preserves_supplied_polynomial_substitution_witness() {
    let family = ExponentialFamily {
        unmarked: expression("z"),
        marked: expression("compose(x*(1+x)^2,z^2+z^3+z^4)"),
    };
    let routes = exponential_row_routes(&family).unwrap();
    assert_eq!(routes[0].route, "polynomial_substitution_paper");
    assert_eq!(routes[0].evidence, ExponentialEvidenceLevel::PaperProof);
    assert_eq!(routes[0].premises[0].route, "nonnegative_quartic_seed");
    assert_eq!(
        family.rows(16).unwrap(),
        ExponentialFamily {
            unmarked: expression("z"),
            marked: expression("(z^2+z^3+z^4)*(1+z^2+z^3+z^4)^2"),
        }
        .rows(16)
        .unwrap()
    );
}
