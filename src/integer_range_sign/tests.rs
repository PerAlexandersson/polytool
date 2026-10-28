//! Unit tests of the exact index sign checker: brute-force comparisons on
//! pseudo-random inputs, edge cases, budgets and certificate replay.

use super::*;

fn q(value: i64) -> Q {
    Q::from_integer(BigInt::from(value))
}

fn qr(numerator: i64, denominator: i64) -> Q {
    Q::new(BigInt::from(numerator), BigInt::from(denominator))
}

fn int(value: i64) -> BigInt {
    BigInt::from(value)
}

const CONDITIONS: [IndexSignCondition; 6] = [
    IndexSignCondition::Positive,
    IndexSignCondition::Nonnegative,
    IndexSignCondition::Negative,
    IndexSignCondition::Nonpositive,
    IndexSignCondition::Nonzero,
    IndexSignCondition::Zero,
];

/// A small deterministic generator (xorshift64*).
struct Generator(u64);

impl Generator {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_f491_4f6c_dd1d)
    }

    fn below(&mut self, bound: u64) -> u64 {
        self.next() % bound
    }

    fn between(&mut self, low: i64, high: i64) -> i64 {
        low + self.below((high - low + 1) as u64) as i64
    }
}

/// `prod (n - root)` times `scale`.
fn from_roots(roots: &[Q], scale: Q) -> Vec<Q> {
    let mut p = vec![scale];
    for root in roots {
        p = vec_poly::mul(&p, &[-root.clone(), Q::one()]);
    }
    p
}

/// A pseudo-random polynomial with an integer `B` such that every real root
/// lies in `(-B, B)`.
fn random_polynomial(generator: &mut Generator) -> (Vec<Q>, i64) {
    if generator.below(2) == 0 {
        let degree = generator.below(6) as usize;
        let p: Vec<Q> = (0..=degree)
            .map(|_| qr(generator.between(-6, 6), generator.between(1, 3)))
            .collect();
        let bound = cauchy_test_bound(&p);
        (p, bound)
    } else {
        // Integer roots, double roots and roots between consecutive integers.
        let count = generator.below(5) as usize;
        let mut roots = Vec::new();
        for _ in 0..count {
            let root = match generator.below(3) {
                0 => q(generator.between(-8, 12)),
                1 => qr(2 * generator.between(-8, 12) + 1, 2),
                _ => qr(generator.between(-30, 40), 3),
            };
            if generator.below(4) == 0 {
                roots.push(root.clone());
            }
            roots.push(root);
        }
        let scale = q(if generator.below(2) == 0 { 1 } else { -2 });
        (from_roots(&roots, scale), 15)
    }
}

/// `1 + sum |p_i / p_d|`, a root bound computed independently of the module.
fn cauchy_test_bound(p: &[Q]) -> i64 {
    let p = trimmed(p);
    if p.len() < 2 {
        return 0;
    }
    let leading = p.last().unwrap().abs();
    let sum = p.iter().fold(Q::zero(), |acc, c| acc + c.abs() / &leading);
    sum.ceil().to_integer().try_into().unwrap()
}

/// The smallest failing index by brute force, `None` when the condition
/// holds; an unbounded range is scanned up to `bound + 2`, beyond which no
/// root of either polynomial lies.
fn brute_force(
    numerator: &[Q],
    denominator: &[Q],
    first: i64,
    last: Option<i64>,
    condition: IndexSignCondition,
    bound: i64,
) -> Option<(BigInt, Option<Q>)> {
    let end = last.unwrap_or_else(|| first.max(bound) + 2);
    for n in first..=end {
        let value_q = vec_poly::evaluate(denominator, &q(n));
        if value_q.is_zero() {
            return Some((int(n), None));
        }
        let value = vec_poly::evaluate(numerator, &q(n)) / value_q;
        if !condition.is_satisfied_by(IndexSign::of(&value)) {
            return Some((int(n), Some(value)));
        }
    }
    None
}

fn range(first: i64, last: Option<i64>) -> IntegerRange {
    IntegerRange {
        first: int(first),
        last: last.map(int),
    }
}

fn check_against_brute_force(
    numerator: &[Q],
    denominator: &[Q],
    first: i64,
    last: Option<i64>,
    condition: IndexSignCondition,
    bound: i64,
) {
    let decision =
        decide_index_rational_function_sign(numerator, denominator, &range(first, last), condition);
    let expected = brute_force(numerator, denominator, first, last, condition, bound);
    match (&decision, expected) {
        (IndexSignDecision::Holds { certificate }, None) => {
            certificate.verify().unwrap_or_else(|error| {
                panic!("{numerator:?}/{denominator:?} on {first}..{last:?} {condition:?}: {error}")
            });
            let json = serde_json::to_string(&decision).unwrap();
            let back: IndexSignDecision = serde_json::from_str(&json).unwrap();
            assert_eq!(back, decision);
        }
        (IndexSignDecision::FailsAt { index, value }, Some((expected_index, expected_value))) => {
            assert_eq!(
                (index, value),
                (&expected_index, &expected_value),
                "{numerator:?}/{denominator:?} on {first}..{last:?} {condition:?}"
            );
        }
        (decision, expected) => panic!(
            "{numerator:?}/{denominator:?} on {first}..{last:?} {condition:?}: got \
             {decision:?}, brute force {expected:?}"
        ),
    }
}

#[test]
fn polynomials_agree_with_brute_force() {
    let mut generator = Generator(0x9e37_79b9_7f4a_7c15);
    for _ in 0..400 {
        let (p, bound) = random_polynomial(&mut generator);
        let first = generator.between(-12, 12);
        let last = match generator.below(3) {
            0 => None,
            1 => Some(first + generator.between(-2, 30)),
            _ => Some(first),
        };
        for condition in CONDITIONS {
            check_against_brute_force(&p, &[q(1)], first, last, condition, bound);
        }
    }
}

#[test]
fn rational_functions_agree_with_brute_force() {
    let mut generator = Generator(0x0123_4567_89ab_cdef);
    for _ in 0..300 {
        let (numerator, numerator_bound) = random_polynomial(&mut generator);
        let (mut denominator, denominator_bound) = random_polynomial(&mut generator);
        if vec_poly::is_zero(&denominator) {
            denominator = vec![q(3)];
        }
        let bound = numerator_bound.max(denominator_bound);
        let first = generator.between(-12, 12);
        let last = match generator.below(2) {
            0 => None,
            _ => Some(first + generator.between(-1, 30)),
        };
        for condition in CONDITIONS {
            check_against_brute_force(&numerator, &denominator, first, last, condition, bound);
        }
    }
}

#[test]
fn zero_polynomials_and_empty_ranges() {
    let zero: Vec<Q> = Vec::new();
    for condition in CONDITIONS {
        let decision = decide_index_polynomial_sign(&zero, &range(3, None), condition);
        let expected = matches!(
            condition,
            IndexSignCondition::Nonnegative
                | IndexSignCondition::Nonpositive
                | IndexSignCondition::Zero
        );
        assert_eq!(decision.holds(), expected, "{condition:?}");
        if !expected {
            assert_eq!(
                decision,
                IndexSignDecision::FailsAt {
                    index: int(3),
                    value: Some(Q::zero())
                }
            );
        }
        // An empty range satisfies everything, even a zero denominator.
        let empty =
            decide_index_rational_function_sign(&[q(-1)], &zero, &range(3, Some(2)), condition);
        assert!(empty.holds());
        // A zero denominator fails at the first index.
        assert_eq!(
            decide_index_rational_function_sign(&[q(1)], &[q(0)], &range(-4, None), condition),
            IndexSignDecision::FailsAt {
                index: int(-4),
                value: None
            }
        );
    }
    // 0 / (n - 5) on [0, 10]: zero except where undefined.
    assert_eq!(
        decide_index_rational_function_sign(
            &zero,
            &[q(-5), q(1)],
            &range(0, Some(10)),
            IndexSignCondition::Zero
        ),
        IndexSignDecision::FailsAt {
            index: int(5),
            value: None
        }
    );
}

#[test]
fn huge_roots_are_not_enumerated() {
    let big = BigInt::from(10).pow(30);
    let big_q = Q::from_integer(big.clone());
    // 10^30 - n.
    let p = vec![big_q.clone(), q(-1)];
    assert_eq!(
        decide_index_polynomial_sign(&p, &range(0, None), IndexSignCondition::Positive),
        IndexSignDecision::FailsAt {
            index: big.clone(),
            value: Some(Q::zero())
        }
    );
    assert_eq!(
        decide_index_polynomial_sign(&p, &range(0, None), IndexSignCondition::Nonnegative),
        IndexSignDecision::FailsAt {
            index: &big + 1,
            value: Some(q(-1))
        }
    );
    // (n - 10^30)^2 >= 0 everywhere, = 0 only at 10^30.
    let square = from_roots(&[big_q.clone(), big_q.clone()], q(1));
    let decision =
        decide_index_polynomial_sign(&square, &range(-5, None), IndexSignCondition::Nonnegative);
    match &decision {
        IndexSignDecision::Holds { certificate } => certificate.verify().unwrap(),
        other => panic!("{other:?}"),
    }
    assert_eq!(
        decide_index_polynomial_sign(&square, &range(-5, None), IndexSignCondition::Positive)
            .failing_index(),
        Some(&big)
    );
    // (n - 10^30 - 1/2)(n + 7) != 0 on [0, +inf), with a root between integers.
    let between = from_roots(&[big_q + qr(1, 2), q(-7)], q(3));
    let decision =
        decide_index_polynomial_sign(&between, &range(0, None), IndexSignCondition::Nonzero);
    match &decision {
        IndexSignDecision::Holds { certificate } => {
            certificate.verify().unwrap();
            let IndexSignMethod::SignSamples {
                numerator_root_cells,
                ..
            } = &certificate.method
            else {
                panic!("expected sign samples");
            };
            assert_eq!(
                numerator_root_cells,
                &vec![IndexRootCell {
                    lower: big.clone(),
                    root_count: 1
                }]
            );
        }
        other => panic!("{other:?}"),
    }
    // A huge negative first index.
    let decision = decide_index_polynomial_sign(
        &[q(1), q(0), q(1)],
        &IntegerRange::starting_at(-&big),
        IndexSignCondition::Positive,
    );
    assert!(decision.holds());
}

#[test]
fn taylor_shift_fast_path_and_its_replay() {
    // (n - 3)(n - 4) + 1 > 0 for n >= 4: Taylor shift at 4 is 1 + u + u^2.
    let p = vec_poly::add(&from_roots(&[q(3), q(4)], q(1)), &[q(1)]);
    let decision = decide_index_polynomial_sign(&p, &range(4, None), IndexSignCondition::Positive);
    let IndexSignDecision::Holds { certificate } = &decision else {
        panic!("{decision:?}");
    };
    assert!(matches!(
        certificate.method,
        IndexSignMethod::TaylorShift { .. }
    ));
    certificate.verify().unwrap();
    let mut tampered = (**certificate).clone();
    tampered.condition = IndexSignCondition::Negative;
    assert!(matches!(
        tampered.verify(),
        Err(IndexSignCertificateError::MethodDoesNotApply(_))
    ));
    let mut tampered = (**certificate).clone();
    tampered.range.first = int(0);
    assert!(tampered.verify().is_err());
}

#[test]
fn tampered_sample_certificates_are_rejected() {
    // (n - 2)(2n - 7) >= 0 fails at 3 only; (n - 2)^2 (n - 9/2)^2 >= 0 holds by samples.
    let p = from_roots(&[q(2), q(2), qr(9, 2), qr(9, 2)], q(1));
    let decision =
        decide_index_polynomial_sign(&p, &range(0, Some(40)), IndexSignCondition::Nonnegative);
    let IndexSignDecision::Holds { certificate } = &decision else {
        panic!("{decision:?}");
    };
    certificate.verify().unwrap();
    let IndexSignMethod::SignSamples {
        samples,
        numerator_root_cells,
        denominator_root_cells,
    } = &certificate.method
    else {
        panic!("expected sign samples");
    };
    assert!(samples
        .iter()
        .any(|s| s.index == int(2) && s.numerator_sign == IndexSign::Zero));
    assert_eq!(numerator_root_cells.len(), 1);
    assert!(denominator_root_cells.is_empty());
    let rebuild = |samples: Vec<IndexSignSample>, cells: Vec<IndexRootCell>| {
        let mut tampered = (**certificate).clone();
        tampered.method = IndexSignMethod::SignSamples {
            samples,
            numerator_root_cells: cells,
            denominator_root_cells: Vec::new(),
        };
        tampered.verify()
    };
    // A wrong sign.
    let mut wrong = samples.clone();
    wrong[0].numerator_sign = IndexSign::Negative;
    assert!(matches!(
        rebuild(wrong, numerator_root_cells.clone()),
        Err(IndexSignCertificateError::SampleSign { .. })
    ));
    // A dropped sample hides a root between samples.
    let kept: Vec<IndexSignSample> = samples
        .iter()
        .filter(|s| s.index == int(0) || s.index == int(40))
        .cloned()
        .collect();
    assert!(matches!(
        rebuild(kept, numerator_root_cells.clone()),
        Err(IndexSignCertificateError::RootBetweenSamples { .. })
    ));
    // A wrong root cell.
    let mut cells = numerator_root_cells.clone();
    cells[0].lower += 1;
    assert!(matches!(
        rebuild(samples.clone(), cells),
        Err(IndexSignCertificateError::RootCell { .. })
    ));
    // A stronger condition fails on replay.
    let mut stronger = (**certificate).clone();
    stronger.condition = IndexSignCondition::Positive;
    assert_eq!(
        stronger.verify(),
        Err(IndexSignCertificateError::ConditionFails { index: int(2) })
    );
    // An unbounded copy needs a nonzero last sample with no roots beyond.
    let mut unbounded = (**certificate).clone();
    unbounded.range.last = None;
    unbounded.numerator = vec_poly::mul(&p, &[q(-50), q(1)]);
    assert!(unbounded.verify().is_err());
}

#[test]
fn budgets_give_undecided() {
    let p = from_roots(&[q(1), q(2), q(3)], q(1));
    let tight_degree = IndexSignBudget {
        maximum_degree: 2,
        ..IndexSignBudget::default()
    };
    assert!(matches!(
        decide_index_rational_function_sign_with_budget(
            &p,
            &[q(1)],
            &range(0, None),
            IndexSignCondition::Nonnegative,
            &tight_degree
        ),
        IndexSignDecision::Undecided { .. }
    ));
    let no_steps = IndexSignBudget {
        maximum_bisection_steps: 0,
        ..IndexSignBudget::default()
    };
    assert!(matches!(
        decide_index_rational_function_sign_with_budget(
            &p,
            &[q(1)],
            &range(0, None),
            IndexSignCondition::Nonnegative,
            &no_steps
        ),
        IndexSignDecision::Undecided { .. }
    ));
    // The Taylor-shift fast path needs no bisection.
    assert!(decide_index_rational_function_sign_with_budget(
        &p,
        &[q(1)],
        &range(3, None),
        IndexSignCondition::Nonnegative,
        &no_steps
    )
    .holds());
}

#[test]
fn json_form_uses_exact_strings() {
    let decision = decide_index_rational_function_sign(
        &[q(1), q(-2)],
        &[q(3)],
        &range(0, None),
        IndexSignCondition::Nonnegative,
    );
    let value = serde_json::to_value(&decision).unwrap();
    assert_eq!(
        value,
        serde_json::json!({"decision": "fails_at", "index": "1", "value": "-1/3"})
    );
    let decision = decide_index_polynomial_sign(
        &[qr(1, 2)],
        &range(-3, Some(7)),
        IndexSignCondition::Positive,
    );
    let value = serde_json::to_value(&decision).unwrap();
    assert_eq!(value["decision"], "holds");
    assert_eq!(
        value["certificate"]["numerator"],
        serde_json::json!(["1/2"])
    );
    assert_eq!(
        value["certificate"]["range"],
        serde_json::json!({"first": "-3", "last": "7"})
    );
    assert_eq!(value["certificate"]["method"]["kind"], "taylor_shift");
    let back: IndexSignDecision = serde_json::from_value(value).unwrap();
    assert_eq!(back, decision);
    assert!(serde_json::from_str::<IntegerRange>(r#"{"first": "x"}"#).is_err());
    assert_eq!(
        serde_json::from_str::<IntegerRange>(r#"{"first": "5"}"#).unwrap(),
        IntegerRange::starting_at(5)
    );
}

#[test]
fn root_domain_helpers_delegate_to_the_exact_checker() {
    use crate::recurrence::{bivariate_leader_nonvanishing, index_polynomial_sign_from, BivarPoly};
    use std::cmp::Ordering;
    // n^2 - 5n + 7 > 0 at every integer (its minimum 3/4 is at n = 5/2), but
    // the Taylor shift at 0 has a negative coefficient.
    assert_eq!(
        index_polynomial_sign_from(&[q(7), q(-5), q(1)], 0),
        Some(Ordering::Greater)
    );
    assert_eq!(
        index_polynomial_sign_from(&[q(-7), q(5), q(-1)], 0),
        Some(Ordering::Less)
    );
    // (n - 2)(n - 3) vanishes at integers; 2n - 5 changes sign between them.
    assert_eq!(index_polynomial_sign_from(&[q(6), q(-5), q(1)], 0), None);
    assert_eq!(index_polynomial_sign_from(&[q(-5), q(2)], 0), None);
    // f(m, x) = (2m - 5) + (2m - 5) x: nonzero at every integer m, with no
    // coefficient of constant sign; (m - 3)(1 + x) vanishes at m = 3.
    let f = BivarPoly::from_coefficients_normalized(vec![vec![q(-5), q(-5)], vec![q(2), q(2)]]);
    assert!(bivariate_leader_nonvanishing(&f, 0));
    let g = BivarPoly::from_coefficients_normalized(vec![vec![q(-3), q(-3)], vec![q(1), q(1)]]);
    assert!(!bivariate_leader_nonvanishing(&g, 0));
    assert!(bivariate_leader_nonvanishing(&g, 4));
}
