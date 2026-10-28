//! End-to-end tests of the exact index sign checker through the public API:
//! brute-force comparisons on pseudo-random rational functions, the
//! certificate replay from JSON, and the documented conventions.

use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};
use polytool::{
    decide_index_polynomial_sign, decide_index_rational_function_sign, BigRational, IndexSign,
    IndexSignCertificate, IndexSignCondition, IndexSignDecision, IndexSignMethod, IntegerRange,
};

fn q(value: i64) -> BigRational {
    BigRational::from_integer(BigInt::from(value))
}

fn qr(numerator: i64, denominator: i64) -> BigRational {
    BigRational::new(BigInt::from(numerator), BigInt::from(denominator))
}

const CONDITIONS: [IndexSignCondition; 6] = [
    IndexSignCondition::Positive,
    IndexSignCondition::Nonnegative,
    IndexSignCondition::Negative,
    IndexSignCondition::Nonpositive,
    IndexSignCondition::Nonzero,
    IndexSignCondition::Zero,
];

/// A linear congruential generator, enough for reproducible test inputs.
struct Generator(u64);

impl Generator {
    fn between(&mut self, low: i64, high: i64) -> i64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        low + ((self.0 >> 33) % (high - low + 1) as u64) as i64
    }
}

fn evaluate(p: &[BigRational], n: i64) -> BigRational {
    p.iter()
        .rev()
        .fold(BigRational::zero(), |acc, c| acc * q(n) + c)
}

/// `scale * prod (n - root)`.
fn from_roots(roots: &[BigRational], scale: BigRational) -> Vec<BigRational> {
    let mut p = vec![scale];
    for root in roots {
        let mut next = vec![BigRational::zero(); p.len() + 1];
        for (i, c) in p.iter().enumerate() {
            next[i + 1] += c;
            next[i] -= c * root;
        }
        p = next;
    }
    p
}

/// A polynomial with all real roots in `[-18, 18]`: built from roots (integers,
/// half-integers, repeated roots) times a linear-free positive or negative
/// quadratic `n^2 + c` with `c > 0` (no real roots).
fn random_polynomial(generator: &mut Generator) -> Vec<BigRational> {
    let mut roots = Vec::new();
    for _ in 0..generator.between(0, 4) {
        let root = qr(generator.between(-18, 18), generator.between(1, 2));
        if generator.between(0, 3) == 0 {
            roots.push(root.clone());
        }
        roots.push(root);
    }
    let scale = qr(
        generator.between(1, 5) * if generator.between(0, 1) == 0 { 1 } else { -1 },
        2,
    );
    let p = from_roots(&roots, scale);
    if generator.between(0, 2) == 0 {
        // Times n^2 + c, c > 0.
        let c = q(generator.between(1, 4));
        let mut product = vec![BigRational::zero(); p.len() + 2];
        for (i, a) in p.iter().enumerate() {
            product[i] += a * &c;
            product[i + 2] += a;
        }
        product
    } else {
        p
    }
}

fn brute_force(
    numerator: &[BigRational],
    denominator: &[BigRational],
    first: i64,
    last: Option<i64>,
    condition: IndexSignCondition,
) -> Option<(BigInt, Option<BigRational>)> {
    // Every root lies in [-18, 18], so scanning to max(first, 18) + 2 suffices.
    let end = last.unwrap_or(first.max(18) + 2);
    (first..=end).find_map(|n| {
        let d = evaluate(denominator, n);
        if d.is_zero() {
            return Some((BigInt::from(n), None));
        }
        let value = evaluate(numerator, n) / d;
        let sign = if value.is_positive() {
            IndexSign::Positive
        } else if value.is_negative() {
            IndexSign::Negative
        } else {
            IndexSign::Zero
        };
        (!condition.is_satisfied_by(sign)).then(|| (BigInt::from(n), Some(value)))
    })
}

#[test]
fn decisions_agree_with_brute_force_and_replay_from_json() {
    let mut generator = Generator(2024);
    for case in 0..250 {
        let numerator = random_polynomial(&mut generator);
        let denominator = if case % 3 == 0 {
            vec![q(1)]
        } else {
            random_polynomial(&mut generator)
        };
        let first = generator.between(-12, 12);
        let last = match generator.between(0, 2) {
            0 => None,
            _ => Some(first + generator.between(-1, 25)),
        };
        let range = IntegerRange {
            first: BigInt::from(first),
            last: last.map(BigInt::from),
        };
        for condition in CONDITIONS {
            let decision =
                decide_index_rational_function_sign(&numerator, &denominator, &range, condition);
            let expected = brute_force(&numerator, &denominator, first, last, condition);
            match (&decision, &expected) {
                (IndexSignDecision::Holds { certificate }, None) => {
                    let json = serde_json::to_string(certificate).unwrap();
                    let replayed: IndexSignCertificate = serde_json::from_str(&json).unwrap();
                    replayed.verify().unwrap();
                }
                (IndexSignDecision::FailsAt { index, value }, Some((n, v))) => {
                    assert_eq!((index, value), (n, v));
                }
                _ => panic!(
                    "case {case}: {numerator:?} / {denominator:?} on {range} {condition:?}: \
                     {decision:?} but brute force gives {expected:?}"
                ),
            }
        }
    }
}

#[test]
fn denominator_roots_inside_and_outside_the_range() {
    // (n + 1) / ((n - 4)(2n - 1)): the denominator vanishes at 4 (and 1/2).
    let numerator = [q(1), q(1)];
    let denominator = from_roots(&[q(4), qr(1, 2)], q(2));
    let inside = IntegerRange::closed(1, 10);
    for condition in CONDITIONS {
        assert_eq!(
            decide_index_rational_function_sign(&numerator, &denominator, &inside, condition),
            if matches!(
                condition,
                IndexSignCondition::Negative
                    | IndexSignCondition::Nonpositive
                    | IndexSignCondition::Nonzero
            ) {
                // Negative on [1, 3], then undefined at 4.
                IndexSignDecision::FailsAt {
                    index: BigInt::from(4),
                    value: None,
                }
            } else {
                IndexSignDecision::FailsAt {
                    index: BigInt::from(1),
                    value: Some(qr(-2, 3)),
                }
            }
        );
    }
    // On [5, +inf) the root 4 is outside: positive throughout.
    let outside = IntegerRange::starting_at(5);
    let decision = decide_index_rational_function_sign(
        &numerator,
        &denominator,
        &outside,
        IndexSignCondition::Positive,
    );
    let IndexSignDecision::Holds { certificate } = decision else {
        panic!("expected a proof");
    };
    certificate.verify().unwrap();
    // On [1, 3] the sign is negative: -2/3, -3/2, -4/5.
    let decision = decide_index_rational_function_sign(
        &numerator,
        &denominator,
        &IntegerRange::closed(1, 3),
        IndexSignCondition::Negative,
    );
    assert!(decision.holds(), "{decision:?}");
}

#[test]
fn huge_and_clustered_roots() {
    let big = BigInt::from(10).pow(30);
    let root = BigRational::from_integer(big.clone());
    // (n - 10^30)(n - 10^30 - 1/3)(n - 10^30 - 2/3) changes sign three times
    // inside (10^30, 10^30 + 1) and once at 10^30.
    let p = from_roots(
        &[root.clone(), &root + qr(1, 3), &root + qr(2, 3)],
        BigRational::one(),
    );
    let decision = decide_index_polynomial_sign(
        &p,
        &IntegerRange::starting_at(0),
        IndexSignCondition::Nonzero,
    );
    assert_eq!(decision.failing_index(), Some(&big));
    let decision = decide_index_polynomial_sign(
        &p,
        &IntegerRange::starting_at(&big + 1),
        IndexSignCondition::Positive,
    );
    let IndexSignDecision::Holds { certificate } = &decision else {
        panic!("{decision:?}");
    };
    certificate.verify().unwrap();
    let decision = decide_index_polynomial_sign(
        &p,
        &IntegerRange::closed(0, &big - 1),
        IndexSignCondition::Negative,
    );
    let IndexSignDecision::Holds { certificate } = &decision else {
        panic!("{decision:?}");
    };
    assert!(matches!(
        certificate.method,
        IndexSignMethod::SignSamples { .. }
    ));
    certificate.verify().unwrap();
    assert!(decision
        .describe("p(n)")
        .contains("p(n) < 0 for every integer n in [0, "));
}

/// The example of `doc/library.md`.
#[test]
fn library_documentation_example() {
    let q = |n: i64| BigRational::from_integer(BigInt::from(n));
    let decision = decide_index_rational_function_sign(
        &[q(7), q(-5), q(1)],
        &[q(-5), q(2)],
        &IntegerRange::starting_at(3),
        IndexSignCondition::Positive,
    );
    let IndexSignDecision::Holds { certificate } = &decision else {
        unreachable!()
    };
    certificate.verify().unwrap();
    let decision = decide_index_rational_function_sign(
        &[q(7), q(-5), q(1)],
        &[q(-5), q(2)],
        &IntegerRange::starting_at(0),
        IndexSignCondition::Positive,
    );
    assert_eq!(decision.failing_index(), Some(&BigInt::from(0)));
}
