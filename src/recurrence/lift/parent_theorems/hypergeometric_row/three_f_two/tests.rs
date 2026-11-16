//! Exact parameter and parity tests for terminating three-over-two recognition.

use super::*;

fn q(n: i64, d: i64) -> Q {
    Q::new(n.into(), d.into())
}

fn ratio() -> Ratio {
    Ratio {
        scale: IndexFunction::constant(q(-4, 1)),
        numerators: vec![(q(-1, 1), q(-1, 1)), (q(1, 1), q(1, 2)), (q(3, 2), q(1, 2))],
        denominators: vec![affine_constant(0), affine_constant(1), affine_constant(2)],
    }
}

#[test]
fn weak_boundary_and_parity_swap_are_recognized() {
    let mut ratio = ratio();
    let matching = find_matching(&ratio, 1).unwrap();
    assert_eq!(matching.lower, [q(2, 1), q(3, 1)]);
    assert_eq!(matching.upper, [q(5, 2), q(3, 1)]);
    assert!(!matching.strict);
    // At degree 2m the pair is {m+5/2,m+2}; at degree 2m+1
    // it is {m+5/2,m+3}: the affine factors exchange theorem roles.
    ratio.numerators.reverse();
    ratio.denominators.reverse();
    assert_eq!(find_matching(&ratio, 1), Some(matching));
    assert!(find_matching(&ratio, 0).is_none());
}

#[test]
fn strict_range_and_either_scale_sign_are_recognized() {
    let mut ratio = ratio();
    ratio.denominators = vec![affine_constant(0), (q(-2, 3), q(0, 1)), (q(-1, 3), q(0, 1))];
    for scale in [q(-4, 27), q(4, 27)] {
        ratio.scale = IndexFunction::constant(scale.clone());
        let matching = find_matching(&ratio, 1).unwrap();
        assert!(matching.strict);
        assert_eq!(matching.scale, scale);
    }
}

#[test]
fn every_weak_inequality_is_required() {
    for (lower, upper) in [
        ([2, 1], [1, 3]), // a1 <= b1
        ([1, 3], [1, 3]), // a2 - 1 <= b1
        ([3, 1], [3, 2]), // a1 <= b2
        ([1, 3], [3, 2]), // a2 <= b2
    ] {
        assert_eq!(
            parameter_range(&lower.map(|x| q(x, 1)), &upper.map(|x| q(x, 1))),
            None
        );
    }
    for i in 0..4 {
        for bad in [0, -1] {
            let mut parameters = [1, 1, 2, 3];
            parameters[i] = bad;
            assert_eq!(
                parameter_range(
                    &[parameters[0], parameters[1]].map(|x| q(x, 1)),
                    &[parameters[2], parameters[3]].map(|x| q(x, 1))
                ),
                None
            );
        }
    }
}

#[test]
fn nonconstant_scale_wrong_slopes_and_inconsistent_parities_are_rejected() {
    let mut r = ratio();
    r.scale = IndexFunction::polynomial(vec![q(-4, 1), q(-1, 1)]);
    assert!(find_matching(&r, 1).is_none());
    r = ratio();
    r.scale = IndexFunction::constant(Q::zero());
    assert!(find_matching(&r, 1).is_none());
    r = ratio();
    r.numerators[1].1 = q(1, 3);
    assert!(find_matching(&r, 1).is_none());
    r = ratio();
    r.numerators[1].0 = q(4, 3);
    assert!(find_matching(&r, 1).is_none());
    r = ratio();
    r.denominators[1].0 = q(100, 1);
    assert!(find_matching(&r, 1).is_none());
    r = ratio();
    r.denominators[1].0 = q(-1, 1);
    assert!(find_matching(&r, 1).is_none());
    r = ratio();
    r.denominators[1].1 = q(1, 1);
    assert!(find_matching(&r, 1).is_none());
}

#[test]
fn degree_shift_and_single_parity_subsequence_match_exactly() {
    let mut r = ratio();
    // M=n-1, and affine upper parameters are M/2+2 and M/2+5/2.
    r.numerators[0].0 = q(0, 1);
    r.numerators[1].0 = q(1, 2);
    r.numerators[2].0 = q(1, 1);
    assert_eq!(
        find_matching(&r, 2).unwrap().termination,
        (q(-1, 1), q(1, 1))
    );
    // Restrict to even degrees M=2n: no exchange of the upper pair.
    r = ratio();
    r.numerators[0].1 = q(-2, 1);
    r.numerators[1].1 = q(1, 1);
    r.numerators[2].1 = q(1, 1);
    assert!(find_matching(&r, 1).is_some());
}

#[test]
fn a_strict_lower_parameter_order_is_preferred_to_a_weak_one() {
    let mut r = ratio();
    r.denominators = vec![affine_constant(0), (q(3, 2), q(0, 1)), affine_constant(0)];
    let matching = find_matching(&r, 1).unwrap();
    assert!(matching.strict);
    assert_eq!(matching.lower, [q(1, 1), q(5, 2)]);
}
