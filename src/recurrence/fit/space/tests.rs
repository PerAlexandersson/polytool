//! Tests for the enclosing module.

use super::*;
#[test]
fn two_singular_identities_combine_to_a_regular_one() {
    let cols = vec![
        RecurrenceMonomial {
            lag: 0,
            derivative_order: 0,
            index_degree: 0,
            variable_degree: 0,
        },
        RecurrenceMonomial {
            lag: 0,
            derivative_order: 0,
            index_degree: 1,
            variable_degree: 0,
        },
        RecurrenceMonomial {
            lag: 1,
            derivative_order: 0,
            index_degree: 0,
            variable_degree: 0,
        },
        RecurrenceMonomial {
            lag: 1,
            derivative_order: 0,
            index_degree: 1,
            variable_degree: 0,
        },
    ];
    let mut result = empty_result(1, cols.clone(), 2);
    for (i, values) in [[0, 1, 0, 1], [-1, 1, -1, 1]].into_iter().enumerate() {
        let v = values
            .into_iter()
            .map(|n| Q::from_integer(n.into()))
            .collect::<Vec<_>>();
        result
            .representatives
            .push(representative(&cols, &v, 1, i, "P"));
    }
    assert_eq!(
        result.representatives[0].singular_ranks,
        Some(vec!["0".into()])
    );
    assert_eq!(
        result.representatives[1].singular_ranks,
        Some(vec!["1".into()])
    );
    let combined = leader_coverage(&result, &Default::default())
        .unwrap()
        .unwrap();
    assert_eq!(combined.common_singular_ranks, Some(vec![]));
    let rec = combined.to_recurrence(1).unwrap();
    assert_eq!(
        rec.denominator.unwrap().coeffs,
        vec![
            vec![Q::from_integer(5.into())],
            vec![Q::from_integer((-6).into())],
            vec![Q::from_integer(2.into())]
        ]
    );
    // The backward target index is m=n+1: n^2+(n-1)^2 becomes 2m^2-6m+5.
}
