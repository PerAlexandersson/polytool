//! Shared-kernel compatibility at the library and CLI boundaries.

use crate::common::run_polytool;

#[test]
fn sparse_rows_retain_the_public_polytool_path_and_all_errors() {
    use exact_algebra::sparse_modular::SparseModRow;
    let rows = [SparseModRow::from_entries([(0, 1)], 3, 101).unwrap()];
    let result: polytool::linalg::SparseModEliminationResult =
        polytool::linalg::sparse_modular_linear_system_consistency(rows, 1, 101).unwrap();
    assert!(result.consistent);
    assert_eq!(result.solution, Some(vec![3]));
    assert!(matches!(
        polytool::linalg::sparse_modular_linear_system_consistency(Vec::new(), 0, 9),
        Err(polytool::linalg::SparseModEliminationError::ModulusNotPrime { modulus: 9 })
    ));
}

#[test]
fn cli_modular_prefilter_and_exact_fallback_keep_the_same_record() {
    let input =
        "[1]\n[1,1]\n[1,2,1]\n[1,3,3,1]\n[1,4,6,4,1]\n[1,5,10,10,5,1]\n[1,6,15,20,15,6,1]\n";
    let modular = run_polytool(&["recurrence", "--json", "--modular-prefilter"], input);
    let exact = run_polytool(&["recurrence", "--json", "--no-modular-prefilter"], input);
    assert!(modular.status.success());
    assert!(exact.status.success());
    assert_eq!(modular.stdout, exact.stdout);
    let malformed = run_polytool(&["recurrence", "--json"], "[1/nope]\n");
    assert_eq!(malformed.status.code(), Some(2));
    assert!(malformed.stdout.is_empty());
}
