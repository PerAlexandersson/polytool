//! The optional test hook is usable by consumers, including nested guards.
#![cfg(feature = "test-hooks")]
use exact_algebra::integer_elimination::ForcedFallback;
use exact_algebra::solve_linear_system;
use num_rational::BigRational;

#[test]
fn public_fallback_hook_preserves_solutions() {
    let a = vec![vec![BigRational::from_integer(3.into())]];
    let b = vec![BigRational::from_integer(2.into())];
    let reference = solve_linear_system(&a, &b);
    let guard = ForcedFallback::engage();
    assert_eq!(solve_linear_system(&a, &b), reference);
    {
        let _nested = ForcedFallback::engage();
        assert_eq!(solve_linear_system(&a, &b), reference);
    }
    assert_eq!(solve_linear_system(&a, &b), reference);
    drop(guard);
    assert_eq!(solve_linear_system(&a, &b), reference);
}
