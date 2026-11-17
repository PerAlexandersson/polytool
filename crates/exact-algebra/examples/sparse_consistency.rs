//! Construct and solve a small sparse prime-field system.

use exact_algebra::sparse_modular::{sparse_modular_linear_system_consistency, SparseModRow};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let rows = [
        SparseModRow::from_entries([(0, 1), (1, 1)], 5, 101)?,
        SparseModRow::from_entries([(0, 1)], 2, 101)?,
    ];
    let result = sparse_modular_linear_system_consistency(rows, 2, 101)?;
    println!(
        "consistent={} rank={} solution={:?}",
        result.consistent, result.rank, result.solution
    );
    Ok(())
}
