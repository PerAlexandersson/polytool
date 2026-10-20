//! cargo run -p polytool --example recurrence_representations
//! The canonical OEIS generator is the only source of the input rows.
use polytool::recurrence::*;

fn main() {
    let rows: Vec<Vec<BigRational>> = polytool::oeis::by_id("A171608")
        .unwrap()
        .generate_rows(40)
        .unwrap()
        .into_iter()
        .map(|row| row.into_iter().map(BigRational::from_integer).collect())
        .collect();
    let search = AdaptiveSearchOptions {
        max_rec_len: 3,
        max_diff_deg: 2,
        max_var_deg: 3,
        max_idx_deg: 1,
        try_denominator: true,
        max_denom_idx_deg: 1,
        max_denom_var_deg: 0,
        ..Default::default()
    };
    let report = find_recurrence_representations_rational(
        &rows,
        &search,
        AdaptiveSearchBudget::limited(2000),
    )
    .unwrap();
    println!("Finite exact fits; {:?} search.", report.status);
    for &i in &report.frontier {
        let r = &report.representations[i];
        println!(
            "length {}, derivative order {}, checked from n={}: {}",
            r.complexity.rec_len, r.complexity.diff_deg, r.first_checked_index, r.recurrence
        );
    }
}
