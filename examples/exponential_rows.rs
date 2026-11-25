//! Sparse labelled-block rows and exact in-process checking.
use polytool::exponential_rows::{
    check_exponential_rows, BlockWeightSpec, ExponentialFamily, ExponentialRowCheck,
    ExponentialRowsRequest,
};
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let request = ExponentialRowsRequest {
        family: ExponentialFamily {
            unmarked: BlockWeightSpec::Sizes { sizes: vec![1, 2] },
            marked: BlockWeightSpec::Tail { minimum: 3 },
        },
        rows: vec![6, 9, 12],
        check: ExponentialRowCheck::RealRooted,
        summary: false,
        finite_rank: vec![],
    };
    let report = check_exponential_rows(&request)?;
    assert!(report.first_failure.is_none());
    for row in report.rows {
        println!("n={}: {:?}", row.index, row.coefficients);
    }
    Ok(())
}
