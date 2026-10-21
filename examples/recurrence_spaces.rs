//! Complete lag/derivative trade-offs for A171608.
//! Source: sqrt-of-2/real-rooted-oeis-proofs, commit
//! 9bb3970e5ba7877b77572c066bd56fbf8f63f05c, ProofsOeis/A171608.lean.
//! Run: cargo run -p polytool --example recurrence_spaces > spaces.json
use polytool::recurrence::{
    find_recurrence_spaces, verify_recurrence_spaces, RecurrenceSpaceRequest,
};
use serde_json::json;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let shapes: Vec<_> = [(3, 0), (2, 1), (1, 2)]
        .into_iter()
        .map(|(length, derivative_order)| {
            json!({
                "kind": "box", "length": length, "derivative_order": derivative_order,
                "index_degree": 1, "variable_degree": 2, "leading_index_degree": 1
            })
        })
        .collect();
    let request: RecurrenceSpaceRequest = serde_json::from_value(json!({
        "source": {
            "coefficients": [["0", "0", "-1"], ["0"], ["0", "2"], ["0"]],
            "seeds": [["1"], ["2"], ["0", "2"], ["0", "3"]]
        }, "shapes": shapes
    }))?;
    let report = find_recurrence_spaces(&request)?;
    verify_recurrence_spaces(&request, &report)?;
    println!("{}", serde_json::to_string_pretty(&report)?);
    Ok(())
}
