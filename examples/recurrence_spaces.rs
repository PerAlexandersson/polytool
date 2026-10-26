//! Complete lag/derivative trade-offs for A171608.
//! Source: the generated A171608 recurrence definition (rows checked against
//! the OEIS data).
//! Run: cargo run -p polytool --example recurrence_spaces > spaces.json
//! Add `-- hermite` for a polynomial-index source, `-- automatic` for n! X^n,
//! or `-- rescaled` for A007318 rows (1+X)^n after dividing coefficient k by k!.
//! All source generation uses ConstantLagSource or IndexedLagSource.
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
    let input = match std::env::args().nth(1).as_deref() {
        Some("hermite") => {
            json!({"indexed_source":{"coefficients":[[["-2"],["-2"]],[["0","2"]]],"seeds":[["1"],["0","2"]]},
            "include_leader_coverage":true,"shapes":[{"kind":"box","length":1,"derivative_order":1,"variable_degree":1}]})
        }
        Some("automatic") => {
            json!({"indexed_source":{"coefficients":[[["0","1"],["0","1"]]],"seeds":[["1"]]},
            "include_leader_coverage":true,"shapes":[{"kind":"automatic_variable_degree","length":1,"index_degree":1}]})
        }
        Some("rescaled") => {
            let derivative = json!({"kind":"box","length":1,"derivative_order":1,
                "index_degree":1,"variable_degree":1,"leading_index_degree":1});
            let three_term = json!({"kind":"box","length":2,"index_degree":1,
                "variable_degree":1,"leading_index_degree":1});
            json!({"source":{"coefficients":[["1","1"]],"seeds":[["1"]],"origin":0},
            "coefficient_rescaling":{"kind":"factorial_power","power":1},
            "shapes":[derivative, three_term]})
        }
        None => json!({
        "source": {
            "coefficients": [["0", "0", "-1"], ["0"], ["0", "2"], ["0"]],
            "seeds": [["1"], ["2"], ["0", "2"], ["0", "3"]]
        }, "shapes": shapes, "include_leader_coverage": true
        }),
        Some(_) => return Err("expected no argument, hermite, automatic, or rescaled".into()),
    };
    let request: RecurrenceSpaceRequest = serde_json::from_value(input)?;
    let report = find_recurrence_spaces(&request)?;
    verify_recurrence_spaces(&request, &report)?;
    println!("{}", serde_json::to_string_pretty(&report)?);
    Ok(())
}
