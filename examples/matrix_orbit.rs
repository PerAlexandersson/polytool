//! A supplied stationary state system with an exact scalar recurrence.

use polytool::recurrence::{check_matrix_orbit, verify_matrix_orbit, MatrixOrbitRequest};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // M = [[1,1],[x,1]], seed (1,0), projection (1,1).
    // Its scalar rows have P0=1, P1=1+x and P(n+2)=2P(n+1)+(x-1)P(n).
    let request: MatrixOrbitRequest = serde_json::from_value(serde_json::json!({
        "recurrence": {
            "schema": "polytool.recurrence.v1", "first_index": 0,
            "initial_polynomials": [["1"], ["1", "1"]],
            "recurrence": {"terms": [
                {"offset":1, "deriv_order":0, "sign":"none",
                 "coeff":{"coeffs":[["2"]]}},
                {"offset":2, "deriv_order":0, "sign":"none",
                 "coeff":{"coeffs":[["-1","1"]]}}
            ]}
        },
        "state_system": {
            "matrix": [[["1"],["1"]], [["0","1"],["1"]]],
            "start_vector": [["1"],["0"]], "projection": ["1","1"]
        }
    }))?;
    let report = check_matrix_orbit(&request)?;
    assert_eq!(report.status, "proved");
    verify_matrix_orbit(&report)?;
    println!("{}", report.reason);
    Ok(())
}
