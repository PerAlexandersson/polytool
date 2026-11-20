//! Exact all-index certificate for the 4-by-n tiling generating function.
//!
//! The JSON output can be replayed with `polytool bkw effective --verify`.

use polytool::bkw::effective::{certify_effective_bkw, EffectiveBkwRequest};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let request: EffectiveBkwRequest = serde_json::from_value(serde_json::json!({
        "generating_function": {
            "schema": "polytool.generating-function.v1", "first_index": 0,
            "numerator": [["1"], ["0", "-1"]],
            "denominator": [["1"], ["-1", "-1"], ["0", "-2", "-1"],
                            ["0", "0", "1", "1"]]
        },
        "target_zero_count": {"slope": 1, "offsets": [0, -1]}
    }))?;
    let report = certify_effective_bkw(&request)?;
    assert_eq!(report.status, "proved", "{}", report.reason);
    report.verify()?;
    println!("{}", serde_json::to_string_pretty(&report)?);
    Ok(())
}
