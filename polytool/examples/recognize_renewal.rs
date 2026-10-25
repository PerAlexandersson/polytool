//! Read existing polytool recurrence JSON from stdin and report conditional recognition.
use polytool::recurrence::RecurrenceJson;
use polytool::renewal::recognize_weighted_renewal;
use std::io::{self, Read};
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut input = String::new();
    io::stdin().read_to_string(&mut input)?;
    let record: RecurrenceJson = serde_json::from_str(&input)?;
    let (recurrence, first_index, rows) = record.to_recurrence_parts()?;
    let result = match recognize_weighted_renewal(&recurrence) {
        Ok(family) => serde_json::json!({
            "recognized": true,
            "ell": family.ell(),
            "seed": family.seed().iter().map(ToString::to_string).collect::<Vec<_>>(),
            "max_interlacing_step": family.max_interlacing_step(),
            "canonical_initial_data": family.matches_initial_rows(first_index, &rows),
            "identity_status": "caller_must_establish_recurrence_for_all_indices",
            "formalization_status": "paper_proof_not_lean_checked"
        }),
        Err(reason) => serde_json::json!({"recognized": false, "reason": reason}),
    };
    println!("{}", serde_json::to_string_pretty(&result)?);
    Ok(())
}
