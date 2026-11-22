//! Exact refinement CLI discovery, replay and tamper rejection.
use crate::common::{run_polytool_snapshot_with_deadline as run, stdout};
use polytool::recurrence::RecurrenceJson;

#[test]
fn discovers_and_replays_the_derivative_baseline() {
    let entry = polytool::oeis::by_id("A124324").unwrap();
    let (rule, seeds) = entry.recurrence_parts().unwrap();
    let record =
        RecurrenceJson::from_recurrence_rational(&rule, entry.recurrence_first_index, &seeds, None);
    let request =
        serde_json::json!({"recurrence":record,"limits":{"maximum_candidates":1}}).to_string();
    let output = run(&["interlacing-refine", "--compact"], &request);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let mut report: serde_json::Value = serde_json::from_str(&stdout(&output)).unwrap();
    assert_eq!(report["status"], "proved");
    let replay = run(&["interlacing-refine", "--verify"], &report.to_string());
    assert!(replay.status.success());
    report["candidates"][0]["identity"] = serde_json::json!("finite_fit");
    let invalid = run(&["interlacing-refine", "--verify"], &report.to_string());
    assert_eq!(invalid.status.code(), Some(2));
}
