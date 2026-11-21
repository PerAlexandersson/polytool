//! Effective-BKW CLI, recurrence-discovery integration and replay binding.

use crate::common::run_polytool as run;
use polytool::recurrence::{
    check_effective_bkw_theorem, discover_recurrence_proofs, prove_parent_real_rootedness,
    ParentRecurrenceRecord, RealRootednessTheoremOptions, RecurrenceJson, VeroneseParentProofRoute,
    VeroneseParentProofStatus,
};
use serde_json::{json, Value};

fn tiling_record() -> RecurrenceJson {
    let terms = [(1, vec![1, 1]), (2, vec![0, 2, 1]), (3, vec![0, 0, -1, -1])];
    serde_json::from_value(json!({"schema":"polytool.recurrence.v1","first_index":0,
        "initial_polynomials":[["1"],["1"],["1","3","1"]],
        "recurrence":{"terms":terms.iter().map(|(offset,c)|json!({"offset":offset,"deriv_order":0,"sign":"none","coeff":{"coeffs":[c.iter().map(ToString::to_string).collect::<Vec<_>>()]}})).collect::<Vec<_>>()}
    })).unwrap()
}
#[test]
fn effective_cli_and_report_replay() {
    let request = json!({"generating_function":{"schema":"polytool.generating-function.v1","first_index":0,"numerator":[["1"]],"denominator":[["1"],["0","-2"],["1"]]},"interval":{"lower":{"kind":"rational","value":"-1"},"upper":{"kind":"rational","value":"1"}},"target_zero_count":{"slope":1,"offsets":[0]}});
    let result = run(&["bkw", "effective", "--json"], &request.to_string());
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let mut report: Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(report["status"], "proved");
    let replay = run(
        &["bkw", "effective", "--verify", "--json"],
        &report.to_string(),
    );
    assert!(replay.status.success());
    assert_eq!(report["certificate"]["cutoff"], json!(0));
    assert_eq!(report["certificate"]["domination_cutoff"], json!(0));
    report["certificate"]["cutoff"] = json!(1);
    assert_eq!(
        run(&["bkw", "effective", "--verify"], &report.to_string())
            .status
            .code(),
        Some(2)
    );
    assert_eq!(
        run(
            &["bkw", "effective", "--request", "-", "--recurrence", "-"],
            ""
        )
        .status
        .code(),
        Some(2)
    );
}

#[test]
fn prefix_shift_failures_are_replayable_not_proved_reports() {
    for (prefix, second, reason) in [
        (
            json!(["1", "0", "1"]),
            json!(["1"]),
            "not simple real-rooted",
        ),
        (
            json!(["1", "2", "1"]),
            json!(["1"]),
            "not simple real-rooted",
        ),
        (
            json!(["0"]),
            json!(["1", "0", "1"]),
            "not simple real-rooted",
        ),
    ] {
        let record = json!({"schema":"polytool.recurrence.v1","first_index":0,
        "initial_polynomials":[prefix,second,["1"],["1"]],
        "recurrence":{"terms":[
            {"offset":1,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["1"]]}},
            {"offset":2,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["0","1"]]}}
        ]}});
        let result = run(
            &["bkw", "effective", "--recurrence", "-", "--json"],
            &record.to_string(),
        );
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        let report: Value = serde_json::from_slice(&result.stdout).unwrap();
        assert_eq!(report["status"], "not_proved");
        assert!(
            report["reason"].as_str().unwrap().contains(reason),
            "{report}"
        );
        assert!(report["certificate"].is_null());
        let replay = run(
            &["bkw", "effective", "--verify", "--json"],
            &report.to_string(),
        );
        assert!(
            replay.status.success(),
            "{}",
            String::from_utf8_lossy(&replay.stderr)
        );
        let mut altered = report;
        altered["source_recurrence"]["initial_polynomials"][0] = json!(["2"]);
        assert!(!run(
            &["bkw", "effective", "--verify", "--json"],
            &altered.to_string()
        )
        .status
        .success());
    }
}
#[test]
fn tiling_route_is_discovered_and_replays_the_same_family() {
    let record = tiling_record();
    let (rule, first, seeds) = record.to_recurrence_parts().unwrap();
    let theorem = check_effective_bkw_theorem(&rule, first, &seeds, &Default::default()).unwrap();
    assert!(theorem.is_proved(), "{}", theorem.normal_form);
    theorem.verify().unwrap();
    let report = discover_recurrence_proofs(
        &record,
        &RealRootednessTheoremOptions {
            maximum_start_delay: 0,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(report.status, "proved");
    assert!(report
        .attempts
        .iter()
        .any(|a| a.theorem.as_str() == "effective_bkw"));
    report.verify().unwrap();
    let parent = prove_parent_real_rootedness(
        &ParentRecurrenceRecord {
            rule: &rule,
            seeds: &seeds,
            first_index: first,
        },
        first,
        &[],
        true,
    )
    .unwrap();
    assert_eq!(parent.status, VeroneseParentProofStatus::Proved);
    assert_eq!(parent.route, Some(VeroneseParentProofRoute::EffectiveBkw));
    let negative = check_effective_bkw_theorem(
        &rule,
        first,
        &seeds,
        &RealRootednessTheoremOptions {
            root_region: Some(polytool::real_rootedness::RealInterval::at_most(
                polytool::BigRational::from_integer(0.into()),
            )),
            maximum_start_delay: 0,
        },
    )
    .unwrap();
    assert!(negative.is_proved());
    negative.verify().unwrap();
    let mut altered = theorem;
    altered.recurrence.first_index = 1;
    assert!(altered.verify().is_err());
    let result = run(
        &["bkw", "effective", "--recurrence", "-", "--json"],
        &serde_json::to_string(&record).unwrap(),
    );
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(
        serde_json::from_slice::<Value>(&result.stdout).unwrap()["status"],
        "proved"
    );
}
