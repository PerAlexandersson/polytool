//! Raw JSON acceptance and serialization corpus, captured before the wire refactor.
//! Intentional certificate additions extend the unknown-tag diagnostic.

use polytool::recurrence::{RealRootednessTheoremReport, RecurrenceError};
use serde::{Deserialize, Serialize};

#[derive(Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Case {
    name: String,
    input: String,
    output: Option<String>,
    pretty_output: Option<String>,
    parse_error: Option<String>,
    verify_error: Option<String>,
}

fn outcome(name: String, input: String) -> Case {
    match serde_json::from_slice::<RealRootednessTheoremReport>(input.as_bytes()) {
        Ok(report) => {
            let from_text: RealRootednessTheoremReport = serde_json::from_str(&input).unwrap();
            assert_eq!(report, from_text);
            Case {
                name,
                input,
                output: Some(serde_json::to_string(&report).unwrap()),
                pretty_output: Some(serde_json::to_string_pretty(&report).unwrap()),
                parse_error: None,
                verify_error: report
                    .verify()
                    .err()
                    .map(|e: RecurrenceError| e.to_string()),
            }
        }
        Err(error) => Case {
            name,
            input,
            output: None,
            pretty_output: None,
            parse_error: Some(error.to_string()),
            verify_error: None,
        },
    }
}

#[test]
fn raw_report_bytes_and_rejection_stages_match_the_pre_wire_corpus() {
    let corpus = include_str!("../fixtures/theorem-report-wire.jsonl");
    let mut count = 0;
    for line in corpus.lines() {
        let expected: Case = serde_json::from_str(line).unwrap();
        let actual = outcome(expected.name.clone(), expected.input.clone());
        assert_eq!(actual, expected, "{}", expected.name);
        count += 1;
    }
    assert_eq!(count, 54, "the complete raw corpus must be present");
}
