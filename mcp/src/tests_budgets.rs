//! Phase-one response budgets and compatibility checks. Remaining report tools
//! acquire summaries and their budget rows in the next phase.

use crate::report_store::tests::TempDirectory;
use crate::*;
use rmcp::{handler::server::wrapper::Parameters, Json};
use serde_json::{json, Value};

fn discover(server: &PolynomialToolsServer, request: Value) -> Value {
    let Json(value) = server
        .discover_recurrence_proofs(Parameters(serde_json::from_value(request).unwrap()))
        .unwrap();
    serde_json::to_value(value).unwrap()
}

fn fingerprint(bytes: &[u8]) -> String {
    let mut hash = 0xcbf2_9ce4_8422_2325u64;
    for &byte in bytes {
        hash ^= u64::from(byte);
        hash = hash.wrapping_mul(0x0100_0000_01b3);
    }
    format!("{hash:016x}:{}", bytes.len())
}

#[test]
fn default_budget_table_and_catalog_proof_verify_chain() {
    let temp = TempDirectory::new();
    let server = PolynomialToolsServer::with_report_store(temp.config());
    let Json(list) = server
        .list_oeis_sequences(Parameters(Default::default()))
        .unwrap();
    let Json(get) = server
        .get_oeis_sequence(Parameters(
            serde_json::from_value(json!({"id":"A001263"})).unwrap(),
        ))
        .unwrap();
    let Json(find) = server.find_recurrence(Parameters(serde_json::from_value(json!({
        "coefficients":[[1],[2],[4],[8],[16]], "max_variable_degree":0,"max_index_degree":0,"max_derivative_order":0,"max_recurrence_length":1
    })).unwrap())).unwrap();
    let proof = discover(&server, json!({"oeis_id":"A001263"}));
    let Json(verify) = server
        .verify_recurrence_proofs(Parameters(
            serde_json::from_value(json!({"report_id":proof["report"]["id"]})).unwrap(),
        ))
        .unwrap();
    let Json(read) = server
        .read_report(Parameters(
            serde_json::from_value(json!({"report_id":proof["report"]["id"]})).unwrap(),
        ))
        .unwrap();
    for (tool, value, budget) in [
        (
            "list_oeis_sequences",
            serde_json::to_value(list).unwrap(),
            8192,
        ),
        (
            "get_oeis_sequence",
            serde_json::to_value(get).unwrap(),
            8192,
        ),
        ("find_recurrence", serde_json::to_value(find).unwrap(), 8192),
        ("discover_recurrence_proofs", proof.clone(), 2048),
        (
            "verify_recurrence_proofs",
            serde_json::to_value(&verify).unwrap(),
            8192,
        ),
        ("read_report", serde_json::to_value(read).unwrap(), 2048),
    ] {
        let bytes = serde_json::to_vec(&value).unwrap().len();
        assert!(bytes <= budget, "{tool}: {bytes} > {budget}");
    }
    let Json(catalog) = server
        .list_oeis_sequences(Parameters(
            serde_json::from_value(json!({"query":"A001263"})).unwrap(),
        ))
        .unwrap();
    let chain = serde_json::to_vec(&catalog).unwrap().len()
        + serde_json::to_vec(&proof).unwrap().len()
        + serde_json::to_vec(&verify).unwrap().len();
    assert!(chain < 1536, "chain uses {chain} bytes");
}

#[test]
fn full_is_unchanged_store_matches_frozen_fingerprint_and_replay_agrees() {
    let temp = TempDirectory::new();
    let server = PolynomialToolsServer::with_report_store(temp.config());
    for id in ["A001263", "A008292", "A003982"] {
        let summary = discover(&server, json!({"oeis_id":id}));
        let full = discover(&server, json!({"oeis_id":id,"detail":"full"}));
        assert!(full.get("report").is_none());
        let stored = std::fs::read(summary["report"]["path"].as_str().unwrap()).unwrap();
        let fixture = include_str!("../../tests/fixtures/discovery-report-fingerprints.tsv");
        let recorded = fixture
            .lines()
            .find(|line| line.starts_with(&format!("{id}\t")))
            .unwrap()
            .split('\t')
            .nth(1)
            .unwrap();
        assert_eq!(fingerprint(&stored), recorded, "{id}");
        let value: Value = serde_json::from_slice(&stored).unwrap();
        assert_eq!(full, value);
        // The old endpoint used exactly this BTreeMap/value serialization.
        assert_eq!(
            serde_json::to_string(&full).unwrap(),
            serde_json::to_string(&value).unwrap()
        );
        let Json(inline) = server
            .verify_recurrence_proofs(Parameters(
                serde_json::from_value(json!({"report_json":full.to_string()})).unwrap(),
            ))
            .unwrap();
        let restarted = PolynomialToolsServer::with_report_store(temp.config());
        let Json(handle) = restarted
            .verify_recurrence_proofs(Parameters(
                serde_json::from_value(json!({"report_id":summary["report"]["id"]})).unwrap(),
            ))
            .unwrap();
        assert_eq!(inline, handle);
        assert_eq!(inline["verified"], true);
    }
}

#[test]
fn projections_source_exclusivity_disabled_store_and_invalid_bounds() {
    let server = PolynomialToolsServer::with_report_store(Default::default());
    let summary = discover(
        &server,
        json!({"oeis_id":"A001263","select":["/attempts/*/status"]}),
    );
    assert!(summary.get("report").is_none());
    assert!(summary["selected"]["/attempts/0/status"].is_string());
    for request in [
        json!({}),
        json!({"oeis_id":"A001263","recurrence_json":"{}"}),
        json!({"oeis_id":"A001263","select":["/*/*"]}),
        json!({"oeis_id":"A001263","detail":"full","select":[]}),
    ] {
        assert!(server
            .discover_recurrence_proofs(Parameters(serde_json::from_value(request).unwrap()))
            .is_err());
    }
    for request in [json!({}), json!({"report_id":"bad","report_json":"{}"})] {
        assert!(server
            .verify_recurrence_proofs(Parameters(serde_json::from_value(request).unwrap()))
            .is_err());
    }
    let temp = TempDirectory::new();
    let server = PolynomialToolsServer::with_report_store(temp.config());
    let summary = discover(&server, json!({"oeis_id":"A001263"}));
    let Json(read) = server.read_report(Parameters(serde_json::from_value(json!({"report_id":summary["report"]["id"],"select":["/attempts/*"],"max_bytes":512})).unwrap())).unwrap();
    assert_eq!(read["truncated"], true);
    assert!(serde_json::to_vec(&read).unwrap().len() <= 512);
}

#[test]
fn exact_refutation_is_opt_in_and_never_shown_as_a_proof() {
    let server = PolynomialToolsServer::with_report_store(Default::default());
    let record = json!({"schema":"polytool.recurrence.v1","first_index":0,
        "initial_polynomials":[["1"]],"recurrence":{"terms":[
        {"offset":1,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["1","1","1"]]}}]}});
    let mut request = json!({"recurrence_json":record.to_string(),"max_representation_candidates":0,
        "max_lift_candidates":0,"max_affine_coordinates":0});
    let proof_only = discover(&server, request.clone());
    assert_ne!(proof_only["status"], "proved");
    assert!(proof_only.get("refuted").is_none());
    request["refute"] = json!(true);
    let refuted = discover(&server, request);
    assert_eq!(refuted["status"], "refuted");
    assert_eq!(refuted["refuted"], true);
    assert!(refuted["route"].is_null());
    assert!(refuted.get("real_rooted_from_index").is_none());
}

#[test]
fn catalog_paging_filters_fields_and_code_opt_in() {
    let server = PolynomialToolsServer::with_report_store(Default::default());
    let Json(verified) = server
        .list_oeis_sequences(Parameters(
            serde_json::from_value(json!({"status":"verified"})).unwrap(),
        ))
        .unwrap();
    assert!(verified.sequences.iter().all(|s| s["status"] == "verified"));
    let mut ids = vec![];
    let mut offset = 0;
    loop {
        let Json(list) = server
            .list_oeis_sequences(Parameters(
                serde_json::from_value(json!({"offset":offset,"limit":50})).unwrap(),
            ))
            .unwrap();
        assert_eq!(list.total, 919);
        assert_eq!(list.count, list.sequences.len());
        for item in list.sequences {
            assert_eq!(item.len(), 3);
            ids.push(item["id"].clone());
        }
        match list.next_offset {
            Some(next) => offset = next,
            None => break,
        }
    }
    assert_eq!(ids.len(), 919);
    assert!(ids.windows(2).all(|w| w[0].as_str() < w[1].as_str()));
    let Json(list) = server.list_oeis_sequences(Parameters(serde_json::from_value(json!({"query":"a001263","status":"validated","layout":"regular_triangle","fields":["id","layout"]})).unwrap())).unwrap();
    assert!(list.sequences.iter().any(|s| s["id"] == "A001263"));
    assert_eq!(list.sequences[0].len(), 2);
    for request in [
        json!({"limit":1001}),
        json!({"limit":0}),
        json!({"fields":["bad"]}),
        json!({"status":"bad"}),
        json!({"layout":"bad"}),
    ] {
        assert!(server
            .list_oeis_sequences(Parameters(serde_json::from_value(request).unwrap()))
            .is_err());
    }
    let Json(compact) = server
        .get_oeis_sequence(Parameters(
            serde_json::from_value(json!({"id":"A001263"})).unwrap(),
        ))
        .unwrap();
    let Json(full) = server
        .get_oeis_sequence(Parameters(
            serde_json::from_value(json!({"id":"A001263","include_code":true})).unwrap(),
        ))
        .unwrap();
    assert!(compact.python.is_none());
    assert!(full.python.is_some());
    let record: polytool::recurrence::RecurrenceJson =
        serde_json::from_str(&compact.recurrence_json).unwrap();
    assert_eq!(record, crate::theorems::catalog_record("A001263").unwrap());
}
