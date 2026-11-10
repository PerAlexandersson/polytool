//! Output fingerprints of the proof discovery and of the parent provers on
//! every bundled OEIS catalog family.
//!
//! For each catalog entry with a runnable recurrence this records a stable
//! hash of the complete JSON of [`discover_recurrence_proofs`] (default
//! options) and of [`prove_parent_real_rootedness`] (first index, no
//! attestations, root-domain route enabled), or of the error text.  The
//! fixture `tests/fixtures/discovery-report-fingerprints.tsv` is the evidence
//! gate for refactorings of the theorem routes: a behaviour-preserving
//! change must leave every fingerprint unchanged.
//!
//! The quick test compares a fixed sample of 48 cheaper entries; the ignored
//! full test compares all of them, with a limit per entry.  Re-bless only for
//! an intentional change of results:
//!
//! ```text
//! POLYTOOL_BLESS_DISCOVERY_FINGERPRINTS=1 cargo test --release \
//!     --test discovery_fingerprints -- --ignored every_
//! ```
//!
//! (`POLYTOOL_DISCOVERY_ENTRY_SECONDS` overrides the per-entry limit of 120 s.)

use polytool::oeis::{by_id, catalog, OeisSequenceDefinition};
use polytool::recurrence::{
    discover_recurrence_proofs, prove_parent_real_rootedness, ParentRecurrenceRecord,
    RealRootednessTheoremOptions, RecurrenceJson,
};
use std::collections::BTreeMap;
use std::sync::mpsc;
use std::time::{Duration, Instant};

const FINGERPRINTS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/fixtures/discovery-report-fingerprints.tsv"
);

/// The quick sample: 48 entries spread over the catalog, chosen among the
/// cheaper ones (about 4 s in release mode; the full sweep takes about 10 min).
const SAMPLE: [&str; 48] = [
    "A001263", "A003982", "A013609", "A013613", "A023531", "A027467", "A029600", "A034870",
    "A038221", "A038255", "A038279", "A038301", "A038319", "A038323", "A038327", "A038331",
    "A038336", "A039991", "A052553", "A056241", "A084534", "A085880", "A093562", "A093644",
    "A096940", "A102537", "A103633", "A111595", "A120909", "A124848", "A127647", "A129184",
    "A130595", "A132200", "A133081", "A134309", "A136158", "A136394", "A147644", "A152828",
    "A159854", "A163936", "A172185", "A180047", "A193726", "A201972", "A249247", "A298668",
];

/// 64-bit FNV-1a of the text (stable across platforms and toolchains).
fn fingerprint(text: &str) -> String {
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    for byte in text.bytes() {
        hash ^= u64::from(byte);
        hash = hash.wrapping_mul(0x0100_0000_01b3);
    }
    format!("{hash:016x}:{}", text.len())
}

/// The two fingerprints (discovery report, parent proof) of one entry.
fn entry_fingerprints(entry: &OeisSequenceDefinition) -> Option<(String, String)> {
    let (rule, seeds) = entry.recurrence_parts().ok()?;
    let first_index = entry.recurrence_first_index;
    let record = RecurrenceJson::from_recurrence_rational(&rule, first_index, &seeds, None);
    let discovery =
        match discover_recurrence_proofs(&record, &RealRootednessTheoremOptions::default()) {
            Ok(report) => serde_json::to_string(&report).unwrap(),
            Err(error) => format!("error: {error}"),
        };
    let parent_record = ParentRecurrenceRecord {
        rule: &rule,
        seeds: &seeds,
        first_index,
    };
    let parent = match prove_parent_real_rootedness(&parent_record, first_index, &[], true) {
        Ok(proof) => serde_json::to_string(&proof).unwrap(),
        Err(error) => format!("error: {error}"),
    };
    Some((fingerprint(&discovery), fingerprint(&parent)))
}

/// `entry_fingerprints` on a worker thread, failing after the per-entry limit.
fn limited_fingerprints(
    entry: &'static OeisSequenceDefinition,
    limit: Duration,
) -> Option<(String, String)> {
    let (sender, receiver) = mpsc::channel();
    std::thread::Builder::new()
        .stack_size(256 << 20)
        .spawn(move || {
            let _ = sender.send(entry_fingerprints(entry));
        })
        .unwrap();
    receiver
        .recv_timeout(limit)
        .unwrap_or_else(|_| panic!("{}: exceeded the per-entry limit {limit:?}", entry.id))
}

fn entry_limit() -> Duration {
    let seconds = std::env::var("POLYTOOL_DISCOVERY_ENTRY_SECONDS")
        .ok()
        .and_then(|text| text.parse().ok())
        .unwrap_or(120);
    Duration::from_secs(seconds)
}

fn recorded() -> BTreeMap<String, String> {
    std::fs::read_to_string(FINGERPRINTS)
        .unwrap_or_else(|error| panic!("failed to read {FINGERPRINTS}: {error}"))
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            let (id, rest) = line.split_once('\t').unwrap();
            (id.to_string(), rest.to_string())
        })
        .collect()
}

fn check_entries(entries: Vec<&'static OeisSequenceDefinition>) {
    let expected = recorded();
    let limit = entry_limit();
    let mut checked = 0;
    let print_times = std::env::var_os("POLYTOOL_DISCOVERY_TIMES").is_some();
    for entry in entries {
        let started = Instant::now();
        let Some((discovery, parent)) = limited_fingerprints(entry, limit) else {
            continue;
        };
        if print_times {
            eprintln!("time {} {:.3}", entry.id, started.elapsed().as_secs_f64());
        }
        let line = expected
            .get(entry.id)
            .unwrap_or_else(|| panic!("{} has no recorded fingerprint", entry.id));
        assert_eq!(&format!("{discovery}\t{parent}"), line, "{}", entry.id);
        checked += 1;
    }
    assert!(checked > 0);
}

#[test]
fn sample_of_catalog_entries_matches_recorded_fingerprints() {
    check_entries(
        SAMPLE
            .iter()
            .map(|id| by_id(id).unwrap_or_else(|| panic!("{id} left the catalog")))
            .collect(),
    );
}

/// All catalog entries (about 10 minutes in release mode).
#[test]
#[ignore = "full catalog sweep; run with --release -- --ignored every_"]
fn every_catalog_entry_matches_recorded_fingerprints() {
    if std::env::var_os("POLYTOOL_BLESS_DISCOVERY_FINGERPRINTS").is_some() {
        let limit = entry_limit();
        let started = Instant::now();
        let mut text = String::new();
        for entry in catalog() {
            if let Some((discovery, parent)) = limited_fingerprints(entry, limit) {
                text.push_str(&format!("{}\t{discovery}\t{parent}\n", entry.id));
            }
        }
        std::fs::write(FINGERPRINTS, text).unwrap();
        eprintln!("blessed in {:?}", started.elapsed());
        return;
    }
    let expected = recorded();
    check_entries(catalog().iter().collect());
    assert!(expected.len() >= 860);
}
