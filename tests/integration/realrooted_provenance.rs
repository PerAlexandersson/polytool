//! Standalone inventory guard for every RealRooted citation in library/CLI sources.

use std::{collections::BTreeSet, path::Path};

fn citations(path: &Path, found: &mut BTreeSet<String>) {
    for entry in std::fs::read_dir(path).unwrap() {
        let path = entry.unwrap().path();
        if path.is_dir() {
            citations(&path, found);
        } else if path.extension().is_some_and(|e| e == "rs") {
            let text = std::fs::read_to_string(path).unwrap();
            for tail in text.split("RealRooted.").skip(1) {
                let name: String = tail
                    .chars()
                    .take_while(|c| c.is_ascii_alphanumeric() || *c == '_' || *c == '.')
                    .collect();
                let name = name.trim_end_matches('.');
                if name.starts_with(|c: char| c.is_ascii_alphabetic() || c == '_') {
                    found.insert(format!("RealRooted.{name}"));
                }
            }
        }
    }
}

#[test]
fn all_realrooted_citations_match_reviewed_inventory() {
    let mut found = BTreeSet::new();
    citations(
        &Path::new(env!("CARGO_MANIFEST_DIR")).join("src"),
        &mut found,
    );
    let inventory = include_str!("../fixtures/realrooted-declarations.txt");
    let revision = inventory
        .lines()
        .next()
        .unwrap()
        .strip_prefix("# RealRooted commit ")
        .expect("pinned upstream revision");
    assert_eq!(revision.len(), 40);
    assert!(revision.bytes().all(|c| c.is_ascii_hexdigit()));
    let listed: BTreeSet<String> = inventory
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let (kind, name) = line.split_once('\t').expect("inventory kind and name");
            assert!(matches!(kind, "declaration" | "namespace" | "module"));
            name.to_string()
        })
        .collect();
    assert!(!listed.is_empty());
    assert_eq!(
        found, listed,
        "Review changed citations against an independent inventory; see scripts/realrooted-declarations.py --help"
    );
}
