//! Recurrence overfit-detection fixtures.

use polytool::recurrence::parse_rational_coeff;
use std::fs;
use std::path::{Path, PathBuf};

const BASE: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/fixtures/recurrence-overfits");

#[derive(Debug)]
struct OverfitFixture {
    id: String,
    rows_file: String,
    recurrence_file: String,
    rows: usize,
}

fn parse_manifest() -> Vec<OverfitFixture> {
    let manifest = fs::read_to_string(Path::new(BASE).join("manifest.tsv")).expect("read manifest");
    manifest
        .lines()
        .skip(1)
        .map(|line| {
            let cols = line.split('\t').collect::<Vec<_>>();
            assert_eq!(cols.len(), 5, "manifest row should have 5 columns: {line}");
            OverfitFixture {
                id: cols[0].to_string(),
                rows_file: cols[1].to_string(),
                recurrence_file: cols[2].to_string(),
                rows: cols[3].parse().expect("row count"),
            }
        })
        .collect()
}

fn fixture_path(relative: &str) -> PathBuf {
    Path::new(BASE).join(relative)
}

#[test]
fn overfit_fixture_manifest_matches_files() {
    for fixture in parse_manifest() {
        let rows_path = fixture_path(&fixture.rows_file);
        let recurrence_path = fixture_path(&fixture.recurrence_file);
        assert!(rows_path.exists(), "{} rows file exists", fixture.id);
        assert!(
            recurrence_path.exists(),
            "{} recurrence file exists",
            fixture.id
        );

        let rows = fs::read_to_string(rows_path).expect("read rows");
        let row_count = rows.lines().filter(|line| !line.trim().is_empty()).count();
        assert_eq!(row_count, fixture.rows, "{} row count", fixture.id);

        for line in rows.lines().filter(|line| !line.trim().is_empty()) {
            for coeff in line.split(',') {
                parse_rational_coeff(coeff.trim())
                    .unwrap_or_else(|err| panic!("{} row coefficient: {err}", fixture.id));
            }
        }
    }
}

#[test]
fn overfit_recurrences_are_intentionally_suspicious() {
    for fixture in parse_manifest() {
        let recurrence =
            fs::read_to_string(fixture_path(&fixture.recurrence_file)).expect("read recurrence");
        assert!(
            recurrence.contains("P(n)") && recurrence.contains("P(n-"),
            "{} stores a recurrence-shaped expression",
            fixture.id
        );

        let longest_number = recurrence
            .split(|ch: char| !(ch.is_ascii_digit()))
            .map(str::len)
            .max()
            .unwrap_or(0);
        let slash_count = recurrence.matches('/').count();

        match fixture.id.as_str() {
            "A181000" => {
                assert!(
                    longest_number > 100,
                    "A181000 should preserve the extreme coefficient-size regression"
                );
                assert!(
                    slash_count > 20,
                    "A181000 should preserve the rational-overfit regression"
                );
            }
            "A177970" => {
                assert!(
                    longest_number > 10,
                    "A177970 should preserve the high-rational-coefficient regression"
                );
                assert!(
                    slash_count > 10,
                    "A177970 should preserve the rational-overfit regression"
                );
            }
            other => panic!("unexpected overfit fixture {other}"),
        }
    }
}
