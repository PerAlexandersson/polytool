//! Recurrence JSON fixtures: the manifest matches the files and the rows replay.

use polytool::recurrence::{
    format_rational_coeff, parse_rational_coeff, BigRational, RecurrenceJson,
};
use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};

const BASE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/fixtures/recurrence-benchmarks"
);

#[derive(Debug)]
struct FixtureManifestRow {
    slug: String,
    rows_file: String,
    json_file: String,
}

fn parse_manifest() -> Vec<FixtureManifestRow> {
    let manifest = fs::read_to_string(Path::new(BASE).join("manifest.tsv"))
        .expect("read recurrence benchmark manifest");
    manifest
        .lines()
        .skip(1)
        .map(|line| {
            let cols: Vec<&str> = line.split('\t').collect();
            assert_eq!(cols.len(), 7, "manifest row should have 7 columns: {line}");
            FixtureManifestRow {
                slug: cols[0].to_string(),
                rows_file: cols[5].to_string(),
                json_file: cols[6].to_string(),
            }
        })
        .collect()
}

fn parse_row(line: &str) -> Vec<BigRational> {
    line.split(',')
        .map(|coeff| parse_rational_coeff(coeff.trim()).expect("parse rational coefficient"))
        .collect()
}

fn read_rows(path: PathBuf) -> Vec<Vec<BigRational>> {
    fs::read_to_string(path)
        .expect("read coefficient rows")
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(parse_row)
        .collect()
}

fn format_rows(rows: &[Vec<BigRational>]) -> Vec<String> {
    let zero = parse_rational_coeff("0").expect("parse zero");
    rows.iter()
        .map(|row| {
            let end = row
                .iter()
                .rposition(|coeff| coeff != &zero)
                .map_or(1, |i| i + 1);
            row[..end]
                .iter()
                .map(format_rational_coeff)
                .collect::<Vec<_>>()
                .join(", ")
        })
        .collect()
}

fn fixture_files(dir: &str, extension: &str) -> BTreeSet<String> {
    fs::read_dir(Path::new(BASE).join(dir))
        .unwrap_or_else(|err| panic!("read fixture directory {dir}: {err}"))
        .map(|entry| entry.expect("read fixture directory entry").path())
        .filter(|path| path.extension().is_some_and(|ext| ext == extension))
        .map(|path| {
            format!(
                "{}/{}",
                dir,
                path.file_name()
                    .expect("fixture path has file name")
                    .to_string_lossy()
            )
        })
        .collect()
}

#[test]
fn recurrence_benchmark_manifest_matches_generated_files() {
    let fixtures = parse_manifest();
    let expected_rows = fixtures
        .iter()
        .map(|fixture| fixture.rows_file.clone())
        .collect::<BTreeSet<_>>();
    let expected_json = fixtures
        .iter()
        .map(|fixture| fixture.json_file.clone())
        .collect::<BTreeSet<_>>();

    assert_eq!(
        expected_rows,
        fixture_files("rows", "txt"),
        "manifest rows_file entries should match generated row files exactly"
    );
    assert_eq!(
        expected_json,
        fixture_files("json", "json"),
        "manifest json_file entries should match generated JSON files exactly"
    );
}

#[test]
fn recurrence_json_fixtures_regenerate_raw_rows() {
    let fixtures = parse_manifest();
    assert!(
        fixtures.len() >= 45,
        "expected the synthetic plus OEIS benchmark fixture suite"
    );

    for fixture in fixtures {
        let expected = read_rows(Path::new(BASE).join(&fixture.rows_file));
        let json_text =
            fs::read_to_string(Path::new(BASE).join(&fixture.json_file)).expect("read JSON file");
        let recurrence_json: RecurrenceJson =
            serde_json::from_str(&json_text).expect("parse recurrence JSON");
        let skip_prefix = recurrence_json
            .search
            .as_ref()
            .map(|search| search.skip_prefix)
            .unwrap_or(0);
        let expected_generated = expected
            .get(skip_prefix..)
            .unwrap_or_else(|| panic!("{} has invalid skip_prefix", fixture.slug));
        let (recurrence, first_index, initial_polys) = recurrence_json
            .to_recurrence_parts()
            .expect("convert recurrence JSON");
        let generated = recurrence
            .generate_rows_rational(&initial_polys, first_index, expected_generated.len())
            .unwrap_or_else(|err| panic!("generate rows for {}: {err}", fixture.slug));
        assert_eq!(
            format_rows(&generated),
            format_rows(expected_generated),
            "JSON recurrence should regenerate {}",
            fixture.slug
        );

        let extended = recurrence
            .generate_rows_rational(&initial_polys, first_index, expected_generated.len() + 5)
            .unwrap_or_else(|err| panic!("extend rows for {}: {err}", fixture.slug));
        assert_eq!(extended.len(), expected_generated.len() + 5);
        assert_eq!(
            format_rows(&extended[..expected_generated.len()]),
            format_rows(expected_generated),
            "extended JSON recurrence should preserve fixture prefix {}",
            fixture.slug
        );
    }
}
