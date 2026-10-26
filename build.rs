//! Build script.
//!
//! - Embeds the abbreviated Git commit (or `POLYTOOL_GIT_COMMIT`) as the
//!   `POLYTOOL_GIT_COMMIT` environment variable for `--version` output.
//! - Generates the static OEIS catalog (`OEIS_CATALOG`, one accessor function
//!   per A-number and the test-only validation rows) from the tracked data
//!   file `data/oeis-catalog.jsonl` into `$OUT_DIR/oeis_catalog_generated.rs`,
//!   which `src/oeis.rs` includes.  The record schema is shared with the
//!   library through `src/oeis_catalog_schema.rs`; an invalid data file fails
//!   the build with the offending line or A-number.

use std::env;
use std::fmt::Write as _;
use std::path::Path;
use std::process::Command;

#[allow(dead_code)]
mod oeis_catalog {
    use serde::{Deserialize, Serialize};
    include!("src/oeis_catalog_schema.rs");
}

use oeis_catalog::{OeisCatalogRecord, OeisLayout, OeisSequenceStatus};

const OEIS_CATALOG_DATA: &str = "data/oeis-catalog.jsonl";
const OEIS_CATALOG_SCHEMA: &str = "src/oeis_catalog_schema.rs";

const ABBREVIATED_COMMIT_LENGTH: usize = 12;

fn normalized_commit(value: &str) -> Option<String> {
    let value = value.trim();
    if value.len() < ABBREVIATED_COMMIT_LENGTH
        || !value.bytes().all(|byte| byte.is_ascii_hexdigit())
    {
        return None;
    }
    Some(value[..ABBREVIATED_COMMIT_LENGTH].to_ascii_lowercase())
}

fn git_output(manifest_dir: &Path, args: &[&str]) -> Option<String> {
    let output = Command::new("git")
        .arg("-C")
        .arg(manifest_dir)
        .args(args)
        .output()
        .ok()?;
    if !output.status.success() {
        return None;
    }
    String::from_utf8(output.stdout).ok()
}

fn git_commit(manifest_dir: &Path) -> Option<String> {
    git_output(manifest_dir, &["ls-files", "--error-unmatch", "Cargo.toml"])?;
    normalized_commit(&git_output(
        manifest_dir,
        &["rev-parse", "--verify", "HEAD"],
    )?)
}

fn track_git_ref(manifest_dir: &Path) {
    let Some(head_path) = git_output(manifest_dir, &["rev-parse", "--git-path", "HEAD"]) else {
        return;
    };
    println!("cargo:rerun-if-changed={}", head_path.trim());

    if let Some(reference) =
        git_output(manifest_dir, &["symbolic-ref", "-q", "HEAD"]).and_then(|reference| {
            git_output(manifest_dir, &["rev-parse", "--git-path", reference.trim()])
        })
    {
        println!("cargo:rerun-if-changed={}", reference.trim());
    }

    if let Some(packed_refs) = git_output(manifest_dir, &["rev-parse", "--git-path", "packed-refs"])
    {
        println!("cargo:rerun-if-changed={}", packed_refs.trim());
    }
}

/// Rust string literal for `value` (`Debug` escaping is valid literal syntax).
fn rust_string(value: &str) -> String {
    format!("{value:?}")
}

/// Rustdoc line for a catalog accessor, with Markdown metacharacters escaped.
fn accessor_doc(id: &str, name: &str) -> String {
    let mut text = String::new();
    for (index, word) in name.split_whitespace().enumerate() {
        if index > 0 {
            text.push(' ');
        }
        for character in word.chars() {
            if matches!(character, '\\' | '`' | '*' | '_' | '[' | ']' | '<' | '>') {
                text.push('\\');
            }
            text.push(character);
        }
    }
    format!("/// {id}: {text}")
}

fn render_rows(rows: &[Vec<String>], indent: &str) -> String {
    if rows.is_empty() {
        return "&[]".to_string();
    }
    let rendered = rows
        .iter()
        .map(|row| {
            let values = row
                .iter()
                .map(|value| rust_string(value))
                .collect::<Vec<_>>()
                .join(", ");
            format!("{indent}    &[{values}]")
        })
        .collect::<Vec<_>>()
        .join(",\n");
    format!("&[\n{rendered},\n{indent}]")
}

fn render_monomials(terms: &[(usize, usize, String)], indent: &str) -> String {
    if terms.is_empty() {
        return "&[]".to_string();
    }
    let rendered = terms
        .iter()
        .map(|(n_degree, t_degree, value)| {
            format!(
                "{indent}    SparseCoefficient::new({n_degree}, {t_degree}, {})",
                rust_string(value)
            )
        })
        .collect::<Vec<_>>()
        .join(",\n");
    format!("&[\n{rendered},\n{indent}]")
}

fn render_optional_monomials(terms: Option<&[(usize, usize, String)]>, indent: &str) -> String {
    match terms {
        None => "None".to_string(),
        Some(terms) => format!("Some({})", render_monomials(terms, indent)),
    }
}

fn layout_variant(layout: OeisLayout) -> &'static str {
    match layout {
        OeisLayout::RegularTriangle => "RegularTriangle",
        OeisLayout::Table => "Table",
    }
}

fn status_variant(status: OeisSequenceStatus) -> &'static str {
    match status {
        OeisSequenceStatus::Verified => "Verified",
        OeisSequenceStatus::Validated => "Validated",
        OeisSequenceStatus::Experimental => "Experimental",
    }
}

/// Render one `OeisSequenceDefinition` literal.
fn render_definition(out: &mut String, record: &OeisCatalogRecord) {
    let provenance = &record.provenance;
    let recurrence = &record.recurrence;
    let fixture_slug = provenance.fixture_slug.as_deref().unwrap_or("");
    let _ = writeln!(out, "    OeisSequenceDefinition {{");
    let _ = writeln!(out, "        id: {},", rust_string(&record.id));
    let _ = writeln!(out, "        name: {},", rust_string(&record.name));
    let _ = writeln!(
        out,
        "        status: OeisSequenceStatus::{},",
        status_variant(provenance.status)
    );
    let _ = writeln!(
        out,
        "        layout: OeisLayout::{},",
        layout_variant(record.layout)
    );
    let _ = writeln!(out, "        first_row: {},", record.first_row);
    let _ = writeln!(
        out,
        "        flattened_offset: {},",
        record.flattened_offset
    );
    let _ = writeln!(
        out,
        "        bfile_prefix_verified: {},",
        record.bfile_prefix_verified
    );
    let _ = writeln!(out, "        fixture_slug: {},", rust_string(fixture_slug));
    let _ = writeln!(
        out,
        "        fixture_row_offset: {},",
        record.fixture_row_offset()
    );
    let _ = writeln!(
        out,
        "        prefix_rows: {},",
        render_rows(&record.prefix_rows, "        ")
    );
    let _ = writeln!(
        out,
        "        recurrence_first_index: {},",
        record.recurrence_first_index
    );
    let _ = writeln!(
        out,
        "        recurrence_first_width: {},",
        record.recurrence_first_width
    );
    let _ = writeln!(
        out,
        "        initial_rows: {},",
        render_rows(&record.initial_rows, "        ")
    );
    let _ = writeln!(out, "        recurrence: SparseRecurrence {{");
    let _ = writeln!(out, "            terms: &[");
    for term in &recurrence.terms {
        let indent = "                    ";
        let _ = writeln!(out, "                SparseRecurrenceTerm {{");
        let _ = writeln!(out, "{indent}offset: {},", term.offset);
        let _ = writeln!(out, "{indent}derivative_order: {},", term.derivative_order);
        let _ = writeln!(out, "{indent}alternating_sign: {},", term.alternating_sign);
        let _ = writeln!(
            out,
            "{indent}coefficient: {},",
            render_monomials(&term.coefficient, indent)
        );
        let _ = writeln!(out, "                }},");
    }
    let _ = writeln!(out, "            ],");
    let _ = writeln!(
        out,
        "            denominator: {},",
        render_optional_monomials(recurrence.denominator.as_deref(), "            ")
    );
    let _ = writeln!(
        out,
        "            inhomogeneous: {},",
        render_optional_monomials(recurrence.inhomogeneous.as_deref(), "            ")
    );
    let _ = writeln!(out, "        }},");
    let _ = writeln!(out, "        source_rows: {},", provenance.source_rows);
    let _ = writeln!(
        out,
        "        verification_rows: {},",
        provenance.verification_rows
    );
    let _ = writeln!(
        out,
        "        rows_sha256: {},",
        rust_string(&provenance.rows_sha256)
    );
    let _ = writeln!(out, "    }},");
}

/// Render the catalog source from parsed (and validated) records.
fn render_oeis_catalog(records: &[OeisCatalogRecord]) -> String {
    let mut out = String::new();
    out.push_str(
        "// @generated by build.rs from data/oeis-catalog.jsonl; do not edit by hand.\n\n",
    );
    out.push_str("/// Every bundled OEIS sequence definition, in catalog order.\n");
    out.push_str("pub static OEIS_CATALOG: &[OeisSequenceDefinition] = &[\n");
    for record in records {
        render_definition(&mut out, record);
    }
    out.push_str("];\n\n");
    out.push_str("#[cfg(test)]\n");
    out.push_str(
        "/// Test-only `(id, row index, expected row)` triples imported from OEIS b-files.\n",
    );
    out.push_str("pub static OEIS_IMPORTED_VALIDATION_ROWS: &[(&str, usize, &[&str])] = &[\n");
    for record in records {
        if let Some(row) = &record.provenance.validation_row {
            let values = row
                .terms
                .iter()
                .map(|value| rust_string(value))
                .collect::<Vec<_>>()
                .join(", ");
            let _ = writeln!(
                out,
                "    ({}, {}, &[{values}]),",
                rust_string(&record.id),
                row.index
            );
        }
    }
    out.push_str("];\n\n");
    for (index, record) in records.iter().enumerate() {
        let _ = writeln!(out, "{}", accessor_doc(&record.id, &record.name));
        let _ = writeln!(out, "#[allow(non_snake_case)]");
        let _ = writeln!(
            out,
            "pub fn {}() -> &'static OeisSequenceDefinition {{",
            record.id
        );
        let _ = writeln!(out, "    &OEIS_CATALOG[{index}]");
        let _ = writeln!(out, "}}\n");
    }
    out
}

fn generate_oeis_catalog(manifest_dir: &Path) {
    println!("cargo:rerun-if-changed={OEIS_CATALOG_DATA}");
    println!("cargo:rerun-if-changed={OEIS_CATALOG_SCHEMA}");
    let data_path = manifest_dir.join(OEIS_CATALOG_DATA);
    let text = std::fs::read_to_string(&data_path)
        .unwrap_or_else(|error| panic!("failed to read {}: {error}", data_path.display()));
    let records = OeisCatalogRecord::parse_jsonl(&text)
        .unwrap_or_else(|error| panic!("invalid {OEIS_CATALOG_DATA}: {error}"));
    let out_dir = env::var_os("OUT_DIR")
        .map(std::path::PathBuf::from)
        .expect("Cargo sets OUT_DIR for build scripts");
    let out_path = out_dir.join("oeis_catalog_generated.rs");
    std::fs::write(&out_path, render_oeis_catalog(&records))
        .unwrap_or_else(|error| panic!("failed to write {}: {error}", out_path.display()));
}

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    println!("cargo:rerun-if-env-changed=POLYTOOL_GIT_COMMIT");
    let manifest_dir = env::var_os("CARGO_MANIFEST_DIR")
        .map(std::path::PathBuf::from)
        .expect("Cargo sets CARGO_MANIFEST_DIR for build scripts");
    generate_oeis_catalog(&manifest_dir);

    let commit = match env::var("POLYTOOL_GIT_COMMIT") {
        Ok(value) => normalized_commit(&value),
        Err(_) => {
            track_git_ref(&manifest_dir);
            git_commit(&manifest_dir)
        }
    };
    println!(
        "cargo:rustc-env=POLYTOOL_GIT_COMMIT={}",
        commit.as_deref().unwrap_or("")
    );
}
