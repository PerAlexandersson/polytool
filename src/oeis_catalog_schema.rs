// Schema of the bundled OEIS catalog data, `data/oeis-catalog.jsonl`.
//
// This file is shared verbatim: `src/oeis.rs` includes it into the public
// `oeis` module, and `build.rs` includes it to turn the data file into the
// static catalog.  It therefore has no inner attributes or `use` items; the
// including scope provides `serde::{Deserialize, Serialize}`.

/// Publication status of one bundled recurrence.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum OeisSequenceStatus {
    /// The recurrence was checked on rows excluded from fitting.
    Verified,
    /// The recurrence matches a current OEIS prefix, but fitting provenance is unavailable.
    Validated,
    /// The recurrence generates known rows but still lacks a genuine holdout.
    Experimental,
}

impl OeisSequenceStatus {
    /// Stable lowercase name used by the CLI and MCP server.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Verified => "verified",
            Self::Validated => "validated",
            Self::Experimental => "experimental",
        }
    }
}

/// Shape of the OEIS data represented by polynomial coefficient rows.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum OeisLayout {
    /// OEIS `tabl`: successive rows have lengths `1, 2, 3, ...`.
    RegularTriangle,
    /// OEIS `tabf` or another explicitly row-grouped table.
    Table,
}

impl OeisLayout {
    /// Stable lowercase name used by structured output.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::RegularTriangle => "regular_triangle",
            Self::Table => "table",
        }
    }
}

/// Where the recurrence of a catalog record came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum OeisProvenanceSource {
    /// Recurrence from `polytool recurrence --json` on the bundled
    /// `fixtures/recurrence-benchmarks` rows, checked on held-out rows.
    BenchmarkFixture,
    /// Recurrence originally fitted by an earlier version of the polytool
    /// recurrence search and transcribed into generated Lean definitions;
    /// validated against an OEIS data prefix (16 generated rows).
    LeanDefinition,
    /// Recurrence fitted by recurrence search on cached rows of a curated
    /// queue; verified when it was checked on held-out rows.
    CuratedQueue,
    /// Recurrence found by `polytool oeis add`: the adaptive recurrence search
    /// on a prefix of the OEIS rows, then an exact check of every remaining
    /// available term.
    PolytoolRecurrenceSearch,
}

impl OeisProvenanceSource {
    /// Stable lowercase name used in the data file and structured output.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::BenchmarkFixture => "benchmark_fixture",
            Self::LeanDefinition => "lean_definition",
            Self::CuratedQueue => "curated_queue",
            Self::PolytoolRecurrenceSearch => "polytool_recurrence_search",
        }
    }
}

/// One complete OEIS row recorded as validation evidence.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OeisValidationRow {
    /// Zero-based row index in catalog order (prefix rows included).
    pub index: usize,
    /// The expected integer coefficients, as decimal strings.
    pub terms: Vec<String>,
}

/// Provenance and validation evidence of one catalog record.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OeisCatalogProvenance {
    /// How the recurrence was obtained.
    pub source: OeisProvenanceSource,
    /// Publication status; becomes [`OeisSequenceDefinition::status`].
    pub status: OeisSequenceStatus,
    /// Number of OEIS-derived rows used when the entry was generated
    /// (fitted plus held out); becomes [`OeisSequenceDefinition::source_rows`].
    pub source_rows: usize,
    /// Number of those rows held out from fitting and checked exactly;
    /// becomes [`OeisSequenceDefinition::verification_rows`].
    pub verification_rows: usize,
    /// SHA-256 of the source row data at generation time (lowercase hex);
    /// becomes [`OeisSequenceDefinition::rows_sha256`].
    pub rows_sha256: String,
    /// Recurrence benchmark slug (only for
    /// [`OeisProvenanceSource::BenchmarkFixture`]).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fixture_slug: Option<String>,
    /// Number of catalog rows preceding the first fixture row; when absent it
    /// equals the number of prefix rows.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fixture_row_offset: Option<usize>,
    /// The last complete OEIS row that was compared, kept as a regression
    /// check.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub validation_row: Option<OeisValidationRow>,
    /// Number of leading rows given to the recurrence search.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fit_rows: Option<usize>,
    /// Number of flattened OEIS terms compared exactly with generated rows.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub checked_terms: Option<usize>,
    /// Kind of OEIS input the terms came from (`oeis_internal_format` for the
    /// `%S/%T/%U` data lines, `bfile` for a b-file).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub terms_source: Option<String>,
    /// Adaptive search options and candidate budget used by the search.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub search_options: Option<serde_json::Value>,
    /// Version of polytool that performed the search.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub polytool_version: Option<String>,
}

/// One recurrence term `sign(n) c(n,t) D^d P_{n-r}(t)` of a catalog record.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OeisCatalogTerm {
    /// Recurrence offset `r`.
    pub offset: usize,
    /// Order `d` of the `t`-derivative.
    pub derivative_order: usize,
    /// Whether the term is multiplied by `(-1)^n`.
    pub alternating_sign: bool,
    /// Nonzero monomials `[n_degree, t_degree, "value"]` of `c(n,t)`.
    pub coefficient: Vec<(usize, usize, String)>,
}

/// Recurrence `f(n,t) P_n = sum of terms + g(n,t)` of a catalog record.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OeisCatalogRecurrence {
    /// Right-hand-side terms.
    pub terms: Vec<OeisCatalogTerm>,
    /// Monomials of `f(n,t)`; absent means 1.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub denominator: Option<Vec<(usize, usize, String)>>,
    /// Monomials of `g(n,t)`; absent means homogeneous.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub inhomogeneous: Option<Vec<(usize, usize, String)>>,
}

/// One line of `data/oeis-catalog.jsonl`: an executable recurrence description
/// (the fields of [`OeisSequenceDefinition`]) plus its provenance.
///
/// The data file holds one record per line, sorted by A-number, each written
/// by [`OeisCatalogRecord::to_jsonl_line`] (so key order is fixed).  The build
/// script turns it into the static [`catalog`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OeisCatalogRecord {
    /// OEIS A-number, upper case.
    pub id: String,
    /// OEIS title of the sequence.
    pub name: String,
    /// Row-shape convention of the OEIS data.
    pub layout: OeisLayout,
    /// Displayed index of the first row.
    pub first_row: i64,
    /// Flattened OEIS index of the first coefficient.
    pub flattened_offset: i64,
    /// Whether the flattened prefix agrees with OEIS data.
    pub bfile_prefix_verified: bool,
    /// Complete rows before the recurrence's first row.
    pub prefix_rows: Vec<Vec<String>>,
    /// Recurrence index `n` of `initial_rows[0]`.
    pub recurrence_first_index: usize,
    /// Structural width of the first row after the prefix.
    pub recurrence_first_width: usize,
    /// Initial rows of the recurrence (exact rational strings).
    pub initial_rows: Vec<Vec<String>>,
    /// The recurrence.
    pub recurrence: OeisCatalogRecurrence,
    /// Where the recurrence came from and how it was checked.
    pub provenance: OeisCatalogProvenance,
}

/// Why a catalog data file or one catalog record was rejected.
///
/// Every variant that concerns one record carries its A-number in `id`.  The
/// [`Display`](std::fmt::Display) text of each variant is the message of the
/// former `String` error.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OeisCatalogError {
    /// A data-file line is not a valid JSON catalog record.
    InvalidJsonLine {
        /// One-based line number in the data file.
        line: usize,
        /// The deserialization error.
        message: String,
    },
    /// The record id is not an upper-case seven-character A-number.
    InvalidId {
        /// The offending id.
        id: String,
    },
    /// The record name is empty.
    EmptyName {
        /// The record id.
        id: String,
    },
    /// A prefix or initial row entry is not an exact rational string.
    InvalidRowCoefficient {
        /// The record id.
        id: String,
        /// The offending entry.
        value: String,
    },
    /// The record has no initial rows.
    NoInitialRows {
        /// The record id.
        id: String,
    },
    /// A recurrence coefficient is not an exact rational string.
    InvalidRecurrenceCoefficient {
        /// The record id.
        id: String,
        /// The offending entry.
        value: String,
    },
    /// A recurrence polynomial lists the same monomial twice.
    DuplicateMonomial {
        /// The record id.
        id: String,
        /// Power of `n` of the repeated monomial.
        n_degree: usize,
        /// Power of `t` of the repeated monomial.
        t_degree: usize,
    },
    /// A recurrence term has offset 0.
    ZeroOffsetTerm {
        /// The record id.
        id: String,
    },
    /// `rows_sha256` is not 64 lowercase hexadecimal digits.
    InvalidRowsSha256 {
        /// The record id.
        id: String,
    },
    /// `verification_rows` exceeds `source_rows`.
    MoreVerificationRowsThanSourceRows {
        /// The record id.
        id: String,
    },
    /// The status is `verified` but no rows were held out.
    VerifiedWithoutHeldOutRows {
        /// The record id.
        id: String,
    },
    /// Rows were held out but the status is not `verified`.
    HeldOutRowsWithoutVerifiedStatus {
        /// The record id.
        id: String,
    },
    /// `fixture_slug` is present unless the source is `benchmark_fixture`, or
    /// missing for it.
    FixtureSlugMismatch {
        /// The record id.
        id: String,
    },
    /// `fixture_slug` is the empty string.
    EmptyFixtureSlug {
        /// The record id.
        id: String,
    },
    /// A `polytool_recurrence_search` record lacks `fit_rows`, `checked_terms`
    /// or `search_options`.
    MissingSearchProvenance {
        /// The record id.
        id: String,
    },
    /// Two consecutive records are not in strictly increasing A-number order.
    UnsortedRecords {
        /// The earlier record's id.
        first: String,
        /// The following record's id.
        second: String,
    },
}

impl std::fmt::Display for OeisCatalogError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidJsonLine { line, message } => write!(f, "line {line}: {message}"),
            Self::InvalidId { id } => write!(f, "`{id}` is not an upper-case A-number"),
            Self::EmptyName { id } => write!(f, "{id}: empty name"),
            Self::InvalidRowCoefficient { id, value } => {
                write!(f, "{id}: invalid row coefficient `{value}`")
            }
            Self::NoInitialRows { id } => write!(f, "{id}: no initial rows"),
            Self::InvalidRecurrenceCoefficient { id, value } => {
                write!(f, "{id}: invalid recurrence coefficient `{value}`")
            }
            Self::DuplicateMonomial {
                id,
                n_degree,
                t_degree,
            } => write!(f, "{id}: duplicate monomial n^{n_degree} t^{t_degree}"),
            Self::ZeroOffsetTerm { id } => write!(f, "{id}: recurrence term with offset 0"),
            Self::InvalidRowsSha256 { id } => {
                write!(f, "{id}: rows_sha256 is not a lowercase SHA-256")
            }
            Self::MoreVerificationRowsThanSourceRows { id } => {
                write!(f, "{id}: more verification rows than source rows")
            }
            Self::VerifiedWithoutHeldOutRows { id } => {
                write!(f, "{id}: verified without held-out rows")
            }
            Self::HeldOutRowsWithoutVerifiedStatus { id } => {
                write!(f, "{id}: held-out rows recorded but status is not verified")
            }
            Self::FixtureSlugMismatch { id } => write!(
                f,
                "{id}: fixture_slug must be present exactly for benchmark_fixture"
            ),
            Self::EmptyFixtureSlug { id } => write!(f, "{id}: empty fixture_slug"),
            Self::MissingSearchProvenance { id } => write!(
                f,
                "{id}: polytool_recurrence_search needs fit_rows, checked_terms and search_options"
            ),
            Self::UnsortedRecords { first, second } => write!(
                f,
                "records are not strictly sorted by A-number: {first} before {second}"
            ),
        }
    }
}

impl std::error::Error for OeisCatalogError {}

impl OeisCatalogRecord {
    /// Serialize as one canonical data-file line (no trailing newline).
    pub fn to_jsonl_line(&self) -> String {
        serde_json::to_string(self).expect("catalog records serialize")
    }

    /// Number of catalog rows preceding the first benchmark-fixture row.
    pub fn fixture_row_offset(&self) -> usize {
        self.provenance
            .fixture_row_offset
            .unwrap_or(self.prefix_rows.len())
    }

    /// Parse a whole data file, then check it with [`Self::validate_all`].
    /// Blank lines are ignored.
    pub fn parse_jsonl(text: &str) -> Result<Vec<Self>, OeisCatalogError> {
        let mut records = Vec::new();
        for (line_index, line) in text.lines().enumerate() {
            if line.trim().is_empty() {
                continue;
            }
            let record: Self = serde_json::from_str(line)
                .map_err(|error| OeisCatalogError::InvalidJsonLine {
                    line: line_index + 1,
                    message: error.to_string(),
                })?;
            records.push(record);
        }
        Self::validate_all(&records)?;
        Ok(records)
    }

    /// Check one record's internal consistency: A-number syntax, exact
    /// rational coefficient strings, and provenance fields matching the
    /// status and source.
    pub fn validate(&self) -> Result<(), OeisCatalogError> {
        let id = &self.id;
        let error_id = || id.clone();
        let bytes = id.as_bytes();
        if bytes.len() != 7 || bytes[0] != b'A' || !bytes[1..].iter().all(u8::is_ascii_digit) {
            return Err(OeisCatalogError::InvalidId { id: error_id() });
        }
        if self.name.trim().is_empty() {
            return Err(OeisCatalogError::EmptyName { id: error_id() });
        }
        let is_rational = |value: &str| {
            let value = value.strip_prefix('-').unwrap_or(value);
            let mut parts = value.splitn(2, '/');
            let numerator = parts.next().unwrap_or("");
            let digits = |part: &str| !part.is_empty() && part.bytes().all(|b| b.is_ascii_digit());
            digits(numerator)
                && parts
                    .next()
                    .is_none_or(|denominator| {
                        digits(denominator) && denominator.bytes().any(|b| b != b'0')
                    })
        };
        let rows = self.prefix_rows.iter().chain(&self.initial_rows);
        for value in rows.flatten() {
            if !is_rational(value) {
                return Err(OeisCatalogError::InvalidRowCoefficient {
                    id: error_id(),
                    value: value.clone(),
                });
            }
        }
        if self.initial_rows.is_empty() {
            return Err(OeisCatalogError::NoInitialRows { id: error_id() });
        }
        let recurrence = &self.recurrence;
        let polynomials = recurrence
            .terms
            .iter()
            .map(|term| &term.coefficient)
            .chain(recurrence.denominator.iter())
            .chain(recurrence.inhomogeneous.iter());
        for polynomial in polynomials {
            let mut seen = Vec::new();
            for (n_degree, t_degree, value) in polynomial {
                if !is_rational(value) {
                    return Err(OeisCatalogError::InvalidRecurrenceCoefficient {
                        id: error_id(),
                        value: value.clone(),
                    });
                }
                if seen.contains(&(n_degree, t_degree)) {
                    return Err(OeisCatalogError::DuplicateMonomial {
                        id: error_id(),
                        n_degree: *n_degree,
                        t_degree: *t_degree,
                    });
                }
                seen.push((n_degree, t_degree));
            }
        }
        if recurrence.terms.iter().any(|term| term.offset == 0) {
            return Err(OeisCatalogError::ZeroOffsetTerm { id: error_id() });
        }
        let provenance = &self.provenance;
        let sha = &provenance.rows_sha256;
        if sha.len() != 64 || !sha.bytes().all(|b| matches!(b, b'0'..=b'9' | b'a'..=b'f')) {
            return Err(OeisCatalogError::InvalidRowsSha256 { id: error_id() });
        }
        if provenance.verification_rows > provenance.source_rows {
            return Err(OeisCatalogError::MoreVerificationRowsThanSourceRows { id: error_id() });
        }
        if provenance.status == OeisSequenceStatus::Verified && provenance.verification_rows == 0 {
            return Err(OeisCatalogError::VerifiedWithoutHeldOutRows { id: error_id() });
        }
        if provenance.status != OeisSequenceStatus::Verified && provenance.verification_rows != 0 {
            return Err(OeisCatalogError::HeldOutRowsWithoutVerifiedStatus { id: error_id() });
        }
        let is_fixture = provenance.source == OeisProvenanceSource::BenchmarkFixture;
        if is_fixture != provenance.fixture_slug.is_some() {
            return Err(OeisCatalogError::FixtureSlugMismatch { id: error_id() });
        }
        if provenance.fixture_slug.as_deref() == Some("") {
            return Err(OeisCatalogError::EmptyFixtureSlug { id: error_id() });
        }
        if provenance.source == OeisProvenanceSource::PolytoolRecurrenceSearch
            && (provenance.fit_rows.is_none()
                || provenance.checked_terms.is_none()
                || provenance.search_options.is_none())
        {
            return Err(OeisCatalogError::MissingSearchProvenance { id: error_id() });
        }
        Ok(())
    }

    /// Validate every record and check that ids are strictly increasing.
    pub fn validate_all(records: &[Self]) -> Result<(), OeisCatalogError> {
        for record in records {
            record.validate()?;
        }
        for pair in records.windows(2) {
            if pair[0].id >= pair[1].id {
                return Err(OeisCatalogError::UnsortedRecords {
                    first: pair[0].id.clone(),
                    second: pair[1].id.clone(),
                });
            }
        }
        Ok(())
    }
}
