//! The stable `polytool.recurrence.v1` JSON schema and its conversions.

use super::*;
use serde::{Deserialize, Serialize};

/// Current JSON schema tag for serialized recurrences.
pub const RECURRENCE_JSON_SCHEMA: &str = "polytool.recurrence.v1";

/// Stable JSON representation of a recurrence together with initial rows.
#[derive(Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceJson {
    /// Schema identifier. Currently `polytool.recurrence.v1`.
    pub schema: String,
    /// Index of the first polynomial in `initial_polynomials`.
    pub first_index: usize,
    /// Initial coefficient rows, in ascending powers of `t`.
    pub initial_polynomials: Vec<Vec<String>>,
    /// The recurrence data.
    pub recurrence: RecurrenceJsonData,
    /// Optional metadata from an adaptive search.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub search: Option<RecurrenceJsonSearch>,
    /// Optional presentation with seeds and a starting index, one `$...$`
    /// formula per line. Ignored when evaluating the recurrence. Use
    /// [`Self::with_latex`] to refresh it after editing a record.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub latex: Option<String>,
}

impl std::fmt::Debug for RecurrenceJson {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut record = f.debug_struct("RecurrenceJson");
        record
            .field("schema", &self.schema)
            .field("first_index", &self.first_index)
            .field("initial_polynomials", &self.initial_polynomials)
            .field("recurrence", &self.recurrence)
            .field("search", &self.search);
        // Preserve existing diagnostic fingerprints when presentation is absent.
        if self.latex.is_some() {
            record.field("latex", &self.latex);
        }
        record.finish()
    }
}

/// Stable JSON representation of a recurrence without initial rows.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceJsonData {
    /// Right-hand-side terms `sign(n) c(n,t) D^d P_{n-r}(t)`.
    pub terms: Vec<RecurrenceTermJson>,
    /// Left-hand factor `f(n,t)`; omitted when it is the constant 1.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub denominator: Option<BivarPolyJson>,
    /// Additive term `g(n,t)`; omitted for a homogeneous recurrence.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub inhomogeneous: Option<BivarPolyJson>,
}

/// Stable JSON representation of one recurrence term.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceTermJson {
    /// Lag `r`: the term refers to `P_{n-r}`.
    pub offset: usize,
    /// Derivative order `d` in `t` applied to `P_{n-r}`.
    pub deriv_order: usize,
    /// Extra sign factor: `none`, or `alternating_n` for a multiplier `(-1)^n`.
    pub sign: String,
    /// Coefficient polynomial `c(n,t)`.
    pub coeff: BivarPolyJson,
}

/// Stable JSON representation of a bivariate polynomial in `(n,t)`.
///
/// `coeffs[i][j]` is the coefficient of `n^i t^j`, stored as an exact rational
/// string.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct BivarPolyJson {
    /// Coefficient array: `coeffs[i][j]` is the coefficient of `n^i t^j`, an integer or
    /// `num/den` string. Rows are padded to equal length on output and may be ragged on
    /// input.
    pub coeffs: Vec<Vec<String>>,
}

/// Optional metadata describing how an adaptive recurrence was found.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceJsonSearch {
    /// The recurrence as text, `f(n,t) P(n) = ...`.
    pub recurrence_text: String,
    /// Number of input rows supplied to the search, before `skip_prefix`.
    pub source_rows: usize,
    /// Number of leading input rows that were ignored.
    pub skip_prefix: usize,
    /// Number of unknown coefficients of the winning linear system.
    pub unknowns: usize,
    /// Weighted unknown count used to order the candidates.
    pub weighted_unknowns: usize,
    /// Number of scalar equations of the winning linear system.
    pub equations: usize,
    /// Number of rows, after `skip_prefix`, used to fit the winning system.
    pub fit_polynomials: usize,
    /// Number of held-out rows verified after fitting.
    pub verification_polynomials: usize,
    /// Held-out coefficient slots of the actual recurrence before cancellation;
    /// not an independent-equation count. Absent in older records.
    #[serde(default)]
    pub verification_equations: usize,
    /// Number of candidates that were actually solved before this one was accepted.
    pub candidates_tried: usize,
    /// Bounds of the winning candidate configuration.
    pub options: RecurrenceOptionsJson,
    /// Deterministic work units charged by the search that found this
    /// recurrence (see [`AdaptiveSearchBudget`]); present only when the search
    /// ran under a work limit (`polytool recurrence --max-work`), absent in
    /// older records.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub work_units: Option<u64>,
}

/// Search-space options for the recurrence that was found.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceOptionsJson {
    /// Maximum degree in `t` of the coefficient polynomials.
    pub var_deg: usize,
    /// Maximum degree in `n` of the coefficient polynomials.
    pub idx_deg: usize,
    /// Maximum derivative order.
    pub diff_deg: usize,
    /// Maximum lag (number of previous rows used).
    pub rec_len: usize,
    /// False when an inhomogeneous term was allowed.
    pub homogeneous: bool,
    /// Maximum degree in `t` of the inhomogeneous term.
    pub inhomo_var_deg: usize,
    /// Maximum degree in `n` of the inhomogeneous term.
    pub inhomo_idx_deg: usize,
    /// Degree bound in `t` of the left-hand factor `f(n,t)`.
    pub denom_var_deg: usize,
    /// Degree bound in `n` of the left-hand factor `f(n,t)`.
    pub denom_idx_deg: usize,
    /// Whether terms multiplied by `(-1)^n` were allowed.
    pub alternating_sign: bool,
}

impl From<&RecurrenceOptions> for RecurrenceOptionsJson {
    fn from(opts: &RecurrenceOptions) -> Self {
        Self {
            var_deg: opts.var_deg,
            idx_deg: opts.idx_deg,
            diff_deg: opts.diff_deg,
            rec_len: opts.rec_len,
            homogeneous: opts.homogeneous,
            inhomo_var_deg: opts.inhomo_var_deg,
            inhomo_idx_deg: opts.inhomo_idx_deg,
            denom_var_deg: opts.denom_var_deg,
            denom_idx_deg: opts.denom_idx_deg,
            alternating_sign: opts.alternating_sign,
        }
    }
}

impl RecurrenceJson {
    /// Render the initial rows and recurrence, one `$...$` formula per line,
    /// with every recurrence coefficient factored over `Q`.
    ///
    /// The recurrence starts immediately after all supplied initial rows.
    /// Search records include a LaTeX comment marking their finite-fit status.
    /// Each coefficient is factored exactly and verified by re-expansion (see
    /// [`BivarPoly::irreducible_factorization`]); a coefficient whose
    /// factorization does not finish within the default budget stays expanded
    /// and a `%` comment line says so.  The stored `latex` field is ignored.
    /// Invalid records, insufficient initial rows and index overflow return
    /// an error.  [`Self::to_latex_expanded`] never factors.
    pub fn to_latex(&self) -> Result<String, RecurrenceError> {
        self.latex_presentation(true)
    }

    /// Like [`Self::to_latex`], with every coefficient fully expanded.
    pub fn to_latex_expanded(&self) -> Result<String, RecurrenceError> {
        self.latex_presentation(false)
    }

    fn latex_presentation(&self, factored: bool) -> Result<String, RecurrenceError> {
        let (rule, first_index, seeds) = self.to_recurrence_parts()?;
        if seeds.len() < rule.max_offset() {
            return Err(RecurrenceError::InvalidInput(
                "not enough initial rows for LaTeX recurrence definition".into(),
            ));
        }
        let start = first_index.checked_add(seeds.len()).ok_or_else(|| {
            RecurrenceError::InvalidInput("LaTeX recurrence index overflows".into())
        })?;
        let mut lines = Vec::new();
        if self.search.is_some() {
            lines.push("% Candidate fitted to finite rows; validity for all n is unproved.".into());
        }
        let (formula, notes) = if factored {
            let report = rule.to_latex_factored_with_budget(
                crate::factorization::DEFAULT_RECOMBINATION_SUBSET_LIMIT,
            );
            (report.formula, report.unfactored_notes)
        } else {
            (rule.to_latex(), Vec::new())
        };
        for note in notes {
            lines.push(format!("% {note}"));
        }
        for (offset, row) in seeds.into_iter().enumerate() {
            let polynomial = BivarPoly { coeffs: vec![row] }.to_latex();
            lines.push(format!("$P_{{{}}} = {polynomial}$", first_index + offset));
        }
        lines.push(format!("${formula} \\quad (n \\geq {start})$"));
        Ok(lines.join("\n"))
    }

    /// Refresh the optional LaTeX presentation from this record's exact data.
    pub fn with_latex(mut self) -> Result<Self, RecurrenceError> {
        self.latex = Some(self.to_latex()?);
        Ok(self)
    }

    /// Relabel the first row, substituting the index in every coefficient.
    ///
    /// Seeds and search counts are unchanged. Alternating signs are shifted
    /// together with the coefficient polynomials, so replay gives the same rows.
    pub fn reindexed(&self, first_index: usize) -> Result<Self, RecurrenceError> {
        let (rule, old_index, seeds) = self.to_recurrence_parts()?;
        let difference = old_index as i128 - first_index as i128;
        let shift = i64::try_from(difference).map_err(|_| {
            RecurrenceError::InvalidInput("recurrence index shift exceeds i64".into())
        })?;
        let rule = rule.shift_index(shift);
        let mut search = self.search.clone();
        if let Some(search) = &mut search {
            search.recurrence_text = rule.to_string();
        }
        let record = Self::from_recurrence_rational(&rule, first_index, &seeds, search);
        if self.latex.is_some() {
            record.with_latex()
        } else {
            Ok(record)
        }
    }

    /// Create a JSON record from a recurrence and initial rows.
    pub fn from_recurrence_rational(
        recurrence: &Recurrence,
        first_index: usize,
        initial_polynomials: &[Vec<BigRational>],
        search: Option<RecurrenceJsonSearch>,
    ) -> Self {
        Self {
            schema: RECURRENCE_JSON_SCHEMA.to_string(),
            first_index,
            initial_polynomials: initial_polynomials
                .iter()
                .map(|row| row.iter().map(format_rational_coeff).collect())
                .collect(),
            recurrence: RecurrenceJsonData::from_recurrence(recurrence),
            search,
            latex: None,
        }
    }

    /// Convert this JSON record back to a recurrence and initial rows.
    pub fn to_recurrence_parts(
        &self,
    ) -> Result<(Recurrence, usize, Vec<Vec<BigRational>>), RecurrenceError> {
        if self.schema != RECURRENCE_JSON_SCHEMA {
            return Err(RecurrenceError::InvalidInput(format!(
                "unsupported recurrence JSON schema `{}`",
                self.schema
            )));
        }
        let recurrence = self.recurrence.to_recurrence()?;
        let initial_polynomials = self
            .initial_polynomials
            .iter()
            .enumerate()
            .map(|(row_idx, row)| {
                row.iter()
                    .enumerate()
                    .map(|(coeff_idx, coeff)| {
                        parse_rational_coeff(coeff).map_err(|e| {
                            RecurrenceError::from(e).with_context(format!(
                                "invalid initial coefficient at row {}, coefficient {}",
                                row_idx + 1,
                                coeff_idx
                            ))
                        })
                    })
                    .collect::<Result<Vec<_>, _>>()
                    .map(trim_poly_rational)
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok((recurrence, self.first_index, initial_polynomials))
    }
}

impl RecurrenceJsonData {
    fn from_recurrence(recurrence: &Recurrence) -> Self {
        Self {
            terms: recurrence
                .terms
                .iter()
                .map(RecurrenceTermJson::from_term)
                .collect(),
            denominator: recurrence
                .denominator
                .as_ref()
                .map(BivarPolyJson::from_bivar),
            inhomogeneous: recurrence
                .inhomogeneous
                .as_ref()
                .map(BivarPolyJson::from_bivar),
        }
    }

    fn to_recurrence(&self) -> Result<Recurrence, RecurrenceError> {
        Ok(Recurrence {
            terms: self
                .terms
                .iter()
                .map(RecurrenceTermJson::to_term)
                .collect::<Result<Vec<_>, _>>()?,
            denominator: self
                .denominator
                .as_ref()
                .map(BivarPolyJson::to_bivar)
                .transpose()?,
            inhomogeneous: self
                .inhomogeneous
                .as_ref()
                .map(BivarPolyJson::to_bivar)
                .transpose()?,
        })
    }
}

impl RecurrenceTermJson {
    fn from_term(term: &RecurrenceTerm) -> Self {
        Self {
            offset: term.offset,
            deriv_order: term.deriv_order,
            sign: match term.sign {
                RecurrenceSign::None => "none",
                RecurrenceSign::AlternatingN => "alternating_n",
            }
            .to_string(),
            coeff: BivarPolyJson::from_bivar(&term.coeff),
        }
    }

    fn to_term(&self) -> Result<RecurrenceTerm, RecurrenceError> {
        let sign = match self.sign.as_str() {
            "none" => RecurrenceSign::None,
            "alternating_n" => RecurrenceSign::AlternatingN,
            other => {
                return Err(RecurrenceError::InvalidInput(format!(
                    "unsupported recurrence sign `{other}`"
                )))
            }
        };
        Ok(RecurrenceTerm {
            offset: self.offset,
            deriv_order: self.deriv_order,
            sign,
            coeff: self.coeff.to_bivar()?,
        })
    }
}

impl BivarPolyJson {
    fn from_bivar(poly: &BivarPoly) -> Self {
        Self {
            coeffs: poly
                .normalized_coefficients()
                .iter()
                .map(|row| row.iter().map(format_rational_coeff).collect())
                .collect(),
        }
    }

    fn to_bivar(&self) -> Result<BivarPoly, RecurrenceError> {
        let coeffs = self
            .coeffs
            .iter()
            .enumerate()
            .map(|(i, row)| {
                row.iter()
                    .enumerate()
                    .map(|(j, coeff)| {
                        parse_rational_coeff(coeff).map_err(|e| {
                            RecurrenceError::from(e).with_context(format!(
                                "invalid bivariate coefficient ({i},{j}) `{coeff}`"
                            ))
                        })
                    })
                    .collect::<Result<Vec<_>, _>>()
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(BivarPoly::from_coefficients_normalized(coeffs))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::recurrence::test_support::*;
    use num_traits::One;

    #[test]
    fn recurrence_json_roundtrip_generates_rows() {
        let polys: Vec<Vec<i64>> = vec![
            vec![1],
            vec![1],
            vec![2],
            vec![3],
            vec![5],
            vec![8],
            vec![13],
            vec![21],
        ];
        let result = find_recurrence_adaptive(&polys, &AdaptiveSearchOptions::default()).unwrap();
        let rational_polys = i64_polys_to_rational(&polys);
        let initial_count = result.recurrence.max_offset();
        let recurrence_json = RecurrenceJson::from_recurrence_rational(
            &result.recurrence,
            1,
            &rational_polys[..initial_count],
            Some(RecurrenceJsonSearch {
                recurrence_text: result.recurrence.to_string(),
                source_rows: polys.len(),
                skip_prefix: 0,
                unknowns: result.num_unknowns,
                weighted_unknowns: result.weighted_unknowns,
                equations: result.num_equations,
                fit_polynomials: result.fit_polynomials,
                verification_polynomials: result.verification_polynomials,
                verification_equations: result.verification_equations,
                candidates_tried: result.candidates_tried,
                work_units: None,
                options: RecurrenceOptionsJson::from(&result.opts),
            }),
        );

        let encoded = serde_json::to_string(&recurrence_json).unwrap();
        let decoded: RecurrenceJson = serde_json::from_str(&encoded).unwrap();
        let (recurrence, first_index, initial_polys) = decoded.to_recurrence_parts().unwrap();
        let generated = recurrence
            .generate_rows_rational(&initial_polys, first_index, polys.len())
            .unwrap();
        assert_eq!(generated, rational_polys);
    }

    #[test]
    fn recurrence_json_normalizes_ragged_bivariate_coefficients() {
        let decoded = BivarPolyJson {
            coeffs: vec![
                vec!["1".to_string()],
                vec!["0".to_string(), "1".to_string()],
            ],
        }
        .to_bivar()
        .unwrap();
        assert_eq!(
            decoded.coeffs.iter().map(Vec::len).collect::<Vec<_>>(),
            [2, 2]
        );
        assert_eq!(decoded.to_string(), "1 + nt");

        let encoded = BivarPolyJson::from_bivar(&BivarPoly {
            coeffs: vec![vec![br(1, 1)], vec![br(0, 1), br(1, 1)]],
        });
        assert_eq!(
            encoded.coeffs.iter().map(Vec::len).collect::<Vec<_>>(),
            [2, 2]
        );
    }

    #[test]
    fn recurrence_json_generates_denominator_rows() {
        let mut polys: Vec<Vec<BigRational>> = vec![vec![BigRational::one()]];
        for n in 2..=8 {
            let prev = polys.last().unwrap()[0].clone();
            polys.push(vec![prev / BigRational::from_integer(BigInt::from(n + 1))]);
        }
        let opts = RecurrenceOptions {
            var_deg: 0,
            idx_deg: 0,
            diff_deg: 0,
            rec_len: 1,
            homogeneous: true,
            denom_var_deg: 0,
            denom_idx_deg: 1,
            ..Default::default()
        };
        let recurrence =
            find_polynomial_recurrence_rational(&polys, &opts).expect("should find denominator");
        let recurrence_json =
            RecurrenceJson::from_recurrence_rational(&recurrence, 1, &polys[..1], None);
        let (recurrence, first_index, initial_polys) =
            recurrence_json.to_recurrence_parts().unwrap();
        let generated = recurrence
            .generate_rows_rational(&initial_polys, first_index, polys.len())
            .unwrap();
        assert_eq!(generated, polys);
    }
}
