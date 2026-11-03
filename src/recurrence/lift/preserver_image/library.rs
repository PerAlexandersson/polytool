//! The library of proved parents of the preserver-image lifts
//! ([`PreserverParentLibrary`]): named families with their defining
//! recurrences and, optionally, the OEIS catalog entries whose parent
//! statement the provers settle; the parent proofs are computed once per
//! library.

use super::interval::interval_parent_report;
use crate::oeis;
use crate::recurrence::lift::hadamard_factor::prepare::{prepare_rule, PreparedRule};
use crate::recurrence::lift::parent_proof::{prove_parent, ParentRecord};
use crate::recurrence::lift::section_certify::x_free_leader;
use crate::recurrence::{
    check_three_term_recurrence_theorem, BivarPoly, GammaLiftParent, RealRootednessTheoremOptions,
    RealRootednessTheoremReport, Recurrence, RecurrenceError, RecurrenceJson, RecurrenceTerm,
    VeroneseParentProof, VeroneseParentProofRoute, VeroneseParentProofStatus,
};
use crate::BigRational;
use num_bigint::BigInt;
use num_traits::Zero;

type Q = BigRational;

/// Largest number of parent rows kept per parent.
const LIBRARY_ROWS: usize = 40;

/// The root region every row of a library parent is proved to have.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PreserverParentRegion {
    /// Every row from the first index on is zero or real-rooted with every
    /// root in `(-inf, 0]` (the shared parent provers).
    Nonpositive,
    /// Every row from the first index on is real-rooted with every root in
    /// the bounded interval of [`PreserverParent::interval_proof`] (the
    /// Ma--Wang interval check or the three-term theorem with a root
    /// region), and the statement above is not proved.
    Interval,
    /// Every row from the first index on is real-rooted (the three-term
    /// theorem without a root region).
    RealLine,
    /// Not proved (kept only on request).
    Unproved,
}

/// One candidate parent with its proofs.
#[derive(Debug, Clone)]
pub struct PreserverParent {
    /// A short name (`eulerian`, an OEIS A-number, ...).
    pub name: String,
    /// The parent rule with its first index and identifying seeds.
    pub recurrence: RecurrenceJson,
    /// The parent statement with roots in `(-inf, 0]` and its attempts.
    pub proof: VeroneseParentProof,
    /// The three-term theorem on the whole line, when it proves the rows
    /// real-rooted and the statement above is not proved.
    pub whole_line_proof: Option<RealRootednessTheoremReport>,
    /// A theorem report proving every row real-rooted with every root in a
    /// bounded interval (the `preserver_image` interval provers), whether or
    /// not the statement above is proved.
    pub interval_proof: Option<RealRootednessTheoremReport>,
    /// The proved region.
    pub region: PreserverParentRegion,
    /// The rows `P_(p0), P_(p0 + 1), ...` (the first rows of the family).
    rows: Vec<Vec<Q>>,
}

impl PreserverParent {
    /// The first index `p0`.
    pub fn first_index(&self) -> usize {
        self.recurrence.first_index
    }

    /// The precomputed rows from `p0` on.
    pub fn rows(&self) -> &[Vec<Q>] {
        &self.rows
    }
}

/// Bounds of [`PreserverParentLibrary::build`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PreserverParentLibraryOptions {
    /// Include the named families ([`PreserverParentLibrary::named_family_records`]).
    /// Default true.
    pub include_named_families: bool,
    /// Include OEIS catalog entries ([`crate::oeis::catalog`], in catalog
    /// order).  Default false.
    pub include_oeis_catalog: bool,
    /// Largest number of catalog entries examined.  Default 1000.
    pub max_catalog_entries: usize,
    /// Largest order of a catalog rule.  Default 2.
    pub max_catalog_order: usize,
    /// Keep parents whose statement is not proved.  Default false.
    pub keep_unproved: bool,
    /// Try the root-domain parent route.  Default true.
    pub try_root_domain: bool,
}

impl Default for PreserverParentLibraryOptions {
    fn default() -> Self {
        Self {
            include_named_families: true,
            include_oeis_catalog: false,
            max_catalog_entries: 1000,
            max_catalog_order: 2,
            keep_unproved: false,
            try_root_domain: true,
        }
    }
}

/// The parents tried by [`crate::recurrence::find_preserver_image_lifts`],
/// in a fixed order, with their proofs computed once.
#[derive(Debug, Clone, Default)]
pub struct PreserverParentLibrary {
    /// The parents.
    pub parents: Vec<PreserverParent>,
    /// Catalog entries examined.
    pub catalog_entries_examined: usize,
    /// Why records were skipped (name and reason).
    pub skipped: Vec<(String, String)>,
}

fn bivariate(monomials: &[(usize, usize, i64)]) -> BivarPoly {
    let n_size = monomials.iter().map(|m| m.0 + 1).max().unwrap_or(1);
    let x_size = monomials.iter().map(|m| m.1 + 1).max().unwrap_or(1);
    let mut coeffs = vec![vec![Q::zero(); x_size]; n_size];
    for &(i, j, c) in monomials {
        coeffs[i][j] += Q::from_integer(BigInt::from(c));
    }
    BivarPoly::from_coefficients_normalized(coeffs)
}

fn integer_rows(rows: &[&[i64]]) -> Vec<Vec<Q>> {
    rows.iter()
        .map(|row| {
            row.iter()
                .map(|&c| Q::from_integer(BigInt::from(c)))
                .collect()
        })
        .collect()
}

fn record(
    terms: Vec<RecurrenceTerm>,
    denominator: Option<BivarPoly>,
    first_index: usize,
    seeds: &[&[i64]],
) -> RecurrenceJson {
    RecurrenceJson::from_recurrence_rational(
        &Recurrence {
            terms,
            denominator,
            inhomogeneous: None,
        },
        first_index,
        &integer_rows(seeds),
        None,
    )
}

impl PreserverParentLibrary {
    /// The named families with their defining recurrences, in the order
    /// tried: Eulerian `A_n` (`n >= 1`), type-B Eulerian, Narayana
    /// (`sum_k N(n, k) x^(k-1)`), the peak polynomials, the binomial rows
    /// `(1 + x)^n`, Delannoy, the symmetric parts `a_n` of the `1/k`-Eulerian
    /// polynomials for `k = 2, 3`, and the whole-line families Chebyshev
    /// `T_n`, `U_n`, Hermite `He_n` and the shifted Legendre polynomials (see
    /// [`crate::sequences`] for each generator).  A library keeps only the
    /// proved ones: the peak polynomials and the `1/k`-Eulerian symmetric
    /// parts are proved by lifts (gamma, symmetric companion), not by the
    /// parent provers, so they are listed in [`Self::skipped`].
    pub fn named_family_records() -> Vec<(String, RecurrenceJson)> {
        let mut out = vec![
            ("eulerian".into(), GammaLiftParent::eulerian().recurrence),
            (
                "type_b_eulerian".into(),
                GammaLiftParent::type_b_eulerian().recurrence,
            ),
            ("narayana".into(), GammaLiftParent::narayana().recurrence),
        ];
        let (rule, seeds, first_index) =
            crate::sequences::literature::peak_polynomials_recurrence();
        out.push((
            "peak".into(),
            RecurrenceJson::from_recurrence_rational(&rule, first_index, &seeds, None),
        ));
        out.push((
            "binomial".into(),
            record(
                vec![RecurrenceTerm::new(
                    1,
                    0,
                    bivariate(&[(0, 0, 1), (0, 1, 1)]),
                )],
                None,
                0,
                &[&[1]],
            ),
        ));
        out.push((
            "delannoy".into(),
            record(
                vec![
                    RecurrenceTerm::new(1, 0, bivariate(&[(0, 0, 1), (0, 1, 1)])),
                    RecurrenceTerm::new(2, 0, bivariate(&[(0, 1, 1)])),
                ],
                None,
                0,
                &[&[1], &[1, 1]],
            ),
        ));
        for k in [2, 3] {
            let parent =
                GammaLiftParent::one_over_k_eulerian_symmetric_part(&Q::from_integer(k.into()));
            out.push((
                format!("one_over_k_eulerian_symmetric_part_a_k{k}"),
                parent.recurrence,
            ));
        }
        out.push((
            "chebyshev_t".into(),
            record(
                vec![
                    RecurrenceTerm::new(1, 0, bivariate(&[(0, 1, 2)])),
                    RecurrenceTerm::new(2, 0, bivariate(&[(0, 0, -1)])),
                ],
                None,
                0,
                &[&[1], &[0, 1]],
            ),
        ));
        out.push((
            "chebyshev_u".into(),
            record(
                vec![
                    RecurrenceTerm::new(1, 0, bivariate(&[(0, 1, 2)])),
                    RecurrenceTerm::new(2, 0, bivariate(&[(0, 0, -1)])),
                ],
                None,
                0,
                &[&[1], &[0, 2]],
            ),
        ));
        out.push((
            "hermite_he".into(),
            record(
                vec![
                    RecurrenceTerm::new(1, 0, bivariate(&[(0, 1, 1)])),
                    RecurrenceTerm::new(2, 0, bivariate(&[(1, 0, -1), (0, 0, 1)])),
                ],
                None,
                0,
                &[&[1], &[0, 1]],
            ),
        ));
        out.push((
            "shifted_legendre".into(),
            record(
                vec![
                    RecurrenceTerm::new(
                        1,
                        0,
                        bivariate(&[(1, 1, 4), (1, 0, -2), (0, 1, -2), (0, 0, 1)]),
                    ),
                    RecurrenceTerm::new(2, 0, bivariate(&[(1, 0, -1), (0, 0, 1)])),
                ],
                Some(bivariate(&[(1, 0, 1)])),
                1,
                &[&[-1, 2], &[1, -6, 6]],
            ),
        ));
        out
    }

    /// The library of the named families
    /// ([`Self::named_family_records`]) with their proofs.
    ///
    /// # Errors
    ///
    /// [`RecurrenceError`] only from an internal failure of the provers.
    pub fn named_families() -> Result<Self, RecurrenceError> {
        Self::build(&PreserverParentLibraryOptions::default())
    }

    /// Build the library: the named families and/or the catalog entries,
    /// each validated (homogeneous, positive offsets, a leader free of `x`
    /// and of derivatives, nonvanishing at every target), proved (the shared
    /// parent provers with the parent's first index, the interval provers,
    /// then the three-term theorem on the whole line) and kept when proved (or always, with
    /// [`PreserverParentLibraryOptions::keep_unproved`]).  Records that fail
    /// validation are listed in [`Self::skipped`].  Deterministic.
    ///
    /// # Errors
    ///
    /// [`RecurrenceError`] only from an internal failure of the provers.
    pub fn build(options: &PreserverParentLibraryOptions) -> Result<Self, RecurrenceError> {
        let mut library = Self::default();
        if options.include_named_families {
            for (name, record) in Self::named_family_records() {
                library.add(name, &record, options)?;
            }
        }
        if options.include_oeis_catalog {
            for entry in oeis::catalog().iter().take(options.max_catalog_entries) {
                library.add_catalog_entry(entry, options)?;
            }
        }
        Ok(library)
    }

    /// A library of the listed catalog entries only (in the given order),
    /// built as by [`Self::build`] (for focused tests).
    #[cfg(test)]
    pub(super) fn from_catalog_entries(
        ids: &[&str],
        options: &PreserverParentLibraryOptions,
    ) -> Result<Self, RecurrenceError> {
        let mut library = Self::default();
        for id in ids {
            let entry = oeis::by_id(id).expect("a catalog entry");
            library.add_catalog_entry(entry, options)?;
        }
        Ok(library)
    }

    /// Validate, prove and (when proved or requested) add one catalog entry.
    fn add_catalog_entry(
        &mut self,
        entry: &oeis::OeisSequenceDefinition,
        options: &PreserverParentLibraryOptions,
    ) -> Result<(), RecurrenceError> {
        self.catalog_entries_examined += 1;
        let (rule, rows) = match entry.recurrence_parts() {
            Ok(parts) => parts,
            Err(error) => {
                self.skipped.push((entry.id.into(), error.to_string()));
                return Ok(());
            }
        };
        if rule.max_offset() > options.max_catalog_order {
            self.skipped.push((
                entry.id.into(),
                format!("order {} is above the bound", rule.max_offset()),
            ));
            return Ok(());
        }
        let record = RecurrenceJson::from_recurrence_rational(
            &rule,
            entry.recurrence_first_index,
            &rows,
            None,
        );
        self.add(entry.id.into(), &record, options)
    }

    /// Validate, prove and (when proved or requested) add one record.
    pub(crate) fn add(
        &mut self,
        name: String,
        record: &RecurrenceJson,
        options: &PreserverParentLibraryOptions,
    ) -> Result<(), RecurrenceError> {
        let prepared = match prepare(record) {
            Ok(prepared) => prepared,
            Err(reason) => {
                self.skipped.push((name, reason));
                return Ok(());
            }
        };
        let canonical = RecurrenceJson::from_recurrence_rational(
            &prepared.rule,
            prepared.first_index,
            &prepared.seeds,
            None,
        );
        if self.parents.iter().any(|p| p.recurrence == canonical) {
            self.skipped
                .push((name, "the same record is already in the library".into()));
            return Ok(());
        }
        let proof = prove_parent(
            &ParentRecord {
                rule: &prepared.rule,
                seeds: &prepared.seeds,
                first_index: prepared.first_index,
            },
            prepared.first_index,
            &[],
            options.try_root_domain,
        )?;
        let proved = proof.status == VeroneseParentProofStatus::Proved
            && proof.route != Some(VeroneseParentProofRoute::Attestation);
        let mut whole_line_proof = None;
        if !proved {
            let report = check_three_term_recurrence_theorem(
                &prepared.rule,
                prepared.first_index,
                &prepared.seeds,
                &RealRootednessTheoremOptions {
                    root_region: None,
                    ..RealRootednessTheoremOptions::default()
                },
            )?;
            if report.is_proved() {
                whole_line_proof = Some(report);
            }
        }
        let interval_proof =
            interval_parent_report(&prepared.rule, prepared.first_index, &prepared.seeds)?;
        let region = if proved {
            PreserverParentRegion::Nonpositive
        } else if interval_proof.is_some() {
            PreserverParentRegion::Interval
        } else if whole_line_proof.is_some() {
            PreserverParentRegion::RealLine
        } else {
            PreserverParentRegion::Unproved
        };
        if region == PreserverParentRegion::Unproved && !options.keep_unproved {
            self.skipped.push((name, "the parent is not proved".into()));
            return Ok(());
        }
        let rows = match prepared.rows(LIBRARY_ROWS) {
            Ok(rows) => rows,
            Err(error) => {
                self.skipped.push((name, error.to_string()));
                return Ok(());
            }
        };
        self.parents.push(PreserverParent {
            name,
            recurrence: canonical,
            proof,
            whole_line_proof,
            interval_proof,
            region,
            rows,
        });
        Ok(())
    }
}

/// Validate a record for use as a parent; `Err` gives the reason.
fn prepare(record: &RecurrenceJson) -> Result<PreparedRule, String> {
    let (rule, first_index, rows) = record.to_recurrence_parts().map_err(|e| e.to_string())?;
    let prepared = prepare_rule(&rule, &rows, first_index, "parent").map_err(|e| e.to_string())?;
    x_free_leader(&prepared.relation).map_err(|e| e.to_string())?;
    Ok(prepared)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequences;

    #[test]
    fn named_records_reproduce_their_generators() {
        let records = PreserverParentLibrary::named_family_records();
        let rows = |name: &str| -> Vec<Vec<Q>> {
            let (_, record) = records.iter().find(|(n, _)| n == name).unwrap();
            prepare(record).unwrap().rows(12).unwrap()
        };
        let rational = |rows: Vec<Vec<BigInt>>| -> Vec<Vec<Q>> {
            rows.into_iter()
                .map(|row| {
                    crate::vec_poly::trim(&row.into_iter().map(Q::from_integer).collect::<Vec<_>>())
                })
                .collect()
        };
        assert_eq!(
            rows("chebyshev_t"),
            rational(sequences::chebyshev_polynomials_t_bigint(11))
        );
        assert_eq!(
            rows("chebyshev_u"),
            rational(sequences::chebyshev_polynomials_u_bigint(11))
        );
        assert_eq!(
            rows("hermite_he"),
            rational(sequences::hermite_polynomials_bigint(11))
        );
        assert_eq!(
            rows("shifted_legendre"),
            rational(sequences::shifted_legendre_polynomials_bigint(12))[1..].to_vec()
        );
        assert_eq!(
            rows("delannoy"),
            rational(sequences::literature::delannoy_polynomials_bigint(11))
        );
        assert_eq!(
            rows("eulerian"),
            rational(sequences::eulerian_polynomials_bigint(12))[1..].to_vec()
        );
    }

    #[test]
    fn the_named_library_proves_its_parents() {
        let library = PreserverParentLibrary::named_families().unwrap();
        let region = |name: &str| {
            library
                .parents
                .iter()
                .find(|p| p.name == name)
                .map(|p| p.region)
        };
        assert_eq!(region("eulerian"), Some(PreserverParentRegion::Nonpositive));
        assert_eq!(region("narayana"), Some(PreserverParentRegion::Nonpositive));
        assert_eq!(region("binomial"), Some(PreserverParentRegion::Nonpositive));
        assert_eq!(region("hermite_he"), Some(PreserverParentRegion::RealLine));
        assert_eq!(region("chebyshev_t"), Some(PreserverParentRegion::Interval));
        assert_eq!(
            region("shifted_legendre"),
            Some(PreserverParentRegion::Interval)
        );
        let interval = |name: &str| {
            library
                .parents
                .iter()
                .find(|p| p.name == name)
                .and_then(|p| p.interval_proof.as_ref())
                .and_then(|report| report.conclusion.as_ref())
                .and_then(|c| c.root_region.as_ref())
                .map(|r| r.to_string())
        };
        assert_eq!(interval("chebyshev_t").as_deref(), Some("(-1, 1)"));
        assert_eq!(interval("shifted_legendre").as_deref(), Some("(0, 1)"));
        assert_eq!(interval("hermite_he"), None);
    }
}
