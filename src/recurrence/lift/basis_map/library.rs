//! [`BasisMapLibrary`]: the basis maps that the lift search tries, each
//! with its proved [`BasisMapTheoremReport`].

use super::map::BasisMap;
use super::theorems::{
    check_basis_map_preserves_pf, BasisMapTheoremOptions, BasisMapTheoremReport,
};
use crate::recurrence::RecurrenceError;
use crate::vec_poly;
use num_traits::{One, Zero};

/// Rows compared to recognise duplicate and identity maps and to require
/// triangular rows.
const LIBRARY_FINGERPRINT_ROWS: usize = 8;

/// One map of a [`BasisMapLibrary`].
#[derive(Debug, Clone, PartialEq)]
pub struct BasisMapLibraryEntry {
    /// The map.
    pub map: BasisMap,
    /// Its theorem report (proved).
    pub theorem: BasisMapTheoremReport,
    /// Names of other admitted families with the same first rows (the same
    /// map up to the checked rows; only the first is tried).
    pub aliases: Vec<String>,
}

/// Options of [`BasisMapLibrary::build`].
#[derive(Debug, Clone, Default)]
pub struct BasisMapLibraryOptions {
    /// Include every bundled OEIS catalog entry whose map is admitted.
    /// Default false ([`BasisMapLibrary::oeis_catalog`] sets it).
    pub include_oeis_catalog: bool,
    /// Further maps, tried first.
    pub maps: Vec<BasisMap>,
    /// Options of the theorem checks.
    pub theorem: BasisMapTheoremOptions,
}

/// The maps tried by [`super::find_basis_map_lifts`].
///
/// A map is admitted when [`check_basis_map_preserves_pf`] proves at least
/// one conclusion, the family has the companion form the transport needs
/// (a rule free of `n` with a constant leader and no exceptional prefix),
/// its first rows are nonzero with strictly increasing degrees (so that the
/// triangular expansion applies), and it is not the identity `u^k -> t^k`.
/// Maps with the same first rows are kept once (the others are recorded as
/// aliases).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct BasisMapLibrary {
    entries: Vec<BasisMapLibraryEntry>,
}

impl BasisMapLibrary {
    /// Build a library (see [`BasisMapLibraryOptions`]).
    ///
    /// # Errors
    ///
    /// [`RecurrenceError`] from a theorem check of an explicitly supplied map
    /// (catalog entries that cannot be checked are skipped).
    pub fn build(options: &BasisMapLibraryOptions) -> Result<Self, RecurrenceError> {
        let mut library = Self::default();
        for map in &options.maps {
            library.admit(map.clone(), &options.theorem)?;
        }
        if options.include_oeis_catalog {
            for definition in crate::oeis::catalog() {
                let Ok(map) = BasisMap::from_oeis(definition.id) else {
                    continue;
                };
                let _ = library.admit(map, &options.theorem);
            }
        }
        Ok(library)
    }

    /// The admitted maps of the bundled OEIS catalog.
    ///
    /// # Errors
    ///
    /// As [`BasisMapLibrary::build`].
    pub fn oeis_catalog() -> Result<Self, RecurrenceError> {
        Self::build(&BasisMapLibraryOptions {
            include_oeis_catalog: true,
            ..BasisMapLibraryOptions::default()
        })
    }

    /// The entries in the order tried.
    pub fn entries(&self) -> &[BasisMapLibraryEntry] {
        &self.entries
    }

    /// Number of entries.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether the library is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Check and admit one map; returns whether it was added as a new entry.
    ///
    /// # Errors
    ///
    /// The errors of [`check_basis_map_preserves_pf`] and of generating the
    /// first rows.
    pub fn admit(
        &mut self,
        map: BasisMap,
        options: &BasisMapTheoremOptions,
    ) -> Result<bool, RecurrenceError> {
        if map.companion_operators().is_err() {
            return Ok(false);
        }
        let rows = map.basis_rows(LIBRARY_FINGERPRINT_ROWS)?;
        let mut previous: Option<usize> = None;
        for row in &rows {
            if vec_poly::is_zero(row) {
                return Ok(false);
            }
            let degree = vec_poly::degree(row).unwrap_or(0);
            if previous.is_some_and(|p| degree <= p) {
                return Ok(false);
            }
            previous = Some(degree);
        }
        let identity = rows.iter().enumerate().all(|(k, row)| {
            row.len() == k + 1 && row[k].is_one() && row[..k].iter().all(Zero::is_zero)
        });
        if identity {
            return Ok(false);
        }
        if let Some(existing) = self.entries.iter_mut().find(|entry| {
            entry
                .map
                .basis_rows(LIBRARY_FINGERPRINT_ROWS)
                .is_ok_and(|other| other == rows)
        }) {
            existing.aliases.push(map.name().to_string());
            return Ok(false);
        }
        let theorem = check_basis_map_preserves_pf(&map, options)?;
        if !theorem.is_proved() {
            return Ok(false);
        }
        self.entries.push(BasisMapLibraryEntry {
            map,
            theorem,
            aliases: Vec::new(),
        });
        Ok(true)
    }
}
