//! [`BasisMapLibrary`]: the basis maps that the lift search tries, each
//! with its proved [`BasisMapTheoremReport`].

use super::map::BasisMap;
use super::published::{fingerprint_candidate, PublishedFamily};
use super::theorems::{
    check_basis_map_preserves_pf, BasisMapShape, BasisMapTheoremOptions, BasisMapTheoremReport,
};
use crate::recurrence::RecurrenceError;
use crate::vec_poly;
use num_traits::{One, Zero};

/// Rows compared to recognise duplicate and identity maps and to require
/// triangular rows.
const LIBRARY_FINGERPRINT_ROWS: usize = 8;

/// The published transformations added by
/// [`BasisMapLibraryOptions::include_published_maps`] (the derangement
/// transformation is not triangular, so it is a forward map only).
const PUBLISHED_LIBRARY_FAMILIES: [PublishedFamily; 4] = [
    PublishedFamily::Narayana(0),
    PublishedFamily::Narayana(1),
    PublishedFamily::Narayana(2),
    PublishedFamily::EulerianTransformation,
];

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
    /// Include every bundled OEIS catalog entry whose map is admitted and
    /// has a companion form (catalog families with index-dependent rules
    /// are not added; the published ones come in through
    /// `include_published_maps`).  Default false ([`BasisMapLibrary::oeis_catalog`] sets it).
    pub include_oeis_catalog: bool,
    /// Further maps, tried first.
    pub maps: Vec<BasisMap>,
    /// Include the published transformations with index-dependent rules
    /// that the backward lift can use: the Narayana transformations
    /// `N_(k,m)`, `m = 0, 1, 2`, and the Eulerian transformation `A°`
    /// (module documentation of `basis_map`, "Published transformations"),
    /// tried after the catalog maps.  Default false
    /// ([`BasisMapLibrary::oeis_catalog`] sets it).
    pub include_published_maps: bool,
    /// Options of the theorem checks.
    pub theorem: BasisMapTheoremOptions,
}

/// The maps tried by [`super::find_basis_map_lifts`].
///
/// A map is admitted when [`check_basis_map_preserves_pf`] proves at least
/// one conclusion, the family has the companion form the transport needs
/// (a rule free of `n` with a constant leader and no exceptional prefix)
/// or is a published transformation with an index-dependent rule (the
/// index-side transport certifies its lifts from Riordan coefficient
/// families),
/// its first rows are nonzero with strictly increasing degrees (so that the
/// triangular expansion applies), and it is not the identity `u^k -> t^k`.
/// Maps with the same first rows are kept once (the others are recorded as
/// aliases).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct BasisMapLibrary {
    entries: Vec<BasisMapLibraryEntry>,
    /// The first [`LIBRARY_FINGERPRINT_ROWS`] basis rows of each entry (the
    /// duplicate test), parallel to `entries`.
    first_rows: Vec<Vec<Vec<crate::BigRational>>>,
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
                // Catalog maps keep the companion-form admission (and so
                // their order); the published index-dependent maps are
                // added below.
                if map.companion_operators().is_err() {
                    continue;
                }
                let _ = library.admit(map, &options.theorem);
            }
        }
        if options.include_published_maps {
            for family in PUBLISHED_LIBRARY_FAMILIES {
                let _ = library.admit(family.canonical_map(), &options.theorem);
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
            include_published_maps: true,
            ..BasisMapLibraryOptions::default()
        })
    }

    /// The entries in the order tried.
    pub fn entries(&self) -> &[BasisMapLibraryEntry] {
        &self.entries
    }

    /// Keep only the entries for which `keep` holds (in order).
    pub(in crate::recurrence) fn retain_entries(
        &mut self,
        mut keep: impl FnMut(&BasisMapLibraryEntry) -> bool,
    ) {
        let kept: Vec<bool> = self.entries.iter().map(&mut keep).collect();
        let mut flags = kept.iter();
        self.entries.retain(|_| *flags.next().unwrap_or(&false));
        let mut flags = kept.iter();
        self.first_rows.retain(|_| *flags.next().unwrap_or(&false));
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
        let companion = map.companion_operators().is_ok();
        if !companion && fingerprint_candidate(&map)?.is_none() {
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
        if let Some(position) = self.first_rows.iter().position(|other| *other == rows) {
            self.entries[position].aliases.push(map.name().to_string());
            return Ok(false);
        }
        let theorem = check_basis_map_preserves_pf(&map, options)?;
        if !theorem.is_proved() {
            return Ok(false);
        }
        let published = matches!(
            theorem.shape,
            BasisMapShape::Narayana { .. }
                | BasisMapShape::EulerianTransformation
                | BasisMapShape::DerangementTransformation
        );
        if !companion && !published {
            return Ok(false);
        }
        self.entries.push(BasisMapLibraryEntry {
            map,
            theorem,
            aliases: Vec::new(),
        });
        self.first_rows.push(rows);
        Ok(true)
    }
}
