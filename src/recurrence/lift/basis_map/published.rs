//! The published basis maps with index-dependent rules (module
//! documentation of `basis_map`, "Published transformations"): their
//! canonical families, constructors, and the all-index recognition of a
//! map as one of them.

use super::map::BasisMap;
use crate::recurrence::{
    certify_recurrence_equality, BivarPoly, Recurrence, RecurrenceEqualityCertificate,
    RecurrenceError, RecurrenceJson, RecurrenceTerm,
};
use crate::vec_poly;
use crate::BigRational;
use num_traits::{One, Zero};

type Q = BigRational;

/// Rows compared before the all-index recognition is attempted.
const FINGERPRINT_ROWS: usize = 6;

/// Largest Narayana parameter `m` recognised (the rules grow with `m`).
pub(super) const MAX_NARAYANA_PARAMETER: usize = 64;

/// A published family `P_k` with a real-rootedness theorem for its basis
/// map `u^k -> P_k`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum PublishedFamily {
    /// `N_(k,m)(t) = sum_j C(k,j) C(k+m,j) / C(m+j,j) t^j` (Mao and Wang).
    Narayana(usize),
    /// `1, t A_1(t), t A_2(t), ...` with the Eulerian polynomials `A_k`
    /// (Athanasiadis' `A°`).
    EulerianTransformation,
    /// The derangement polynomials `d_k(t)` by excedances, `d_0 = 1`
    /// (Brändén and Solus).
    DerangementTransformation,
}

fn q(value: i64) -> Q {
    Q::from_integer(value.into())
}

fn bivariate(rows: Vec<Vec<Q>>) -> BivarPoly {
    BivarPoly::from_coefficients_normalized(rows)
}

impl PublishedFamily {
    /// The canonical rule (absolute index `n`, the basis element of `u^n`
    /// at `n`) and seeds, as proved in RealRooted (module documentation).
    pub(super) fn canonical_parts(self) -> (Recurrence, Vec<Vec<Q>>) {
        match self {
            // (n + 2m) P_n = ((n + 2m) + (3n + 2m - 2) t) P_(n-1)
            //               + 2 t (1 - t) P'_(n-1),  n >= 1, P_0 = 1.
            Self::Narayana(m) => {
                let m = q(m as i64);
                let two_m = &m + &m;
                let rule = Recurrence {
                    terms: vec![
                        RecurrenceTerm::new(
                            1,
                            0,
                            bivariate(vec![vec![two_m.clone(), &two_m - q(2)], vec![q(1), q(3)]]),
                        ),
                        RecurrenceTerm::new(1, 1, bivariate(vec![vec![q(0), q(2), q(-2)]])),
                    ],
                    denominator: Some(bivariate(vec![vec![two_m], vec![q(1)]])),
                    inhomogeneous: None,
                };
                (rule, vec![vec![q(1)]])
            }
            // P_n = n t P_(n-1) + t (1 - t) P'_(n-1), n >= 1, P_0 = 1.
            Self::EulerianTransformation => {
                let rule = Recurrence {
                    terms: vec![
                        RecurrenceTerm::new(1, 0, bivariate(vec![vec![q(0)], vec![q(0), q(1)]])),
                        RecurrenceTerm::new(1, 1, bivariate(vec![vec![q(0), q(1), q(-1)]])),
                    ],
                    denominator: None,
                    inhomogeneous: None,
                };
                (rule, vec![vec![q(1)]])
            }
            // P_n = (n - 1) t (P_(n-1) + P_(n-2)) + t (1 - t) P'_(n-1),
            // n >= 2, P_0 = 1, P_1 = 0.
            Self::DerangementTransformation => {
                let index_times_t = || bivariate(vec![vec![q(0), q(-1)], vec![q(0), q(1)]]);
                let rule = Recurrence {
                    terms: vec![
                        RecurrenceTerm::new(1, 0, index_times_t()),
                        RecurrenceTerm::new(1, 1, bivariate(vec![vec![q(0), q(1), q(-1)]])),
                        RecurrenceTerm::new(2, 0, index_times_t()),
                    ],
                    denominator: None,
                    inhomogeneous: None,
                };
                (rule, vec![vec![q(1)], vec![q(0)]])
            }
        }
    }

    /// The canonical map `u^k -> P_k`.
    pub(super) fn canonical_map(self) -> BasisMap {
        let (rule, seeds) = self.canonical_parts();
        BasisMap::new(&self.name(), rule, 0, seeds)
            .expect("the canonical published families replay")
    }

    /// A short name.
    pub(super) fn name(self) -> String {
        match self {
            Self::Narayana(m) => format!("Narayana transformation N_(k,{m})"),
            Self::EulerianTransformation => "Eulerian transformation A°".into(),
            Self::DerangementTransformation => "derangement transformation".into(),
        }
    }

    /// The canonical record with first index `first_index` (the basis
    /// element of `u^k` at `first_index + k`).
    fn canonical_record(self, first_index: usize) -> RecurrenceJson {
        let (rule, seeds) = self.canonical_parts();
        let rule = rule.shift_index(-(first_index as i64));
        RecurrenceJson::from_recurrence_rational(&rule, first_index, &seeds, None)
    }
}

impl BasisMap {
    /// The Narayana transformation `u^k -> N_(k,m)(t)`, `N_(k,m)(t) =
    /// sum_j C(k,j) C(k+m,j) / C(m+j,j) t^j` (Mao and Wang,
    /// arXiv:2607.01572, Eq. (2); `m = 0, 1` are the types B and A), given
    /// by the RealRooted recurrence `narayanaPolynomial_deriv_lag_rec`
    /// (module documentation of `basis_map`, "Published transformations").
    ///
    /// # Panics
    ///
    /// When `m` exceeds 64 (the recognised range).
    pub fn narayana_transformation(m: usize) -> Self {
        assert!(
            m <= MAX_NARAYANA_PARAMETER,
            "the Narayana parameter is at most {MAX_NARAYANA_PARAMETER}"
        );
        PublishedFamily::Narayana(m).canonical_map()
    }

    /// The Eulerian transformation `A°`: `u^0 -> 1`, `u^k -> t A_k(t)` for
    /// `k >= 1` with the Eulerian polynomials `A_k` (Brändén and Jochemko;
    /// Athanasiadis, arXiv:2302.00754).
    pub fn eulerian_transformation() -> Self {
        PublishedFamily::EulerianTransformation.canonical_map()
    }

    /// The derangement transformation `u^k -> d_k(t)`, the derangement
    /// polynomials by excedances with `d_0 = 1`, `d_1 = 0` (Brändén and
    /// Solus, arXiv:1808.04141, Section 3.2).  Not injective (`d_1 = 0`): a
    /// forward map only.
    pub fn derangement_transformation() -> Self {
        PublishedFamily::DerangementTransformation.canonical_map()
    }
}

/// The candidate published family whose first rows equal those of `map`.
pub(super) fn fingerprint_candidate(
    map: &BasisMap,
) -> Result<Option<PublishedFamily>, RecurrenceError> {
    let rows = map.basis_rows(FINGERPRINT_ROWS)?;
    let mut candidates = vec![
        PublishedFamily::EulerianTransformation,
        PublishedFamily::DerangementTransformation,
    ];
    // N_(2,m) = 1 + c t + t^2 with c = 2 (m + 2) / (m + 1), so m = (4 - c) / (c - 2).
    if let Some(row) = rows.get(2) {
        let row = vec_poly::trim(row);
        if row.len() == 3 && row[0].is_one() && row[2].is_one() {
            let c = &row[1];
            let two = q(2);
            if c > &two {
                let m = (q(4) - c) / (c - &two);
                if m.is_integer() && m >= Q::zero() {
                    if let Ok(m) = usize::try_from(m.to_integer()) {
                        if m <= MAX_NARAYANA_PARAMETER {
                            candidates.push(PublishedFamily::Narayana(m));
                        }
                    }
                }
            }
        }
    }
    for candidate in candidates {
        if candidate.canonical_map().basis_rows(FINGERPRINT_ROWS)? == rows {
            return Ok(Some(candidate));
        }
    }
    Ok(None)
}

/// Recognise `map` as a published family for every index: the first rows
/// agree, and the all-index equality certificate of the canonical record
/// (with the map's first index) and the map's record is computed.  `None`
/// when no published family has the same first rows.
///
/// # Errors
///
/// [`RecurrenceError`] when the rows cannot be generated or a record is
/// malformed.
pub(super) fn recognize_published_family(
    map: &BasisMap,
) -> Result<Option<(PublishedFamily, RecurrenceEqualityCertificate)>, RecurrenceError> {
    let Some(family) = fingerprint_candidate(map)? else {
        return Ok(None);
    };
    let certificate = certify_recurrence_equality(
        &family.canonical_record(map.first_index()),
        &map.recurrence_json(),
    )?;
    Ok(Some((family, certificate)))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn narayana_row(m: usize, k: usize) -> Vec<Q> {
        fn binomial(n: usize, k: usize) -> Q {
            let mut value = Q::one();
            for i in 0..k {
                value = value * q((n - i) as i64) / q((i + 1) as i64);
            }
            value
        }
        (0..=k)
            .map(|j| binomial(k, j) * binomial(k + m, j) / binomial(m + j, j))
            .collect()
    }

    #[test]
    fn canonical_narayana_rows_are_the_closed_form() {
        for m in 0..4 {
            let rows = BasisMap::narayana_transformation(m).basis_rows(9).unwrap();
            for (k, row) in rows.iter().enumerate() {
                assert_eq!(row, &narayana_row(m, k), "m = {m}, k = {k}");
            }
        }
    }

    #[test]
    fn canonical_eulerian_and_derangement_rows() {
        let e = BasisMap::eulerian_transformation().basis_rows(5).unwrap();
        assert_eq!(e[0], vec![q(1)]);
        assert_eq!(e[1], vec![q(0), q(1)]);
        assert_eq!(e[3], vec![q(0), q(1), q(4), q(1)]);
        assert_eq!(e[4], vec![q(0), q(1), q(11), q(11), q(1)]);
        let d = BasisMap::derangement_transformation()
            .basis_rows(6)
            .unwrap();
        assert_eq!(d[1], vec![q(0)]);
        assert_eq!(d[2], vec![q(0), q(1)]);
        assert_eq!(d[3], vec![q(0), q(1), q(1)]);
        assert_eq!(d[4], vec![q(0), q(1), q(7), q(1)]);
        // d_5 = x + 21 x^2 + 21 x^3 + x^4 (44 derangements of 5).
        assert_eq!(d[5], vec![q(0), q(1), q(21), q(21), q(1)]);
    }

    #[test]
    fn recognition_is_certified_for_all_indices() {
        for family in [
            PublishedFamily::Narayana(0),
            PublishedFamily::Narayana(3),
            PublishedFamily::EulerianTransformation,
            PublishedFamily::DerangementTransformation,
        ] {
            let (found, certificate) = recognize_published_family(&family.canonical_map())
                .unwrap()
                .unwrap();
            assert_eq!(found, family);
            assert!(certificate.is_certified(), "{family:?}");
        }
        // The binomial map is none of them.
        let binomial = BasisMap::from_oeis("A007318").unwrap();
        assert!(recognize_published_family(&binomial).unwrap().is_none());
    }
}
