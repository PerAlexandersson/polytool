//! Real-rootedness proofs of named families from `polytool::sequences`.
//!
//! For each family, the test builds a recurrence with seed rows, checks that
//! the generator's first rows equal the rows the recurrence produces
//! (exactly, over the rationals), and runs the proving route that settles
//! the family, asserting that it reports `proved` with the expected theorem:
//!
//! - [`check_three_term_recurrence_theorem`] (orthogonal-polynomial
//!   recurrences, also in reciprocal coordinates as for the Narayana rows);
//! - [`check_liu_wang_theorem`] (first-order derivative recurrences);
//! - [`transform_recurrence`] with a root-domain request (a completed Rolle
//!   induction transported back through the `f`-polynomial map).
//!
//! What a passing test establishes depends on where the recurrence comes
//! from, which each test states:
//!
//! - **defining recurrence**: the generator is computed by exactly this
//!   recurrence (Eulerian, type-B Eulerian, Chebyshev, Hermite, 1/k-Eulerian),
//!   so the theorem applies to every row of the generator; the row comparison
//!   only guards against the two drifting apart;
//! - **OEIS catalog record** or **documented / fitted recurrence**: the
//!   generator is a closed form or a different algorithm, and its identity
//!   with the recurrence is checked on finitely many rows only.  The proof
//!   is then a proof for the recurrence family, which agrees with the
//!   generator on the checked rows.
//!
//! Where the induction cannot start at the first generator row (the Eulerian
//! row `A_0 = A_1 = 1` does not raise the degree), the theorem checks move
//! their induction start themselves and record it in `start_delay`; the
//! skipped rows enter as seed rows, checked exactly.  The root-domain tests
//! still start the transform one row later by hand.
//!
//! The Baxter triangle is proved through a Hadamard-factor lift
//! ([`find_hadamard_factor_lifts`]): the identity with its parent is
//! certified exactly and the parent is proved by Liu--Wang.
//!
//! Families whose only proof is the crate-internal periodic-product parent
//! route (for example the Hoggatt rows with `m = 1, q = 1`, which are
//! `(1 + x)^(n-1)`) are not listed here.

use num_bigint::BigInt;
use polytool::oeis;
use polytool::recurrence::{
    check_liu_wang_theorem, check_three_term_recurrence_theorem, find_hadamard_factor_lifts,
    transform_recurrence, verify_recurrence_transform, BigRational, BivarPoly, DegreeSchedule,
    HadamardFactorIdentityStatus, HadamardFactorLiftOptions, LiftVerdict,
    RealRootednessTheoremOptions, RealRootednessTheoremReport, Recurrence, RecurrenceJson,
    RecurrenceTerm, RecurrenceTransformRequest, RecurrenceTransformation, RootDomainRequest,
    RootEvidence, RootProperty, RootSubject, VeroneseParentProofRoute,
};
use polytool::sequences::{self, literature};

type Q = BigRational;

/// `(i, j, c)`: the monomial `c n^i x^j`.
type Monomial = (usize, usize, Q);

fn q(value: i64) -> Q {
    Q::from_integer(BigInt::from(value))
}

fn rational_rows(rows: Vec<Vec<BigInt>>) -> Vec<Vec<Q>> {
    rows.into_iter()
        .map(|row| row.into_iter().map(Q::from_integer).collect())
        .collect()
}

/// Trailing zeros removed (the zero polynomial is empty).
fn trimmed(row: &[Q]) -> Vec<Q> {
    let mut row = row.to_vec();
    while row.last().is_some_and(|c| *c == q(0)) {
        row.pop();
    }
    row
}

/// `sum c n^i x^j` from `(i, j, c)`.
fn bivariate(monomials: &[Monomial]) -> BivarPoly {
    let n_size = monomials.iter().map(|m| m.0 + 1).max().unwrap_or(1);
    let x_size = monomials.iter().map(|m| m.1 + 1).max().unwrap_or(1);
    let mut coeffs = vec![vec![q(0); x_size]; n_size];
    for (i, j, c) in monomials {
        coeffs[*i][*j] += c.clone();
    }
    BivarPoly::from_coefficients_normalized(coeffs)
}

/// `f(n) P_n = sum c_k(n, x) D^d P_(n-k)` from `(k, d, monomials of c_k)`.
fn rule(terms: &[(usize, usize, &[Monomial])], denominator: Option<&[Monomial]>) -> Recurrence {
    Recurrence {
        terms: terms
            .iter()
            .map(|(offset, order, monomials)| {
                RecurrenceTerm::new(*offset, *order, bivariate(monomials))
            })
            .collect(),
        denominator: denominator.map(bivariate),
        inhomogeneous: None,
    }
}

/// `P_n = (1 + (c (n - 1) + e) x) P_(n-1) + c x (1 - x) P'_(n-1)`.
fn derivative_rule(c: Q, e: Q) -> Recurrence {
    rule(
        &[
            (
                1,
                0,
                &[(0, 0, q(1)), (1, 1, c.clone()), (0, 1, e - c.clone())],
            ),
            (1, 1, &[(0, 1, c.clone()), (0, 2, -c)]),
        ],
        None,
    )
}

/// A family `P_(first_index), P_(first_index + 1), ...` given by a rule and
/// seed rows, with generator rows `generator[i]` = `P_(first_index + i)`.
struct Family {
    rule: Recurrence,
    first_index: usize,
    seeds: Vec<Vec<Q>>,
}

impl Family {
    fn new(rule: Recurrence, first_index: usize, generator: &[Vec<Q>], seed_count: usize) -> Self {
        let family = Self {
            rule,
            first_index,
            seeds: generator[..seed_count].to_vec(),
        };
        family.assert_matches(generator);
        family
    }

    /// The recurrence rows equal the generator rows.
    fn assert_matches(&self, generator: &[Vec<Q>]) {
        let rows = self
            .rule
            .generate_rows_rational(&self.seeds, self.first_index, generator.len())
            .expect("the recurrence replays");
        for (i, (left, right)) in rows.iter().zip(generator).enumerate() {
            assert_eq!(
                trimmed(left),
                trimmed(right),
                "recurrence and generator differ at index {}",
                self.first_index + i
            );
        }
    }

    fn three_term(&self) -> RealRootednessTheoremReport {
        let report = check_three_term_recurrence_theorem(
            &self.rule,
            self.first_index,
            &self.seeds,
            &RealRootednessTheoremOptions::default(),
        )
        .expect("the check runs");
        assert_proved(report)
    }

    fn liu_wang(&self) -> RealRootednessTheoremReport {
        let report = check_liu_wang_theorem(
            &self.rule,
            self.first_index,
            &self.seeds,
            &RealRootednessTheoremOptions::default(),
        )
        .expect("the check runs");
        let report = assert_proved(report);
        let region = report.conclusion.as_ref().unwrap().root_region.as_ref();
        assert_eq!(
            region.map(ToString::to_string).as_deref(),
            Some("(-inf, 0)")
        );
        report
    }

    /// The root-domain route through the `f`-polynomial map with the degree
    /// schedule `d_n = first_degree + (n - first_index)`: a completed Rolle
    /// induction certifies every source row real-rooted in `(-inf, 0]`.
    fn root_domain_f_polynomial(&self, first_degree: usize) {
        let source = RecurrenceJson::from_recurrence_rational(
            &self.rule,
            self.first_index,
            &self.seeds,
            None,
        );
        let transformation = RecurrenceTransformation::FPolynomial {
            schedule: DegreeSchedule::new(self.first_index, first_degree, 1),
        };
        let mut request = RecurrenceTransformRequest::from_source(transformation, source, 8);
        request.root_domain = Some(RootDomainRequest::default());
        let report = transform_recurrence(&request).expect("the transform runs");
        verify_recurrence_transform(&report).expect("the report replays");
        let claims = &report
            .root_domain
            .as_ref()
            .expect("root-domain report")
            .claims;
        let claim = claims
            .iter()
            .find(|claim| {
                claim.subject == RootSubject::Source
                    && claim.component == 0
                    && claim.property == RootProperty::RealRootedIn
                    && matches!(
                        claim.evidence,
                        RootEvidence::CompletedRecurrenceInduction
                            | RootEvidence::CertifiedFamilyInvariant
                    )
                    && claim.conditional_on.is_empty()
                    && claim.last_index.is_none()
            })
            .unwrap_or_else(|| panic!("no unconditional source claim: {claims:#?}"));
        assert!(claim.first_index <= self.first_index);
        let interval = claim.interval.as_ref().expect("an interval");
        assert!(interval.is_index_independent());
        assert_eq!(interval.to_string(), "(-inf, 0]");
    }
}

fn assert_proved(report: RealRootednessTheoremReport) -> RealRootednessTheoremReport {
    assert!(report.is_proved(), "{report:#?}");
    report.verify().expect("the report replays");
    report
}

const ROWS: usize = 30;

// ---------------------------------------------------------------------------
// Defining recurrences: the proof applies to the generator itself.
// ---------------------------------------------------------------------------

#[test]
fn chebyshev_first_kind_three_term() {
    // T_n = 2x T_(n-1) - T_(n-2), the loop of `chebyshev_polynomials_t_bigint`.
    let rows = rational_rows(sequences::chebyshev_polynomials_t_bigint(ROWS - 1));
    let rule = rule(&[(1, 0, &[(0, 1, q(2))]), (2, 0, &[(0, 0, q(-1))])], None);
    Family::new(rule, 0, &rows, 2).three_term();
}

#[test]
fn chebyshev_second_kind_three_term() {
    let rows = rational_rows(sequences::chebyshev_polynomials_u_bigint(ROWS - 1));
    let rule = rule(&[(1, 0, &[(0, 1, q(2))]), (2, 0, &[(0, 0, q(-1))])], None);
    Family::new(rule, 0, &rows, 2).three_term();
}

#[test]
fn hermite_three_term() {
    // He_n = x He_(n-1) - (n - 1) He_(n-2).
    let rows = rational_rows(sequences::hermite_polynomials_bigint(ROWS - 1));
    let rule = rule(
        &[
            (1, 0, &[(0, 1, q(1))]),
            (2, 0, &[(1, 0, q(-1)), (0, 0, q(1))]),
        ],
        None,
    );
    Family::new(rule, 0, &rows, 2).three_term();
}

#[test]
fn type_b_eulerian_liu_wang() {
    // B_n = (1 + (2n - 1) x) B_(n-1) + 2x(1 - x) B'_(n-1).
    let rows = rational_rows(sequences::type_b_eulerian_polynomials_bigint(ROWS - 1));
    Family::new(derivative_rule(q(2), q(1)), 0, &rows, 1).liu_wang();
}

#[test]
fn one_over_k_eulerian_liu_wang() {
    // A_n = (1 + k (n - 1) x) A_(n-1) + k x (1 - x) A'_(n-1).  A_0 = A_1 = 1,
    // so the degree step fails at n = 1 and the check moves the induction
    // start one row later by itself (A_0 enters as a checked seed row).
    for k in [q(2), q(3), Q::new(1.into(), 2.into())] {
        let rows = literature::one_over_k_eulerian_polynomials_rational(ROWS - 1, &k);
        let report = Family::new(derivative_rule(k, q(0)), 0, &rows, 1).liu_wang();
        assert_eq!(report.start_delay, 1);
        assert_eq!(report.conclusion.as_ref().unwrap().first_index, 0);
    }
}

#[test]
fn eulerian_liu_wang_from_the_first_row() {
    // A_n = (1 + (n - 1) x) A_(n-1) + x (1 - x) A'_(n-1) from A_0 = 1: the
    // induction start moves to A_1, and A_0 = A_1 = 1 also interlace (exact),
    // so the interlacing claim covers every n >= 1.
    let rows = rational_rows(sequences::eulerian_polynomials_bigint(ROWS - 1));
    let report = Family::new(derivative_rule(q(1), q(0)), 0, &rows, 1).liu_wang();
    assert_eq!(report.start_delay, 1);
    let conclusion = report.conclusion.as_ref().unwrap();
    assert_eq!(conclusion.first_index, 0);
    assert_eq!(conclusion.interlacing_from, Some(1));
}

#[test]
fn eulerian_root_domain() {
    // A_n = (1 + (n - 1) x) A_(n-1) + x (1 - x) A'_(n-1).  A_0 = A_1 = 1, so
    // the induction starts at A_1 (degree 0); A_0 = 1 needs no check.
    let rows = rational_rows(sequences::eulerian_polynomials_bigint(ROWS - 1));
    Family::new(derivative_rule(q(1), q(0)), 0, &rows, 1);
    Family::new(derivative_rule(q(1), q(0)), 1, &rows[1..], 1).root_domain_f_polynomial(0);
}

// ---------------------------------------------------------------------------
// Catalog or documented recurrences: identity checked on finitely many rows.
// ---------------------------------------------------------------------------

/// The OEIS catalog family `id`, its seeds advanced by `advance` rows, so that
/// `generator[0]` is catalog recurrence row `advance`.
fn catalog_family(id: &str, advance: usize, generator: &[Vec<Q>]) -> Family {
    let entry = oeis::by_id(id).expect("catalog entry");
    assert!(entry.prefix_rows.is_empty());
    let (rule, seeds) = entry.recurrence_parts().expect("catalog recurrence");
    let first_index = entry.recurrence_first_index;
    let rows = rule
        .generate_rows_rational(&seeds, first_index, advance + seeds.len())
        .expect("the catalog recurrence replays");
    let family = Family {
        rule,
        first_index: first_index + advance,
        seeds: rows[advance..].to_vec(),
    };
    family.assert_matches(generator);
    family
}

/// `s`-Eulerian polynomials for `s = (s_of(1), ..., s_of(n))`, `n = 1..=count`.
fn s_eulerian_rows(count: usize, s_of: fn(usize) -> usize) -> Vec<Vec<Q>> {
    (1..=count)
        .map(|n| {
            let s: Vec<usize> = (1..=n).map(s_of).collect();
            sequences::s_eulerian_polynomial_bigint(&s)
                .into_iter()
                .map(Q::from_integer)
                .collect()
        })
        .collect()
}

#[test]
fn narayana_three_term_via_catalog() {
    // A001263 row m is N_(m+1); the generator is the closed form
    // C(n,k) C(n,k-1) / n, compared on ROWS rows.
    let rows = rational_rows(sequences::narayana_polynomials_bigint(ROWS));
    catalog_family("A001263", 0, &rows[1..]).three_term();
    // Hoggatt polynomials with m = 2, q = 1 are the Narayana rows.
    let hoggatt = rational_rows(literature::hoggatt_polynomials_bigint(ROWS, 2, 1));
    assert_eq!(hoggatt[1..], rows[1..]);
}

#[test]
fn shifted_legendre_three_term() {
    // The documented n P~_n = (2n - 1)(2x - 1) P~_(n-1) - (n - 1) P~_(n-2);
    // the generator is the closed form, compared on ROWS rows.
    let rows = rational_rows(sequences::shifted_legendre_polynomials_bigint(ROWS - 1));
    let rule = rule(
        &[
            (
                1,
                0,
                &[(0, 0, q(1)), (0, 1, q(-2)), (1, 0, q(-2)), (1, 1, q(4))],
            ),
            (2, 0, &[(0, 0, q(1)), (1, 0, q(-1))]),
        ],
        Some(&[(1, 0, q(1))]),
    );
    Family::new(rule, 0, &rows, 2).three_term();
}

#[test]
fn s_eulerian_families_via_catalog() {
    // Compared with the generator on 18 rows (n = 1..=18).
    // s = (2, 4, ..., 2n): the type-B Eulerian rows A060187 (catalog row 0 is
    // the empty s); Liu--Wang.
    let rows = s_eulerian_rows(18, |i| 2 * i);
    catalog_family("A060187", 1, &rows).liu_wang();
    // s = (1, 3, ..., 2n - 1): A185410 (the 1/2-Eulerian rows); Liu--Wang.
    let rows = s_eulerian_rows(18, |i| 2 * i - 1);
    catalog_family("A185410", 1, &rows).liu_wang();
    // s = (2, 3, ..., n + 1): the Eulerian rows A008292 from A_2 on; the
    // root-domain route (degree schedule d_n = n - 1 from catalog index 1).
    let rows = s_eulerian_rows(18, |i| i + 1);
    catalog_family("A008292", 1, &rows).root_domain_f_polynomial(1);
}

#[test]
fn s_eulerian_multiples_of_three_fitted() {
    // s = (3, 6, ..., 3n): fitted P_n = (1 + (3n - 1) x) P_(n-1)
    // + 3 x (1 - x) P'_(n-1), P_1 = 1 + 2x, compared on 18 rows; Liu--Wang.
    let rows = s_eulerian_rows(18, |i| 3 * i);
    Family::new(derivative_rule(q(3), q(2)), 1, &rows, 1).liu_wang();
}

#[test]
fn refined_eulerian_second_letter_fitted() {
    // rows[n][1] (permutations of [n] with first letter 2), n >= 3: fitted
    // R_n = (1 + (n - 2) x) R_(n-1) + x (1 - x) R'_(n-1) with R_3 = 2x,
    // compared on 25 rows (R_2 = x is checked by the comparison too);
    // root-domain route with d_n = n - 2.
    let rows: Vec<Vec<Q>> = sequences::refined_eulerian_polynomials_bigint(27)
        .into_iter()
        .skip(2)
        .map(|mut row| {
            row.swap_remove(1)
                .into_iter()
                .map(Q::from_integer)
                .collect()
        })
        .collect();
    // In the rule index m = n - 1: P_m = (1 + (m - 1) x) P_(m-1) + x (1 - x) P'_(m-1).
    Family::new(derivative_rule(q(1), q(0)), 1, &rows, 1);
    Family::new(derivative_rule(q(1), q(0)), 2, &rows[1..], 1).root_domain_f_polynomial(1);
}

#[test]
fn baxter_hadamard_lift_via_catalog() {
    // A056939 (Baxter triangle): T(n, k) = binom(n+2, k+1) p(n, k), the
    // quotient rows satisfying a first-order derivative recurrence proved by
    // Liu--Wang; the identity is certified with an invariant row annihilator.
    // The Hoggatt rows with m = 3, q = 1 (the MacMahon box product, a closed
    // form) are catalog row n - 1, compared on ROWS rows.
    let entry = oeis::by_id("A056939").expect("catalog entry");
    let (rule, seeds) = entry.recurrence_parts().expect("catalog recurrence");
    let first_index = entry.recurrence_first_index;
    let catalog_rows = rule
        .generate_rows_rational(&seeds, first_index, ROWS)
        .expect("the catalog recurrence replays");
    let hoggatt = rational_rows(literature::hoggatt_polynomials_bigint(ROWS, 3, 1));
    for (offset, row) in catalog_rows.iter().enumerate() {
        assert_eq!(trimmed(row), trimmed(&hoggatt[offset + 1]), "row {offset}");
    }
    let report = find_hadamard_factor_lifts(
        &rule,
        &seeds,
        first_index,
        &HadamardFactorLiftOptions::default(),
    )
    .expect("the search runs");
    report.verify().expect("the report replays");
    let lift = report.proved_lifts().next().expect("a proved lift");
    assert_eq!(lift.factor_text, "binom(n+2, k+1)");
    assert!(matches!(
        lift.identity,
        HadamardFactorIdentityStatus::Certified { .. }
    ));
    assert_eq!(
        lift.parent_proof.route,
        Some(VeroneseParentProofRoute::LiuWang)
    );
    assert_eq!(lift.lift_summary.verdict, LiftVerdict::Proved);
}

// ---------------------------------------------------------------------------
// Gamma lifts and library recurrences of fitted families.
// ---------------------------------------------------------------------------

mod gamma_and_library_recurrences {
    use super::{q, rational_rows, Family, Q, ROWS};
    use polytool::oeis;
    use polytool::recurrence::{
        certify_gamma_lift, find_symmetric_companion_lifts, GammaLiftCandidate, GammaLiftOptions,
        GammaLiftParent, LiftVerdict, SymmetricCompanionLiftOptions, VeroneseSectionGauge,
    };
    use polytool::sequences::literature;

    /// `W_n = 2^(n-1) Gamma(A_n)(z/4)` from `n = n0`, with the Eulerian parent
    /// at index `n + shift`; asserts the proved verdict and the root region.
    fn peak_gamma_lift(family: &Family, shift: i64) {
        let candidate = GammaLiftCandidate {
            parent: GammaLiftParent::eulerian(),
            index_shift: shift,
            scale: "1/4".into(),
            gauge: VeroneseSectionGauge::Ratio {
                numerator: vec!["2".into()],
                denominator: vec!["1".into()],
                initial_value: "1".into(),
            },
        };
        let certification = certify_gamma_lift(
            &family.rule,
            &family.seeds,
            family.first_index,
            &candidate,
            &GammaLiftOptions::default(),
        )
        .expect("the certification runs");
        assert!(certification.is_proved(), "{certification:#?}");
        certification.verify().expect("the certification replays");
        assert_eq!(
            certification
                .lift_summary
                .root_transfer
                .input_root_region
                .to_string(),
            "(-inf, 0]"
        );
    }

    #[test]
    fn peak_polynomials_gamma_lift_from_eulerian() {
        // Defining recurrence: `peak_polynomials_bigint` is computed by the
        // rule of `peak_polynomials_recurrence` from W_1 on (W_0 = 1 is
        // constant).  The same rows are the valley polynomials.
        let rows = rational_rows(literature::peak_polynomials_bigint(ROWS - 1));
        let (rule, seeds, first_index) = literature::peak_polynomials_recurrence();
        assert_eq!(seeds[..], rows[1..2]);
        let family = Family::new(rule, first_index, &rows[1..], 1);
        peak_gamma_lift(&family, 0);
    }

    #[test]
    fn peak_polynomials_gamma_lift_via_catalog() {
        // A008303 row m is W_(m+1); compared with the generator on ROWS rows.
        let rows = rational_rows(literature::peak_polynomials_bigint(ROWS));
        let entry = oeis::by_id("A008303").expect("catalog entry");
        let (rule, seeds) = entry.recurrence_parts().expect("catalog recurrence");
        let family = Family::new(rule, entry.recurrence_first_index, &rows[1..], seeds.len());
        assert_eq!(family.seeds, seeds);
        peak_gamma_lift(&family, 1);
    }

    #[test]
    fn one_over_k_eulerian_symmetric_parts_symmetric_companion() {
        // The library second-order recurrence of the symmetric parts a_n
        // (proved in its rustdoc; it could define the generator), compared
        // with the Stapledon decompositions on ROWS rows.  The
        // symmetric-companion lift proves a_n, n >= 2, real-rooted for
        // k = 2, 3 (a_1 = 1).  For k = 1/2 the companion is certified but the
        // sign hypotheses fail (k - 1 < 0), and the parts b_n are open.
        for (k, proved) in [
            (q(2), true),
            (q(3), true),
            (Q::new(1.into(), 2.into()), false),
        ] {
            let parts: Vec<Vec<Q>> =
                literature::one_over_k_eulerian_symmetric_parts_rational(ROWS, &k)
                    .into_iter()
                    .map(|(a, _)| a)
                    .collect();
            let (rule, seeds, first_index) =
                literature::one_over_k_eulerian_symmetric_part_a_recurrence(&k);
            assert_eq!(seeds[..], parts[1..3]);
            Family::new(rule.clone(), first_index, &parts[1..], 2);
            let family = Family::new(rule, 2, &parts[2..], 2);
            let report = find_symmetric_companion_lifts(
                &family.rule,
                &family.seeds,
                family.first_index,
                &SymmetricCompanionLiftOptions::default(),
            )
            .expect("the search runs");
            match report.proving_lift() {
                Some(lift) => {
                    assert!(proved, "k = {k}");
                    let summary = lift.lift_summary.as_ref().expect("a summary");
                    assert_eq!(summary.verdict, LiftVerdict::Proved);
                    assert_eq!(
                        summary.root_transfer.input_root_region.to_string(),
                        "(-inf, 0]"
                    );
                }
                None => {
                    assert!(!proved, "k = {k}: {report:#?}");
                    assert!(report.certified_lifts().count() > 0);
                }
            }
        }
    }
}
