//! The route certificate of a [`RealRootednessTheoremReport`]: at most one
//! [`RouteCertificate`] per report, of the report's own theorem, replayed by
//! [`RealRootednessTheoremReport::verify`] through the certificate's own
//! module.

use super::*;

// Keep the typed certificate, its theorem, its replay, its delay bound and
// its existing public accessor together. Adding a certificate needs one row.
macro_rules! route_certificates {
    ($( $(#[$doc:meta])* $variant:ident($payload:ty) => $module:ident,
        $accessor:ident: $view:ty, $bound:expr; )*) => {
        /// A route certificate, tagged in JSON by its theorem's snake-case name.
        /// Each report carries at most one; replay checks that it is its own theorem.
        #[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
        #[serde(rename_all = "snake_case")]
        pub enum RouteCertificate {
            $($(#[$doc])* $variant($payload),)*
        }

        impl RouteCertificate {
            /// The theorem whose reports carry this certificate.
            pub fn theorem(&self) -> RealRootednessTheorem {
                match self { $(Self::$variant(_) => RealRootednessTheorem::$variant,)* }
            }

            /// The bound passed to sub-checks, checked against its route's limit.
            pub(super) fn sub_check_start_delay_bound(&self) -> Result<Option<usize>, RecurrenceError> {
                match self { $(Self::$variant(certificate) => ($bound)(certificate),)* }
            }

            /// Replay a certificate on a recomputed report of its own theorem.
            pub(super) fn verify_in_report(&self, report: &RealRootednessTheoremReport)
                -> Result<(), RecurrenceError> {
                match self {
                    $(Self::$variant(certificate) =>
                        $module::verify_certificate_in_report(certificate, report),)*
                }
            }
        }

        impl RealRootednessTheoremReport {
            $(#[doc = concat!("The [`", stringify!($view), "`] of this route, when recorded.")]
            pub fn $accessor(&self) -> Option<&$view> {
                match &self.certificate {
                    Some(RouteCertificate::$variant(certificate)) => Some(certificate),
                    _ => None,
                }
            })*
        }
    };
}

route_certificates! {
    /// [`check_iterated_preserver_theorem`]: the recognized step, the
    /// induction start, the degree step and the region (JSON key
    /// `iterated_preserver`).
    IteratedPreserver(IteratedPreserverCertificate) => iterated_preserver,
        iterated_preserver_certificate: IteratedPreserverCertificate,
    |_: &IteratedPreserverCertificate| Ok(None);
    /// [`check_rational_pf_riordan_theorem`]: exact generating functions and PF evidence.
    RationalPfRiordan(Box<RationalPfRiordanCertificate>) => rational_pf_riordan,
        rational_pf_riordan_certificate: RationalPfRiordanCertificate,
    |_: &RationalPfRiordanCertificate| Ok(None);
    /// Shared reversal transform and a proved theorem on the transformed family.
    DegreeLawReversal(Box<DegreeLawReversalCertificate>) => degree_law_reversal,
        degree_law_reversal_certificate: DegreeLawReversalCertificate,
    |c: &DegreeLawReversalCertificate| {
        if c.maximum_start_delay > index_section::MAXIMUM_PARENT_START_DELAY {
            return Err(RecurrenceError::VerificationFailed(
                "reversal start-delay bound exceeds 16".into(),
            ));
        }
        Ok(Some(c.maximum_start_delay))
    };
    /// [`check_variable_delay_theorem`]: the delay and the seed rows as exact
    /// steps of the window (JSON key `variable_delay`).
    VariableDelay(VariableDelayCertificate) => variable_delay,
        variable_delay_certificate: VariableDelayCertificate,
    |_: &VariableDelayCertificate| Ok(None);
    /// [`check_residue_class_subsequence_theorem`]: the period and the class
    /// families with their theorem reports (JSON key `residue_class_subsequences`).
    ResidueClassSubsequences(ResidueClassCertificate) => residue_class,
        residue_class_certificate: ResidueClassCertificate,
    |c: &ResidueClassCertificate| {
        if c.class_maximum_start_delay > residue_class::MAXIMUM_CLASS_START_DELAY {
            return Err(RecurrenceError::VerificationFailed(format!(
                "the class start-delay bound exceeds {}",
                residue_class::MAXIMUM_CLASS_START_DELAY
            )));
        }
        Ok(Some(c.class_maximum_start_delay))
    };
    /// [`check_index_section_lift_theorem`]: the period, the transform, the
    /// parent family and its theorem reports (JSON key `index_section_lift`).
    IndexSectionLift(IndexSectionCertificate) => index_section,
        index_section_certificate: IndexSectionCertificate,
    |c: &IndexSectionCertificate| {
        if c.parent_maximum_start_delay > index_section::MAXIMUM_PARENT_START_DELAY {
            return Err(RecurrenceError::VerificationFailed(format!(
                "the parent start-delay bound exceeds {}",
                index_section::MAXIMUM_PARENT_START_DELAY
            )));
        }
        Ok(Some(c.parent_maximum_start_delay))
    };
    /// [`check_factor_gauge_theorem`]: the gauge factor, the scalar
    /// multipliers, the common factor, the residual window and its
    /// coefficient polynomials (JSON key `factor_gauge`).
    FactorGauge(FactorGaugeCertificate) => factor_gauge,
        factor_gauge_certificate: FactorGaugeCertificate,
    |_: &FactorGaugeCertificate| Ok(None);
    /// [`check_lag_two_euler_theorem`]: normalized operators and exact seed steps.
    LagTwoEuler(LagTwoEulerCertificate) => lag_two_euler,
        lag_two_euler_certificate: LagTwoEulerCertificate,
    |_: &LagTwoEulerCertificate| Ok(None);
    /// [`check_binomial_perturbed_euler_orbit_theorem`]: parameters, exact
    /// coordinate transport and canonical seeds.
    BinomialPerturbedEulerOrbit(BinomialPerturbedEulerOrbitCertificate) => binomial_perturbed_euler_orbit,
        binomial_perturbed_euler_orbit_certificate: BinomialPerturbedEulerOrbitCertificate,
    |_: &BinomialPerturbedEulerOrbitCertificate| Ok(None);
}
