//! The route certificate of a [`RealRootednessTheoremReport`]: at most one
//! [`RouteCertificate`] per report, of the report's own theorem, replayed by
//! [`RealRootednessTheoremReport::verify`] through the certificate's own
//! module.

use super::*;

/// The certificate of the theorem routes whose conclusion rests on more
/// than the hypotheses of the report, replayed by
/// [`RealRootednessTheoremReport::verify`].  Each variant is named after its
/// theorem; in JSON the certificate is an object with that theorem's name as
/// its single key, for example `{"factor_gauge": {...}}`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum RouteCertificate {
    /// [`check_iterated_preserver_theorem`]: the recognized step, the
    /// induction start, the degree step and the region (JSON key
    /// `iterated_preserver`).
    IteratedPreserver(IteratedPreserverCertificate),
    /// [`check_rational_pf_riordan_theorem`]: exact generating functions and PF evidence.
    RationalPfRiordan(Box<RationalPfRiordanCertificate>),
    /// [`check_variable_delay_theorem`]: the delay and the seed rows as exact
    /// steps of the window (JSON key `variable_delay`).
    VariableDelay(VariableDelayCertificate),
    /// [`check_residue_class_subsequence_theorem`]: the period and the class
    /// families with their theorem reports (JSON key `residue_class_subsequences`).
    ResidueClassSubsequences(ResidueClassCertificate),
    /// [`check_index_section_lift_theorem`]: the period, the transform, the
    /// parent family and its theorem reports (JSON key `index_section_lift`).
    IndexSectionLift(IndexSectionCertificate),
    /// [`check_factor_gauge_theorem`]: the gauge factor, the scalar
    /// multipliers, the common factor, the residual window and its
    /// coefficient polynomials (JSON key `factor_gauge`).
    FactorGauge(FactorGaugeCertificate),
    /// [`check_lag_two_euler_theorem`]: normalized operators and exact seed steps.
    LagTwoEuler(LagTwoEulerCertificate),
    /// [`check_binomial_perturbed_euler_orbit_theorem`]: parameters, exact
    /// coordinate transport and canonical seeds.
    BinomialPerturbedEulerOrbit(BinomialPerturbedEulerOrbitCertificate),
}

impl RouteCertificate {
    /// The theorem whose reports carry this certificate.
    pub fn theorem(&self) -> RealRootednessTheorem {
        match self {
            RouteCertificate::RationalPfRiordan(_) => RealRootednessTheorem::RationalPfRiordan,
            RouteCertificate::IteratedPreserver(_) => RealRootednessTheorem::IteratedPreserver,
            RouteCertificate::VariableDelay(_) => RealRootednessTheorem::VariableDelay,
            RouteCertificate::LagTwoEuler(_) => RealRootednessTheorem::LagTwoEuler,
            RouteCertificate::BinomialPerturbedEulerOrbit(_) => {
                RealRootednessTheorem::BinomialPerturbedEulerOrbit
            }
            RouteCertificate::ResidueClassSubsequences(_) => {
                RealRootednessTheorem::ResidueClassSubsequences
            }
            RouteCertificate::IndexSectionLift(_) => RealRootednessTheorem::IndexSectionLift,
            RouteCertificate::FactorGauge(_) => RealRootednessTheorem::FactorGauge,
        }
    }

    /// The start-delay bound that a check keeping its own start passed to
    /// its sub-checks (the residue classes, the parent family of an index
    /// section), checked against the largest allowed bound; `None` for the
    /// other certificates.
    pub(super) fn sub_check_start_delay_bound(&self) -> Result<Option<usize>, RecurrenceError> {
        match self {
            RouteCertificate::ResidueClassSubsequences(certificate) => {
                if certificate.class_maximum_start_delay > residue_class::MAXIMUM_CLASS_START_DELAY
                {
                    return Err(RecurrenceError::VerificationFailed(format!(
                        "the class start-delay bound exceeds {}",
                        residue_class::MAXIMUM_CLASS_START_DELAY
                    )));
                }
                Ok(Some(certificate.class_maximum_start_delay))
            }
            RouteCertificate::IndexSectionLift(certificate) => {
                if certificate.parent_maximum_start_delay
                    > index_section::MAXIMUM_PARENT_START_DELAY
                {
                    return Err(RecurrenceError::VerificationFailed(format!(
                        "the parent start-delay bound exceeds {}",
                        index_section::MAXIMUM_PARENT_START_DELAY
                    )));
                }
                Ok(Some(certificate.parent_maximum_start_delay))
            }
            RouteCertificate::IteratedPreserver(_)
            | RouteCertificate::LagTwoEuler(_)
            | RouteCertificate::BinomialPerturbedEulerOrbit(_)
            | RouteCertificate::VariableDelay(_)
            | RouteCertificate::RationalPfRiordan(_)
            | RouteCertificate::FactorGauge(_) => Ok(None),
        }
    }

    /// Replay the certificate on `report`, a report of its own theorem that
    /// equals its recomputation.
    pub(super) fn verify_in_report(
        &self,
        report: &RealRootednessTheoremReport,
    ) -> Result<(), RecurrenceError> {
        match self {
            RouteCertificate::IteratedPreserver(certificate) => {
                iterated_preserver::verify_certificate_in_report(certificate, report)
            }
            RouteCertificate::VariableDelay(certificate) => {
                variable_delay::verify_certificate_in_report(certificate, report)
            }
            RouteCertificate::ResidueClassSubsequences(certificate) => {
                residue_class::verify_certificate_in_report(certificate, report)
            }
            RouteCertificate::IndexSectionLift(certificate) => {
                index_section::verify_certificate_in_report(certificate, report)
            }
            RouteCertificate::LagTwoEuler(certificate) => {
                lag_two_euler::verify_certificate_in_report(certificate, report)
            }
            RouteCertificate::BinomialPerturbedEulerOrbit(certificate) => {
                binomial_perturbed_euler_orbit::verify_certificate_in_report(certificate, report)
            }
            RouteCertificate::FactorGauge(certificate) => {
                factor_gauge::verify_certificate_in_report(certificate, report)
            }
            RouteCertificate::RationalPfRiordan(certificate) => {
                rational_pf_riordan::verify_certificate_in_report(certificate, report)
            }
        }
    }
}

impl RealRootednessTheoremReport {
    /// The rational PF Riordan certificate, if recorded.
    pub fn rational_pf_riordan_certificate(&self) -> Option<&RationalPfRiordanCertificate> {
        match &self.certificate {
            Some(RouteCertificate::RationalPfRiordan(certificate)) => Some(certificate),
            _ => None,
        }
    }

    /// The binomial-perturbed Euler-orbit certificate, if recorded.
    pub fn binomial_perturbed_euler_orbit_certificate(
        &self,
    ) -> Option<&BinomialPerturbedEulerOrbitCertificate> {
        match &self.certificate {
            Some(RouteCertificate::BinomialPerturbedEulerOrbit(certificate)) => Some(certificate),
            _ => None,
        }
    }

    /// The positive lag-two Euler route certificate, if recorded.
    pub fn lag_two_euler_certificate(&self) -> Option<&LagTwoEulerCertificate> {
        match &self.certificate {
            Some(RouteCertificate::LagTwoEuler(certificate)) => Some(certificate),
            _ => None,
        }
    }

    /// The route certificate of [`check_iterated_preserver_theorem`], if
    /// recorded.
    pub fn iterated_preserver_certificate(&self) -> Option<&IteratedPreserverCertificate> {
        match &self.certificate {
            Some(RouteCertificate::IteratedPreserver(certificate)) => Some(certificate),
            _ => None,
        }
    }

    /// The route certificate of [`check_variable_delay_theorem`], if
    /// recorded.
    pub fn variable_delay_certificate(&self) -> Option<&VariableDelayCertificate> {
        match &self.certificate {
            Some(RouteCertificate::VariableDelay(certificate)) => Some(certificate),
            _ => None,
        }
    }

    /// The route certificate of [`check_residue_class_subsequence_theorem`],
    /// if recorded.
    pub fn residue_class_certificate(&self) -> Option<&ResidueClassCertificate> {
        match &self.certificate {
            Some(RouteCertificate::ResidueClassSubsequences(certificate)) => Some(certificate),
            _ => None,
        }
    }

    /// The route certificate of [`check_index_section_lift_theorem`], if
    /// recorded.
    pub fn index_section_certificate(&self) -> Option<&IndexSectionCertificate> {
        match &self.certificate {
            Some(RouteCertificate::IndexSectionLift(certificate)) => Some(certificate),
            _ => None,
        }
    }

    /// The route certificate of [`check_factor_gauge_theorem`], if recorded.
    pub fn factor_gauge_certificate(&self) -> Option<&FactorGaugeCertificate> {
        match &self.certificate {
            Some(RouteCertificate::FactorGauge(certificate)) => Some(certificate),
            _ => None,
        }
    }
}
