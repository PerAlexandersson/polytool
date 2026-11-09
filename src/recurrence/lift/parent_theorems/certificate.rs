//! The route certificate of a [`RealRootednessTheoremReport`] and the JSON
//! form of the report.
//!
//! A report carries at most one [`RouteCertificate`], of its own theorem.
//! In JSON each certificate keeps the key it had as a report field
//! (`iterated_preserver`, `variable_delay`, `residue_classes`,
//! `index_section`, `factor_gauge`), so reports serialize exactly as before
//! the enum existed: stored certificates are identified by a hash of their
//! serialization.  The private wire structs below spell out that form field
//! for field, in order and with the same omission rules; [`ReportWire`]
//! rejects unknown keys and a report with two certificate keys.

use super::*;

/// The certificate of the theorem routes whose conclusion rests on more
/// than the hypotheses of the report, replayed by
/// [`RealRootednessTheoremReport::verify`].  In the JSON of a report it is
/// stored under the route's own key (the snake-case name of the variant's
/// field accessor, for example `iterated_preserver`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RouteCertificate {
    /// [`check_iterated_preserver_theorem`]: the recognized step, the
    /// induction start, the degree step and the region (JSON key
    /// `iterated_preserver`).
    IteratedPreserver(IteratedPreserverCertificate),
    /// [`check_variable_delay_theorem`]: the delay and the seed rows as exact
    /// steps of the window (JSON key `variable_delay`).
    VariableDelay(VariableDelayCertificate),
    /// [`check_residue_class_subsequence_theorem`]: the period and the class
    /// families with their theorem reports (JSON key `residue_classes`).
    ResidueClasses(ResidueClassCertificate),
    /// [`check_index_section_lift_theorem`]: the period, the transform, the
    /// parent family and its theorem reports (JSON key `index_section`).
    IndexSection(IndexSectionCertificate),
    /// [`check_factor_gauge_theorem`]: the gauge factor, the scalar
    /// multipliers, the common factor, the residual window and its
    /// coefficient polynomials (JSON key `factor_gauge`).
    FactorGauge(FactorGaugeCertificate),
}

impl RouteCertificate {
    /// The theorem whose reports carry this certificate.
    pub fn theorem(&self) -> RealRootednessTheorem {
        match self {
            RouteCertificate::IteratedPreserver(_) => RealRootednessTheorem::IteratedPreserver,
            RouteCertificate::VariableDelay(_) => RealRootednessTheorem::VariableDelay,
            RouteCertificate::ResidueClasses(_) => RealRootednessTheorem::ResidueClassSubsequences,
            RouteCertificate::IndexSection(_) => RealRootednessTheorem::IndexSectionLift,
            RouteCertificate::FactorGauge(_) => RealRootednessTheorem::FactorGauge,
        }
    }

    /// The start-delay bound that a check keeping its own start passed to
    /// its sub-checks (the residue classes, the parent family of an index
    /// section), checked against the largest allowed bound; `None` for the
    /// other certificates.
    pub(super) fn sub_check_start_delay_bound(&self) -> Result<Option<usize>, RecurrenceError> {
        match self {
            RouteCertificate::ResidueClasses(certificate) => {
                if certificate.class_maximum_start_delay > residue_class::MAXIMUM_CLASS_START_DELAY
                {
                    return Err(RecurrenceError::VerificationFailed(format!(
                        "the class start-delay bound exceeds {}",
                        residue_class::MAXIMUM_CLASS_START_DELAY
                    )));
                }
                Ok(Some(certificate.class_maximum_start_delay))
            }
            RouteCertificate::IndexSection(certificate) => {
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
            | RouteCertificate::VariableDelay(_)
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
            RouteCertificate::ResidueClasses(certificate) => {
                residue_class::verify_certificate_in_report(certificate, report)
            }
            RouteCertificate::IndexSection(certificate) => {
                index_section::verify_certificate_in_report(certificate, report)
            }
            RouteCertificate::FactorGauge(certificate) => {
                factor_gauge::verify_certificate_in_report(certificate, report)
            }
        }
    }
}

impl RealRootednessTheoremReport {
    /// The route certificate of [`check_iterated_preserver_theorem`], if
    /// recorded.
    pub fn iterated_preserver(&self) -> Option<&IteratedPreserverCertificate> {
        match &self.certificate {
            Some(RouteCertificate::IteratedPreserver(certificate)) => Some(certificate),
            _ => None,
        }
    }

    /// The route certificate of [`check_variable_delay_theorem`], if
    /// recorded.
    pub fn variable_delay(&self) -> Option<&VariableDelayCertificate> {
        match &self.certificate {
            Some(RouteCertificate::VariableDelay(certificate)) => Some(certificate),
            _ => None,
        }
    }

    /// The route certificate of [`check_residue_class_subsequence_theorem`],
    /// if recorded.
    pub fn residue_classes(&self) -> Option<&ResidueClassCertificate> {
        match &self.certificate {
            Some(RouteCertificate::ResidueClasses(certificate)) => Some(certificate),
            _ => None,
        }
    }

    /// The route certificate of [`check_index_section_lift_theorem`], if
    /// recorded.
    pub fn index_section(&self) -> Option<&IndexSectionCertificate> {
        match &self.certificate {
            Some(RouteCertificate::IndexSection(certificate)) => Some(certificate),
            _ => None,
        }
    }

    /// The route certificate of [`check_factor_gauge_theorem`], if recorded.
    pub fn factor_gauge(&self) -> Option<&FactorGaugeCertificate> {
        match &self.certificate {
            Some(RouteCertificate::FactorGauge(certificate)) => Some(certificate),
            _ => None,
        }
    }
}

// ---------------------------------------------------------------------------
// JSON form
// ---------------------------------------------------------------------------

/// The JSON form of a report, read: every field of the schema, unknown keys
/// rejected.
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct ReportWire {
    schema: String,
    theorem: RealRootednessTheorem,
    reference: String,
    recurrence: RecurrenceJson,
    #[serde(default)]
    start_delay: usize,
    #[serde(default)]
    requested_root_region: Option<RootInterval>,
    status: RealRootednessTheoremStatus,
    normal_form: String,
    hypotheses: Vec<LiftHypothesis>,
    #[serde(default)]
    conclusion: Option<RealRootednessTheoremConclusion>,
    #[serde(default)]
    trivially_real_rooted: bool,
    lean_replay: String,
    #[serde(default)]
    lean_route: Option<String>,
    #[serde(default)]
    lean_declarations: Vec<String>,
    #[serde(default)]
    left_hand_factor_normalization: Option<LeftHandFactorNormalization>,
    #[serde(default)]
    iterated_preserver: Option<IteratedPreserverCertificate>,
    #[serde(default)]
    variable_delay: Option<VariableDelayCertificate>,
    #[serde(default)]
    residue_classes: Option<ResidueClassCertificate>,
    #[serde(default)]
    index_section: Option<IndexSectionCertificate>,
    #[serde(default)]
    factor_gauge: Option<FactorGaugeCertificate>,
}

impl TryFrom<ReportWire> for RealRootednessTheoremReport {
    type Error = String;

    fn try_from(wire: ReportWire) -> Result<Self, Self::Error> {
        let certificates = [
            wire.iterated_preserver
                .map(RouteCertificate::IteratedPreserver),
            wire.variable_delay.map(RouteCertificate::VariableDelay),
            wire.residue_classes.map(RouteCertificate::ResidueClasses),
            wire.index_section.map(RouteCertificate::IndexSection),
            wire.factor_gauge.map(RouteCertificate::FactorGauge),
        ];
        let mut certificates = certificates.into_iter().flatten();
        let certificate = certificates.next();
        if certificates.next().is_some() {
            return Err("a theorem report has more than one route certificate".into());
        }
        Ok(Self {
            schema: wire.schema,
            theorem: wire.theorem,
            reference: wire.reference,
            recurrence: wire.recurrence,
            start_delay: wire.start_delay,
            requested_root_region: wire.requested_root_region,
            status: wire.status,
            normal_form: wire.normal_form,
            hypotheses: wire.hypotheses,
            conclusion: wire.conclusion,
            trivially_real_rooted: wire.trivially_real_rooted,
            lean_replay: wire.lean_replay,
            lean_route: wire.lean_route,
            lean_declarations: wire.lean_declarations,
            left_hand_factor_normalization: wire.left_hand_factor_normalization,
            certificate,
        })
    }
}

impl<'de> Deserialize<'de> for RealRootednessTheoremReport {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        ReportWire::deserialize(deserializer)?
            .try_into()
            .map_err(serde::de::Error::custom)
    }
}

/// The JSON form of a report, written: the fields of [`ReportWire`] in the
/// same order with the same omission rules, borrowed from the report.
#[derive(Serialize)]
struct ReportWireRef<'a> {
    schema: &'a str,
    theorem: RealRootednessTheorem,
    reference: &'a str,
    recurrence: &'a RecurrenceJson,
    #[serde(skip_serializing_if = "is_zero")]
    start_delay: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    requested_root_region: Option<&'a RootInterval>,
    status: RealRootednessTheoremStatus,
    normal_form: &'a str,
    hypotheses: &'a [LiftHypothesis],
    #[serde(skip_serializing_if = "Option::is_none")]
    conclusion: Option<&'a RealRootednessTheoremConclusion>,
    #[serde(skip_serializing_if = "std::ops::Not::not")]
    trivially_real_rooted: bool,
    lean_replay: &'a str,
    #[serde(skip_serializing_if = "Option::is_none")]
    lean_route: Option<&'a str>,
    #[serde(skip_serializing_if = "is_empty_list")]
    lean_declarations: &'a [String],
    #[serde(skip_serializing_if = "Option::is_none")]
    left_hand_factor_normalization: Option<&'a LeftHandFactorNormalization>,
    #[serde(skip_serializing_if = "Option::is_none")]
    iterated_preserver: Option<&'a IteratedPreserverCertificate>,
    #[serde(skip_serializing_if = "Option::is_none")]
    variable_delay: Option<&'a VariableDelayCertificate>,
    #[serde(skip_serializing_if = "Option::is_none")]
    residue_classes: Option<&'a ResidueClassCertificate>,
    #[serde(skip_serializing_if = "Option::is_none")]
    index_section: Option<&'a IndexSectionCertificate>,
    #[serde(skip_serializing_if = "Option::is_none")]
    factor_gauge: Option<&'a FactorGaugeCertificate>,
}

fn is_empty_list(list: &&[String]) -> bool {
    list.is_empty()
}

impl Serialize for RealRootednessTheoremReport {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        ReportWireRef {
            schema: &self.schema,
            theorem: self.theorem,
            reference: &self.reference,
            recurrence: &self.recurrence,
            start_delay: self.start_delay,
            requested_root_region: self.requested_root_region.as_ref(),
            status: self.status,
            normal_form: &self.normal_form,
            hypotheses: &self.hypotheses,
            conclusion: self.conclusion.as_ref(),
            trivially_real_rooted: self.trivially_real_rooted,
            lean_replay: &self.lean_replay,
            lean_route: self.lean_route.as_deref(),
            lean_declarations: &self.lean_declarations,
            left_hand_factor_normalization: self.left_hand_factor_normalization.as_ref(),
            iterated_preserver: self.iterated_preserver(),
            variable_delay: self.variable_delay(),
            residue_classes: self.residue_classes(),
            index_section: self.index_section(),
            factor_gauge: self.factor_gauge(),
        }
        .serialize(serializer)
    }
}
