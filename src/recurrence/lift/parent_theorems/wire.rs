//! Private, current-format Serde contract of theorem reports.
//!
//! Field order, defaults, skips and the diagnostic struct name mirror the
//! original derive. Parsing remains structural; theorem/certificate agreement
//! is still checked only by report verification. The remote derive builds
//! the public report directly and preserves Serde's original struct and
//! sequence-length diagnostics as well as its structural parsing rules.

use super::*;

#[derive(Deserialize)]
#[serde(
    deny_unknown_fields,
    rename = "RealRootednessTheoremReport",
    remote = "RealRootednessTheoremReport"
)]
struct ReportWire {
    schema: String,
    theorem: RealRootednessTheorem,
    reference: String,
    recurrence: RecurrenceJson,
    #[serde(default, skip_serializing_if = "is_zero")]
    start_delay: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    requested_root_region: Option<RootInterval>,
    status: RealRootednessTheoremStatus,
    normal_form: String,
    hypotheses: Vec<LiftHypothesis>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    conclusion: Option<RealRootednessTheoremConclusion>,
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    trivially_real_rooted: bool,
    lean_replay: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    lean_route: Option<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    lean_declarations: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    left_hand_factor_normalization: Option<LeftHandFactorNormalization>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    certificate: Option<RouteCertificate>,
}

#[derive(Serialize)]
#[serde(rename = "RealRootednessTheoremReport")]
struct ReportWireRef<'a> {
    schema: &'a String,
    theorem: RealRootednessTheorem,
    reference: &'a String,
    recurrence: &'a RecurrenceJson,
    #[serde(skip_serializing_if = "is_zero")]
    start_delay: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    requested_root_region: &'a Option<RootInterval>,
    status: RealRootednessTheoremStatus,
    normal_form: &'a String,
    hypotheses: &'a Vec<LiftHypothesis>,
    #[serde(skip_serializing_if = "Option::is_none")]
    conclusion: &'a Option<RealRootednessTheoremConclusion>,
    #[serde(skip_serializing_if = "std::ops::Not::not")]
    trivially_real_rooted: bool,
    lean_replay: &'a String,
    #[serde(skip_serializing_if = "Option::is_none")]
    lean_route: &'a Option<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    lean_declarations: &'a Vec<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    left_hand_factor_normalization: &'a Option<LeftHandFactorNormalization>,
    #[serde(skip_serializing_if = "Option::is_none")]
    certificate: &'a Option<RouteCertificate>,
}

impl<'a> From<&'a RealRootednessTheoremReport> for ReportWireRef<'a> {
    fn from(report: &'a RealRootednessTheoremReport) -> Self {
        Self {
            schema: &report.schema,
            theorem: report.theorem,
            reference: &report.reference,
            recurrence: &report.recurrence,
            start_delay: report.start_delay,
            requested_root_region: &report.requested_root_region,
            status: report.status,
            normal_form: &report.normal_form,
            hypotheses: &report.hypotheses,
            conclusion: &report.conclusion,
            trivially_real_rooted: report.trivially_real_rooted,
            lean_replay: &report.lean_replay,
            lean_route: &report.lean_route,
            lean_declarations: &report.lean_declarations,
            left_hand_factor_normalization: &report.left_hand_factor_normalization,
            certificate: &report.certificate,
        }
    }
}

impl Serialize for RealRootednessTheoremReport {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        ReportWireRef::from(self).serialize(serializer)
    }
}

impl<'de> Deserialize<'de> for RealRootednessTheoremReport {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        ReportWire::deserialize(deserializer)
    }
}
