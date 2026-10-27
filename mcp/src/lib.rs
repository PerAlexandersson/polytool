//! Model Context Protocol server exposing the exact polytool library.
//!
//! [`PolynomialToolsServer`] implements the MCP tools; the request and
//! response types below are their JSON schemas.  Exact integers and rationals
//! are serialized as strings.

mod enumeration;
mod family;
mod limits;
mod parsing;
mod planar_network;
mod recurrence_options;
mod report_items;
mod requests;
mod responses;
mod server;
#[cfg(test)]
mod tests;
#[cfg(test)]
mod tests_validation;
mod theorems;

pub use requests::EnumerateRealRootedRequest;
pub use requests::RecurrenceTheoremsRequest;
pub use requests::{
    AssessGeneratingDescriptionToolRequest, RecurrenceImplicitToolRequest,
    VerifyRecurrenceImplicitToolRequest,
};
pub use requests::{
    BicyclicSievingRequest, BigIntCoefficientInput, BigIntPolynomialBatchInput,
    BigIntPolynomialInput, BigIntPolynomialPairRequest, BrentiRecognitionRequest,
    BrentiWeightsInput, CheckPolynomialFamilyRequest, CyclicSievingFixedCountsInput,
    CyclicSievingRequest, CyclicSievingSequenceRequest, EhrhartHstarMode, EhrhartHstarRequest,
    FamilyCheckOptions, FindRecurrenceRequest, GenerateOeisRowsRequest,
    GenerateRecurrenceRowsRequest, GenerateSequenceRequest, GetOeisSequenceRequest,
    HStarInequalitiesRequest, LaceCheckRequest, ListOeisSequencesRequest, PolynomialBatchInput,
    PolynomialInput, RationalCoefficientInput, RecurrenceSearchOptionsInput,
    RecurrenceSpacesToolRequest, SequenceKind, TotallyNonnegativeNetworkRequest,
    VerifyRecurrenceSpacesToolRequest,
};
pub use requests::{LeftProductionRequest, PfNetworkRequest};
pub use requests::{
    RecurrenceRepresentationSearchToolRequest, VerifyRecurrenceRepresentationSearchToolRequest,
};
pub use requests::{
    RecurrenceTransformToolRequest, RootEnclosuresToolRequest, VerifyRecurrenceTransformToolRequest,
};
pub use responses::{
    BicyclicSievingResponse, BigIntNormalizedPolynomial, BigIntParsePolynomialItem,
    BrentiRecognitionOptionsItem, BrentiRecognitionResponse, ChainInterlacingItem,
    ChainPolynomialItem, ChainPolynomialsItem, CheckPolynomialFamilyResponse,
    CoefficientCriterionItem, CoefficientInequalityCheckItem, CoefficientTestsItem,
    CoefficientTestsResponse, CyclicSievingPowerCheckItem, CyclicSievingReportItem,
    CyclicSievingResponse, CyclicSievingSequenceItemResponse, CyclicSievingSequenceResponse,
    DecompositionItem, DecompositionReport, DecompositionResponse, DiscriminantItem,
    DiscriminantResponse, DisplayedPolynomial, EhrhartHstarResponse, FamilyPolynomialReport,
    FindRecurrenceResponse, GenerateOeisRowsResponse, GenerateRecurrenceRowsResponse,
    GenerateSequenceResponse, GetOeisSequenceResponse, HStarInequalitiesItem,
    HStarInequalitiesResponse, HStarInequalityCheckItem, InterlacingPairItem,
    InterlacingPairResult, InterlacingProfileItem, InterlacingProfileResponse,
    InterlacingSequenceResponse, LaceCheckResponse, ListOeisSequencesResponse, MagicBasisReport,
    NegativeMinorItem, NormalizedPolynomial, OeisPolynomialRow, OeisSequenceSummary,
    OrbitDecompositionItem, OrbitSizeEntryItem, OrbitSizesCheckItem, ParsePolynomialItem,
    ParsePolynomialsResponse, PolynomialPropertiesItem, PolynomialPropertiesResponse,
    RealRootsItem, RealRootsResponse, RecurrenceSearchStatus, ResidueClassItem, ResultantResponse,
    RootOfUnityEvaluationItem, StripDigraphEdgeItem, StripDigraphItem,
    TotallyNonnegativeFailureItem, TotallyNonnegativeNetworkResponse,
};
pub use responses::{
    EnumerateRealRootedLimits, EnumerateRealRootedResponse, EnumerateRealRootedStats,
};
pub use responses::{
    LeftProductionResponse, PfNetworkItem, PfNetworkOptionsItem, PfNetworkResponse,
    PfNetworkStepItem, ProductionNetworkItem, RationalNegativeMinorItem,
};
pub use server::PolynomialToolsServer;

/// Deprecated name of [`TotallyNonnegativeNetworkRequest`].
#[deprecated(
    since = "0.6.0",
    note = "renamed to `TotallyNonnegativeNetworkRequest`"
)]
pub type TnnNetworkRequest = TotallyNonnegativeNetworkRequest;

/// Deprecated name of [`TotallyNonnegativeNetworkResponse`].
#[deprecated(
    since = "0.6.0",
    note = "renamed to `TotallyNonnegativeNetworkResponse`"
)]
pub type TnnNetworkResponse = TotallyNonnegativeNetworkResponse;

/// Deprecated name of [`TotallyNonnegativeFailureItem`].
#[deprecated(since = "0.6.0", note = "renamed to `TotallyNonnegativeFailureItem`")]
pub type TnnFailureItem = TotallyNonnegativeFailureItem;
