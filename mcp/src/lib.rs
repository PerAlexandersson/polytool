//! Model Context Protocol server exposing the exact polytool library.
//!
//! [`PolynomialToolsServer`] implements the MCP tools; the request and
//! response types below are their JSON schemas.  Exact integers and rationals
//! are serialized as strings.

mod family;
mod limits;
mod parsing;
mod recurrence_options;
mod report_items;
mod requests;
mod responses;
mod server;
#[cfg(test)]
mod tests;

pub use requests::{
    BicyclicSievingRequest, BigIntCoefficientInput, BigIntPolynomialBatchInput,
    BigIntPolynomialInput, BigIntPolynomialPairRequest, CheckPolynomialFamilyRequest,
    CyclicSievingFixedCountsInput, CyclicSievingRequest, CyclicSievingSequenceRequest,
    EhrhartHstarMode, EhrhartHstarRequest, FamilyCheckOptions, FindRecurrenceRequest,
    GenerateOeisRowsRequest, GenerateRecurrenceRowsRequest, GenerateSequenceRequest,
    GetOeisSequenceRequest, HStarInequalitiesRequest, LaceCheckRequest, ListOeisSequencesRequest,
    PolynomialBatchInput, PolynomialInput, RationalCoefficientInput, RecurrenceSearchOptionsInput,
    RecurrenceSpacesToolRequest, SequenceKind, VerifyRecurrenceSpacesToolRequest,
};
pub use responses::{
    BicyclicSievingResponse, BigIntNormalizedPolynomial, BigIntParsePolynomialItem,
    CheckPolynomialFamilyResponse, CoefficientCriterionItem, CoefficientInequalityCheckItem,
    CoefficientTestsItem, CoefficientTestsResponse, CyclicSievingPowerCheckItem,
    CyclicSievingReportItem, CyclicSievingResponse, CyclicSievingSequenceItemResponse,
    CyclicSievingSequenceResponse, DecompositionItem, DecompositionReport, DecompositionResponse,
    DiscriminantItem, DiscriminantResponse, DisplayedPolynomial, EhrhartHstarResponse,
    FamilyPolynomialReport, FindRecurrenceResponse, GenerateOeisRowsResponse,
    GenerateRecurrenceRowsResponse, GenerateSequenceResponse, GetOeisSequenceResponse,
    HStarInequalitiesItem, HStarInequalitiesResponse, HStarInequalityCheckItem,
    InterlacingPairItem, InterlacingPairResult, InterlacingProfileItem, InterlacingProfileResponse,
    InterlacingSequenceResponse, LaceCheckResponse, ListOeisSequencesResponse, MagicBasisReport,
    NormalizedPolynomial, OeisPolynomialRow, OeisSequenceSummary, OrbitDecompositionItem,
    OrbitSizeEntryItem, OrbitSizesCheckItem, ParsePolynomialItem, ParsePolynomialsResponse,
    PolynomialPropertiesItem, PolynomialPropertiesResponse, RealRootsItem, RealRootsResponse,
    RecurrenceSearchStatus, ResidueClassItem, ResultantResponse, RootOfUnityEvaluationItem,
};
pub use server::PolynomialToolsServer;
