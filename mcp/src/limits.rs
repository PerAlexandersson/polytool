//! Input-size limits enforced by the MCP tools before any exact computation.

pub(crate) const MCP_MAX_TEXT_BYTES: usize = 1_048_576;
pub(crate) const MCP_MAX_BATCH_POLYNOMIALS: usize = 512;
pub(crate) const MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL: usize = 4_096;
pub(crate) const MCP_MAX_TOTAL_COEFFICIENTS: usize = 65_536;
pub(crate) const MCP_MAX_COEFFICIENT_TEXT_BYTES: usize = 16_384;
pub(crate) const MCP_MAX_CYCLIC_ORDER: usize = 10_000;
pub(crate) const MCP_MAX_SEQUENCE_N: usize = 200;
pub(crate) const MCP_MAX_RECURRENCE_ROWS: usize = 200;
pub(crate) const MCP_MAX_RECURRENCE_INPUT_ROWS: usize = 256;
pub(crate) const MCP_MAX_RECURRENCE_DEGREE: usize = 16;
pub(crate) const MCP_MAX_RECURRENCE_LENGTH: usize = 32;
pub(crate) const MCP_MAX_RECURRENCE_CANDIDATES: usize = 50_000;
pub(crate) const MCP_MAX_RECURRENCE_UNKNOWNS: usize = 20_000;
pub(crate) const MCP_MAX_LACE_MATRIX_CELLS: usize = 16_384;
pub(crate) const MCP_MAX_LACE_MINORS: usize = 250_000;
pub(crate) const MCP_MAX_TNN_NETWORK_ROWS: usize = 128;
pub(crate) const MCP_MAX_BRENTI_ROWS: usize = 64;
pub(crate) const MCP_MAX_BRENTI_DEGREE: usize = 64;
pub(crate) const MCP_MAX_BRENTI_SEARCH_NODES: usize = 1_000_000;
pub(crate) const MCP_MAX_BRENTI_INTEGER_NODES: usize = 100_000;
pub(crate) const MCP_MAX_PF_NETWORK_ROWS: usize = 64;
pub(crate) const MCP_MAX_PF_NETWORK_DEGREE: usize = 64;
pub(crate) const MCP_MAX_PF_NETWORK_STEP: usize = 8;
pub(crate) const MCP_MAX_LEFT_PRODUCTION_ROWS: usize = 32;
/// `enumerate_real_rooted`: largest degree (the search is exponential in the
/// number of free coefficients, so degree alone does not bound the work).
pub(crate) const MCP_MAX_ENUMERATE_DEGREE: usize = 32;
/// `enumerate_real_rooted`: largest bit length of a prescribed leading coefficient.
pub(crate) const MCP_MAX_ENUMERATE_PREFIX_BITS: u64 = 256;
/// `enumerate_real_rooted`: default and largest number of tower nodes
/// (about 5000 nodes per second at degree 32 in release builds).
pub(crate) const MCP_DEFAULT_ENUMERATE_NODES: u64 = 5_000;
pub(crate) const MCP_MAX_ENUMERATE_NODES: u64 = 25_000;
/// `enumerate_real_rooted`: default and largest number of returned polynomials.
pub(crate) const MCP_DEFAULT_ENUMERATE_RESULTS: usize = 500;
pub(crate) const MCP_MAX_ENUMERATE_RESULTS: usize = 2_000;
/// `enumerate_real_rooted`: default and largest coefficient bit length of a
/// derivative-tower polynomial.
pub(crate) const MCP_DEFAULT_ENUMERATE_COEFFICIENT_BITS: u64 = 1_024;
pub(crate) const MCP_MAX_ENUMERATE_COEFFICIENT_BITS: u64 = 2_048;
