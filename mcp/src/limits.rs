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
