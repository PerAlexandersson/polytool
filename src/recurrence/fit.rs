//! Recurrence fitting: fixed-shape and vector search, candidate enumeration,
//! adaptive search, representation menus, and complete bounded recurrence
//! spaces with replayable certificates (also after factorial coefficient
//! rescaling).

use super::*;

mod adaptive;
pub use adaptive::*;
mod candidates;
use candidates::*;
mod rescaling;
pub use rescaling::*;
mod representations;
pub use representations::*;
mod search;
pub use search::*;
mod space;
pub use space::*;
mod vector_search;
pub use vector_search::*;
