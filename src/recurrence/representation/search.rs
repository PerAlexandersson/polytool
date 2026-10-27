//! Bounded queue search over equivalent polynomial representations of one
//! sequence: alternative recurrences, changes of coordinates, searches in the
//! new coordinates, and certified relations fed into later steps.
//!
//! # Queue
//!
//! [`search_representations`] starts from a source: a defining
//! `polytool.recurrence.v1` rule with seeds and origin, a bundled OEIS family,
//! or rows only (exploratory fits). The source rows are representation 0.
//! The queue is breadth first. For every representation below
//! `limits.max_depth`, each transformation of the request-specified set is
//! applied to its rows:
//!
//! 1. the composite map from the original source is reduced to a certified
//!    normal form (module `representation::map`); a map certified equal to an
//!    existing representation's map is a duplicate (`duplicate_map`) and is
//!    not searched again; rows are never compared for deduplication;
//! 2. the rows are mapped exactly and the source rows are reconstructed by
//!    inverting the whole chain, so the original source stays the reference;
//! 3. known relations are transported: every certified scalar description of
//!    the parent (at most `limits.max_sources_per_edge`, simplest first)
//!    becomes the defining rule of a replayable step,
//!    [`transform_recurrence`] (Möbius class, gamma, symmetric decomposition)
//!    or [`implicit_transform_recurrence`] (diagonal maps, derivative
//!    antecedents, with a reconstruction-certified successor fit); without a
//!    certified parent relation the step maps rows only. A parent with
//!    certified implicit relations but no certified scalar rule feeds them to
//!    diagonal maps, derivative antecedents and the lossy `derivative` map
//!    through [`RecurrenceTransformation::transport_implicit`] (an
//!    [`ImplicitRelationStep`]), carrying its coefficient conditions where
//!    the map allows; the `derivative` map always takes this route;
//! 4. for a factorial diagonal map on a constant-index monic certified source,
//!    complete recurrence spaces are certified ([`find_recurrence_spaces`]
//!    with `coefficient_rescaling`) on a bounded menu of supports; after a
//!    Möbius-class map, a polynomial diagonal weight or `derivative`, the
//!    direct implicit search certifies them ([`search_implicit_relations`],
//!    transported-source and weight theorems; `implicit_spaces`, one
//!    successor support per report);
//! 5. the new coordinates are searched by the bounded fit menu
//!    ([`find_recurrence_representations_rational`]).
//!
//! # Vector representations
//!
//! The symmetric decomposition gives two-component states `(a_n, b_n)`.
//! They are expanded like scalar representations: every scalar bijective map
//! acts componentwise ([`RecurrenceTransformation::map_state`]; a schedule
//! acts on `b` at `d_n - 1`), the chain is inverted componentwise, and normal
//! forms keep one componentwise suffix per component. A certified vector
//! system of the parent is transported exactly through Möbius-class maps
//! ([`transport_vector_componentwise`], a [`VectorTransportStep`]); its
//! evidence is `transported_identity` when the parent system and every degree
//! condition are certified (degree bounds of the components are inherited
//! from a certified decomposition: `deg a_n <= d_n`, `deg b_n <= d_n - 1`).
//! Gamma coordinates of the components of a certified decomposition at the
//! same schedule have certified symmetry by construction. States are fitted
//! by a first-order vector system (`vector_fit`, finite evidence). A vector
//! system with its initial state carries a [`VectorSystemAssessment`]: it
//! determines the states when its shared left factor is a nonzero polynomial
//! at every index, certified by one `x`-coefficient without integer roots in
//! the domain; it is then an `equivalent_generating_description` when its
//! annihilation is certified and it starts at the row origin. Two vector
//! representations with identical equivalent systems and initial states are
//! the same sequence. Lossy and vector-producing maps are rejected on vector
//! states; the component count enters the ranking.
//!
//! # The rank-two gamma case
//!
//! When the even and odd parts of a gamma-transported rule are independent,
//! the step reports the pair (`transported.system`) and the representation
//! gets a `transported_gamma_system` description with both relations. They
//! are exact consequences of the source rule where the rows are symmetric,
//! so the description is certified only when that symmetry is: by its own
//! step, or by another step of the same edge whose source rule is invariant
//! under reversal (the symmetry is a property of the rows). Otherwise its
//! evidence is `checked_supplied_rows` and its label `finite_fit`.
//!
//! Every certified scalar description of the new representation is then
//! available to its children. A representation whose equivalent generating
//! description coincides with one of another representation (same relations,
//! seeds, first index and conditions) is the same sequence, certified by the
//! uniqueness arguments of both descriptions; it is kept (`equal_sequence_as`)
//! but not expanded.
//!
//! # Lossy maps
//!
//! A `derivative` map `Q_n = D^m P_n` forgets `[x^0..x^(m-1)] P_n`. The node
//! records these coefficients (`lost_data`), and the chain is inverted with
//! them. Its implicit-relation steps also record a [`LossyReconstruction`]:
//! the description of the parent through `Q` (pulled-back relations and
//! conditions) is assessed with the conditions `[x^k] (L P)_n = 0` that pin
//! the lost coefficients, and again without every condition that refers to
//! them; the first can be equivalent, the second is a consequence only.
//!
//! # Labels and ranking
//!
//! Each [`RepresentationDescription`] carries one of three labels:
//! `finite_fit` (fitted to finitely many rows, or depending on a condition
//! checked on rows only), `certified_annihilating_relation` (annihilation
//! proved for every index of the domain, uniqueness not certified) and
//! `equivalent_generating_description` (certified annihilation together with
//! sufficient initial or boundary data and a certified uniqueness argument,
//! [`assess_generating_description`]). Descriptions are ranked by shift length,
//! derivative order, coefficient (`x`) degree, index degree, component count
//! and reconstruction cost (maps to invert plus coefficient conditions), with
//! the label as the final tie break.
//!
//! # Limits and replay
//!
//! Depth, representation count, backend steps, space supports, rows, degrees,
//! fit candidates and elimination work are bounded; every limit that cuts work
//! is reported as a [`RepresentationBudgetEvent`] and makes the status
//! `budget_exhausted`. [`verify_representation_search`] replays the whole
//! search from the embedded request, replays every embedded step and space
//! certificate with its own verifier, re-derives every representation from
//! the original source rows, and re-checks every description on its rows.
//!
//! Related modules: `representation::report` (request and report types),
//! `transform::registry` (registry and transport), `implicit::report`
//! (implicit transport and reconstruction), `fit::space` and `fit::rescaling`
//! (certified complete spaces), `fit::representations` (fit menus).

use super::map::{
    apply_chain, chain_lost_data, chain_offsets, invert_chain_with_lost_data, MapNormalForm,
};
use super::report::{conventions, label_rank};
use super::*;
use crate::polynomial::Polynomial;
use crate::recurrence::implicit::algebra::{binomial, roots_in_domain};
use crate::recurrence::implicit::{
    biv_from_json, coefficient_conditions_of, decode_poly, encode_poly, finite_theorem_applies,
    implicit_map_condition_statuses, valuation_exceptional_indices, weight_regularity,
};
use crate::vec_poly;
use num_traits::Zero;
use std::collections::VecDeque;

type Q = BigRational;

mod describe;
mod implicit_step;
mod queue;
mod replay;
mod scalar_steps;
mod source;
#[cfg(test)]
mod tests;
mod vector;

use describe::*;
use implicit_step::*;
use queue::*;
use source::*;
use vector::*;

pub use queue::search_representations;
pub use replay::verify_representation_search;
