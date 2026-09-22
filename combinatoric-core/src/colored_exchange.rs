//! Finite fixed-content colored families, symmetric cycle exchange witnesses,
//! exact weighted color marginals, and independently verified decomposition
//! certificates.
//!
//! A *word* is a map from positions `0..n` to colors `0..k`, stored as a
//! `Vec<usize>`. A [`FixedContentFamily`] is a nonempty, sorted, deduplicated
//! set of words of equal length and equal color content (the *quota*
//! `r_0, ..., r_{k-1}`). Zero quotas are allowed, and so is the singleton
//! family consisting of the empty word.
//!
//! # Symmetric cycle exchange
//!
//! The family satisfies *symmetric cycle exchange* (SCE) when, for every two
//! members `u, v` and every position `p` with `u(p) != v(p)`, there are
//! distinct positions `p_1 = p, p_2, ..., p_l` with distinct colors
//! `c_i = u(p_i)` such that `v(p_i) = c_{i+1}` cyclically, and both offspring
//! (`u` overwritten by `v` on these positions, and `v` overwritten by `u`)
//! belong to the family. Such a position list is a [`CycleWitness`]. It is a
//! simple directed cycle through the arrow `u(p) -> v(p)` in the *color
//! transition multigraph* of the pair, whose arrows are the disagreement
//! positions. The two-color case is exactly symmetric basis exchange of an
//! ordinary matroid.
//!
//! Witness search enumerates simple cycles of that multigraph by depth-first
//! traversal in increasing position order, never subsets of positions. Its
//! worst case is exponential in the number of colors. Family witness and
//! axiom checks have variants accepting an explicit step budget. Exhaustion
//! is reported as a distinct
//! inconclusive status and is never treated as a failed axiom.
//!
//! # Marginals and certificates
//!
//! [`FixedContentFamily::weighted_color_marginal`] computes exact `BigUint`
//! coefficients of the weighted one-color projection. The decomposition
//! certificates [`ProductCertificate`], [`BinaryPortCertificate`], and
//! [`PrivateBlockCertificate`] store their factor families explicitly. Their
//! `verify` methods reconstruct the original family exactly, check every
//! claimed quota, palette, level, and allowed-position condition, and check
//! SCE of each factor independently. They deliberately do not infer any
//! Lorentzian or log-concavity conclusion; those implications remain
//! research statements outside this module.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use num_bigint::BigUint;
use num_traits::{One, Zero};

use crate::permutation::next_permutation;

/// Validation failure when constructing a [`FixedContentFamily`].
///
/// Word indices refer to the caller's input order, before sorting and
/// deduplication.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FamilyError {
    /// No words were supplied, or a restriction removed every word.
    Empty,
    /// A word has a different length from the first word.
    RaggedWord {
        index: usize,
        expected: usize,
        actual: usize,
    },
    /// A word uses a color `>= num_colors`.
    ColorOutOfRange {
        index: usize,
        position: usize,
        color: usize,
        num_colors: usize,
    },
    /// A word has a different color content from the first word.
    ContentMismatch {
        index: usize,
        expected: Vec<usize>,
        actual: Vec<usize>,
    },
    /// A per-position list has the wrong number of entries.
    WrongLength {
        what: &'static str,
        expected: usize,
        actual: usize,
    },
}

impl fmt::Display for FamilyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FamilyError::Empty => write!(f, "family must contain at least one word"),
            FamilyError::RaggedWord {
                index,
                expected,
                actual,
            } => write!(f, "word {index} has length {actual}, expected {expected}"),
            FamilyError::ColorOutOfRange {
                index,
                position,
                color,
                num_colors,
            } => write!(
                f,
                "word {index} position {position} uses color {color} >= {num_colors}"
            ),
            FamilyError::ContentMismatch {
                index,
                expected,
                actual,
            } => write!(
                f,
                "word {index} has content {actual:?}, expected {expected:?}"
            ),
            FamilyError::WrongLength {
                what,
                expected,
                actual,
            } => write!(f, "{what} has length {actual}, expected {expected}"),
        }
    }
}

impl std::error::Error for FamilyError {}

/// Rejected exchange query: bad indices, a bad position, or two words that
/// agree at the prescribed position (including `first == second`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum InvalidQuery {
    WordIndexOutOfRange {
        index: usize,
        len: usize,
    },
    PositionOutOfRange {
        position: usize,
        num_positions: usize,
    },
    /// The two words have different lengths (only for the free-standing
    /// cycle enumeration; family members always have equal length).
    LengthMismatch {
        first: usize,
        second: usize,
    },
    /// The words carry the same color at the prescribed position, so there is
    /// no exchange requirement there.
    WordsAgree {
        first: usize,
        second: usize,
        position: usize,
    },
}

impl fmt::Display for InvalidQuery {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            InvalidQuery::WordIndexOutOfRange { index, len } => {
                write!(f, "word index {index} out of range for {len} words")
            }
            InvalidQuery::PositionOutOfRange {
                position,
                num_positions,
            } => write!(
                f,
                "position {position} out of range for {num_positions} positions"
            ),
            InvalidQuery::LengthMismatch { first, second } => {
                write!(f, "words have lengths {first} and {second}")
            }
            InvalidQuery::WordsAgree {
                first,
                second,
                position,
            } => write!(f, "words {first} and {second} agree at position {position}"),
        }
    }
}

impl std::error::Error for InvalidQuery {}

/// An exchange witness: the positions `p_1, ..., p_l` of a simple directed
/// cycle in the color transition multigraph of an ordered word pair.
///
/// Positions are zero-based. Searches always return the prescribed position
/// first; the verifier accepts any rotation. The field is public so callers
/// and tests can construct or corrupt witnesses deliberately.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CycleWitness {
    pub positions: Vec<usize>,
}

impl CycleWitness {
    /// Number of exchanged positions, between two and the number of colors
    /// for a valid witness.
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    /// True when the witness lists no positions (never valid).
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }

    /// The two offspring `(first', second')`: `first` overwritten by `second`
    /// on the witness positions, and vice versa.
    ///
    /// Panics if a witness position is out of range for either word.
    pub fn offspring(&self, first: &[usize], second: &[usize]) -> (Vec<usize>, Vec<usize>) {
        let mut a = first.to_vec();
        let mut b = second.to_vec();
        for &p in &self.positions {
            a[p] = second[p];
            b[p] = first[p];
        }
        (a, b)
    }
}

/// Why a submitted [`CycleWitness`] fails for a given query.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WitnessError {
    Query(InvalidQuery),
    /// Fewer than two positions.
    TooShort {
        len: usize,
    },
    PositionOutOfRange {
        position: usize,
        num_positions: usize,
    },
    RepeatedPosition {
        position: usize,
    },
    /// Two witness positions carry the same color in the first word.
    RepeatedColor {
        color: usize,
    },
    /// `second(p_i) != first(p_{i+1})` for the cyclically next position.
    BrokenCycle {
        index: usize,
    },
    /// The prescribed disagreement position is not on the cycle.
    MissingPrescribedPosition {
        position: usize,
    },
    /// An offspring word is not a member of the family.
    OffspringMissing {
        word: Vec<usize>,
    },
}

impl fmt::Display for WitnessError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            WitnessError::Query(q) => write!(f, "invalid query: {q}"),
            WitnessError::TooShort { len } => write!(f, "witness has {len} positions, need >= 2"),
            WitnessError::PositionOutOfRange {
                position,
                num_positions,
            } => write!(
                f,
                "witness position {position} out of range for {num_positions} positions"
            ),
            WitnessError::RepeatedPosition { position } => {
                write!(f, "witness repeats position {position}")
            }
            WitnessError::RepeatedColor { color } => {
                write!(f, "witness repeats starting color {color}")
            }
            WitnessError::BrokenCycle { index } => {
                write!(f, "cycle identity fails after witness entry {index}")
            }
            WitnessError::MissingPrescribedPosition { position } => {
                write!(f, "witness omits prescribed position {position}")
            }
            WitnessError::OffspringMissing { word } => {
                write!(f, "offspring {word:?} is not in the family")
            }
        }
    }
}

impl std::error::Error for WitnessError {}

impl From<InvalidQuery> for WitnessError {
    fn from(q: InvalidQuery) -> Self {
        WitnessError::Query(q)
    }
}

/// A proved failure of symmetric cycle exchange: no simple cycle through the
/// prescribed disagreement has both offspring in the family. Word indices
/// refer to the family's canonical sorted order; the words are included.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExchangeObstruction {
    pub first: usize,
    pub second: usize,
    pub first_word: Vec<usize>,
    pub second_word: Vec<usize>,
    pub position: usize,
}

impl fmt::Display for ExchangeObstruction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "no symmetric exchange for words {} {:?} and {} {:?} at position {}",
            self.first, self.first_word, self.second, self.second_word, self.position
        )
    }
}

impl std::error::Error for ExchangeObstruction {}

/// Outcome of a budgeted single-pair witness search.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WitnessSearch {
    /// The first witness in traversal order whose offspring both belong to
    /// the family.
    Found(CycleWitness),
    /// Every simple cycle through the prescribed arrow was enumerated and
    /// none has both offspring in the family.
    Absent,
    /// The step budget ran out before the enumeration finished. This is
    /// inconclusive, not a failure.
    BudgetExhausted { steps: u64 },
}

/// Outcome of a budgeted whole-family exchange check.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ExchangeStatus {
    /// Every pair and disagreement position has a witness.
    Satisfied,
    /// A pair and position without any witness.
    Obstructed(ExchangeObstruction),
    /// The total step budget ran out while examining the recorded pair and
    /// position. Inconclusive.
    BudgetExhausted {
        first: usize,
        second: usize,
        position: usize,
        steps: u64,
    },
}

/// Rejected marginal request.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MarginalError {
    ColorOutOfRange {
        color: usize,
        num_colors: usize,
    },
    /// `weights.len()` differs from the number of positions.
    WeightRows {
        expected: usize,
        actual: usize,
    },
    /// `weights[row].len()` differs from the number of colors.
    WeightRowLength {
        row: usize,
        expected: usize,
        actual: usize,
    },
    /// The quadratic matrix form requires the chosen color to have quota two.
    QuotaNotTwo {
        color: usize,
        quota: usize,
    },
}

impl fmt::Display for MarginalError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MarginalError::ColorOutOfRange { color, num_colors } => {
                write!(f, "color {color} out of range for {num_colors} colors")
            }
            MarginalError::WeightRows { expected, actual } => {
                write!(f, "weights have {actual} rows, expected {expected}")
            }
            MarginalError::WeightRowLength {
                row,
                expected,
                actual,
            } => write!(
                f,
                "weight row {row} has {actual} entries, expected {expected}"
            ),
            MarginalError::QuotaNotTwo { color, quota } => {
                write!(f, "color {color} has quota {quota}, expected 2")
            }
        }
    }
}

impl std::error::Error for MarginalError {}

/// Convert an uppercase-letter word such as `"ABCA"` to colors `A=0, B=1, ...`.
/// Returns `None` if any character is not an ASCII uppercase letter.
pub fn letters_to_word(letters: &str) -> Option<Vec<usize>> {
    letters
        .chars()
        .map(|ch| ch.is_ascii_uppercase().then(|| (ch as u8 - b'A') as usize))
        .collect()
}

/// Render a word with letters `A..Z` for colors below 26 and `[c]` otherwise.
pub fn word_to_letters(word: &[usize]) -> String {
    word.iter()
        .map(|&c| {
            if c < 26 {
                ((b'A' + c as u8) as char).to_string()
            } else {
                format!("[{c}]")
            }
        })
        .collect()
}

/// A nonempty, sorted, deduplicated set of equal-length words with equal
/// color content.
///
/// Word indices used throughout this module refer to the canonical
/// lexicographically sorted order returned by [`FixedContentFamily::words`].
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct FixedContentFamily {
    num_colors: usize,
    quota: Vec<usize>,
    words: Vec<Vec<usize>>,
}

impl FixedContentFamily {
    /// Validate, sort, and deduplicate `words` over colors `0..num_colors`.
    ///
    /// Errors: [`FamilyError::Empty`] for no words, and ragged, out-of-range,
    /// or content-mismatch errors reported with the caller's input index.
    /// Colors with zero quota are permitted; so is `try_new(k, vec![vec![]])`.
    pub fn try_new(num_colors: usize, mut words: Vec<Vec<usize>>) -> Result<Self, FamilyError> {
        let n = words.first().ok_or(FamilyError::Empty)?.len();
        let mut quota = vec![0usize; num_colors];
        for (index, word) in words.iter().enumerate() {
            if word.len() != n {
                return Err(FamilyError::RaggedWord {
                    index,
                    expected: n,
                    actual: word.len(),
                });
            }
            let mut content = vec![0usize; num_colors];
            for (position, &color) in word.iter().enumerate() {
                if color >= num_colors {
                    return Err(FamilyError::ColorOutOfRange {
                        index,
                        position,
                        color,
                        num_colors,
                    });
                }
                content[color] += 1;
            }
            if index == 0 {
                quota = content;
            } else if content != quota {
                return Err(FamilyError::ContentMismatch {
                    index,
                    expected: quota,
                    actual: content,
                });
            }
        }
        words.sort();
        words.dedup();
        Ok(Self {
            num_colors,
            quota,
            words,
        })
    }

    /// Build a family from uppercase-letter words (`A=0, B=1, ...`).
    ///
    /// A non-uppercase character is reported as a color out of range at that
    /// word and position.
    pub fn from_letter_words(num_colors: usize, letters: &[&str]) -> Result<Self, FamilyError> {
        let mut words = Vec::with_capacity(letters.len());
        for (index, text) in letters.iter().enumerate() {
            match letters_to_word(text) {
                Some(word) => words.push(word),
                None => {
                    let position = text
                        .chars()
                        .position(|ch| !ch.is_ascii_uppercase())
                        .unwrap_or(0);
                    return Err(FamilyError::ColorOutOfRange {
                        index,
                        position,
                        color: usize::MAX,
                        num_colors,
                    });
                }
            }
        }
        Self::try_new(num_colors, words)
    }

    /// All words with the given quota, using `quota.len()` colors.
    ///
    /// Enumerates the multinomial coefficient `n! / prod r_c!` words with
    /// [`next_permutation`]; this is exponential in `n`. The all-zero quota
    /// gives the singleton empty word. Every such family satisfies SCE.
    pub fn all_assignments(quota: &[usize]) -> Self {
        let mut word: Vec<usize> = quota
            .iter()
            .enumerate()
            .flat_map(|(c, &r)| std::iter::repeat_n(c, r))
            .collect();
        let mut words = Vec::new();
        loop {
            words.push(word.clone());
            if !next_permutation(&mut word) {
                break;
            }
        }
        Self {
            num_colors: quota.len(),
            quota: quota.to_vec(),
            words,
        }
    }

    /// Keep only the words whose color at every position `p` lies in
    /// `allowed[p]`. Colors outside `0..num_colors` in a list are ignored.
    ///
    /// Errors: [`FamilyError::WrongLength`] if `allowed.len()` is not the
    /// number of positions; [`FamilyError::Empty`] if nothing survives.
    /// Allowed-cell restriction preserves SCE when the result is nonempty.
    pub fn restricted_to_allowed_colors(
        &self,
        allowed: &[Vec<usize>],
    ) -> Result<Self, FamilyError> {
        let n = self.num_positions();
        if allowed.len() != n {
            return Err(FamilyError::WrongLength {
                what: "allowed color lists",
                expected: n,
                actual: allowed.len(),
            });
        }
        let words: Vec<Vec<usize>> = self
            .words
            .iter()
            .filter(|w| w.iter().zip(allowed).all(|(c, list)| list.contains(c)))
            .cloned()
            .collect();
        if words.is_empty() {
            return Err(FamilyError::Empty);
        }
        Ok(Self {
            num_colors: self.num_colors,
            quota: self.quota.clone(),
            words,
        })
    }

    /// For each position, the sorted list of colors occurring there.
    pub fn occurring_colors_by_position(&self) -> Vec<Vec<usize>> {
        (0..self.num_positions())
            .map(|p| {
                let set: BTreeSet<usize> = self.words.iter().map(|w| w[p]).collect();
                set.into_iter().collect()
            })
            .collect()
    }

    /// True when the family consists of *all* quota-compatible words whose
    /// colors lie in the occurring color lists. Such list/quota families are
    /// proved SCE. Cost: enumerates [`Self::all_assignments`] for the quota.
    pub fn is_list_quota_saturated(&self) -> bool {
        let lists = self.occurring_colors_by_position();
        match Self::all_assignments(&self.quota).restricted_to_allowed_colors(&lists) {
            Ok(full) => full.words.len() == self.words.len(),
            Err(_) => false,
        }
    }

    /// Canonically sorted, deduplicated words.
    pub fn words(&self) -> &[Vec<usize>] {
        &self.words
    }

    /// Word at canonical index `i`; panics if out of range.
    pub fn word(&self, i: usize) -> &[usize] {
        &self.words[i]
    }

    /// Canonical index of a word, if it is a member.
    pub fn index_of(&self, word: &[usize]) -> Option<usize> {
        self.words.binary_search_by(|w| w.as_slice().cmp(word)).ok()
    }

    /// Membership test.
    pub fn contains(&self, word: &[usize]) -> bool {
        self.index_of(word).is_some()
    }

    /// Number of distinct words.
    pub fn len(&self) -> usize {
        self.words.len()
    }

    /// Always false: families are nonempty by construction.
    pub fn is_empty(&self) -> bool {
        self.words.is_empty()
    }

    /// Color content `r_0, ..., r_{k-1}` shared by all words.
    pub fn quota(&self) -> &[usize] {
        &self.quota
    }

    /// Number of colors `k` (including zero-quota colors).
    pub fn num_colors(&self) -> usize {
        self.num_colors
    }

    /// Word length `n`.
    pub fn num_positions(&self) -> usize {
        self.words[0].len()
    }

    fn query(
        &self,
        first: usize,
        second: usize,
        position: usize,
    ) -> Result<(&[usize], &[usize]), InvalidQuery> {
        let len = self.words.len();
        for index in [first, second] {
            if index >= len {
                return Err(InvalidQuery::WordIndexOutOfRange { index, len });
            }
        }
        let num_positions = self.num_positions();
        if position >= num_positions {
            return Err(InvalidQuery::PositionOutOfRange {
                position,
                num_positions,
            });
        }
        let (u, v) = (&self.words[first], &self.words[second]);
        if u[position] == v[position] {
            return Err(InvalidQuery::WordsAgree {
                first,
                second,
                position,
            });
        }
        Ok((u, v))
    }

    fn search_pair(
        &self,
        u: &[usize],
        v: &[usize],
        position: usize,
        max_steps: Option<u64>,
    ) -> (WitnessSearch, u64) {
        let mut found = None;
        let outcome = traverse_simple_cycles(u, v, position, max_steps, &mut |path| {
            let witness = CycleWitness {
                positions: path.to_vec(),
            };
            let (a, b) = witness.offspring(u, v);
            if self.contains(&a) && self.contains(&b) {
                found = Some(witness);
                true
            } else {
                false
            }
        });
        let steps = match outcome {
            Traversal::Stopped { steps }
            | Traversal::Completed { steps }
            | Traversal::Exhausted { steps } => steps,
        };
        let result = match (found, outcome) {
            (Some(w), _) => WitnessSearch::Found(w),
            (None, Traversal::Exhausted { steps }) => WitnessSearch::BudgetExhausted { steps },
            (None, _) => WitnessSearch::Absent,
        };
        (result, steps)
    }

    /// First exchange witness (in traversal order) for words `first`,
    /// `second` at `prescribed_position`, or `None` when no simple cycle
    /// through that arrow has both offspring in the family.
    ///
    /// Unbounded: worst case exponential in the number of colors. Use
    /// [`Self::exchange_witness_bounded`] to impose a step budget. Errors on
    /// bad indices, a bad position, or agreement at the position.
    pub fn exchange_witness(
        &self,
        first: usize,
        second: usize,
        prescribed_position: usize,
    ) -> Result<Option<CycleWitness>, InvalidQuery> {
        let (u, v) = self.query(first, second, prescribed_position)?;
        Ok(match self.search_pair(u, v, prescribed_position, None).0 {
            WitnessSearch::Found(w) => Some(w),
            _ => None,
        })
    }

    /// Like [`Self::exchange_witness`] with at most `max_steps` traversal
    /// steps (one step per arrow extension). Exhaustion yields
    /// [`WitnessSearch::BudgetExhausted`], never `Absent`.
    pub fn exchange_witness_bounded(
        &self,
        first: usize,
        second: usize,
        prescribed_position: usize,
        max_steps: u64,
    ) -> Result<WitnessSearch, InvalidQuery> {
        let (u, v) = self.query(first, second, prescribed_position)?;
        Ok(self
            .search_pair(u, v, prescribed_position, Some(max_steps))
            .0)
    }

    /// Validate a submitted witness for the query: distinct positions in
    /// range, distinct starting colors, the directed cycle identities,
    /// inclusion of the prescribed position (any rotation), and membership
    /// of both offspring. For word length `n`, palette size `k`, and family
    /// size `m`, uses `O(n + k)` temporary space and `O(k + n log(m + 1))`
    /// time, including full-word membership searches.
    pub fn verify_exchange(
        &self,
        first: usize,
        second: usize,
        prescribed_position: usize,
        witness: &CycleWitness,
    ) -> Result<(), WitnessError> {
        let (u, v) = self.query(first, second, prescribed_position)?;
        let ps = &witness.positions;
        if ps.len() < 2 {
            return Err(WitnessError::TooShort { len: ps.len() });
        }
        let n = self.num_positions();
        let mut seen_position = vec![false; n];
        let mut seen_color = vec![false; self.num_colors];
        for &p in ps {
            if p >= n {
                return Err(WitnessError::PositionOutOfRange {
                    position: p,
                    num_positions: n,
                });
            }
            if seen_position[p] {
                return Err(WitnessError::RepeatedPosition { position: p });
            }
            seen_position[p] = true;
            if seen_color[u[p]] {
                return Err(WitnessError::RepeatedColor { color: u[p] });
            }
            seen_color[u[p]] = true;
        }
        for (i, &p) in ps.iter().enumerate() {
            let next = ps[(i + 1) % ps.len()];
            if v[p] != u[next] {
                return Err(WitnessError::BrokenCycle { index: i });
            }
        }
        if !seen_position[prescribed_position] {
            return Err(WitnessError::MissingPrescribedPosition {
                position: prescribed_position,
            });
        }
        let (a, b) = witness.offspring(u, v);
        for word in [a, b] {
            if !self.contains(&word) {
                return Err(WitnessError::OffspringMissing { word });
            }
        }
        Ok(())
    }

    /// Check symmetric cycle exchange for every unordered pair and every
    /// disagreement position, without a budget.
    ///
    /// Cost: `O(|F|^2)` pairs, each with a cycle traversal that is exponential
    /// in the number of colors in the worst case (for three colors only two-
    /// and three-cycles exist). Returns the first obstruction in canonical
    /// pair order `(first < second)` and increasing position.
    pub fn check_symmetric_cycle_exchange(&self) -> Result<(), ExchangeObstruction> {
        match self.check_exchange_with_budget(None) {
            ExchangeStatus::Satisfied => Ok(()),
            ExchangeStatus::Obstructed(o) => Err(o),
            ExchangeStatus::BudgetExhausted { .. } => {
                unreachable!("unbounded exchange check cannot exhaust its budget")
            }
        }
    }

    /// Budgeted variant of [`Self::check_symmetric_cycle_exchange`] with a
    /// total budget of `max_steps` traversal steps across all pairs.
    /// Exhaustion returns [`ExchangeStatus::BudgetExhausted`] naming the pair
    /// and position under examination; it is not an obstruction.
    pub fn check_symmetric_cycle_exchange_bounded(&self, max_steps: u64) -> ExchangeStatus {
        self.check_exchange_with_budget(Some(max_steps))
    }

    fn check_exchange_with_budget(&self, max_steps: Option<u64>) -> ExchangeStatus {
        let n = self.num_positions();
        let mut used: u64 = 0;
        for i in 0..self.words.len() {
            for j in i + 1..self.words.len() {
                let (u, v) = (&self.words[i], &self.words[j]);
                let mut satisfied = vec![false; n];
                for p in 0..n {
                    if u[p] == v[p] || satisfied[p] {
                        continue;
                    }
                    let (result, steps) = self.search_pair(
                        u,
                        v,
                        p,
                        max_steps.map(|limit| limit.saturating_sub(used)),
                    );
                    used = used.saturating_add(steps);
                    match result {
                        WitnessSearch::Found(w) => {
                            for &q in &w.positions {
                                satisfied[q] = true;
                            }
                        }
                        WitnessSearch::Absent => {
                            return ExchangeStatus::Obstructed(ExchangeObstruction {
                                first: i,
                                second: j,
                                first_word: u.clone(),
                                second_word: v.clone(),
                                position: p,
                            });
                        }
                        WitnessSearch::BudgetExhausted { .. } => {
                            return ExchangeStatus::BudgetExhausted {
                                first: i,
                                second: j,
                                position: p,
                                steps: used,
                            };
                        }
                    }
                }
            }
        }
        ExchangeStatus::Satisfied
    }

    fn validate_marginal_weights(
        &self,
        color: usize,
        weights: &[Vec<BigUint>],
    ) -> Result<(), MarginalError> {
        if color >= self.num_colors {
            return Err(MarginalError::ColorOutOfRange {
                color,
                num_colors: self.num_colors,
            });
        }
        let n = self.num_positions();
        if weights.len() != n {
            return Err(MarginalError::WeightRows {
                expected: n,
                actual: weights.len(),
            });
        }
        for (row, w) in weights.iter().enumerate() {
            if w.len() != self.num_colors {
                return Err(MarginalError::WeightRowLength {
                    row,
                    expected: self.num_colors,
                    actual: w.len(),
                });
            }
        }
        Ok(())
    }

    fn fold_color_marginal<W>(&self, color: usize, mut weight: W) -> BTreeMap<Vec<usize>, BigUint>
    where
        W: FnMut(&[usize]) -> BigUint,
    {
        let mut marginal: BTreeMap<Vec<usize>, BigUint> = BTreeMap::new();
        for w in &self.words {
            let value = weight(w);
            if value.is_zero() {
                continue;
            }
            let key: Vec<usize> = (0..w.len()).filter(|&p| w[p] == color).collect();
            *marginal.entry(key).or_insert_with(BigUint::zero) += value;
        }
        marginal
    }

    /// Exact weighted color marginal
    /// `B_c(x; w) = sum_u prod_{u(e)=c} x_e prod_{u(e)!=c} w[e][u(e)]`
    /// as a sparse map from the sorted position set `{e : u(e) = c}` to its
    /// coefficient.
    ///
    /// `weights` must have exactly `num_positions` rows of exactly
    /// `num_colors` entries; the entries in column `color` are ignored.
    /// Words of total weight zero contribute nothing, and keys with a zero
    /// coefficient are omitted, so the map is exactly the support of the
    /// polynomial. A zero-quota color yields at most the empty key. Cost is
    /// `O(|F| n)` big-integer multiplications.
    pub fn weighted_color_marginal(
        &self,
        color: usize,
        weights: &[Vec<BigUint>],
    ) -> Result<BTreeMap<Vec<usize>, BigUint>, MarginalError> {
        self.validate_marginal_weights(color, weights)?;
        Ok(self.fold_color_marginal(color, |w| {
            let mut value = BigUint::one();
            for (e, &d) in w.iter().enumerate() {
                if d != color {
                    if weights[e][d].is_zero() {
                        return BigUint::zero();
                    }
                    value *= &weights[e][d];
                }
            }
            value
        }))
    }

    /// Unweighted fiber counts: the number of words with each color-`color`
    /// position set. Same traversal as the weighted marginal with all
    /// weights one.
    pub fn color_fiber_counts(
        &self,
        color: usize,
    ) -> Result<BTreeMap<Vec<usize>, BigUint>, MarginalError> {
        if color >= self.num_colors {
            return Err(MarginalError::ColorOutOfRange {
                color,
                num_colors: self.num_colors,
            });
        }
        Ok(self.fold_color_marginal(color, |_| BigUint::one()))
    }

    /// Support of the color-`color` projection: the distinct sorted position
    /// sets `{e : u(e) = color}`. For an SCE family this is a matroid basis
    /// family (structural consequence, not re-verified here).
    pub fn color_support(&self, color: usize) -> Result<BTreeSet<Vec<usize>>, MarginalError> {
        Ok(self.color_fiber_counts(color)?.into_keys().collect())
    }

    /// For a color of quota two, the symmetric `n x n` matrix `H` with
    /// `H[i][j]` the marginal coefficient of `x_i x_j` (zero diagonal). This
    /// is the Hessian of the quadratic marginal; no spectral conclusion is
    /// drawn here. Errors as in [`Self::weighted_color_marginal`] plus
    /// [`MarginalError::QuotaNotTwo`].
    pub fn weighted_quadratic_marginal_matrix(
        &self,
        color: usize,
        weights: &[Vec<BigUint>],
    ) -> Result<Vec<Vec<BigUint>>, MarginalError> {
        let marginal = self.weighted_color_marginal(color, weights)?;
        if self.quota[color] != 2 {
            return Err(MarginalError::QuotaNotTwo {
                color,
                quota: self.quota[color],
            });
        }
        let n = self.num_positions();
        let mut matrix = vec![vec![BigUint::zero(); n]; n];
        for (key, value) in marginal {
            matrix[key[0]][key[1]] = value.clone();
            matrix[key[1]][key[0]] = value;
        }
        Ok(matrix)
    }

    /// Sorted colors occurring at any of `positions` in any word.
    fn palette(&self, positions: &[usize]) -> Vec<usize> {
        let set: BTreeSet<usize> = self
            .words
            .iter()
            .flat_map(|w| positions.iter().map(move |&p| w[p]))
            .collect();
        set.into_iter().collect()
    }

    fn check_position_subset(&self, positions: &[usize]) -> Result<(), CertificateError> {
        let n = self.num_positions();
        let strictly_increasing = positions.windows(2).all(|w| w[0] < w[1]);
        let in_range = positions.iter().all(|&p| p < n);
        if positions.is_empty() || positions.len() >= n || !strictly_increasing || !in_range {
            return Err(CertificateError::InvalidPositions {
                positions: positions.to_vec(),
                num_positions: n,
            });
        }
        Ok(())
    }

    fn compare_reconstruction(
        &self,
        reconstructed: BTreeSet<Vec<usize>>,
        pairs: usize,
    ) -> Result<(), CertificateError> {
        let original: BTreeSet<Vec<usize>> = self.words.iter().cloned().collect();
        let missing: Vec<Vec<usize>> = original.difference(&reconstructed).cloned().collect();
        let extra: Vec<Vec<usize>> = reconstructed.difference(&original).cloned().collect();
        if !missing.is_empty() || !extra.is_empty() || pairs != self.words.len() {
            return Err(CertificateError::ReconstructionMismatch {
                pairs,
                family_size: self.words.len(),
                missing,
                extra,
            });
        }
        Ok(())
    }

    /// Discover a [`ProductCertificate`] across the explicit `cut` (sorted
    /// left positions), then verify it. The left factor lives on the cut
    /// positions in increasing order, the right factor on the complement.
    ///
    /// Errors: invalid cut; a restriction side without fixed content
    /// ([`CertificateError::FactorInvalid`]); the family is not the product
    /// of its restrictions ([`CertificateError::ReconstructionMismatch`]); or
    /// a factor fails SCE. No search over cuts is performed.
    pub fn product_certificate(
        &self,
        cut: &[usize],
    ) -> Result<ProductCertificate, CertificateError> {
        self.check_position_subset(cut)?;
        let right_positions = complement(cut, self.num_positions());
        let left = build_factor(
            Factor::Left,
            self.num_colors,
            self.words.iter().map(|w| restrict(w, cut)).collect(),
        )?;
        let right = build_factor(
            Factor::Right,
            self.num_colors,
            self.words
                .iter()
                .map(|w| restrict(w, &right_positions))
                .collect(),
        )?;
        let certificate = ProductCertificate {
            cut: cut.to_vec(),
            left,
            right,
        };
        certificate.verify(self)?;
        Ok(certificate)
    }

    /// Discover a [`BinaryPortCertificate`] across the explicit `cut`, then
    /// verify it. The side palettes must share exactly two colors `a < b`;
    /// the left `a`-count must take exactly the values `m, m+1`. Ports are
    /// appended as the last position of each factor: the left port is `a`
    /// at count `m` and `b` at count `m+1`, the right port is the opposite.
    ///
    /// Errors: invalid cut, [`CertificateError::SharedPaletteMismatch`],
    /// [`CertificateError::CountLevels`] (a product when the width is zero,
    /// higher connectivity otherwise), factor construction, reconstruction,
    /// or factor SCE failures.
    pub fn binary_port_certificate(
        &self,
        cut: &[usize],
    ) -> Result<BinaryPortCertificate, CertificateError> {
        self.check_position_subset(cut)?;
        let right_positions = complement(cut, self.num_positions());
        let left_palette = self.palette(cut);
        let right_palette = self.palette(&right_positions);
        let shared: Vec<usize> = left_palette
            .iter()
            .filter(|c| right_palette.contains(c))
            .copied()
            .collect();
        if shared.len() != 2 {
            return Err(CertificateError::SharedPaletteMismatch { shared });
        }
        let (a, b) = (shared[0], shared[1]);
        let counts: Vec<usize> = self
            .words
            .iter()
            .map(|w| cut.iter().filter(|&&p| w[p] == a).count())
            .collect();
        let low = *counts.iter().min().expect("nonempty family");
        let high = *counts.iter().max().expect("nonempty family");
        if high != low + 1 {
            let values: BTreeSet<usize> = counts.iter().copied().collect();
            return Err(CertificateError::CountLevels {
                color: a,
                values: values.into_iter().collect(),
            });
        }
        let mut left_words = Vec::with_capacity(self.words.len());
        let mut right_words = Vec::with_capacity(self.words.len());
        for (w, &count) in self.words.iter().zip(&counts) {
            let mut l = restrict(w, cut);
            let mut r = restrict(w, &right_positions);
            l.push(if count == low { a } else { b });
            r.push(if count == low { b } else { a });
            left_words.push(l);
            right_words.push(r);
        }
        let certificate = BinaryPortCertificate {
            cut: cut.to_vec(),
            counted_color: a,
            other_color: b,
            low_count: low,
            left: build_factor(Factor::Left, self.num_colors, left_words)?,
            right: build_factor(Factor::Right, self.num_colors, right_words)?,
        };
        certificate.verify(self)?;
        Ok(certificate)
    }

    /// Discover a [`PrivateBlockCertificate`] for the explicit `block` of
    /// positions and sorted `private_colors`, then verify it.
    ///
    /// Every occurrence of a private color must lie in the block, and every
    /// word must have exactly one block position with a non-private
    /// (*transmitted*) color. The compressed factor keeps the outside
    /// positions in increasing order followed by a port carrying the
    /// transmitted color; the internal factor keeps the block positions in
    /// increasing order with the active position recolored by the fresh
    /// color `alpha = num_colors`. Allowed active positions are recorded per
    /// transmitted color as block-relative indices.
    ///
    /// Errors: invalid block or colors,
    /// [`CertificateError::PrivateColorOutsideBlock`],
    /// [`CertificateError::ActiveCountNotOne`], factor construction,
    /// reconstruction, or factor SCE failures.
    pub fn private_block_certificate(
        &self,
        block: &[usize],
        private_colors: &[usize],
    ) -> Result<PrivateBlockCertificate, CertificateError> {
        self.check_position_subset(block)?;
        let is_private = self.private_color_mask(private_colors)?;
        let outside = complement(block, self.num_positions());
        let alpha = self.num_colors;
        let mut compressed_words = Vec::with_capacity(self.words.len());
        let mut internal_words = Vec::with_capacity(self.words.len());
        let mut allowed: BTreeMap<usize, BTreeSet<usize>> = BTreeMap::new();
        for (index, w) in self.words.iter().enumerate() {
            check_private_block_word(index, w, block, &outside, &is_private)?;
            let active = block
                .iter()
                .position(|&p| !is_private[w[p]])
                .expect("exactly one active position");
            let transmitted = w[block[active]];
            let mut g = restrict(w, &outside);
            g.push(transmitted);
            compressed_words.push(g);
            let mut h = restrict(w, block);
            h[active] = alpha;
            internal_words.push(h);
            allowed.entry(transmitted).or_default().insert(active);
        }
        let certificate = PrivateBlockCertificate {
            block: block.to_vec(),
            private_colors: private_colors.to_vec(),
            compressed: build_factor(Factor::Compressed, self.num_colors, compressed_words)?,
            internal: build_factor(Factor::Internal, alpha + 1, internal_words)?,
            allowed_positions: allowed
                .into_iter()
                .map(|(d, set)| (d, set.into_iter().collect()))
                .collect(),
        };
        certificate.verify(self)?;
        Ok(certificate)
    }

    fn private_color_mask(&self, private_colors: &[usize]) -> Result<Vec<bool>, CertificateError> {
        let increasing = private_colors.windows(2).all(|w| w[0] < w[1]);
        let in_range = private_colors.iter().all(|&c| c < self.num_colors);
        if private_colors.is_empty() || !increasing || !in_range {
            return Err(CertificateError::InvalidColors {
                colors: private_colors.to_vec(),
                num_colors: self.num_colors,
            });
        }
        let mut mask = vec![false; self.num_colors];
        for &c in private_colors {
            mask[c] = true;
        }
        Ok(mask)
    }
}

/// Outcome of a raw cycle traversal.
#[derive(Debug, Clone, PartialEq, Eq)]
enum Traversal {
    /// The visitor asked to stop.
    Stopped { steps: u64 },
    /// Every simple cycle through the prescribed arrow was visited.
    Completed { steps: u64 },
    /// The step budget ran out.
    Exhausted { steps: u64 },
}

enum Control {
    Continue,
    Stop,
    Exhausted,
}

struct CycleContext<'a> {
    second: &'a [usize],
    by_color: Vec<Vec<usize>>,
    target: usize,
    path: Vec<usize>,
    used: Vec<bool>,
    steps: u64,
    max_steps: Option<u64>,
}

impl CycleContext<'_> {
    fn extend<F: FnMut(&[usize]) -> bool>(&mut self, current: usize, visit: &mut F) -> Control {
        if current == self.target {
            return if visit(&self.path) {
                Control::Stop
            } else {
                Control::Continue
            };
        }
        for idx in 0..self.by_color[current].len() {
            let q = self.by_color[current][idx];
            let next = self.second[q];
            if next != self.target && self.used[next] {
                continue;
            }
            if self.max_steps.is_some_and(|limit| self.steps >= limit) {
                return Control::Exhausted;
            }
            self.steps = self.steps.saturating_add(1);
            self.path.push(q);
            self.used[current] = true;
            let control = self.extend(next, visit);
            self.used[current] = false;
            self.path.pop();
            if !matches!(control, Control::Continue) {
                return control;
            }
        }
        Control::Continue
    }
}

/// Depth-first enumeration of simple directed cycles through the arrow
/// `first[position] -> second[position]`. Requires equal lengths and a
/// disagreement at `position`. Visits paths in lexicographic position order.
fn traverse_simple_cycles<F: FnMut(&[usize]) -> bool>(
    first: &[usize],
    second: &[usize],
    position: usize,
    max_steps: Option<u64>,
    visit: &mut F,
) -> Traversal {
    // Raw callers may use sparse labels, including usize::MAX. Compress
    // the observed palette so allocations depend only on word length.
    let labels: BTreeSet<usize> = first.iter().chain(second).copied().collect();
    let dense: BTreeMap<usize, usize> = labels
        .into_iter()
        .enumerate()
        .map(|(index, label)| (label, index))
        .collect();
    let first: Vec<usize> = first.iter().map(|c| dense[c]).collect();
    let second: Vec<usize> = second.iter().map(|c| dense[c]).collect();
    let num_colors = dense.len();
    let mut by_color = vec![Vec::new(); num_colors];
    for (q, (&a, &b)) in first.iter().zip(&second).enumerate() {
        if a != b {
            by_color[a].push(q);
        }
    }
    let target = first[position];
    let mut used = vec![false; num_colors];
    used[target] = true;
    let mut ctx = CycleContext {
        second: &second,
        by_color,
        target,
        path: vec![position],
        used,
        steps: 0,
        max_steps,
    };
    match ctx.extend(second[position], visit) {
        Control::Stop => Traversal::Stopped { steps: ctx.steps },
        Control::Continue => Traversal::Completed { steps: ctx.steps },
        Control::Exhausted => Traversal::Exhausted { steps: ctx.steps },
    }
}

/// All simple exchange cycles through `position` for an arbitrary word pair,
/// independent of any family, in deterministic traversal order (each starts
/// with `position`).
///
/// Errors on unequal lengths, an out-of-range position, or agreement at the
/// position. Unbounded and exponential in the number of colors in the worst
/// case; for three colors at most two- and three-cycles exist.
pub fn simple_exchange_cycles(
    first: &[usize],
    second: &[usize],
    position: usize,
) -> Result<Vec<CycleWitness>, InvalidQuery> {
    if first.len() != second.len() {
        return Err(InvalidQuery::LengthMismatch {
            first: first.len(),
            second: second.len(),
        });
    }
    if position >= first.len() {
        return Err(InvalidQuery::PositionOutOfRange {
            position,
            num_positions: first.len(),
        });
    }
    if first[position] == second[position] {
        return Err(InvalidQuery::WordsAgree {
            first: 0,
            second: 1,
            position,
        });
    }
    let mut cycles = Vec::new();
    traverse_simple_cycles(first, second, position, None, &mut |path| {
        cycles.push(CycleWitness {
            positions: path.to_vec(),
        });
        false
    });
    Ok(cycles)
}

/// Role of a factor family inside a certificate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Factor {
    Left,
    Right,
    Compressed,
    Internal,
}

impl fmt::Display for Factor {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Factor::Left => "left",
            Factor::Right => "right",
            Factor::Compressed => "compressed",
            Factor::Internal => "internal",
        })
    }
}

/// Why a decomposition certificate is rejected.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CertificateError {
    /// Cut or block is not a nonempty, proper, strictly increasing, in-range
    /// position list.
    InvalidPositions {
        positions: Vec<usize>,
        num_positions: usize,
    },
    /// Color list is empty, not strictly increasing, or out of range; or a
    /// pair of shared colors coincides or is out of range.
    InvalidColors {
        colors: Vec<usize>,
        num_colors: usize,
    },
    /// A factor could not be built as a fixed-content family.
    FactorInvalid { factor: Factor, error: FamilyError },
    /// A factor has the wrong number of positions or colors.
    FactorShape {
        factor: Factor,
        expected_positions: usize,
        actual_positions: usize,
        expected_colors: usize,
        actual_colors: usize,
    },
    /// A factor's quota disagrees with the value forced by the family.
    QuotaMismatch {
        factor: Factor,
        expected: Vec<usize>,
        actual: Vec<usize>,
    },
    /// The two side palettes do not intersect in exactly the claimed pair.
    SharedPaletteMismatch { shared: Vec<usize> },
    /// The left counts of the counted color are not exactly two consecutive
    /// values `m, m+1` matching the claimed `low_count`.
    CountLevels { color: usize, values: Vec<usize> },
    /// A factor word's port color is not one of the two shared colors.
    PortColorInvalid { factor: Factor, word: Vec<usize> },
    /// A factor word's non-port count of the counted color does not match
    /// the level announced by its port color.
    PortLevelMismatch {
        factor: Factor,
        word: Vec<usize>,
        count: usize,
    },
    /// One of the two port colors never occurs in a factor.
    PortLevelMissing { factor: Factor, color: usize },
    /// A private color occurs outside the block.
    PrivateColorOutsideBlock {
        word: usize,
        position: usize,
        color: usize,
    },
    /// A word does not have exactly one non-private position in the block.
    ActiveCountNotOne { word: usize, count: usize },
    /// The allowed-position keys differ from the colors occurring at the
    /// compressed port.
    AllowedPositionsKeys {
        expected: Vec<usize>,
        actual: Vec<usize>,
    },
    /// An allowed-position list is empty, not strictly increasing, or not a
    /// block-relative index.
    AllowedPositionsInvalid { color: usize, positions: Vec<usize> },
    /// The reconstructed word set or pair count differs from the family.
    ReconstructionMismatch {
        pairs: usize,
        family_size: usize,
        missing: Vec<Vec<usize>>,
        extra: Vec<Vec<usize>>,
    },
    /// A submitted private-block factor or allowed list contains data that
    /// occurs in no compatible substitution pair.
    UnusedPrivateBlockData,
    /// A factor fails symmetric cycle exchange.
    FactorNotExchange {
        factor: Factor,
        obstruction: ExchangeObstruction,
    },
}

impl fmt::Display for CertificateError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CertificateError::InvalidPositions {
                positions,
                num_positions,
            } => write!(
                f,
                "positions {positions:?} are not a nonempty proper increasing subset of 0..{num_positions}"
            ),
            CertificateError::InvalidColors { colors, num_colors } => {
                write!(f, "colors {colors:?} are not valid distinct colors in 0..{num_colors}")
            }
            CertificateError::FactorInvalid { factor, error } => {
                write!(f, "{factor} factor is not a fixed-content family: {error}")
            }
            CertificateError::FactorShape {
                factor,
                expected_positions,
                actual_positions,
                expected_colors,
                actual_colors,
            } => write!(
                f,
                "{factor} factor has {actual_positions} positions and {actual_colors} colors, expected {expected_positions} and {expected_colors}"
            ),
            CertificateError::QuotaMismatch {
                factor,
                expected,
                actual,
            } => write!(f, "{factor} factor quota {actual:?}, expected {expected:?}"),
            CertificateError::SharedPaletteMismatch { shared } => {
                write!(f, "side palettes share colors {shared:?}, expected exactly the claimed pair")
            }
            CertificateError::CountLevels { color, values } => write!(
                f,
                "left counts of color {color} take values {values:?}, expected exactly the claimed two consecutive levels"
            ),
            CertificateError::PortColorInvalid { factor, word } => {
                write!(f, "{factor} factor word {word:?} has a port color outside the shared pair")
            }
            CertificateError::PortLevelMismatch {
                factor,
                word,
                count,
            } => write!(
                f,
                "{factor} factor word {word:?} has count {count}, inconsistent with its port color"
            ),
            CertificateError::PortLevelMissing { factor, color } => {
                write!(f, "{factor} factor has no word with port color {color}")
            }
            CertificateError::PrivateColorOutsideBlock {
                word,
                position,
                color,
            } => write!(f, "word {word} has private color {color} at position {position} outside the block"),
            CertificateError::ActiveCountNotOne { word, count } => {
                write!(f, "word {word} has {count} non-private block positions, expected 1")
            }
            CertificateError::AllowedPositionsKeys { expected, actual } => write!(
                f,
                "allowed-position colors {actual:?} differ from port colors {expected:?}"
            ),
            CertificateError::AllowedPositionsInvalid { color, positions } => {
                write!(f, "allowed positions {positions:?} for color {color} are invalid")
            }
            CertificateError::ReconstructionMismatch {
                pairs,
                family_size,
                missing,
                extra,
            } => write!(
                f,
                "reconstruction produced {pairs} pairs for {family_size} words; missing {missing:?}, extra {extra:?}"
            ),
            CertificateError::UnusedPrivateBlockData => {
                write!(f, "private-block certificate contains unused factor states or allowed positions")
            }
            CertificateError::FactorNotExchange { factor, obstruction } => {
                write!(f, "{factor} factor fails exchange: {obstruction}")
            }
        }
    }
}

impl std::error::Error for CertificateError {}

fn restrict(word: &[usize], positions: &[usize]) -> Vec<usize> {
    positions.iter().map(|&p| word[p]).collect()
}

fn complement(positions: &[usize], n: usize) -> Vec<usize> {
    (0..n).filter(|p| !positions.contains(p)).collect()
}

fn build_factor(
    factor: Factor,
    num_colors: usize,
    words: Vec<Vec<usize>>,
) -> Result<FixedContentFamily, CertificateError> {
    FixedContentFamily::try_new(num_colors, words)
        .map_err(|error| CertificateError::FactorInvalid { factor, error })
}

fn check_factor_shape(
    factor: Factor,
    family: &FixedContentFamily,
    expected_positions: usize,
    expected_colors: usize,
) -> Result<(), CertificateError> {
    if family.num_positions() != expected_positions || family.num_colors() != expected_colors {
        return Err(CertificateError::FactorShape {
            factor,
            expected_positions,
            actual_positions: family.num_positions(),
            expected_colors,
            actual_colors: family.num_colors(),
        });
    }
    Ok(())
}

fn check_factor_quota(
    factor: Factor,
    family: &FixedContentFamily,
    expected: Vec<usize>,
) -> Result<(), CertificateError> {
    if family.quota() != expected.as_slice() {
        return Err(CertificateError::QuotaMismatch {
            factor,
            expected,
            actual: family.quota().to_vec(),
        });
    }
    Ok(())
}

fn check_factor_exchange(
    factor: Factor,
    family: &FixedContentFamily,
) -> Result<(), CertificateError> {
    family
        .check_symmetric_cycle_exchange()
        .map_err(|obstruction| CertificateError::FactorNotExchange {
            factor,
            obstruction,
        })
}

fn content(word: &[usize], num_colors: usize) -> Vec<usize> {
    let mut counts = vec![0usize; num_colors];
    for &c in word {
        counts[c] += 1;
    }
    counts
}

/// Assemble a word of length `n` from values at `left_positions` and at
/// `right_positions`, which must partition `0..n`.
fn interleave(
    left_positions: &[usize],
    left: &[usize],
    right_positions: &[usize],
    right: &[usize],
    n: usize,
) -> Vec<usize> {
    let mut word = vec![0usize; n];
    for (&p, &c) in left_positions.iter().zip(left) {
        word[p] = c;
    }
    for (&p, &c) in right_positions.iter().zip(right) {
        word[p] = c;
    }
    word
}

fn check_private_block_word(
    index: usize,
    word: &[usize],
    block: &[usize],
    outside: &[usize],
    is_private: &[bool],
) -> Result<(), CertificateError> {
    for &q in outside {
        if is_private[word[q]] {
            return Err(CertificateError::PrivateColorOutsideBlock {
                word: index,
                position: q,
                color: word[q],
            });
        }
    }
    let count = block.iter().filter(|&&p| !is_private[word[p]]).count();
    if count != 1 {
        return Err(CertificateError::ActiveCountNotOne { word: index, count });
    }
    Ok(())
}

/// Certificate that a family is the direct product of its restrictions to a
/// position cut and its complement.
///
/// `left` lives on `cut` (increasing order), `right` on the complement. All
/// fields are public so tests can corrupt them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProductCertificate {
    pub cut: Vec<usize>,
    pub left: FixedContentFamily,
    pub right: FixedContentFamily,
}

impl ProductCertificate {
    /// Verify against `family`: valid cut, factor shapes, factor quotas equal
    /// to the restricted content, exact reconstruction (`|L|*|R| = |F|` and
    /// equal word sets), and SCE of both factors independently.
    ///
    /// Direct products of SCE families are SCE (proved), so success implies
    /// SCE of `family`; the family itself is not re-checked. Cost is the
    /// reconstruction plus two factor exchange checks.
    pub fn verify(&self, family: &FixedContentFamily) -> Result<(), CertificateError> {
        family.check_position_subset(&self.cut)?;
        let n = family.num_positions();
        let k = family.num_colors();
        let right_positions = complement(&self.cut, n);
        check_factor_shape(Factor::Left, &self.left, self.cut.len(), k)?;
        check_factor_shape(Factor::Right, &self.right, right_positions.len(), k)?;
        check_factor_quota(
            Factor::Left,
            &self.left,
            content(&restrict(family.word(0), &self.cut), k),
        )?;
        check_factor_quota(
            Factor::Right,
            &self.right,
            content(&restrict(family.word(0), &right_positions), k),
        )?;
        let mut reconstructed = BTreeSet::new();
        let mut pairs = 0usize;
        for l in self.left.words() {
            for r in self.right.words() {
                pairs += 1;
                reconstructed.insert(interleave(&self.cut, l, &right_positions, r, n));
            }
        }
        family.compare_reconstruction(reconstructed, pairs)?;
        check_factor_exchange(Factor::Left, &self.left)?;
        check_factor_exchange(Factor::Right, &self.right)
    }
}

/// Certificate that a family is the binary-port gluing of two factors across
/// a cut whose side palettes share exactly two colors and whose left count
/// of `counted_color` takes exactly the values `low_count, low_count + 1`.
///
/// Each factor has its port as its last position. A left word with left
/// count `low_count` carries port `counted_color`, one with count
/// `low_count + 1` carries `other_color`; right ports are the opposite.
/// Gluing pairs words with different port colors and drops the ports.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BinaryPortCertificate {
    pub cut: Vec<usize>,
    pub counted_color: usize,
    pub other_color: usize,
    pub low_count: usize,
    pub left: FixedContentFamily,
    pub right: FixedContentFamily,
}

impl BinaryPortCertificate {
    /// Verify against `family`: valid cut and colors, side palettes sharing
    /// exactly the two colors, factor shapes, combined quota
    /// (`left + right = family + e_a + e_b`), port colors and count levels
    /// for every factor word, both levels present in both factors, exact
    /// reconstruction (compatible pair count `= |F|` and equal word sets),
    /// and SCE of both factors independently.
    ///
    /// Gluing two SCE factors sharing exactly two colors gives an SCE family
    /// (proved), so success implies SCE of `family`; the family itself is
    /// not re-checked.
    pub fn verify(&self, family: &FixedContentFamily) -> Result<(), CertificateError> {
        family.check_position_subset(&self.cut)?;
        let n = family.num_positions();
        let k = family.num_colors();
        let (a, b) = (self.counted_color, self.other_color);
        if a == b || a >= k || b >= k {
            return Err(CertificateError::InvalidColors {
                colors: vec![a, b],
                num_colors: k,
            });
        }
        let right_positions = complement(&self.cut, n);
        let left_palette = family.palette(&self.cut);
        let right_palette = family.palette(&right_positions);
        let shared: Vec<usize> = left_palette
            .iter()
            .filter(|c| right_palette.contains(c))
            .copied()
            .collect();
        if shared != [a.min(b), a.max(b)] {
            return Err(CertificateError::SharedPaletteMismatch { shared });
        }
        check_factor_shape(Factor::Left, &self.left, self.cut.len() + 1, k)?;
        check_factor_shape(Factor::Right, &self.right, right_positions.len() + 1, k)?;
        let mut expected_total = family.quota().to_vec();
        expected_total[a] += 1;
        expected_total[b] += 1;
        let actual_total: Vec<usize> = self
            .left
            .quota()
            .iter()
            .zip(self.right.quota())
            .map(|(x, y)| x + y)
            .collect();
        if actual_total != expected_total {
            return Err(CertificateError::QuotaMismatch {
                factor: Factor::Left,
                expected: expected_total,
                actual: actual_total,
            });
        }
        let quota_a = family.quota()[a];
        let high_count =
            self.low_count
                .checked_add(1)
                .ok_or_else(|| CertificateError::CountLevels {
                    color: a,
                    values: vec![self.low_count],
                })?;
        if high_count > quota_a {
            return Err(CertificateError::CountLevels {
                color: a,
                values: vec![self.low_count, high_count],
            });
        }
        self.check_levels(Factor::Left, &self.left, a, self.low_count, high_count)?;
        self.check_levels(
            Factor::Right,
            &self.right,
            a,
            quota_a - self.low_count - 1,
            quota_a - self.low_count,
        )?;
        let mut reconstructed = BTreeSet::new();
        let mut pairs = 0usize;
        for l in self.left.words() {
            let (l_body, l_port) = l.split_at(self.cut.len());
            for r in self.right.words() {
                let (r_body, r_port) = r.split_at(right_positions.len());
                if l_port[0] != r_port[0] {
                    pairs += 1;
                    reconstructed.insert(interleave(
                        &self.cut,
                        l_body,
                        &right_positions,
                        r_body,
                        n,
                    ));
                }
            }
        }
        family.compare_reconstruction(reconstructed, pairs)?;
        check_factor_exchange(Factor::Left, &self.left)?;
        check_factor_exchange(Factor::Right, &self.right)
    }

    /// Every word's non-port `a`-count must be `count_at_a` when its port is
    /// `a` and `count_at_b` when its port is `b`; both ports must occur.
    fn check_levels(
        &self,
        factor: Factor,
        side: &FixedContentFamily,
        a: usize,
        count_at_a: usize,
        count_at_b: usize,
    ) -> Result<(), CertificateError> {
        let b = self.other_color;
        let body_len = side.num_positions() - 1;
        let (mut saw_a, mut saw_b) = (false, false);
        for w in side.words() {
            let count = w[..body_len].iter().filter(|&&c| c == a).count();
            let port = w[body_len];
            let expected = if port == a {
                saw_a = true;
                count_at_a
            } else if port == b {
                saw_b = true;
                count_at_b
            } else {
                return Err(CertificateError::PortColorInvalid {
                    factor,
                    word: w.clone(),
                });
            };
            if count != expected {
                return Err(CertificateError::PortLevelMismatch {
                    factor,
                    word: w.clone(),
                    count,
                });
            }
        }
        for (seen, color) in [(saw_a, a), (saw_b, b)] {
            if !seen {
                return Err(CertificateError::PortLevelMissing { factor, color });
            }
        }
        Ok(())
    }
}

/// Certificate that a family is a singleton-color substitution with allowed
/// active positions: a block `P` of positions carries all occurrences of the
/// `private_colors`, and every word has exactly one non-private (active)
/// position in `P`.
///
/// `compressed` is the family outside the block plus a final port position
/// carrying the transmitted color (`num_colors` colors, private colors of
/// zero quota). `internal` is the family of block restrictions with the
/// active position recolored by the fresh color `alpha = num_colors`
/// (`num_colors + 1` colors). `allowed_positions` maps each transmitted
/// color `d` occurring at the port to the sorted block-relative active
/// positions permitted with `d`. The family equals
/// `{ g[h; port] : g in compressed, h in internal, pos_alpha(h) in J_{g(port)} }`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PrivateBlockCertificate {
    pub block: Vec<usize>,
    pub private_colors: Vec<usize>,
    pub compressed: FixedContentFamily,
    pub internal: FixedContentFamily,
    pub allowed_positions: BTreeMap<usize, Vec<usize>>,
}

impl PrivateBlockCertificate {
    /// Fresh color index used for the active position inside `internal`.
    pub fn alpha(&self) -> usize {
        self.compressed.num_colors()
    }

    /// Verify against `family`: valid block and private colors, private
    /// colors confined to the block, exactly one active block position per
    /// word, factor shapes and quotas (compressed keeps the non-private
    /// quotas with private quotas zero; internal keeps the private quotas
    /// with `alpha` of quota one), allowed-position keys equal to the port
    /// colors with valid block-relative lists, exact reconstruction (triple
    /// count `= |F|` and equal word sets), no unused factor states or allowed
    /// positions, and SCE of both factors independently.
    ///
    /// Singleton-color substitution followed by allowed-cell restriction of
    /// SCE factors gives an SCE family (proved), so success implies SCE of
    /// `family`; the family itself is not re-checked.
    pub fn verify(&self, family: &FixedContentFamily) -> Result<(), CertificateError> {
        family.check_position_subset(&self.block)?;
        let is_private = family.private_color_mask(&self.private_colors)?;
        let n = family.num_positions();
        let k = family.num_colors();
        let outside = complement(&self.block, n);
        for (index, w) in family.words().iter().enumerate() {
            check_private_block_word(index, w, &self.block, &outside, &is_private)?;
        }
        let alpha = k;
        check_factor_shape(Factor::Compressed, &self.compressed, outside.len() + 1, k)?;
        check_factor_shape(Factor::Internal, &self.internal, self.block.len(), k + 1)?;
        let mut compressed_quota = family.quota().to_vec();
        let mut internal_quota = vec![0usize; k + 1];
        for c in 0..k {
            if is_private[c] {
                internal_quota[c] = compressed_quota[c];
                compressed_quota[c] = 0;
            }
        }
        internal_quota[alpha] = 1;
        check_factor_quota(Factor::Compressed, &self.compressed, compressed_quota)?;
        check_factor_quota(Factor::Internal, &self.internal, internal_quota)?;
        let port_colors: BTreeSet<usize> = self
            .compressed
            .words()
            .iter()
            .map(|g| g[outside.len()])
            .collect();
        let keys: Vec<usize> = self.allowed_positions.keys().copied().collect();
        let expected_keys: Vec<usize> = port_colors.into_iter().collect();
        if keys != expected_keys {
            return Err(CertificateError::AllowedPositionsKeys {
                expected: expected_keys,
                actual: keys,
            });
        }
        for (&color, positions) in &self.allowed_positions {
            let increasing = positions.windows(2).all(|w| w[0] < w[1]);
            if positions.is_empty()
                || !increasing
                || positions.iter().any(|&j| j >= self.block.len())
            {
                return Err(CertificateError::AllowedPositionsInvalid {
                    color,
                    positions: positions.clone(),
                });
            }
        }
        let mut reconstructed = BTreeSet::new();
        let mut pairs = 0usize;
        for g in self.compressed.words() {
            let (g_body, g_port) = g.split_at(outside.len());
            let d = g_port[0];
            let allowed = &self.allowed_positions[&d];
            for h in self.internal.words() {
                let active = h
                    .iter()
                    .position(|&c| c == alpha)
                    .expect("alpha has quota one");
                if allowed.binary_search(&active).is_err() {
                    continue;
                }
                pairs += 1;
                let mut block_word = h.clone();
                block_word[active] = d;
                reconstructed.insert(interleave(&outside, g_body, &self.block, &block_word, n));
            }
        }
        family.compare_reconstruction(reconstructed, pairs)?;
        let active_positions: BTreeSet<usize> = self
            .internal
            .words()
            .iter()
            .map(|h| {
                h.iter()
                    .position(|&c| c == alpha)
                    .expect("alpha has quota one")
            })
            .collect();
        let listed_positions: BTreeSet<usize> =
            self.allowed_positions.values().flatten().copied().collect();
        if active_positions != listed_positions {
            return Err(CertificateError::UnusedPrivateBlockData);
        }
        // Every nonempty J_d meets the internal support now, and every
        // internal state is used for some port color. Reconstruction then
        // forces both factors and all J_d to be the actual restrictions.
        check_factor_exchange(Factor::Compressed, &self.compressed)?;
        check_factor_exchange(Factor::Internal, &self.internal)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn family(letters: &[&str]) -> FixedContentFamily {
        let k = letters
            .iter()
            .flat_map(|w| w.bytes())
            .map(|b| (b - b'A') as usize + 1)
            .max()
            .unwrap_or(0);
        FixedContentFamily::from_letter_words(k, letters).unwrap()
    }

    #[test]
    fn validation_rejects_bad_inputs() {
        assert_eq!(
            FixedContentFamily::try_new(2, vec![]).unwrap_err(),
            FamilyError::Empty
        );
        assert!(matches!(
            FixedContentFamily::try_new(2, vec![vec![0, 1], vec![0]]).unwrap_err(),
            FamilyError::RaggedWord { index: 1, .. }
        ));
        assert!(matches!(
            FixedContentFamily::try_new(2, vec![vec![0, 2]]).unwrap_err(),
            FamilyError::ColorOutOfRange {
                index: 0,
                position: 1,
                color: 2,
                ..
            }
        ));
        assert!(matches!(
            FixedContentFamily::try_new(2, vec![vec![0, 1], vec![1, 1]]).unwrap_err(),
            FamilyError::ContentMismatch { index: 1, .. }
        ));
        let dup = FixedContentFamily::try_new(2, vec![vec![1, 0], vec![0, 1], vec![1, 0]]).unwrap();
        assert_eq!(dup.words(), &[vec![0, 1], vec![1, 0]]);
        assert_eq!(dup.index_of(&[1, 0]), Some(1));
    }

    #[test]
    fn zero_quota_and_empty_word_are_allowed() {
        let f = FixedContentFamily::try_new(3, vec![vec![]]).unwrap();
        assert_eq!(f.quota(), &[0, 0, 0]);
        assert_eq!(f.num_positions(), 0);
        assert_eq!(f.check_symmetric_cycle_exchange(), Ok(()));
        let marginal = f.color_fiber_counts(1).unwrap();
        assert_eq!(marginal.len(), 1);
        assert_eq!(marginal[&Vec::<usize>::new()], BigUint::one());
        let g = FixedContentFamily::all_assignments(&[0, 2, 0]);
        assert_eq!(g.words(), &[vec![1, 1]]);
        assert_eq!(
            FixedContentFamily::all_assignments(&[]).words(),
            &[Vec::<usize>::new()]
        );
    }

    #[test]
    fn cyclic_three_word_family_exchanges_by_full_cycles() {
        let f = family(&["ABC", "BCA", "CAB"]);
        assert_eq!(f.check_symmetric_cycle_exchange(), Ok(()));
        let w = f.exchange_witness(0, 1, 0).unwrap().unwrap();
        assert_eq!(w.positions.len(), 3);
        assert_eq!(f.verify_exchange(0, 1, 0, &w), Ok(()));
        assert!(!f.is_list_quota_saturated());
        assert!(FixedContentFamily::all_assignments(&[1, 1, 1]).is_list_quota_saturated());
    }

    #[test]
    fn two_color_case_is_matroid_exchange() {
        let bad = family(&["ABAB", "BABA"]);
        let obstruction = bad.check_symmetric_cycle_exchange().unwrap_err();
        assert_eq!(
            (obstruction.first, obstruction.second, obstruction.position),
            (0, 1, 0)
        );
        let good = family(&["ABAB", "ABBA", "BAAB", "BABA"]);
        assert_eq!(good.check_symmetric_cycle_exchange(), Ok(()));
    }

    #[test]
    fn query_errors_are_reported() {
        let f = family(&["AB", "BA"]);
        assert_eq!(
            f.exchange_witness(0, 2, 0).unwrap_err(),
            InvalidQuery::WordIndexOutOfRange { index: 2, len: 2 }
        );
        assert_eq!(
            f.exchange_witness(0, 1, 2).unwrap_err(),
            InvalidQuery::PositionOutOfRange {
                position: 2,
                num_positions: 2
            }
        );
        assert_eq!(
            f.exchange_witness(1, 1, 0).unwrap_err(),
            InvalidQuery::WordsAgree {
                first: 1,
                second: 1,
                position: 0
            }
        );
        assert_eq!(
            simple_exchange_cycles(&[0, 1], &[0, 1, 1], 0).unwrap_err(),
            InvalidQuery::LengthMismatch {
                first: 2,
                second: 3
            }
        );
    }

    #[test]
    fn budget_exhaustion_is_reported_not_failure() {
        let f = family(&["ABC", "BCA", "CAB"]);
        assert_eq!(
            f.exchange_witness_bounded(0, 1, 0, 1).unwrap(),
            WitnessSearch::BudgetExhausted { steps: 1 }
        );
        assert!(matches!(
            f.check_symmetric_cycle_exchange_bounded(1),
            ExchangeStatus::BudgetExhausted {
                first: 0,
                second: 1,
                position: 0,
                ..
            }
        ));
        assert_eq!(
            f.check_symmetric_cycle_exchange_bounded(1_000),
            ExchangeStatus::Satisfied
        );
    }

    #[test]
    fn simple_cycles_are_enumerated_deterministically() {
        // A->B at 0, B->A at 1 and 2: two 2-cycles through position 0.
        let cycles = simple_exchange_cycles(&[0, 1, 1], &[1, 0, 0], 0).unwrap();
        assert_eq!(
            cycles,
            vec![
                CycleWitness {
                    positions: vec![0, 1]
                },
                CycleWitness {
                    positions: vec![0, 2]
                }
            ]
        );
    }

    #[test]
    fn marginal_matches_direct_formula() {
        let f = family(&["ABC", "BCA", "CAB"]);
        let w = |p: usize, c: usize| -> BigUint { BigUint::from(10usize.pow((3 * p + c) as u32)) };
        let weights: Vec<Vec<BigUint>> =
            (0..3).map(|p| (0..3).map(|c| w(p, c)).collect()).collect();
        let marginal = f.weighted_color_marginal(0, &weights).unwrap();
        // ABC: A at 0, weight w(1,B) w(2,C); BCA: A at 2, weight w(0,B) w(1,C); CAB: A at 1, w(0,C) w(2,B).
        assert_eq!(marginal[&vec![0]], w(1, 1) * w(2, 2));
        assert_eq!(marginal[&vec![2]], w(0, 1) * w(1, 2));
        assert_eq!(marginal[&vec![1]], w(0, 2) * w(2, 1));
        assert_eq!(
            f.weighted_color_marginal(0, &weights[..2]).unwrap_err(),
            MarginalError::WeightRows {
                expected: 3,
                actual: 2
            }
        );
        assert_eq!(
            f.weighted_quadratic_marginal_matrix(0, &weights)
                .unwrap_err(),
            MarginalError::QuotaNotTwo { color: 0, quota: 1 }
        );
    }
}
