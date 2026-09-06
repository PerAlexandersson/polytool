//! Traced inverse column insertion for `P`-Robinson--Schensted algorithms.
//!
//! This module implements the inverse column step (Algorithm 2 in the
//! Kim--Pylyavskyy `P`-RS correspondence) against an abstract ordered
//! alphabet.  The caller supplies both the strict `P`-order and the ladder
//! predicate; consequently this module contains no graph- or
//! unit-interval-order-specific code.
//!
//! Words and columns use the convention of the algorithm: a word is stored
//! as `(a_m, ..., a_1)` and is scanned from right to left, while a column is
//! stored bottom-to-top and is likewise reversed into active scan order.

use std::collections::BTreeSet;
use std::error::Error;
use std::fmt;

/// Ordered-alphabet data needed by `P`-column insertion.
pub trait PInsertionOrder {
    /// Return whether `left` is strictly below `right` in the `P`-order.
    fn less(&self, left: u32, right: u32) -> bool;

    /// Return whether the distinct values form one `P`-ladder.
    fn is_ladder(&self, values: &[u32]) -> bool;
}

/// The local action taken by one scan block of inverse column insertion.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum InverseColumnEventKind {
    Empty,
    Drag,
    Append,
    OrdinaryBump,
    LadderCopy,
    LadderMove,
}

impl InverseColumnEventKind {
    /// Whether this event is one of the two ladder cases.
    pub fn is_ladder(self) -> bool {
        matches!(self, Self::LadderCopy | Self::LadderMove)
    }
}

/// One maximal scan block in a traced inverse column step.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InverseColumnEvent {
    /// First one-based scan position in the block.
    pub first_position: usize,
    /// Last one-based scan position in the block.
    pub last_position: usize,
    pub kind: InverseColumnEventKind,
    /// Nonempty source letters read by this block, in scan order.
    pub source: Vec<u32>,
    /// Active chain in scan order immediately before the block.
    pub active_before: Vec<u32>,
    /// Active chain in scan order immediately after the block.
    pub active_after: Vec<u32>,
    /// Nonempty remainder letters emitted by the block, in scan order.
    pub emitted: Vec<u32>,
}

impl InverseColumnEvent {
    /// Whether this block contains the given one-based scan position.
    pub fn contains_position(&self, position: usize) -> bool {
        self.first_position <= position && position <= self.last_position
    }
}

/// Boundary behavior of a distinguished active letter during one scan block.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BoundaryCertificate {
    pub boundary_position: usize,
    pub event_index: usize,
    pub event_kind: InverseColumnEventKind,
    pub event_first_position: usize,
    pub event_last_position: usize,
    pub active_before: bool,
    pub active_after: bool,
    pub emitted: bool,
}

/// Full result of one traced inverse column step.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InverseColumnStep {
    /// Remaining active column, stored bottom-to-top.
    pub remaining_chain: Vec<u32>,
    /// Reconstructed previous word, stored in `(a_m, ..., a_1)` order.
    pub previous_word: Vec<Option<u32>>,
    /// Maximal scan events, ordered by increasing scan position.
    pub events: Vec<InverseColumnEvent>,
}

impl InverseColumnStep {
    /// Return the unique scan event containing `position`.
    pub fn event_containing(&self, position: usize) -> Option<(usize, &InverseColumnEvent)> {
        self.events
            .iter()
            .enumerate()
            .find(|(_, event)| event.contains_position(position))
    }

    /// Describe how `distinguished` behaves in the block containing a boundary.
    pub fn boundary_certificate(
        &self,
        boundary_position: usize,
        distinguished: u32,
    ) -> Option<BoundaryCertificate> {
        let (event_index, event) = self.event_containing(boundary_position)?;
        Some(BoundaryCertificate {
            boundary_position,
            event_index,
            event_kind: event.kind,
            event_first_position: event.first_position,
            event_last_position: event.last_position,
            active_before: event.active_before.contains(&distinguished),
            active_after: event.active_after.contains(&distinguished),
            emitted: event.emitted.contains(&distinguished),
        })
    }
}

/// Invalid input or an impossible local state in inverse column insertion.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PInsertionError {
    InvalidDragPosition {
        position: usize,
        word_length: usize,
    },
    NoLadder {
        position: usize,
        active: u32,
        source: u32,
    },
    InvalidLadderMove {
        position: usize,
    },
    InvalidComplementValue {
        value: u32,
        maximum: u32,
    },
}

impl fmt::Display for PInsertionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidDragPosition {
                position,
                word_length,
            } => write!(
                formatter,
                "drag position {position} is outside 1..={word_length}"
            ),
            Self::NoLadder {
                position,
                active,
                source,
            } => write!(
                formatter,
                "no ladder at scan position {position} for active {active} and source {source}"
            ),
            Self::InvalidLadderMove { position } => write!(
                formatter,
                "ladder move at scan position {position} has no valid source cut"
            ),
            Self::InvalidComplementValue { value, maximum } => write!(
                formatter,
                "value {value} is outside the finite alphabet 1..={maximum}"
            ),
        }
    }
}

impl Error for PInsertionError {}

/// Reverse and complement a partial word over the finite alphabet `1..=maximum`.
///
/// Empty positions are fixed.  This is the `hat` involution used for a
/// self-dual finite ordered alphabet.
pub fn reverse_complement(
    values: &[Option<u32>],
    maximum: u32,
) -> Result<Vec<Option<u32>>, PInsertionError> {
    values
        .iter()
        .rev()
        .map(|value| match value {
            None => Ok(None),
            Some(value) if (1..=maximum).contains(value) => Ok(Some(maximum + 1 - value)),
            Some(value) => Err(PInsertionError::InvalidComplementValue {
                value: *value,
                maximum,
            }),
        })
        .collect()
}

/// Run one inverse column step and retain every maximal scan event.
///
/// `drag_positions` are one-based positions in scan order.  The returned
/// column and word are changed back to the stored bottom-to-top/right-to-left
/// convention.
pub fn inverse_column_step_with_trace<O: PInsertionOrder>(
    word: &[Option<u32>],
    chain: &[u32],
    drag_positions: &BTreeSet<usize>,
    order: &O,
) -> Result<InverseColumnStep, PInsertionError> {
    let length = word.len();
    if let Some(&position) = drag_positions
        .iter()
        .find(|&&position| position == 0 || position > length)
    {
        return Err(PInsertionError::InvalidDragPosition {
            position,
            word_length: length,
        });
    }

    let source: Vec<Option<u32>> = word.iter().rev().copied().collect();
    let mut remainder = vec![None; length];
    let mut active: Vec<u32> = chain.iter().rev().copied().collect();
    let mut events = Vec::new();
    let mut offset = 0;

    while offset < length {
        let position = offset + 1;
        let active_before = active.clone();
        let Some(entry) = source[offset] else {
            let kind = if drag_positions.contains(&position) && !active.is_empty() {
                remainder[offset] = Some(active.remove(0));
                InverseColumnEventKind::Drag
            } else {
                InverseColumnEventKind::Empty
            };
            events.push(InverseColumnEvent {
                first_position: position,
                last_position: position,
                kind,
                source: Vec::new(),
                active_before,
                active_after: active.clone(),
                emitted: remainder[offset].into_iter().collect(),
            });
            offset += 1;
            continue;
        };

        let predecessor = active.iter().rposition(|&old| old < entry);
        if predecessor.is_none() || order.less(active[predecessor.unwrap()], entry) {
            let replacement = predecessor.map_or(0, |index| index + 1);
            let kind = if replacement == active.len() {
                active.push(entry);
                InverseColumnEventKind::Append
            } else {
                remainder[offset] = Some(active[replacement]);
                active[replacement] = entry;
                InverseColumnEventKind::OrdinaryBump
            };
            events.push(InverseColumnEvent {
                first_position: position,
                last_position: position,
                kind,
                source: vec![entry],
                active_before,
                active_after: active.clone(),
                emitted: remainder[offset].into_iter().collect(),
            });
            offset += 1;
            continue;
        }

        let row = predecessor.expect("the ordinary case handles no predecessor");
        let mut choices = Vec::new();
        for height in 0..active.len() - row {
            for width in 0..length - offset {
                let source_block = &source[offset..=offset + width];
                let Some(block) = source_block.iter().copied().collect::<Option<Vec<_>>>() else {
                    continue;
                };
                if block.windows(2).any(|pair| pair[0] >= pair[1]) {
                    continue;
                }
                let mut ladder = active[row..=row + height].to_vec();
                ladder.extend(block);
                if order.is_ladder(&ladder) {
                    choices.push((height, width));
                }
            }
        }
        let Some(max_height) = choices.iter().map(|choice| choice.0).max() else {
            return Err(PInsertionError::NoLadder {
                position,
                active: active[row],
                source: entry,
            });
        };
        let width = choices
            .iter()
            .filter(|choice| choice.0 == max_height)
            .map(|choice| choice.1)
            .max()
            .expect("a maximal-height ladder exists");
        let height = max_height;
        let source_block: Vec<u32> = source[offset..=offset + width]
            .iter()
            .map(|value| value.expect("a chosen ladder source block is nonempty"))
            .collect();

        let kind = if source_block[width] < active[row + height] {
            for (block_offset, &value) in source_block.iter().enumerate() {
                remainder[offset + block_offset] = Some(value);
            }
            InverseColumnEventKind::LadderCopy
        } else {
            let old_active = active.clone();
            for index in row..=row + height {
                let Some(first) = source_block
                    .iter()
                    .position(|&value| value > old_active[index])
                else {
                    return Err(PInsertionError::InvalidLadderMove { position });
                };
                let last = if index == row + height {
                    width
                } else {
                    let Some(last) = (0..width)
                        .rev()
                        .find(|&source_index| source_block[source_index] < old_active[index + 1])
                    else {
                        return Err(PInsertionError::InvalidLadderMove { position });
                    };
                    last
                };
                remainder[offset + first] = Some(old_active[index]);
                for source_index in first..last {
                    remainder[offset + source_index + 1] = Some(source_block[source_index]);
                }
                active[index] = source_block[last];
            }
            InverseColumnEventKind::LadderMove
        };

        events.push(InverseColumnEvent {
            first_position: position,
            last_position: position + width,
            kind,
            source: source_block,
            active_before,
            active_after: active.clone(),
            emitted: remainder[offset..=offset + width]
                .iter()
                .filter_map(|&value| value)
                .collect(),
        });
        offset += width + 1;
    }

    active.reverse();
    remainder.reverse();
    Ok(InverseColumnStep {
        remaining_chain: active,
        previous_word: remainder,
        events,
    })
}

/// Run the exact traced implementation and discard its event certificate.
pub fn inverse_column_step<O: PInsertionOrder>(
    word: &[Option<u32>],
    chain: &[u32],
    drag_positions: &BTreeSet<usize>,
    order: &O,
) -> Result<(Vec<u32>, Vec<Option<u32>>), PInsertionError> {
    let result = inverse_column_step_with_trace(word, chain, drag_positions, order)?;
    Ok((result.remaining_chain, result.previous_word))
}

#[cfg(test)]
mod tests {
    use super::*;

    struct PathOrder;

    impl PInsertionOrder for PathOrder {
        fn less(&self, left: u32, right: u32) -> bool {
            left < right && right - left > 1
        }

        fn is_ladder(&self, values: &[u32]) -> bool {
            let mut sorted = values.to_vec();
            sorted.sort_unstable();
            sorted.dedup();
            sorted.len() == values.len() && sorted.windows(2).all(|pair| pair[1] == pair[0] + 1)
        }
    }

    #[test]
    fn reverse_complement_is_an_involution() {
        let word = vec![Some(1), None, Some(3), Some(6)];
        let first = reverse_complement(&word, 6).unwrap();
        assert_eq!(first, vec![Some(1), Some(4), None, Some(6)]);
        assert_eq!(reverse_complement(&first, 6).unwrap(), word);
    }

    #[test]
    fn ordinary_bump_and_append_are_traced() {
        let drag_positions = BTreeSet::new();
        let bump = inverse_column_step_with_trace(&[Some(3)], &[6, 1], &drag_positions, &PathOrder)
            .unwrap();
        assert_eq!(bump.previous_word, vec![Some(6)]);
        assert_eq!(bump.remaining_chain, vec![3, 1]);
        assert_eq!(bump.events[0].kind, InverseColumnEventKind::OrdinaryBump);

        let append =
            inverse_column_step_with_trace(&[Some(5)], &[3, 1], &drag_positions, &PathOrder)
                .unwrap();
        assert_eq!(append.previous_word, vec![None]);
        assert_eq!(append.remaining_chain, vec![5, 3, 1]);
        assert_eq!(append.events[0].kind, InverseColumnEventKind::Append);
    }

    #[test]
    fn drag_removes_the_first_active_letter() {
        let drag_positions = BTreeSet::from([1]);
        let result =
            inverse_column_step_with_trace(&[None], &[3, 1], &drag_positions, &PathOrder).unwrap();
        assert_eq!(result.previous_word, vec![Some(1)]);
        assert_eq!(result.remaining_chain, vec![3]);
        assert_eq!(result.events[0].kind, InverseColumnEventKind::Drag);
    }

    #[test]
    fn ladder_move_preserves_the_active_maximum() {
        let result = inverse_column_step_with_trace(
            &[Some(4), Some(2)],
            &[8, 3, 1],
            &BTreeSet::new(),
            &PathOrder,
        )
        .unwrap();
        assert_eq!(result.previous_word, vec![Some(3), Some(1)]);
        assert_eq!(result.remaining_chain, vec![8, 4, 2]);
        assert_eq!(result.events.len(), 1);
        assert_eq!(result.events[0].kind, InverseColumnEventKind::LadderMove);
        let certificate = result.boundary_certificate(1, 8).unwrap();
        assert!(certificate.active_before);
        assert!(certificate.active_after);
        assert!(!certificate.emitted);
    }

    #[test]
    fn ladder_copy_leaves_the_active_chain_fixed() {
        let result =
            inverse_column_step_with_trace(&[Some(4)], &[10, 5, 3], &BTreeSet::new(), &PathOrder)
                .unwrap();
        assert_eq!(result.previous_word, vec![Some(4)]);
        assert_eq!(result.remaining_chain, vec![10, 5, 3]);
        assert_eq!(result.events[0].kind, InverseColumnEventKind::LadderCopy);
    }

    #[test]
    fn untraced_projection_uses_the_same_implementation() {
        let word = vec![None, Some(4), Some(2)];
        let chain = vec![8, 3, 1];
        let drags = BTreeSet::from([3]);
        let traced = inverse_column_step_with_trace(&word, &chain, &drags, &PathOrder).unwrap();
        let untraced = inverse_column_step(&word, &chain, &drags, &PathOrder).unwrap();
        assert_eq!(untraced, (traced.remaining_chain, traced.previous_word));
    }
}
