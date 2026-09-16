use crate::partition::Partition;
use std::collections::BTreeSet;
use std::fmt;
use std::hash::Hash;

// =========================================================================
// WeakComposition
// =========================================================================

/// A weak integer composition: an ordered tuple of non-negative integers
/// with a fixed number of parts (including trailing zeros).
///
/// Unlike [`Composition`], trailing zeros are **not** stripped — the length
/// encodes the number of variables, which is structural for polynomial bases
/// (key, atom, slide, etc.).
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct WeakComposition(Vec<u32>);

impl WeakComposition {
    /// Create a weak composition from parts. Length is preserved exactly.
    pub fn new(parts: Vec<u32>) -> Self {
        WeakComposition(parts)
    }

    /// Create from a slice.
    pub fn from_slice(s: &[u32]) -> Self {
        WeakComposition(s.to_vec())
    }

    /// The zero weak composition with `k` parts (all zeros).
    pub fn zero(k: usize) -> Self {
        WeakComposition(vec![0; k])
    }

    /// The parts as a slice.
    pub fn parts(&self) -> &[u32] {
        &self.0
    }

    /// Number of parts (= number of variables).
    pub fn num_vars(&self) -> usize {
        self.0.len()
    }

    /// Sum of all parts.
    pub fn degree(&self) -> u32 {
        checked_part_sum(&self.0)
    }

    /// Whether all parts are zero.
    pub fn is_zero(&self) -> bool {
        self.0.iter().all(|&p| p == 0)
    }

    /// Whether the parts are weakly decreasing (dominant composition).
    pub fn is_dominant(&self) -> bool {
        self.0.windows(2).all(|w| w[0] >= w[1])
    }

    /// Remove zero parts to get a (strong) composition.
    pub fn flat(&self) -> Composition {
        let parts: Vec<u32> = self.0.iter().copied().filter(|&p| p > 0).collect();
        Composition::new(parts)
    }

    /// Sort non-zero parts into a partition (the "content" of the composition).
    pub fn content(&self) -> Partition {
        let mut parts: Vec<u32> = self.0.iter().copied().filter(|&p| p > 0).collect();
        parts.sort_unstable_by(|a, b| b.cmp(a));
        Partition::from_sorted(parts)
    }

    /// Convert to a `Composition` by stripping trailing zeros.
    pub fn to_composition(&self) -> Composition {
        Composition::new(self.0.clone())
    }

    /// Create from a `Composition` by padding with zeros to `num_vars` parts.
    /// Panics if `num_vars < c.num_parts()`.
    pub fn from_composition(c: &Composition, num_vars: usize) -> Self {
        assert!(
            num_vars >= c.num_parts(),
            "num_vars ({}) < composition length ({})",
            num_vars,
            c.num_parts()
        );
        let mut parts = c.parts().to_vec();
        parts.resize(num_vars, 0);
        WeakComposition(parts)
    }

    /// All weak compositions of `d` into exactly `k` parts (parts >= 0).
    pub fn all_weak_compositions(d: u32, k: usize) -> Vec<Self> {
        if k == 0 {
            return if d == 0 {
                vec![WeakComposition(vec![])]
            } else {
                vec![]
            };
        }
        let mut results = Vec::new();
        Self::enumerate_helper(d, k, &mut Vec::with_capacity(k), &mut results);
        results
    }

    fn enumerate_helper(
        remaining: u32,
        parts_left: usize,
        current: &mut Vec<u32>,
        results: &mut Vec<WeakComposition>,
    ) {
        if parts_left == 1 {
            let mut c = current.clone();
            c.push(remaining);
            results.push(WeakComposition(c));
            return;
        }
        for part in 0..=remaining {
            current.push(part);
            Self::enumerate_helper(remaining - part, parts_left - 1, current, results);
            current.pop();
        }
    }
}

impl fmt::Display for WeakComposition {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s: Vec<String> = self.0.iter().map(|p| p.to_string()).collect();
        write!(f, "({})", s.join(","))
    }
}

/// An integer composition: an ordered tuple of positive integers.
///
/// Unlike partitions, the order of parts matters. The legacy [`Self::new`]
/// constructor and [`Self::weak_integer_compositions`] can retain zero parts
/// for compatibility with older weak-composition consumers. Use
/// [`Self::try_new`] when a mathematically strong composition is required.
/// Matches `IntegerCompositions` etc. from CombinatoricTools.m.
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Composition(Vec<u32>);

/// An invalid strong-composition or descent-set operation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CompositionError {
    /// A strong composition cannot contain a zero part.
    ZeroPart { index: usize },
    /// A descent must lie strictly between zero and the total size.
    DescentOutOfRange { descent: u32, size: u32 },
    /// The sum of parts cannot be represented by `u32`.
    SizeOverflow,
    /// A ribbon width cannot be represented by `u32`.
    RibbonWidthOverflow,
}

impl Composition {
    /// Create a legacy composition from parts, stripping trailing zeros.
    ///
    /// This constructor preserves interior zero parts for compatibility with
    /// existing weak-composition consumers. Prefer [`Self::try_new`] for an
    /// all-positive composition.
    pub fn new(parts: Vec<u32>) -> Self {
        let mut p = parts;
        while p.last() == Some(&0) {
            p.pop();
        }
        Composition(p)
    }

    /// Create a mathematically strong composition, rejecting every zero part.
    pub fn try_new(parts: Vec<u32>) -> Result<Self, CompositionError> {
        validate_positive_parts(&parts)?;
        Ok(Composition(parts))
    }

    /// The empty composition.
    pub fn empty() -> Self {
        Composition(vec![])
    }

    /// The parts as a slice, in the given order.
    pub fn parts(&self) -> &[u32] {
        &self.0
    }

    /// Number of parts.
    pub fn num_parts(&self) -> usize {
        self.0.len()
    }

    /// Sum of all parts.
    pub fn size(&self) -> u32 {
        checked_part_sum(&self.0)
    }

    /// Whether every part is positive.
    pub fn is_strong(&self) -> bool {
        self.0.iter().all(|&part| part > 0)
    }

    /// Whether this is the empty composition.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Sort into a partition.
    pub fn to_partition(&self) -> Partition {
        Partition::new(self.0.clone())
    }

    // -----------------------------------------------------------------------
    // Enumeration — IntegerCompositions from CombinatoricTools.m
    // -----------------------------------------------------------------------

    /// All compositions of n (with all positive parts).
    /// `IntegerCompositions[n]` in Mathematica.
    pub fn integer_compositions(n: u32) -> Vec<Composition> {
        if n == 0 {
            return vec![Composition::empty()];
        }
        let mut results = Vec::new();
        Self::compositions_helper(n, &mut vec![], &mut results);
        results
    }

    fn compositions_helper(remaining: u32, current: &mut Vec<u32>, results: &mut Vec<Composition>) {
        if remaining == 0 {
            results.push(Composition(current.clone()));
            return;
        }
        for part in 1..=remaining {
            current.push(part);
            Self::compositions_helper(remaining - part, current, results);
            current.pop();
        }
    }

    /// All compositions of n with exactly k parts.
    /// `IntegerCompositions[n, k]` in Mathematica.
    pub fn integer_compositions_k(n: u32, k: usize) -> Vec<Composition> {
        if k == 0 {
            return if n == 0 {
                vec![Composition::empty()]
            } else {
                vec![]
            };
        }
        let mut results = Vec::new();
        Self::compositions_k_helper(n, k, &mut vec![], &mut results);
        results
    }

    fn compositions_k_helper(
        remaining: u32,
        parts_left: usize,
        current: &mut Vec<u32>,
        results: &mut Vec<Composition>,
    ) {
        if parts_left == 0 {
            if remaining == 0 {
                results.push(Composition(current.clone()));
            }
            return;
        }
        let max_part = if parts_left > 1 {
            let required_after_first =
                u32::try_from(parts_left - 1).expect("composition length overflow");
            if remaining < required_after_first {
                return;
            }
            remaining - required_after_first
        } else {
            remaining
        };
        for part in 1..=max_part {
            current.push(part);
            Self::compositions_k_helper(remaining - part, parts_left - 1, current, results);
            current.pop();
        }
    }

    /// All weak compositions of n with exactly k parts (parts >= 0).
    ///
    /// This is a compatibility API and returns `Composition` values that can
    /// contain zero parts. New code that needs a fixed-length weak composition
    /// should use [`WeakComposition::all_weak_compositions`].
    /// `WeakIntegerCompositions[n, k]` in Mathematica.
    pub fn weak_integer_compositions(n: u32, k: usize) -> Vec<Composition> {
        if k == 0 {
            return if n == 0 {
                vec![Composition::empty()]
            } else {
                vec![]
            };
        }
        let mut results = Vec::new();
        Self::weak_compositions_helper(n, k, &mut vec![], &mut results);
        results
    }

    fn weak_compositions_helper(
        remaining: u32,
        parts_left: usize,
        current: &mut Vec<u32>,
        results: &mut Vec<Composition>,
    ) {
        if parts_left == 1 {
            let mut c = current.clone();
            c.push(remaining);
            results.push(Composition(c));
            return;
        }
        for part in 0..=remaining {
            current.push(part);
            Self::weak_compositions_helper(remaining - part, parts_left - 1, current, results);
            current.pop();
        }
    }

    // -----------------------------------------------------------------------
    // Descent set / ribbon conversions
    // -----------------------------------------------------------------------

    /// Convert a strong composition to its descent set, excluding the total.
    /// `CompositionToDescentSet` in Mathematica.
    pub fn composition_to_descent_set(&self) -> BTreeSet<u32> {
        self.try_composition_to_descent_set()
            .expect("descent sets require a strong composition")
    }

    /// Convert a strong composition to its descent set, excluding the total.
    pub fn try_composition_to_descent_set(&self) -> Result<BTreeSet<u32>, CompositionError> {
        validate_positive_parts(&self.0)?;
        let mut set = BTreeSet::new();
        let mut sum = 0u32;
        for (i, &part) in self.0.iter().enumerate() {
            sum = sum
                .checked_add(part)
                .ok_or(CompositionError::SizeOverflow)?;
            if i < self.0.len() - 1 {
                set.insert(sum);
            }
        }
        Ok(set)
    }

    /// Reconstruct composition from descent set and total n.
    /// `DescentSetToComposition` in Mathematica.
    pub fn from_descent_set(descent_set: &BTreeSet<u32>, n: u32) -> Composition {
        Self::try_from_descent_set(descent_set, n)
            .expect("descent set entries must lie strictly between zero and the total")
    }

    /// Reconstruct a strong composition from descents strictly between `0` and `n`.
    pub fn try_from_descent_set(
        descent_set: &BTreeSet<u32>,
        n: u32,
    ) -> Result<Composition, CompositionError> {
        if n == 0 {
            return if descent_set.is_empty() {
                Ok(Composition::empty())
            } else {
                Err(CompositionError::DescentOutOfRange {
                    descent: *descent_set.first().expect("set is nonempty"),
                    size: n,
                })
            };
        }
        let mut parts = Vec::new();
        let mut prev = 0u32;
        for &d in descent_set {
            if d == 0 || d >= n {
                return Err(CompositionError::DescentOutOfRange {
                    descent: d,
                    size: n,
                });
            }
            parts.push(d - prev);
            prev = d;
        }
        parts.push(n - prev);
        Self::try_new(parts)
    }

    /// Convert a strong composition to a ribbon (border strip) as a skew shape.
    /// `CompositionToRibbon` in Mathematica.
    pub fn composition_to_ribbon(&self) -> (Partition, Partition) {
        self.try_composition_to_ribbon()
            .expect("ribbons require a strong composition with representable widths")
    }

    /// Convert a strong composition to a ribbon, reporting invalid parts or overflow.
    pub fn try_composition_to_ribbon(&self) -> Result<(Partition, Partition), CompositionError> {
        validate_positive_parts(&self.0)?;
        if self.is_empty() {
            return Ok((Partition::empty(), Partition::empty()));
        }
        let k = self.num_parts();
        let mut outer_parts = Vec::new();
        let mut inner_parts = Vec::new();

        let reversed: Vec<u32> = self.0.iter().rev().copied().collect();
        let mut offset = 0u32;

        for i in 0..k {
            let width = reversed[i];
            let right = offset
                .checked_add(width)
                .ok_or(CompositionError::RibbonWidthOverflow)?;
            outer_parts.push(right);
            inner_parts.push(offset);
            if i < k - 1 {
                offset = right
                    .checked_sub(1)
                    .ok_or(CompositionError::RibbonWidthOverflow)?;
            }
        }

        outer_parts.reverse();
        inner_parts.reverse();

        Ok((
            Partition::from_sorted(outer_parts),
            Partition::from_sorted(inner_parts),
        ))
    }

    /// All refinements of a composition (split each part into a composition).
    /// `CompositionRefinements` in Mathematica.
    pub fn composition_refinements(&self) -> Vec<Composition> {
        if self.is_empty() {
            return vec![Composition::empty()];
        }
        let part_compositions: Vec<Vec<Vec<u32>>> = self
            .0
            .iter()
            .map(|&p| {
                Composition::integer_compositions(p)
                    .into_iter()
                    .map(|c| c.0)
                    .collect()
            })
            .collect();

        let mut results = vec![vec![]];
        for choices in &part_compositions {
            let mut new_results = Vec::new();
            for existing in &results {
                for choice in choices {
                    let mut combined = existing.clone();
                    combined.extend(choice);
                    new_results.push(combined);
                }
            }
            results = new_results;
        }

        results.into_iter().map(Composition).collect()
    }

    /// Binary word encoding of a strong composition: 1 at part boundaries,
    /// 0 elsewhere.
    /// `CompositionWord` in Mathematica.
    pub fn composition_word(&self) -> Vec<u8> {
        self.try_composition_word()
            .expect("composition words require a strong composition")
    }

    /// Return the binary-word encoding of a strong composition.
    pub fn try_composition_word(&self) -> Result<Vec<u8>, CompositionError> {
        validate_positive_parts(&self.0)?;
        if self.is_empty() {
            return Ok(vec![]);
        }
        let n = self
            .0
            .iter()
            .try_fold(0u32, |total, &part| total.checked_add(part))
            .ok_or(CompositionError::SizeOverflow)? as usize;
        let mut word = vec![0u8; n - 1];
        let mut pos = 0usize;
        for (i, &part) in self.0.iter().enumerate() {
            pos += part as usize;
            if i < self.0.len() - 1 && pos > 0 && pos <= n - 1 {
                word[pos - 1] = 1;
            }
        }
        Ok(word)
    }

    /// Reconstruct composition from binary word.
    /// `WordComposition` in Mathematica.
    pub fn from_composition_word(word: &[u8]) -> Composition {
        if word.is_empty() {
            return Composition(vec![1]);
        }
        let mut parts = Vec::new();
        let mut current = 1u32;
        for &bit in word {
            if bit == 1 {
                parts.push(current);
                current = 1;
            } else {
                current = current
                    .checked_add(1)
                    .expect("composition word length overflow");
            }
        }
        parts.push(current);
        Composition(parts)
    }

    /// Apply the slinky rule (Egge-Loehr-Warrington) to straighten s_α into ±s_λ.
    ///
    /// Returns `Some((partition, sign))` where sign is ±1, or `None` if s_α = 0.
    pub fn composition_slinky(&self) -> Option<(Partition, i64)> {
        if self.is_empty() {
            return Some((Partition::empty(), 1));
        }

        let mut parts: Vec<i64> = self.0.iter().map(|&p| p as i64).collect();
        let mut sign: i64 = 1;

        loop {
            let bad = parts.windows(2).position(|w| w[0] < w[1]);
            match bad {
                None => {
                    let p: Vec<u32> = parts.iter().map(|&x| x as u32).collect();
                    return Some((Partition::from_sorted(p), sign));
                }
                Some(i) => {
                    let a = parts[i];
                    let b = parts[i + 1];
                    if b == a + 1 {
                        return None;
                    }
                    parts[i] = b - 1;
                    parts[i + 1] = a + 1;
                    sign = -sign;

                    if parts[i] <= 0 || parts[i + 1] <= 0 {
                        return None;
                    }
                }
            }
        }
    }
}

fn checked_part_sum(parts: &[u32]) -> u32 {
    parts
        .iter()
        .try_fold(0u32, |total, &part| total.checked_add(part))
        .expect("composition size overflow")
}

fn validate_positive_parts(parts: &[u32]) -> Result<(), CompositionError> {
    parts
        .iter()
        .position(|&part| part == 0)
        .map_or(Ok(()), |index| Err(CompositionError::ZeroPart { index }))
}

impl fmt::Display for Composition {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_empty() {
            write!(f, "∅")
        } else {
            let s: Vec<String> = self.0.iter().map(|p| p.to_string()).collect();
            write!(f, "({})", s.join(","))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compositions_of_4() {
        let comps = Composition::integer_compositions(4);
        assert_eq!(comps.len(), 8);
    }

    #[test]
    fn test_compositions_k() {
        let comps = Composition::integer_compositions_k(5, 3);
        assert_eq!(comps.len(), 6);
        for c in &comps {
            assert_eq!(c.size(), 5);
            assert_eq!(c.num_parts(), 3);
        }
    }

    #[test]
    fn test_compositions_k_impossible_length_is_empty() {
        assert!(Composition::integer_compositions_k(2, 4).is_empty());
    }

    #[test]
    fn test_weak_compositions() {
        let comps = Composition::weak_integer_compositions(3, 2);
        assert_eq!(comps.len(), 4);
    }

    #[test]
    fn test_descent_set_roundtrip() {
        let c = Composition::new(vec![2, 3, 1]);
        let des = c.composition_to_descent_set();
        let n = c.size();
        let c2 = Composition::from_descent_set(&des, n);
        assert_eq!(c, c2);
    }

    #[test]
    fn test_checked_strong_composition_and_descent_contracts() {
        assert_eq!(
            Composition::try_new(vec![1, 0, 2]),
            Err(CompositionError::ZeroPart { index: 1 })
        );
        let legacy_weak = Composition::new(vec![1, 0, 2, 0]);
        assert_eq!(legacy_weak.parts(), &[1, 0, 2]);
        assert!(!legacy_weak.is_strong());
        assert_eq!(
            legacy_weak.try_composition_to_descent_set(),
            Err(CompositionError::ZeroPart { index: 1 })
        );
        assert_eq!(
            legacy_weak.try_composition_word(),
            Err(CompositionError::ZeroPart { index: 1 })
        );

        assert_eq!(
            Composition::try_from_descent_set(&BTreeSet::from([0]), 3),
            Err(CompositionError::DescentOutOfRange {
                descent: 0,
                size: 3
            })
        );
        assert_eq!(
            Composition::try_from_descent_set(&BTreeSet::from([3]), 3),
            Err(CompositionError::DescentOutOfRange {
                descent: 3,
                size: 3
            })
        );
        assert_eq!(
            Composition::try_from_descent_set(&BTreeSet::from([1]), 0),
            Err(CompositionError::DescentOutOfRange {
                descent: 1,
                size: 0
            })
        );
    }

    #[test]
    #[should_panic(expected = "composition size overflow")]
    fn test_composition_size_rejects_overflow() {
        let c = Composition::new(vec![u32::MAX, 1]);

        let _ = c.size();
    }

    #[test]
    fn test_checked_descent_set_rejects_overflow() {
        let c = Composition::new(vec![u32::MAX, 1]);

        assert_eq!(
            c.try_composition_to_descent_set(),
            Err(CompositionError::SizeOverflow)
        );
    }

    #[test]
    fn test_ribbon() {
        let c = Composition::new(vec![2, 3]);
        let (outer, inner) = c.composition_to_ribbon();
        assert_eq!(outer.size() - inner.size(), c.size());
    }

    #[test]
    fn test_checked_ribbon_rejects_overflow() {
        let c = Composition::new(vec![u32::MAX, 2]);

        assert_eq!(
            c.try_composition_to_ribbon(),
            Err(CompositionError::RibbonWidthOverflow)
        );
    }

    #[test]
    fn test_refinements() {
        let c = Composition::new(vec![3]);
        let refs = c.composition_refinements();
        assert_eq!(refs.len(), 4);
    }

    #[test]
    fn test_composition_word_roundtrip() {
        let c = Composition::new(vec![2, 1, 3]);
        let word = c.composition_word();
        let c2 = Composition::from_composition_word(&word);
        assert_eq!(c, c2);
    }

    #[test]
    fn test_slinky_partition() {
        let c = Composition::new(vec![3, 2, 1]);
        let (lam, sign) = c.composition_slinky().unwrap();
        assert_eq!(lam, Partition::new(vec![3, 2, 1]));
        assert_eq!(sign, 1);
    }

    #[test]
    fn test_slinky_zero_12() {
        let c = Composition::new(vec![1, 2]);
        assert!(c.composition_slinky().is_none());
    }

    // WeakComposition tests

    #[test]
    fn test_weak_composition_preserves_trailing_zeros() {
        let wc = WeakComposition::new(vec![2, 0, 1, 0]);
        assert_eq!(wc.num_vars(), 4);
        assert_eq!(wc.parts(), &[2, 0, 1, 0]);
    }

    #[test]
    fn test_weak_composition_degree() {
        let wc = WeakComposition::new(vec![2, 0, 3]);
        assert_eq!(wc.degree(), 5);
    }

    #[test]
    #[should_panic(expected = "composition size overflow")]
    fn test_weak_composition_degree_rejects_overflow() {
        let wc = WeakComposition::new(vec![u32::MAX, 1]);

        let _ = wc.degree();
    }

    #[test]
    fn test_weak_composition_is_dominant() {
        assert!(WeakComposition::new(vec![3, 2, 0]).is_dominant());
        assert!(!WeakComposition::new(vec![1, 2, 0]).is_dominant());
    }

    #[test]
    fn test_weak_composition_flat() {
        let wc = WeakComposition::new(vec![2, 0, 3, 0, 1]);
        assert_eq!(wc.flat(), Composition::new(vec![2, 3, 1]));
    }

    #[test]
    fn test_weak_composition_content() {
        let wc = WeakComposition::new(vec![1, 0, 3, 2]);
        assert_eq!(wc.content(), Partition::new(vec![3, 2, 1]));
    }

    #[test]
    fn test_weak_composition_enumeration_count() {
        // C(d+k-1, k-1) = C(3+3-1, 3-1) = C(5,2) = 10
        assert_eq!(WeakComposition::all_weak_compositions(3, 3).len(), 10);
        // C(0+3-1, 3-1) = C(2,2) = 1
        assert_eq!(WeakComposition::all_weak_compositions(0, 3).len(), 1);
        // C(2+1-1, 1-1) = C(2,0) = 1
        assert_eq!(WeakComposition::all_weak_compositions(2, 1).len(), 1);
        // k=0, d=0 -> 1
        assert_eq!(WeakComposition::all_weak_compositions(0, 0).len(), 1);
        // k=0, d>0 -> 0
        assert_eq!(WeakComposition::all_weak_compositions(3, 0).len(), 0);
    }

    #[test]
    fn test_weak_composition_enumeration_length() {
        for wc in WeakComposition::all_weak_compositions(4, 3) {
            assert_eq!(wc.num_vars(), 3);
            assert_eq!(wc.degree(), 4);
        }
    }

    #[test]
    fn test_weak_composition_from_composition() {
        let c = Composition::new(vec![2, 3, 1]);
        let wc = WeakComposition::from_composition(&c, 5);
        assert_eq!(wc.parts(), &[2, 3, 1, 0, 0]);
    }

    #[test]
    fn test_weak_composition_ordering() {
        let a = WeakComposition::new(vec![0, 3]);
        let b = WeakComposition::new(vec![1, 2]);
        assert!(a < b); // lexicographic
    }
}
