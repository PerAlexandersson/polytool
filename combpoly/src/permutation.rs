//! Core permutation operations: enumeration, pattern avoidance, classical maps,
//! decompositions, and special families.

// ============================================================
// Enumeration
// ============================================================

/// Generate all permutations of \[1..n\] in lexicographic order.
///
/// Uses `u8` elements (1-indexed) for compact storage in permutation statistics.
/// For `usize` elements (0-indexed), see `combinatoric_core::all_permutations_zero_indexed`.
pub fn all_permutations(n: u8) -> Vec<Vec<u8>> {
    if n == 0 {
        return vec![vec![]];
    }
    let mut perms = Vec::new();
    let mut current: Vec<u8> = (1..=n).collect();
    loop {
        perms.push(current.clone());
        if !next_permutation(&mut current) {
            break;
        }
    }
    perms
}

fn next_permutation(perm: &mut [u8]) -> bool {
    let n = perm.len();
    if n <= 1 {
        return false;
    }
    let mut i = n - 2;
    loop {
        if perm[i] < perm[i + 1] {
            break;
        }
        if i == 0 {
            return false;
        }
        i -= 1;
    }
    let mut j = n - 1;
    while perm[j] <= perm[i] {
        j -= 1;
    }
    perm.swap(i, j);
    perm[i + 1..].reverse();
    true
}

/// Parse a sequence string like `"312"` or `"3,1,2"` or `"10,1,2"` into a `Vec<u8>`.
pub fn parse_sequence(s: &str) -> Vec<u8> {
    if s.contains(',') {
        s.split(',')
            .map(|x| x.trim().parse().expect("invalid number in sequence"))
            .collect()
    } else {
        s.chars()
            .map(|c| c.to_digit(10).expect("invalid digit in sequence") as u8)
            .collect()
    }
}

// ============================================================
// Pattern avoidance
// ============================================================

/// Check if `perm` contains `pattern` as a classical permutation pattern.
pub fn contains_pattern(perm: &[u8], pattern: &[u8]) -> bool {
    let n = perm.len();
    let k = pattern.len();
    if k > n {
        return false;
    }
    if k == 0 {
        return true;
    }
    search_pattern(perm, pattern, 0, &mut Vec::with_capacity(k))
}

fn search_pattern(perm: &[u8], pattern: &[u8], start: usize, chosen: &mut Vec<u8>) -> bool {
    if chosen.len() == pattern.len() {
        return true;
    }
    let remaining = pattern.len() - chosen.len();
    for i in start..=perm.len() - remaining {
        let pos = chosen.len();
        let mut valid = true;
        for j in 0..pos {
            let pat_cmp = pattern[j].cmp(&pattern[pos]);
            let perm_cmp = chosen[j].cmp(&perm[i]);
            if pat_cmp != perm_cmp {
                valid = false;
                break;
            }
        }
        if valid {
            chosen.push(perm[i]);
            if search_pattern(perm, pattern, i + 1, chosen) {
                chosen.pop();
                return true;
            }
            chosen.pop();
        }
    }
    false
}

/// Constraints for filtered permutation generation.
///
/// All constraints are checked incrementally during backtracking,
/// so impossible branches are pruned early.
#[derive(Default, Clone)]
pub struct PermConstraints {
    /// Patterns to avoid (classical pattern avoidance).
    pub avoiding: Vec<Vec<u8>>,
    /// Arrow patterns to avoid.
    ///
    /// These are checked when a full permutation has been built, since the
    /// arrow condition depends on the inverse Foata cycle-word map.
    pub avoiding_arrow: Vec<ArrowPattern>,
    /// Only derangements (no fixed points).
    pub derangement: bool,
    /// Only involutions (σ² = id).
    pub involution: bool,
    /// Only alternating (up-down) permutations.
    pub alternating: bool,
    /// First element must equal this value.
    pub starts_with: Option<u8>,
    /// Last element must equal this value.
    pub ends_with: Option<u8>,
}

/// Generate all permutations of \[1..n\] satisfying the given constraints,
/// using incremental backtracking with early pruning.
///
/// This is much faster than generating all of S_n and filtering,
/// especially for pattern-avoiding classes (Catalan-sized output vs n!).
pub fn filtered_permutations(n: u8, constraints: &PermConstraints) -> Vec<Vec<u8>> {
    if n == 0 {
        return vec![vec![]];
    }
    let mut result = Vec::new();
    let available: Vec<u8> = (1..=n).collect();
    build_filtered(
        n,
        constraints,
        &mut Vec::with_capacity(n as usize),
        &available,
        &mut result,
    );
    result
}

/// Convenience wrapper: generate all τ-avoiding permutations of \[1..n\].
///
/// For a single length-3 pattern, uses a dedicated O(n × C_n) recursive
/// generator instead of the general backtracker. For longer patterns or
/// multiple patterns, falls back to incremental backtracking with pruning.
pub fn avoiding_permutations(n: u8, patterns: &[Vec<u8>]) -> Vec<Vec<u8>> {
    // Fast path: single length-3 pattern → dedicated generator
    if patterns.len() == 1 && patterns[0].len() == 3 {
        return gen_avoiding_length3(n, &patterns[0]);
    }
    let constraints = PermConstraints {
        avoiding: patterns.to_vec(),
        ..Default::default()
    };
    filtered_permutations(n, &constraints)
}

/// Generate all permutations of \[1..n\] avoiding the given arrow patterns.
pub fn arrow_avoiding_permutations(n: u8, patterns: &[ArrowPattern]) -> Vec<Vec<u8>> {
    let constraints = PermConstraints {
        avoiding_arrow: patterns.to_vec(),
        ..Default::default()
    };
    filtered_permutations(n, &constraints)
}

/// An arrow pattern `(nu; H)`.
///
/// The word `word` is the classical pattern part `nu`.  The arrows are pairs
/// `(source_label, target_label)`.  Labels are ranks in a selected value set
/// `x_1 < ... < x_size`.  A permutation `sigma` contains the pattern if
/// `sigma` contains the word values as a subsequence and the inverse Foata
/// cycle-word preimage `hat_sigma` satisfies
/// `hat_sigma(x_source_label) = x_target_label` for every arrow.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ArrowPattern {
    pub size: u8,
    pub word: Vec<u8>,
    pub arrows: Vec<(u8, u8)>,
}

impl ArrowPattern {
    /// Construct an arrow pattern, panicking if the labels do not form `[size]`.
    pub fn new(size: u8, word: Vec<u8>, arrows: Vec<(u8, u8)>) -> Self {
        let pattern = Self { size, word, arrows };
        pattern.assert_valid();
        pattern
    }

    fn assert_valid(&self) {
        let mut seen = vec![false; self.size as usize + 1];
        for &label in &self.word {
            assert!(
                (1..=self.size).contains(&label),
                "arrow-pattern word label {} is outside 1..={}",
                label,
                self.size
            );
            seen[label as usize] = true;
        }
        for &(source, target) in &self.arrows {
            assert!(
                (1..=self.size).contains(&source),
                "arrow-pattern source label {} is outside 1..={}",
                source,
                self.size
            );
            assert!(
                (1..=self.size).contains(&target),
                "arrow-pattern target label {} is outside 1..={}",
                target,
                self.size
            );
            seen[source as usize] = true;
            seen[target as usize] = true;
        }
        assert!(
            (1..=self.size).all(|label| seen[label as usize]),
            "arrow-pattern labels must be exactly 1..={}",
            self.size
        );
    }
}

/// Parse an arrow pattern such as `"12;1->3"` or `"231;1->4"`.
///
/// The semicolon and arrow list may be omitted, in which case this is just a
/// classical pattern in arrow-pattern form.  For labels above 9, write the word
/// with commas, for example `"10,2,1;1->10"`.
pub fn parse_arrow_pattern(s: &str) -> ArrowPattern {
    let trimmed = s
        .trim()
        .trim_start_matches('(')
        .trim_end_matches(')')
        .trim();
    let (word_part, arrow_part) = trimmed.split_once(';').unwrap_or((trimmed, ""));
    let word = parse_sequence(word_part.trim());
    let arrows = parse_arrow_list(arrow_part);
    let size = word
        .iter()
        .copied()
        .chain(arrows.iter().flat_map(|&(source, target)| [source, target]))
        .max()
        .unwrap_or(0);
    ArrowPattern::new(size, word, arrows)
}

fn parse_arrow_list(s: &str) -> Vec<(u8, u8)> {
    let s = s.trim();
    if s.is_empty() {
        return Vec::new();
    }
    s.split(',')
        .map(|arrow| {
            let compact = arrow.trim().replace(' ', "").replace('→', "->");
            let (source, target) = compact
                .split_once("->")
                .unwrap_or_else(|| panic!("invalid arrow pattern arrow: {}", arrow));
            (
                source.parse().expect("invalid arrow source label"),
                target.parse().expect("invalid arrow target label"),
            )
        })
        .collect()
}

/// Count occurrences of an arrow pattern in a permutation.
///
/// Occurrences are counted by selected value sets `x_1 < ... < x_k` satisfying
/// both the classical word condition and all arrow conditions.
pub fn count_arrow_pattern_occurrences(perm: &[u8], pattern: &ArrowPattern) -> usize {
    pattern.assert_valid();
    let n = perm.len();
    let k = pattern.size as usize;
    if k > n {
        return 0;
    }
    if k == 0 {
        return usize::from(pattern.word.is_empty() && pattern.arrows.is_empty());
    }

    let positions_by_value = inverse(perm);
    let foata_preimage = inverse_foata_cycle_word_map(perm);
    let mut selected_values = Vec::with_capacity(k);
    count_arrow_pattern_occurrences_rec(
        perm,
        &foata_preimage,
        &positions_by_value,
        pattern,
        1,
        &mut selected_values,
    )
}

fn count_arrow_pattern_occurrences_rec(
    perm: &[u8],
    foata_preimage: &[u8],
    positions_by_value: &[u8],
    pattern: &ArrowPattern,
    next_value: u8,
    selected_values: &mut Vec<u8>,
) -> usize {
    if selected_values.len() == pattern.size as usize {
        return usize::from(selected_values_support_arrow_pattern(
            foata_preimage,
            positions_by_value,
            pattern,
            selected_values,
        ));
    }

    let n = perm.len() as u8;
    let remaining_after_this = pattern.size as usize - selected_values.len() - 1;
    let max_value = n - remaining_after_this as u8;
    let mut count = 0;
    for value in next_value..=max_value {
        selected_values.push(value);
        count += count_arrow_pattern_occurrences_rec(
            perm,
            foata_preimage,
            positions_by_value,
            pattern,
            value + 1,
            selected_values,
        );
        selected_values.pop();
    }
    count
}

fn selected_values_support_arrow_pattern(
    foata_preimage: &[u8],
    positions_by_value: &[u8],
    pattern: &ArrowPattern,
    selected_values: &[u8],
) -> bool {
    let mut previous_position = 0u8;
    for &label in &pattern.word {
        let value = selected_values[(label - 1) as usize];
        let position = positions_by_value[(value - 1) as usize];
        if position <= previous_position {
            return false;
        }
        previous_position = position;
    }

    pattern.arrows.iter().all(|&(source, target)| {
        let source_value = selected_values[(source - 1) as usize];
        let target_value = selected_values[(target - 1) as usize];
        foata_preimage[(source_value - 1) as usize] == target_value
    })
}

/// Check whether `perm` contains `pattern` as an arrow pattern.
pub fn contains_arrow_pattern(perm: &[u8], pattern: &ArrowPattern) -> bool {
    count_arrow_pattern_occurrences(perm, pattern) > 0
}

/// Check whether `perm` avoids `pattern` as an arrow pattern.
pub fn avoids_arrow_pattern(perm: &[u8], pattern: &ArrowPattern) -> bool {
    !contains_arrow_pattern(perm, pattern)
}

// ============================================================
// Dedicated generators for length-3 pattern avoidance
// ============================================================

/// Generate Av_τ(n) for a length-3 pattern using dedicated recursions.
///
/// Uses two base generators:
/// - Av_213: recursive decomposition by position of minimum value
/// - Av_321: incremental insertion of maximum value
///
/// The remaining four patterns are obtained by complement and reverse:
/// - Av_231 = complement(Av_213)
/// - Av_312 = reverse(Av_213)
/// - Av_132 = reverse(complement(Av_213))
/// - Av_123 = complement(Av_321)
fn gen_avoiding_length3(n: u8, pattern: &[u8]) -> Vec<Vec<u8>> {
    match (pattern[0], pattern[1], pattern[2]) {
        (2, 1, 3) => gen_av213(n),
        (2, 3, 1) => apply_complement(n, &gen_av213(n)),
        (3, 1, 2) => apply_reverse(&gen_av213(n)),
        (1, 3, 2) => apply_reverse(&apply_complement(n, &gen_av213(n))),
        (3, 2, 1) => gen_av321(n),
        (1, 2, 3) => apply_complement(n, &gen_av321(n)),
        _ => unreachable!("not a valid length-3 pattern"),
    }
}

/// Generate all 213-avoiding permutations of \[1..n\].
///
/// Recursive decomposition: place value 1 at position p.
/// Left block (positions 1..p-1) gets the largest p-1 values;
/// right block (positions p+1..n) gets the next n-p values.
/// Both blocks independently avoid 213.
fn gen_av213(n: u8) -> Vec<Vec<u8>> {
    // Build Av_213(k) for k = 0, 1, ..., n using dynamic programming.
    let mut av: Vec<Vec<Vec<u8>>> = Vec::with_capacity(n as usize + 1);
    av.push(vec![vec![]]); // Av_213(0)

    for k in 1..=n {
        let mut perms = Vec::new();
        for p in 1..=k {
            // Position of value 1 at position p (1-indexed)
            let left_size = (p - 1) as usize;
            let right_size = (k - p) as usize;
            let shift_left = k - p + 1; // map {1..p-1} → {k-p+2..k}
            let shift_right = 1u8; // map {1..n-p} → {2..n-p+1}

            for left in &av[left_size] {
                for right in &av[right_size] {
                    let mut perm = Vec::with_capacity(k as usize);
                    for &v in left {
                        perm.push(v + shift_left);
                    }
                    perm.push(1);
                    for &v in right {
                        perm.push(v + shift_right);
                    }
                    perms.push(perm);
                }
            }
        }
        av.push(perms);
    }
    av.into_iter().last().unwrap()
}

/// Generate all 321-avoiding permutations of \[1..n\].
///
/// Incremental construction: to go from Av_321(k-1) to Av_321(k),
/// insert the value k at any position where the suffix remains increasing
/// (since k is the new maximum, it can only participate in a 321 pattern
/// as the largest value, requiring a descent in the suffix).
fn gen_av321(n: u8) -> Vec<Vec<u8>> {
    let mut prev: Vec<Vec<u8>> = vec![vec![]]; // Av_321(0)

    for k in 1..=n {
        let mut next = Vec::new();
        for sigma in &prev {
            let m = sigma.len();
            // Find valid insertion positions: those where the suffix is increasing.
            // Start from the end (always valid) and scan left while suffix is increasing.
            let mut pos = m; // insert at end
            loop {
                let mut new_perm = Vec::with_capacity(m + 1);
                new_perm.extend_from_slice(&sigma[..pos]);
                new_perm.push(k);
                new_perm.extend_from_slice(&sigma[pos..]);
                next.push(new_perm);

                if pos == 0 {
                    break;
                }
                pos -= 1;
                // Check: is sigma[pos] < sigma[pos+1]? (suffix still increasing?)
                if pos + 1 < m && sigma[pos] > sigma[pos + 1] {
                    break;
                }
            }
        }
        prev = next;
    }
    prev
}

/// Apply complement: replace each value v by n+1-v.
fn apply_complement(n: u8, perms: &[Vec<u8>]) -> Vec<Vec<u8>> {
    let n1 = n + 1;
    perms
        .iter()
        .map(|p| p.iter().map(|&v| n1 - v).collect())
        .collect()
}

/// Apply reverse: reverse each permutation.
fn apply_reverse(perms: &[Vec<u8>]) -> Vec<Vec<u8>> {
    perms
        .iter()
        .map(|p| p.iter().rev().copied().collect())
        .collect()
}

fn build_filtered(
    n: u8,
    constraints: &PermConstraints,
    current: &mut Vec<u8>,
    available: &[u8],
    result: &mut Vec<Vec<u8>>,
) {
    let pos = current.len(); // 0-indexed position being filled

    if pos == n as usize {
        if constraints
            .avoiding_arrow
            .iter()
            .all(|pat| avoids_arrow_pattern(current, pat))
        {
            result.push(current.clone());
        }
        return;
    }

    for i in 0..available.len() {
        let v = available[i];

        // --- Early pruning checks (before pushing) ---

        // starts_with: position 0 must have the specified value
        if pos == 0 {
            if let Some(sw) = constraints.starts_with {
                if v != sw {
                    continue;
                }
            }
        }

        // ends_with: if this is the last position, must match
        if pos == (n as usize) - 1 {
            if let Some(ew) = constraints.ends_with {
                if v != ew {
                    continue;
                }
            }
        }

        // derangement: v must not equal pos+1 (1-indexed)
        if constraints.derangement && v == (pos as u8) + 1 {
            continue;
        }

        // involution: if we place v at position pos+1,
        // then σ(pos+1) = v, so we need σ(v) = pos+1.
        // If v < pos+1, then position v was already filled — check it.
        // If v > pos+1, we'll check later when position v is filled.
        // If v == pos+1, it's a fixed point (always ok for involutions).
        if constraints.involution && (v as usize) < pos + 1 {
            // Position v (0-indexed: v-1) must have value pos+1
            if current[(v as usize) - 1] != (pos as u8) + 1 {
                continue;
            }
        }

        // alternating: check the up-down condition at the new position
        if constraints.alternating && pos >= 1 {
            let prev = *current.last().unwrap();
            if pos % 2 == 1 {
                // odd position (1-indexed: even): should be an ascent (prev < v)
                if prev >= v {
                    continue;
                }
            } else {
                // even position (1-indexed: odd): should be a descent (prev > v)
                if prev <= v {
                    continue;
                }
            }
        }

        current.push(v);

        // Pattern avoidance: check if appending v creates any forbidden pattern
        let ok = constraints
            .avoiding
            .iter()
            .all(|pat| !has_pattern_ending_here(current, pat));

        if ok {
            // For involutions: if v > pos+1, position v must map back to pos+1.
            // Restrict available values at position v when we get there.
            // (This is handled by the check above when we reach position v.)
            let mut remaining = available.to_vec();
            remaining.remove(i);
            build_filtered(n, constraints, current, &remaining, result);
        }
        current.pop();
    }
}

/// Check if `perm` contains `pattern` using the last element of `perm`.
/// This is faster than full containment check during incremental construction.
fn has_pattern_ending_here(perm: &[u8], pattern: &[u8]) -> bool {
    let k = pattern.len();
    let n = perm.len();
    if k > n {
        return false;
    }
    if k <= 1 {
        return k == 1; // single-element pattern always matches
    }
    // The last element of perm must play the role of the last element of pattern.
    // Find subsequences of length k-1 in perm[..n-1] matching pattern[..k-1],
    // where the last element perm[n-1] has the right relative order with pattern[k-1].
    search_pattern_ending_here(perm, pattern, 0, &mut Vec::with_capacity(k))
}

fn search_pattern_ending_here(
    perm: &[u8],
    pattern: &[u8],
    start: usize,
    chosen: &mut Vec<u8>,
) -> bool {
    let k = pattern.len();
    let n = perm.len();
    let last_val = perm[n - 1];

    if chosen.len() == k - 1 {
        // Check if last_val fits as pattern[k-1]
        for j in 0..chosen.len() {
            let pat_cmp = pattern[j].cmp(&pattern[k - 1]);
            let perm_cmp = chosen[j].cmp(&last_val);
            if pat_cmp != perm_cmp {
                return false;
            }
        }
        return true;
    }

    let remaining = (k - 1) - chosen.len();
    // Only search in perm[start..n-1] (exclude last element)
    if n < 2 {
        return false;
    }
    let end = (n - 1).saturating_sub(remaining);
    for i in start..=end {
        let pos = chosen.len();
        let mut valid = true;
        for j in 0..pos {
            let pat_cmp = pattern[j].cmp(&pattern[pos]);
            let perm_cmp = chosen[j].cmp(&perm[i]);
            if pat_cmp != perm_cmp {
                valid = false;
                break;
            }
        }
        if valid {
            chosen.push(perm[i]);
            if search_pattern_ending_here(perm, pattern, i + 1, chosen) {
                chosen.pop();
                return true;
            }
            chosen.pop();
        }
    }
    false
}

// ============================================================
// Codes and representations
// ============================================================

/// Compute the Lehmer code (inversion table) of a permutation.
///
/// Entry `i` counts the number of positions `j > i` with `σ(j) < σ(i)`.
/// The sum of the code equals the inversion number.
pub fn lehmer_code(perm: &[u8]) -> Vec<u8> {
    let n = perm.len();
    let mut code = vec![0u8; n];
    for i in 0..n {
        for j in i + 1..n {
            if perm[j] < perm[i] {
                code[i] += 1;
            }
        }
    }
    code
}

/// Reconstruct a permutation from its Lehmer code.
///
/// Inverse of [`lehmer_code`]: for each entry `c[i]`, select the `c[i]`-th
/// smallest remaining value.
pub fn from_lehmer_code(code: &[u8]) -> Vec<u8> {
    let n = code.len();
    let mut available: Vec<u8> = (1..=n as u8).collect();
    let mut perm = Vec::with_capacity(n);
    for &c in code {
        perm.push(available.remove(c as usize));
    }
    perm
}

/// Compute a reduced word for a permutation.
///
/// Returns a sequence `[i₁, i₂, …]` of 1-indexed simple transposition
/// indices such that `σ = s_{i₁} · s_{i₂} · …` (leftmost-descent
/// bubble-sort algorithm). The length equals the inversion number.
pub fn reduced_word(perm: &[u8]) -> Vec<u8> {
    let mut pi = perm.to_vec();
    let mut word = Vec::new();
    loop {
        let mut found = false;
        for i in 0..pi.len().saturating_sub(1) {
            if pi[i] > pi[i + 1] {
                word.push((i + 1) as u8);
                pi.swap(i, i + 1);
                found = true;
                break;
            }
        }
        if !found {
            break;
        }
    }
    word
}

/// Reconstruct a permutation from a word of simple transpositions.
///
/// Applies transpositions right-to-left to the identity:
/// `σ = s_{w[0]} · s_{w[1]} · … · s_{w[k−1]}`.
pub fn from_reduced_word(word: &[u8], n: u8) -> Vec<u8> {
    let mut perm: Vec<u8> = (1..=n).collect();
    for &s in word.iter().rev() {
        let i = (s - 1) as usize;
        perm.swap(i, i + 1);
    }
    perm
}

// ============================================================
// Cycle structure
// ============================================================

/// Compute the full cycle decomposition of a permutation, including fixed points.
///
/// Each cycle lists elements starting from the smallest. Cycles are sorted
/// by their smallest element (i.e., by starting position).
pub fn cycle_decomposition(perm: &[u8]) -> Vec<Vec<u8>> {
    let n = perm.len();
    let mut visited = vec![false; n];
    let mut cycles = Vec::new();
    for start in 0..n {
        if visited[start] {
            continue;
        }
        let mut cycle = Vec::new();
        let mut j = start;
        loop {
            visited[j] = true;
            cycle.push((j + 1) as u8);
            j = (perm[j] as usize) - 1;
            if j == start {
                break;
            }
        }
        cycles.push(cycle);
    }
    cycles
}

/// Return the cycle type: a partition (sorted in decreasing order)
/// of cycle lengths.
pub fn cycle_type(perm: &[u8]) -> Vec<usize> {
    let mut lengths: Vec<usize> = cycle_decomposition(perm).iter().map(|c| c.len()).collect();
    lengths.sort_unstable_by(|a, b| b.cmp(a));
    lengths
}

/// Compute the inverse permutation σ⁻¹.
pub fn inverse(perm: &[u8]) -> Vec<u8> {
    let n = perm.len();
    let mut inv = vec![0u8; n];
    for i in 0..n {
        inv[(perm[i] as usize) - 1] = (i + 1) as u8;
    }
    inv
}

// ============================================================
// Classical bijections
// ============================================================

/// Apply the cycle-word variant of Foata's first fundamental transformation.
///
/// Write each cycle of `perm` with its largest element first, order the cycles
/// by increasing largest element, and erase the parentheses.
pub fn foata_cycle_word_map(perm: &[u8]) -> Vec<u8> {
    let mut cycles = cycle_decomposition(perm);
    for cycle in &mut cycles {
        let max_pos = cycle
            .iter()
            .enumerate()
            .max_by_key(|&(_, value)| value)
            .map(|(idx, _)| idx)
            .unwrap();
        cycle.rotate_left(max_pos);
    }
    cycles.sort_unstable_by_key(|cycle| cycle[0]);
    cycles.into_iter().flatten().collect()
}

/// Invert [`foata_cycle_word_map`].
///
/// The cycle starts are exactly the left-to-right maxima of the word.
pub fn inverse_foata_cycle_word_map(word: &[u8]) -> Vec<u8> {
    let n = word.len();
    if n == 0 {
        return Vec::new();
    }
    let mut perm = vec![0u8; n];
    let mut starts = Vec::new();
    let mut max_so_far = 0u8;
    for (idx, &value) in word.iter().enumerate() {
        if value > max_so_far {
            starts.push(idx);
            max_so_far = value;
        }
    }
    starts.push(n);

    for window in starts.windows(2) {
        let cycle = &word[window[0]..window[1]];
        for idx in 0..cycle.len() {
            let source = cycle[idx];
            let target = cycle[(idx + 1) % cycle.len()];
            perm[(source - 1) as usize] = target;
        }
    }
    perm
}

/// Apply Foata's fundamental transformation (insertion-based).
///
/// This is a bijection S_n → S_n satisfying `inv(foata(σ)) = maj(σ)`.
/// The algorithm inserts elements of σ one at a time from left to right.
pub fn foata_map(perm: &[u8]) -> Vec<u8> {
    if perm.len() <= 1 {
        return perm.to_vec();
    }
    let mut q = vec![perm[0]];
    for &s in &perm[1..] {
        foata_insert(&mut q, s);
    }
    q
}

/// Insert element `s` into the current Foata word `q`.
///
/// Marks elements in `q` according to comparison with `s` (relative to the
/// last element of `q`), inserts separators after marked elements, appends `s`,
/// then rotates each block of consecutive non-separator elements to the right.
fn foata_insert(q: &mut Vec<u8>, s: u8) {
    let last = *q.last().unwrap();

    let mut seq: Vec<Option<u8>> = Vec::with_capacity(q.len() * 2 + 1);
    for &x in q.iter() {
        seq.push(Some(x));
        let mark = if last < s { x < s } else { x > s };
        if mark {
            seq.push(None); // separator
        }
    }
    seq.push(Some(s));

    // Group consecutive non-separator elements and rotate each group right
    let mut result = Vec::with_capacity(q.len() + 1);
    let mut group: Vec<u8> = Vec::new();
    for item in seq {
        match item {
            Some(x) => group.push(x),
            None => {
                if !group.is_empty() {
                    let tail = group.pop().unwrap();
                    result.push(tail);
                    result.append(&mut group);
                }
            }
        }
    }
    if !group.is_empty() {
        let tail = group.pop().unwrap();
        result.push(tail);
        result.append(&mut group);
    }

    *q = result;
}

// ============================================================
// Direct sum and skew sum
// ============================================================

/// Compute the direct sum π ⊕ τ.
///
/// Concatenates π and a shifted copy of τ, placing the permutation
/// matrices along the main diagonal of a block matrix.
pub fn direct_sum(pi: &[u8], tau: &[u8]) -> Vec<u8> {
    let n = pi.len() as u8;
    let mut result = pi.to_vec();
    result.extend(tau.iter().map(|&x| x + n));
    result
}

/// Compute the skew sum π ⊖ τ.
///
/// Concatenates a shifted copy of π and τ, placing the permutation
/// matrices along the anti-diagonal of a block matrix.
pub fn skew_sum(pi: &[u8], tau: &[u8]) -> Vec<u8> {
    let m = tau.len() as u8;
    let mut result: Vec<u8> = pi.iter().map(|&x| x + m).collect();
    result.extend_from_slice(tau);
    result
}

/// Decompose a permutation into direct-sum and skew-sum blocks.
///
/// Returns the list of indecomposable blocks (preserving original values).
/// A separable permutation decomposes all the way to singletons.
pub fn split_separable(perm: &[u8]) -> Vec<Vec<u8>> {
    if perm.len() <= 1 {
        return vec![perm.to_vec()];
    }
    let n = perm.len();

    // Find all split points: positions k where perm[..k] and perm[k..]
    // form a direct-sum or skew-sum decomposition.
    let mut splits = Vec::new();
    for k in 1..n {
        let (left, right) = perm.split_at(k);
        let left_max = *left.iter().max().unwrap();
        let left_min = *left.iter().min().unwrap();
        let right_max = *right.iter().max().unwrap();
        let right_min = *right.iter().min().unwrap();
        if left_max < right_min || left_min > right_max {
            splits.push(k);
        }
    }

    if splits.is_empty() {
        return vec![perm.to_vec()];
    }

    let mut boundaries = vec![0];
    boundaries.extend(&splits);
    boundaries.push(n);

    let mut result = Vec::new();
    for w in boundaries.windows(2) {
        result.extend(split_separable(&perm[w[0]..w[1]]));
    }
    result
}

// ============================================================
// Permutation predicates
// ============================================================

/// Check if permutation is alternating (up-down: p₁ < p₂ > p₃ < p₄ > …).
pub fn is_alternating(perm: &[u8]) -> bool {
    for i in 0..perm.len().saturating_sub(1) {
        if i % 2 == 0 {
            if perm[i] >= perm[i + 1] {
                return false;
            }
        } else if perm[i] <= perm[i + 1] {
            return false;
        }
    }
    true
}

/// Check if permutation is a derangement (no fixed points).
pub fn is_derangement(perm: &[u8]) -> bool {
    (0..perm.len()).all(|i| perm[i] as usize != i + 1)
}

/// Check if permutation is an involution (σ² = id).
pub fn is_involution(perm: &[u8]) -> bool {
    (0..perm.len()).all(|i| perm[perm[i] as usize - 1] as usize == i + 1)
}

/// Check if permutation is vexillary (avoids 2143).
pub fn is_vexillary(perm: &[u8]) -> bool {
    !contains_pattern(perm, &[2, 1, 4, 3])
}

/// Check if permutation is Grassmann (at most one descent).
pub fn is_grassmann(perm: &[u8]) -> bool {
    let mut desc_count = 0;
    for i in 0..perm.len().saturating_sub(1) {
        if perm[i] > perm[i + 1] {
            desc_count += 1;
            if desc_count > 1 {
                return false;
            }
        }
    }
    true
}

/// Check if permutation is separable (avoids both 2413 and 3142).
pub fn is_separable(perm: &[u8]) -> bool {
    !contains_pattern(perm, &[2, 4, 1, 3]) && !contains_pattern(perm, &[3, 1, 4, 2])
}

/// Check if permutation is simsun.
///
/// A permutation is simsun if for every k ∈ {3,…,n}, the subsequence of
/// entries ≤ k has no double descent (no three consecutive a > b > c).
/// Counted by the tangent numbers E_n (OEIS A000111).
pub fn is_simsun(perm: &[u8]) -> bool {
    let n = perm.len();
    if n <= 2 {
        return true;
    }
    for k in 3..=n as u8 {
        let sub: Vec<u8> = perm.iter().filter(|&&x| x <= k).copied().collect();
        for i in 0..sub.len().saturating_sub(2) {
            if sub[i] > sub[i + 1] && sub[i + 1] > sub[i + 2] {
                return false;
            }
        }
    }
    true
}

// ============================================================
// Special families
// ============================================================

/// Generate all Stirling permutations of order n.
///
/// A Stirling permutation of the multiset {1,1,2,2,…,n,n} has the property
/// that for each i, all values between the two copies of i are greater than i.
/// Counted by (2n−1)!! (OEIS A001147).
pub fn all_stirling_permutations(n: u8) -> Vec<Vec<u8>> {
    if n == 0 {
        return vec![vec![]];
    }
    if n == 1 {
        return vec![vec![1, 1]];
    }
    let prev = all_stirling_permutations(n - 1);
    let mut result = Vec::new();
    for sp in &prev {
        // Insert two adjacent copies of n at each gap
        for j in 0..=sp.len() {
            let mut new_sp = Vec::with_capacity(sp.len() + 2);
            new_sp.extend_from_slice(&sp[..j]);
            new_sp.push(n);
            new_sp.push(n);
            new_sp.extend_from_slice(&sp[j..]);
            result.push(new_sp);
        }
    }
    result
}

// ============================================================
// Backtracking
// ============================================================

/// Compute the backtrack image of π with respect to a set S.
///
/// Returns the unique σ ∈ S that minimizes (σ(π(1)), σ(π(2)), …, σ(π(n)))
/// lexicographically.
pub fn backtrack_image(pi: &[u8], s_set: &[Vec<u8>]) -> Vec<u8> {
    let mut candidates: Vec<usize> = (0..s_set.len()).collect();

    for &pi_k in pi.iter() {
        let pos = pi_k as usize - 1;
        let min_val = candidates.iter().map(|&idx| s_set[idx][pos]).min().unwrap();
        candidates.retain(|&idx| s_set[idx][pos] == min_val);
    }

    s_set[candidates[0]].clone()
}

/// Compute the backtrack image of a search order with respect to derangements.
///
/// This uses the greedy completion criterion: after assigning the current
/// position, the only forbidden remaining state is a common singleton position
/// and value.
pub fn derangement_backtrack_image(pi: &[u8]) -> Vec<u8> {
    let n = pi.len();
    assert!(n != 1, "D_1 has no derangements");
    let mut image = vec![0u8; n];
    let mut used_values = vec![false; n + 1];

    for (step, &p) in pi.iter().enumerate() {
        let p_idx = p as usize;
        let one_position_will_remain = step + 2 == n;
        let remaining_position = one_position_will_remain.then(|| pi[step + 1]);
        let mut chosen = None;

        for v in 1..=n as u8 {
            if used_values[v as usize] || v == p {
                continue;
            }

            if let Some(r) = remaining_position {
                let remaining_value = (1..=n as u8)
                    .find(|&w| !used_values[w as usize] && w != v)
                    .expect("one value should remain");
                if remaining_value == r {
                    continue;
                }
            }

            chosen = Some(v);
            break;
        }

        let v = chosen.expect("derangement completion should exist");
        image[p_idx - 1] = v;
        used_values[v as usize] = true;
    }

    image
}

/// Compute the backtrack polynomial: Σ_{π ∈ S_n} t^{stat(backtrack_π(S))}.
pub fn backtrack_stat_poly(n: u8, s_set: &[Vec<u8>], stat: crate::statistics::Stat) -> Vec<i64> {
    let all_pi = all_permutations(n);
    let mut coeffs: Vec<i64> = Vec::new();

    for pi in &all_pi {
        let sigma = backtrack_image(pi, s_set);
        let val = crate::statistics::compute(&sigma, stat);
        while coeffs.len() <= val {
            coeffs.push(0);
        }
        coeffs[val] += 1;
    }

    if coeffs.is_empty() {
        coeffs.push(0);
    }
    coeffs
}

/// Compute backtrack polynomial over a given subset of search orders.
pub fn backtrack_stat_poly_subset(
    pi_set: &[Vec<u8>],
    s_set: &[Vec<u8>],
    stat: crate::statistics::Stat,
) -> Vec<i64> {
    let mut coeffs: Vec<i64> = Vec::new();

    for pi in pi_set {
        let sigma = backtrack_image(pi, s_set);
        let val = crate::statistics::compute(&sigma, stat);
        while coeffs.len() <= val {
            coeffs.push(0);
        }
        coeffs[val] += 1;
    }

    if coeffs.is_empty() {
        coeffs.push(0);
    }
    coeffs
}

// ============================================================
// Special families: derangements, involutions
// ============================================================

/// Generate all derangements of \[1..n\] (permutations with no fixed points).
///
/// A derangement σ satisfies σ(i) ≠ i for all i.
/// Counted by subfactorials: D(n) = n! Σ_{k=0}^{n} (-1)^k / k!.
pub fn all_derangements(n: u8) -> Vec<Vec<u8>> {
    if n == 0 {
        return vec![vec![]]; // empty permutation is vacuously a derangement
    }
    all_permutations(n)
        .into_iter()
        .filter(|sigma| sigma.iter().enumerate().all(|(i, &v)| v != (i as u8 + 1)))
        .collect()
}

/// Generate all involutions of \[1..n\] (permutations σ with σ² = id).
///
/// An involution consists only of fixed points and 2-cycles.
/// Counted by A000085: 1, 1, 2, 4, 10, 26, 76, 232, 764, ...
pub fn all_involutions(n: u8) -> Vec<Vec<u8>> {
    let mut result = Vec::new();
    let mut current = vec![0u8; n as usize];
    involution_rec(n, 0, &mut current, &mut result);
    result
}

fn involution_rec(n: u8, pos: u8, current: &mut Vec<u8>, result: &mut Vec<Vec<u8>>) {
    if pos == n {
        result.push(current.clone());
        return;
    }
    if current[pos as usize] != 0 {
        // Already paired, skip
        involution_rec(n, pos + 1, current, result);
        return;
    }
    // Option 1: fixed point
    current[pos as usize] = pos + 1;
    involution_rec(n, pos + 1, current, result);
    current[pos as usize] = 0;

    // Option 2: pair with a later unpaired element
    for j in (pos + 1)..n {
        if current[j as usize] == 0 {
            current[pos as usize] = j + 1;
            current[j as usize] = pos + 1;
            involution_rec(n, pos + 1, current, result);
            current[pos as usize] = 0;
            current[j as usize] = 0;
        }
    }
}

/// Count the number of cycles of a permutation.
pub fn num_cycles(sigma: &[u8]) -> usize {
    let n = sigma.len();
    let mut visited = vec![false; n];
    let mut cycles = 0;
    for i in 0..n {
        if visited[i] {
            continue;
        }
        cycles += 1;
        let mut j = i;
        while !visited[j] {
            visited[j] = true;
            j = (sigma[j] - 1) as usize;
        }
    }
    cycles
}

// ============================================================
// Special families: alternating permutations
// ============================================================

/// Generate all alternating (up-down) permutations of \[1..n\].
///
/// An alternating permutation satisfies w(1) < w(2) > w(3) < w(4) > ...
/// (starting with an ascent). These are counted by the tangent/secant
/// numbers (OEIS A000111).
///
/// Uses a recursive decomposition: place the maximum value n at an
/// odd position (1-indexed), split the remaining elements into left
/// and right subsequences, and recurse.
pub fn alternating_permutations(n: u8) -> Vec<Vec<u8>> {
    alt_perms_cached(n, &mut std::collections::HashMap::new())
}

fn alt_perms_cached(
    n: u8,
    cache: &mut std::collections::HashMap<u8, Vec<Vec<u8>>>,
) -> Vec<Vec<u8>> {
    if let Some(cached) = cache.get(&n) {
        return cached.clone();
    }
    let result = if n == 0 {
        vec![vec![]]
    } else if n == 1 {
        vec![vec![1]]
    } else {
        let mut perms = Vec::new();
        // Place n at position k (0-indexed), where k is odd (1-indexed: even position)
        // k ranges over 1, 3, 5, ... (0-indexed) so that n sits at a "peak"
        let mut k = 1usize; // 0-indexed position
        while k < n as usize {
            let left_size = k as u8;
            let right_size = n - 1 - left_size;
            let alt_left = alt_perms_cached(left_size, cache);
            let alt_right = alt_perms_cached(right_size, cache);
            // Choose which elements go left
            let elems: Vec<u8> = (1..n).collect();
            for subset in subsets_of_size(&elems, left_size as usize) {
                let complement: Vec<u8> = elems
                    .iter()
                    .copied()
                    .filter(|x| !subset.contains(x))
                    .collect();
                for al in &alt_left {
                    for ar in &alt_right {
                        // Map left alternating perm to the chosen subset
                        let mut w = Vec::with_capacity(n as usize);
                        for &i in al {
                            w.push(subset[(i - 1) as usize]);
                        }
                        w.push(n);
                        for &i in ar {
                            w.push(complement[(i - 1) as usize]);
                        }
                        perms.push(w);
                    }
                }
            }
            k += 2;
        }
        perms
    };
    cache.insert(n, result.clone());
    result
}

/// Generate all reverse-alternating (down-up) permutations of \[1..n\].
///
/// A reverse-alternating permutation satisfies w(1) > w(2) < w(3) > w(4) < ...
/// (starting with a descent).
pub fn reverse_alternating_permutations(n: u8) -> Vec<Vec<u8>> {
    // Complement of alternating: replace each entry i with n+1-i
    let alt = alternating_permutations(n);
    alt.into_iter()
        .map(|w| w.iter().map(|&v| n + 1 - v).collect())
        .collect()
}

/// All subsets of `elems` with exactly `k` elements, in lexicographic order.
fn subsets_of_size(elems: &[u8], k: usize) -> Vec<Vec<u8>> {
    let mut result = Vec::new();
    let mut current = Vec::with_capacity(k);
    subsets_rec(elems, k, 0, &mut current, &mut result);
    result
}

fn subsets_rec(
    elems: &[u8],
    k: usize,
    start: usize,
    current: &mut Vec<u8>,
    result: &mut Vec<Vec<u8>>,
) {
    if current.len() == k {
        result.push(current.clone());
        return;
    }
    let remaining = k - current.len();
    for i in start..=(elems.len() - remaining) {
        current.push(elems[i]);
        subsets_rec(elems, k, i + 1, current, result);
        current.pop();
    }
}

// ============================================================
// Special families: k-alternating permutations
// ============================================================

/// Generate all k-alternating permutations of \[1..n\].
///
/// A k-alternating permutation has descents exactly at positions
/// divisible by k (1-indexed). That is:
/// - descent at position j (σ(j) > σ(j+1)) iff k | j,
/// - ascent at position j (σ(j) < σ(j+1)) otherwise.
///
/// Special cases:
/// - k=1: only the decreasing permutation [n, n-1, ..., 1]
/// - k=2: standard alternating (up-down) permutations
/// - n=0 and k>0: the unique empty permutation
///
/// # Panics
///
/// Panics if `k=0`, since divisibility by zero does not define a descent set.
pub fn k_alternating_permutations(n: u8, k: u8) -> Vec<Vec<u8>> {
    assert!(k > 0, "k must be positive for k-alternating permutations");
    if n == 0 {
        return vec![vec![]];
    }
    if k > n {
        // No positions are divisible by k, so the permutation is increasing.
        return vec![(1..=n).collect()];
    }

    // Build the required descent set
    let descent_set: Vec<bool> = (0..n as usize)
        .map(|i| {
            let pos = i + 1; // 1-indexed
            pos < n as usize && pos % k as usize == 0
        })
        .collect();

    // Filter all permutations by descent set match
    // For small n, brute force is fine
    let all = all_permutations(n);
    all.into_iter()
        .filter(|perm| {
            for i in 0..(n as usize - 1) {
                let is_descent = perm[i] > perm[i + 1];
                if is_descent != descent_set[i] {
                    return false;
                }
            }
            true
        })
        .collect()
}

// ============================================================
// Tests
// ============================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::statistics::{compute, Stat};

    // --- Derangements and involutions ---

    #[test]
    fn test_derangement_counts() {
        // Subfactorials: 1, 0, 1, 2, 9, 44, 265, 1854
        let expected = [1, 0, 1, 2, 9, 44, 265, 1854];
        for (n, &e) in expected.iter().enumerate() {
            assert_eq!(all_derangements(n as u8).len(), e, "n={}", n);
        }
    }

    #[test]
    fn test_derangements_no_fixed_points() {
        for n in 1..=6u8 {
            for d in all_derangements(n) {
                assert!(
                    d.iter().enumerate().all(|(i, &v)| v != (i as u8 + 1)),
                    "n={} d={:?} has fixed point",
                    n,
                    d
                );
            }
        }
    }

    #[test]
    fn test_involution_counts() {
        // A000085: 1, 1, 2, 4, 10, 26, 76, 232, 764
        let expected = [1, 1, 2, 4, 10, 26, 76, 232, 764];
        for (n, &e) in expected.iter().enumerate() {
            assert_eq!(all_involutions(n as u8).len(), e, "n={}", n);
        }
    }

    #[test]
    fn test_involutions_are_self_inverse() {
        for n in 1..=6u8 {
            for inv in all_involutions(n) {
                for (i, &v) in inv.iter().enumerate() {
                    assert_eq!(
                        inv[(v - 1) as usize],
                        (i as u8 + 1),
                        "n={} inv={:?} not self-inverse",
                        n,
                        inv
                    );
                }
            }
        }
    }

    #[test]
    fn test_num_cycles() {
        assert_eq!(num_cycles(&[1, 2, 3]), 3); // id: 3 fixed points
        assert_eq!(num_cycles(&[2, 1, 3]), 2); // one 2-cycle + one fixed
        assert_eq!(num_cycles(&[2, 3, 1]), 1); // one 3-cycle
        assert_eq!(num_cycles(&[3, 2, 1]), 2); // (1 3) + (2)
    }

    // --- Enumeration ---

    #[test]
    fn test_all_permutations_counts() {
        assert_eq!(all_permutations(0).len(), 1);
        assert_eq!(all_permutations(1).len(), 1);
        assert_eq!(all_permutations(3).len(), 6);
        assert_eq!(all_permutations(4).len(), 24);
        assert_eq!(all_permutations(5).len(), 120);
    }

    #[test]
    fn test_lexicographic_order() {
        let perms = all_permutations(3);
        assert_eq!(perms[0], vec![1, 2, 3]);
        assert_eq!(perms[5], vec![3, 2, 1]);
    }

    // --- Pattern avoidance ---

    #[test]
    fn test_contains_pattern() {
        assert!(contains_pattern(&[3, 1, 2], &[2, 1]));
        assert!(!contains_pattern(&[1, 2, 3], &[2, 1]));
        assert!(contains_pattern(&[1, 3, 2], &[1, 3, 2]));
        assert!(!contains_pattern(&[1, 2, 3], &[3, 1, 2]));
    }

    #[test]
    fn test_312_avoiding_catalan() {
        // 312-avoiding permutations of S_n are counted by Catalan numbers.
        let pat = parse_sequence("312");
        for (n, catalan) in [(3, 5), (4, 14), (5, 42)] {
            let count = all_permutations(n)
                .iter()
                .filter(|p| !contains_pattern(p, &pat))
                .count();
            assert_eq!(count, catalan, "C_{} should be {}", n, catalan);
        }
    }

    #[test]
    fn test_avoiding_permutations_catalan() {
        // Direct generation should give same results as filter approach.
        for (pat_str, n, catalan) in [
            ("132", 5, 42),
            ("231", 5, 42),
            ("321", 5, 42),
            ("213", 6, 132),
        ] {
            let pat = parse_sequence(pat_str);
            let direct = avoiding_permutations(n, &[pat.clone()]);
            assert_eq!(
                direct.len(),
                catalan,
                "Av_{}({}) should have {} elements",
                pat_str,
                n,
                catalan
            );
            // Verify each element actually avoids the pattern
            for p in &direct {
                assert!(
                    !contains_pattern(p, &pat),
                    "{:?} should avoid {}",
                    p,
                    pat_str
                );
            }
        }
    }

    #[test]
    fn test_avoiding_permutations_matches_filter() {
        // Check that direct generation matches generate-then-filter for n=5.
        for pat_str in &["132", "213", "231", "312", "321", "123"] {
            let pat = parse_sequence(pat_str);
            let direct = avoiding_permutations(5, &[pat.clone()]);
            let filtered: Vec<Vec<u8>> = all_permutations(5)
                .into_iter()
                .filter(|p| !contains_pattern(p, &pat))
                .collect();
            assert_eq!(
                direct.len(),
                filtered.len(),
                "Av_{}(5) count mismatch",
                pat_str
            );
            // Check same elements (both should be sorted since we build left-to-right)
            for p in &direct {
                assert!(
                    filtered.contains(p),
                    "{:?} in direct but not filtered for Av_{}",
                    p,
                    pat_str
                );
            }
        }
    }

    #[test]
    fn test_avoiding_multiple_patterns() {
        // Av_{132,321}(5) = permutations avoiding both 132 and 321
        let pats = vec![parse_sequence("132"), parse_sequence("321")];
        let direct = avoiding_permutations(5, &pats);
        let filtered: Vec<Vec<u8>> = all_permutations(5)
            .into_iter()
            .filter(|p| pats.iter().all(|pat| !contains_pattern(p, pat)))
            .collect();
        assert_eq!(direct.len(), filtered.len());
    }

    #[test]
    fn test_foata_cycle_word_map_paper_example() {
        let preimage = parse_sequence("5637421");
        let image = parse_sequence("3627154");
        assert_eq!(foata_cycle_word_map(&preimage), image);
        assert_eq!(inverse_foata_cycle_word_map(&image), preimage);
    }

    #[test]
    fn test_arrow_pattern_parser() {
        assert_eq!(
            parse_arrow_pattern("(12; 1 -> 3)"),
            ArrowPattern::new(3, vec![1, 2], vec![(1, 3)])
        );
        assert_eq!(
            parse_arrow_pattern("10,2,1,3,4,5,6,7,8,9;1->10"),
            ArrowPattern::new(10, vec![10, 2, 1, 3, 4, 5, 6, 7, 8, 9], vec![(1, 10)])
        );
    }

    #[test]
    fn test_arrow_pattern_contains_paper_example() {
        let pi = parse_sequence("3627154");
        let alpha = parse_arrow_pattern("12;1->3");
        let beta = parse_arrow_pattern("231;1->4");
        assert!(contains_arrow_pattern(&pi, &alpha));
        assert!(!contains_arrow_pattern(&pi, &beta));
    }

    #[test]
    fn test_arrow_pattern_count_agrees_with_contains() {
        let alpha = parse_arrow_pattern("12;1->3");
        for n in 1..=5 {
            for perm in all_permutations(n) {
                assert_eq!(
                    contains_arrow_pattern(&perm, &alpha),
                    count_arrow_pattern_occurrences(&perm, &alpha) > 0,
                    "perm={:?}",
                    perm
                );
            }
        }
    }

    #[test]
    fn test_arrow_avoiding_schroder_counts() {
        // Archer--Laudone: avoiding (12;1->3) is counted by large Schroder
        // numbers S_{n-1}: 1, 2, 6, 22, 90.
        let alpha = parse_arrow_pattern("12;1->3");
        for (n, expected) in [(1, 1), (2, 2), (3, 6), (4, 22), (5, 90)] {
            let direct = arrow_avoiding_permutations(n, std::slice::from_ref(&alpha));
            assert_eq!(direct.len(), expected, "n={}", n);
        }
    }

    #[test]
    fn test_arrow_and_classical_avoidance_counts() {
        // Fu--Yang: a_n(123, (12;1->3)) = 2^n - n.
        let alpha = parse_arrow_pattern("12;1->3");
        for n in 2..=6 {
            let constraints = PermConstraints {
                avoiding: vec![parse_sequence("123")],
                avoiding_arrow: vec![alpha.clone()],
                ..Default::default()
            };
            assert_eq!(
                filtered_permutations(n, &constraints).len(),
                (1usize << n) - n as usize,
                "n={}",
                n
            );
        }

        // Fu--Yang: a_n(321, (12;1->2)) = M_n, the Motzkin numbers.
        let beta = parse_arrow_pattern("12;1->2");
        for (n, expected) in [(1, 1), (2, 2), (3, 4), (4, 9), (5, 21), (6, 51)] {
            let constraints = PermConstraints {
                avoiding: vec![parse_sequence("321")],
                avoiding_arrow: vec![beta.clone()],
                ..Default::default()
            };
            assert_eq!(
                filtered_permutations(n, &constraints).len(),
                expected,
                "n={}",
                n
            );
        }
    }

    #[test]
    fn test_filtered_derangements() {
        let c = PermConstraints {
            derangement: true,
            ..Default::default()
        };
        for n in 1..=6 {
            let direct = filtered_permutations(n, &c);
            let filtered: Vec<Vec<u8>> = all_permutations(n)
                .into_iter()
                .filter(|p| is_derangement(p))
                .collect();
            assert_eq!(direct.len(), filtered.len(), "derangements of S_{}", n);
        }
    }

    #[test]
    fn test_filtered_alternating() {
        let c = PermConstraints {
            alternating: true,
            ..Default::default()
        };
        for n in 1..=7 {
            let direct = filtered_permutations(n, &c);
            let filtered: Vec<Vec<u8>> = all_permutations(n)
                .into_iter()
                .filter(|p| is_alternating(p))
                .collect();
            assert_eq!(direct.len(), filtered.len(), "alternating perms of S_{}", n);
        }
    }

    #[test]
    fn test_filtered_involutions() {
        let c = PermConstraints {
            involution: true,
            ..Default::default()
        };
        // Involutions of S_n: 1, 1, 2, 4, 10, 26, 76
        for (n, expected) in [(1, 1), (2, 2), (3, 4), (4, 10), (5, 26), (6, 76)] {
            let direct = filtered_permutations(n, &c);
            assert_eq!(direct.len(), expected, "involutions of S_{}", n);
            for p in &direct {
                assert!(is_involution(p), "{:?} should be an involution", p);
            }
        }
    }

    #[test]
    fn test_filtered_combined() {
        // 132-avoiding derangements
        let c = PermConstraints {
            avoiding: vec![parse_sequence("132")],
            derangement: true,
            ..Default::default()
        };
        for n in 1..=7 {
            let direct = filtered_permutations(n, &c);
            let filtered: Vec<Vec<u8>> = all_permutations(n)
                .into_iter()
                .filter(|p| !contains_pattern(p, &[1, 3, 2]) && is_derangement(p))
                .collect();
            assert_eq!(
                direct.len(),
                filtered.len(),
                "132-avoiding derangements of S_{}",
                n
            );
        }
    }

    // --- Lehmer code ---

    #[test]
    fn test_lehmer_code_examples() {
        assert_eq!(lehmer_code(&[1, 2, 3]), vec![0, 0, 0]);
        assert_eq!(lehmer_code(&[3, 2, 1]), vec![2, 1, 0]);
        assert_eq!(lehmer_code(&[3, 1, 4, 2]), vec![2, 0, 1, 0]);
    }

    #[test]
    fn test_lehmer_code_empty() {
        assert_eq!(lehmer_code(&[]), Vec::<u8>::new());
        assert_eq!(from_lehmer_code(&[]), Vec::<u8>::new());
    }

    #[test]
    fn test_lehmer_roundtrip() {
        for perm in &all_permutations(5) {
            let code = lehmer_code(perm);
            assert_eq!(from_lehmer_code(&code), *perm);
        }
    }

    #[test]
    fn test_lehmer_sum_equals_inversions() {
        for perm in &all_permutations(5) {
            let code = lehmer_code(perm);
            let code_sum: usize = code.iter().map(|&c| c as usize).sum();
            assert_eq!(code_sum, compute(perm, Stat::Inv));
        }
    }

    #[test]
    fn test_lehmer_bounds() {
        // Entry i of the Lehmer code satisfies 0 ≤ code[i] ≤ n-1-i
        for perm in &all_permutations(5) {
            let code = lehmer_code(perm);
            let n = perm.len();
            for (i, &c) in code.iter().enumerate() {
                assert!(
                    (c as usize) <= n - 1 - i,
                    "Lehmer code[{}] = {} out of range for {:?}",
                    i,
                    c,
                    perm
                );
            }
        }
    }

    // --- Reduced words ---

    #[test]
    fn test_reduced_word_identity() {
        assert_eq!(reduced_word(&[1, 2, 3]), Vec::<u8>::new());
    }

    #[test]
    fn test_reduced_word_w0() {
        // w₀ = [3,2,1] has inv = 3, so reduced word has length 3
        let w = reduced_word(&[3, 2, 1]);
        assert_eq!(w.len(), 3);
        assert_eq!(from_reduced_word(&w, 3), vec![3, 2, 1]);
    }

    #[test]
    fn test_reduced_word_length_equals_inversions() {
        for perm in &all_permutations(5) {
            assert_eq!(reduced_word(perm).len(), compute(perm, Stat::Inv));
        }
    }

    #[test]
    fn test_reduced_word_roundtrip() {
        for perm in &all_permutations(5) {
            let n = perm.len() as u8;
            let word = reduced_word(perm);
            assert_eq!(from_reduced_word(&word, n), *perm);
        }
    }

    #[test]
    fn test_from_reduced_word_generators() {
        // Single transposition s_1 = [2,1,3]
        assert_eq!(from_reduced_word(&[1], 3), vec![2, 1, 3]);
        // s_2 = [1,3,2]
        assert_eq!(from_reduced_word(&[2], 3), vec![1, 3, 2]);
        // Word [1,2] = s_1 then s_2 (left-multiply right-to-left): s_1 ∘ s_2 = [3,1,2]
        assert_eq!(from_reduced_word(&[1, 2], 3), vec![3, 1, 2]);
        // Word [2,1] = s_2 then s_1: s_2 ∘ s_1 = [2,3,1]
        assert_eq!(from_reduced_word(&[2, 1], 3), vec![2, 3, 1]);
    }

    // --- Cycle decomposition ---

    #[test]
    fn test_cycle_decomposition_examples() {
        assert_eq!(
            cycle_decomposition(&[1, 2, 3]),
            vec![vec![1], vec![2], vec![3]]
        );
        assert_eq!(cycle_decomposition(&[2, 3, 1]), vec![vec![1, 2, 3]]);
        assert_eq!(cycle_decomposition(&[2, 1, 3]), vec![vec![1, 2], vec![3]]);
        assert_eq!(
            cycle_decomposition(&[2, 1, 4, 3]),
            vec![vec![1, 2], vec![3, 4]]
        );
    }

    #[test]
    fn test_cycle_decomposition_covers_all_elements() {
        for perm in &all_permutations(5) {
            let cycles = cycle_decomposition(perm);
            let mut all_elems: Vec<u8> = cycles.into_iter().flatten().collect();
            all_elems.sort();
            assert_eq!(all_elems, vec![1, 2, 3, 4, 5]);
        }
    }

    #[test]
    fn test_cycle_type_examples() {
        assert_eq!(cycle_type(&[1, 2, 3]), vec![1, 1, 1]);
        assert_eq!(cycle_type(&[2, 3, 1]), vec![3]);
        assert_eq!(cycle_type(&[2, 1, 4, 3]), vec![2, 2]);
        assert_eq!(cycle_type(&[2, 1, 3]), vec![2, 1]);
    }

    #[test]
    fn test_cycle_count_matches_stat() {
        for perm in &all_permutations(5) {
            assert_eq!(cycle_decomposition(perm).len(), compute(perm, Stat::Cyc));
        }
    }

    #[test]
    fn test_cycle_type_is_partition() {
        // Cycle type should be a partition: sorted decreasingly, summing to n
        for perm in &all_permutations(5) {
            let ct = cycle_type(perm);
            assert_eq!(ct.iter().sum::<usize>(), 5);
            for w in ct.windows(2) {
                assert!(w[0] >= w[1]);
            }
        }
    }

    // --- Inverse ---

    #[test]
    fn test_inverse_examples() {
        assert_eq!(inverse(&[1, 2, 3]), vec![1, 2, 3]);
        assert_eq!(inverse(&[2, 3, 1]), vec![3, 1, 2]);
        assert_eq!(inverse(&[3, 1, 2]), vec![2, 3, 1]);
    }

    #[test]
    fn test_inverse_is_involution() {
        for perm in &all_permutations(5) {
            assert_eq!(inverse(&inverse(perm)), *perm);
        }
    }

    #[test]
    fn test_inverse_composition_is_identity() {
        // σ⁻¹(σ(i)) = i for all i
        for perm in &all_permutations(5) {
            let inv = inverse(perm);
            let n = perm.len();
            for i in 0..n {
                assert_eq!(inv[(perm[i] as usize) - 1], (i + 1) as u8);
            }
        }
    }

    // --- Foata map ---

    #[test]
    fn test_foata_known_value() {
        assert_eq!(foata_map(&[4, 1, 3, 7, 5, 6, 2]), vec![7, 1, 4, 3, 5, 6, 2]);
    }

    #[test]
    fn test_foata_identity() {
        assert_eq!(foata_map(&[1, 2, 3]), vec![1, 2, 3]);
    }

    #[test]
    fn test_foata_singleton() {
        assert_eq!(foata_map(&[1]), vec![1]);
    }

    #[test]
    fn test_foata_maj_to_inv() {
        for perm in &all_permutations(6) {
            let f = foata_map(perm);
            assert_eq!(
                compute(&f, Stat::Inv),
                compute(perm, Stat::Maj),
                "inv(foata(σ)) ≠ maj(σ) for σ = {:?}",
                perm
            );
        }
    }

    #[test]
    fn test_foata_is_bijection() {
        let perms = all_permutations(5);
        let mut images: Vec<Vec<u8>> = perms.iter().map(|p| foata_map(p)).collect();
        images.sort();
        let mut sorted_perms = perms.clone();
        sorted_perms.sort();
        assert_eq!(images, sorted_perms);
    }

    #[test]
    fn test_foata_preserves_des_exc_equidistribution() {
        // Foata sends des to exc (or preserves their equidistribution)
        // Specifically: des(σ) = exc(foata(σ)) would NOT hold in general,
        // but the generating functions over S_n are equidistributed.
        for n in 1..=6u8 {
            let perms = all_permutations(n);
            let mut des_dist = vec![0usize; n as usize];
            let mut foata_exc_dist = vec![0usize; n as usize];
            for p in &perms {
                let d = compute(p, Stat::Des);
                let e = compute(&foata_map(p), Stat::Exc);
                if d < des_dist.len() {
                    des_dist[d] += 1;
                }
                if e < foata_exc_dist.len() {
                    foata_exc_dist[e] += 1;
                }
            }
            // des and exc are equidistributed over S_n
            assert_eq!(des_dist, foata_exc_dist, "Failed for n={}", n);
        }
    }

    // --- Direct sum / skew sum ---

    #[test]
    fn test_direct_sum() {
        assert_eq!(direct_sum(&[2, 1], &[1, 2]), vec![2, 1, 3, 4]);
        assert_eq!(direct_sum(&[1], &[2, 1]), vec![1, 3, 2]);
        assert_eq!(direct_sum(&[], &[1, 2]), vec![1, 2]);
        assert_eq!(direct_sum(&[1, 2], &[]), vec![1, 2]);
    }

    #[test]
    fn test_skew_sum() {
        assert_eq!(skew_sum(&[2, 1], &[1, 2]), vec![4, 3, 1, 2]);
        assert_eq!(skew_sum(&[1], &[1, 2]), vec![3, 1, 2]);
        assert_eq!(skew_sum(&[], &[1, 2]), vec![1, 2]);
        assert_eq!(skew_sum(&[1, 2], &[]), vec![1, 2]);
    }

    #[test]
    fn test_direct_sum_associativity() {
        let a = &[2, 1];
        let b = &[1];
        let c = &[2, 1];
        assert_eq!(
            direct_sum(&direct_sum(a, b), c),
            direct_sum(a, &direct_sum(b, c))
        );
    }

    #[test]
    fn test_direct_sum_inversions() {
        // inv(π ⊕ τ) = inv(π) + inv(τ)
        for pi in &all_permutations(3) {
            for tau in &all_permutations(3) {
                let ds = direct_sum(pi, tau);
                assert_eq!(
                    compute(&ds, Stat::Inv),
                    compute(pi, Stat::Inv) + compute(tau, Stat::Inv)
                );
            }
        }
    }

    // --- Separable permutations ---

    #[test]
    fn test_is_separable() {
        assert!(is_separable(&[1, 2, 3]));
        assert!(is_separable(&[3, 2, 1]));
        assert!(is_separable(&[2, 1, 4, 3]));
        assert!(!is_separable(&[2, 4, 1, 3]));
        assert!(!is_separable(&[3, 1, 4, 2]));
    }

    #[test]
    fn test_separable_counts() {
        // Large Schröder numbers: 1, 2, 6, 22, 90 (OEIS A006318)
        for (n, expected) in [(1, 1), (2, 2), (3, 6), (4, 22), (5, 90)] {
            let count = all_permutations(n)
                .iter()
                .filter(|p| is_separable(p))
                .count();
            assert_eq!(count, expected, "Separable(S_{}) = {}", n, expected);
        }
    }

    #[test]
    fn test_split_separable_identity() {
        assert_eq!(split_separable(&[1, 2, 3]), vec![vec![1], vec![2], vec![3]]);
    }

    #[test]
    fn test_split_separable_w0() {
        assert_eq!(split_separable(&[3, 2, 1]), vec![vec![3], vec![2], vec![1]]);
    }

    #[test]
    fn test_split_separable_mixed() {
        assert_eq!(
            split_separable(&[2, 1, 4, 3]),
            vec![vec![2], vec![1], vec![4], vec![3]]
        );
    }

    #[test]
    fn test_split_separable_singleton_for_indecomposable() {
        // 2413 is not separable and not decomposable
        assert_eq!(split_separable(&[2, 4, 1, 3]), vec![vec![2, 4, 1, 3]]);
    }

    #[test]
    fn test_split_separable_all_singletons_iff_separable() {
        for perm in &all_permutations(5) {
            let blocks = split_separable(perm);
            let all_singletons = blocks.iter().all(|b| b.len() == 1);
            assert_eq!(
                all_singletons,
                is_separable(perm),
                "Mismatch for {:?}",
                perm
            );
        }
    }

    // --- Simsun permutations ---

    #[test]
    fn test_is_simsun() {
        assert!(is_simsun(&[1, 2, 3]));
        assert!(is_simsun(&[1, 3, 2]));
        assert!(!is_simsun(&[3, 2, 1])); // double descent 3 > 2 > 1
        assert!(is_simsun(&[2, 3, 1]));
    }

    #[test]
    fn test_simsun_counts() {
        // Tangent/secant numbers E_n (OEIS A000111): 1, 2, 5, 16, 61
        for (n, expected) in [(1, 1), (2, 2), (3, 5), (4, 16), (5, 61)] {
            let count = all_permutations(n).iter().filter(|p| is_simsun(p)).count();
            assert_eq!(count, expected, "Simsun(S_{}) = {}", n, expected);
        }
    }

    #[test]
    fn test_simsun_small() {
        assert!(is_simsun(&[1]));
        assert!(is_simsun(&[1, 2]));
        assert!(is_simsun(&[2, 1]));
    }

    // --- Stirling permutations ---

    #[test]
    fn test_stirling_counts() {
        // (2n−1)!! = 1, 1, 3, 15, 105 (OEIS A001147)
        assert_eq!(all_stirling_permutations(0).len(), 1);
        assert_eq!(all_stirling_permutations(1).len(), 1);
        assert_eq!(all_stirling_permutations(2).len(), 3);
        assert_eq!(all_stirling_permutations(3).len(), 15);
        assert_eq!(all_stirling_permutations(4).len(), 105);
    }

    #[test]
    fn test_stirling_property() {
        // Between the two copies of i, all values are > i
        for sp in &all_stirling_permutations(4) {
            for i in 1..=4u8 {
                let positions: Vec<usize> = sp
                    .iter()
                    .enumerate()
                    .filter(|&(_, &v)| v == i)
                    .map(|(pos, _)| pos)
                    .collect();
                assert_eq!(positions.len(), 2);
                for &v in &sp[positions[0] + 1..positions[1]] {
                    assert!(v > i, "Violated in {:?} for i={}", sp, i);
                }
            }
        }
    }

    #[test]
    fn test_stirling_no_duplicates() {
        let sps = all_stirling_permutations(4);
        let mut sorted = sps.clone();
        sorted.sort();
        sorted.dedup();
        assert_eq!(sorted.len(), sps.len());
    }

    #[test]
    fn test_stirling_content() {
        // Each Stirling permutation of order n contains exactly two copies of each 1..n
        for sp in &all_stirling_permutations(3) {
            assert_eq!(sp.len(), 6);
            for i in 1..=3u8 {
                assert_eq!(sp.iter().filter(|&&v| v == i).count(), 2);
            }
        }
    }

    // --- Existing predicates ---

    #[test]
    fn test_alternating() {
        assert!(is_alternating(&[1, 3, 2]));
        assert!(is_alternating(&[1, 3, 2, 4]));
        assert!(!is_alternating(&[3, 1, 2]));
        // E_5 = 16
        assert_eq!(
            all_permutations(5)
                .iter()
                .filter(|p| is_alternating(p))
                .count(),
            16
        );
    }

    #[test]
    fn test_derangement() {
        assert!(!is_derangement(&[1, 2, 3]));
        assert!(is_derangement(&[2, 3, 1]));
        assert!(is_derangement(&[2, 1]));
        // D_4 = 9
        assert_eq!(
            all_permutations(4)
                .iter()
                .filter(|p| is_derangement(p))
                .count(),
            9
        );
    }

    #[test]
    fn test_involution() {
        assert!(is_involution(&[1, 2, 3]));
        assert!(is_involution(&[2, 1, 3]));
        assert!(!is_involution(&[2, 3, 1]));
        assert!(is_involution(&[3, 2, 1]));
        // I_4 = 10, I_5 = 26
        assert_eq!(
            all_permutations(4)
                .iter()
                .filter(|p| is_involution(p))
                .count(),
            10
        );
        assert_eq!(
            all_permutations(5)
                .iter()
                .filter(|p| is_involution(p))
                .count(),
            26
        );
    }

    #[test]
    fn test_vexillary() {
        assert_eq!(
            all_permutations(3)
                .iter()
                .filter(|p| is_vexillary(p))
                .count(),
            6
        );
        assert_eq!(
            all_permutations(4)
                .iter()
                .filter(|p| is_vexillary(p))
                .count(),
            23
        );
    }

    #[test]
    fn test_grassmann() {
        assert!(is_grassmann(&[1, 2, 3]));
        assert!(is_grassmann(&[1, 3, 2]));
        assert!(is_grassmann(&[3, 1, 2]));
        assert!(!is_grassmann(&[3, 2, 1]));
        // A(3,0) + A(3,1) = 1 + 4 = 5
        assert_eq!(
            all_permutations(3)
                .iter()
                .filter(|p| is_grassmann(p))
                .count(),
            5
        );
        // A(4,0) + A(4,1) = 1 + 11 = 12
        assert_eq!(
            all_permutations(4)
                .iter()
                .filter(|p| is_grassmann(p))
                .count(),
            12
        );
    }

    // --- Backtracking ---

    #[test]
    fn test_backtrack_image_example() {
        let s_set = vec![vec![2, 3, 1], vec![3, 1, 2]];
        assert_eq!(backtrack_image(&[1, 2, 3], &s_set), vec![2, 3, 1]);
        assert_eq!(backtrack_image(&[2, 1, 3], &s_set), vec![3, 1, 2]);
    }

    #[test]
    fn test_derangement_backtrack_image_matches_generic() {
        for n in 2..=6u8 {
            let derangements: Vec<Vec<u8>> = all_permutations(n)
                .into_iter()
                .filter(|p| is_derangement(p))
                .collect();
            for pi in all_permutations(n) {
                assert_eq!(
                    derangement_backtrack_image(&pi),
                    backtrack_image(&pi, &derangements),
                    "n={} pi={:?}",
                    n,
                    pi
                );
            }
        }
    }

    #[test]
    fn test_backtrack_des_av321_n3() {
        let pat = parse_sequence("321");
        let mut s_set: Vec<Vec<u8>> = all_permutations(3);
        s_set.retain(|p| !contains_pattern(p, &pat));
        let poly = backtrack_stat_poly(3, &s_set, Stat::Des);
        assert_eq!(poly.iter().sum::<i64>(), 6);
    }

    // --- Alternating permutations ---

    #[test]
    fn test_alternating_counts() {
        // A000111: 1, 1, 1, 2, 5, 16, 61, 272, 1385
        let expected = [1, 1, 1, 2, 5, 16, 61, 272, 1385];
        for (n, &e) in expected.iter().enumerate() {
            let alt = alternating_permutations(n as u8);
            assert_eq!(alt.len(), e, "n={}", n);
        }
    }

    #[test]
    fn test_alternating_are_updown() {
        for n in 1..=7u8 {
            for w in alternating_permutations(n) {
                for i in 0..(w.len().saturating_sub(1)) {
                    if i % 2 == 0 {
                        assert!(w[i] < w[i + 1], "n={} w={:?} not ascending at {}", n, w, i);
                    } else {
                        assert!(w[i] > w[i + 1], "n={} w={:?} not descending at {}", n, w, i);
                    }
                }
            }
        }
    }

    #[test]
    fn test_reverse_alternating_counts() {
        // Same counts as alternating
        let expected = [1, 1, 1, 2, 5, 16, 61, 272];
        for (n, &e) in expected.iter().enumerate() {
            let ralt = reverse_alternating_permutations(n as u8);
            assert_eq!(ralt.len(), e, "n={}", n);
        }
    }

    #[test]
    fn test_reverse_alternating_are_downup() {
        for n in 2..=7u8 {
            for w in reverse_alternating_permutations(n) {
                for i in 0..(w.len().saturating_sub(1)) {
                    if i % 2 == 0 {
                        assert!(w[i] > w[i + 1], "n={} w={:?} not descending at {}", n, w, i);
                    } else {
                        assert!(w[i] < w[i + 1], "n={} w={:?} not ascending at {}", n, w, i);
                    }
                }
            }
        }
    }

    // --- k-alternating permutations ---

    #[test]
    fn test_k_alt_k2_equals_alternating() {
        // k=2 should give same as alternating_permutations
        for n in 1..=8u8 {
            let k2 = k_alternating_permutations(n, 2);
            let alt = alternating_permutations(n);
            assert_eq!(k2.len(), alt.len(), "n={}", n);
        }
    }

    #[test]
    fn test_k_alt_k1() {
        // k=1: all descents, so only the decreasing permutation
        for n in 1..=6u8 {
            let k1 = k_alternating_permutations(n, 1);
            assert_eq!(k1.len(), 1, "n={}", n);
            let expected: Vec<u8> = (1..=n).rev().collect();
            assert_eq!(k1[0], expected);
        }
    }

    #[test]
    fn test_k_alt_empty_permutation() {
        assert_eq!(k_alternating_permutations(0, 1), vec![Vec::<u8>::new()]);
        assert_eq!(
            k_alternating_permutations(0, u8::MAX),
            vec![Vec::<u8>::new()]
        );
    }

    #[test]
    #[should_panic(expected = "k must be positive")]
    fn test_k_alt_rejects_zero_k() {
        k_alternating_permutations(0, 0);
    }

    #[test]
    fn test_k_alt_k3_counts() {
        // k=3: descents at positions 3, 6, 9, ...
        // n=3: asc asc desc -> σ(1)<σ(2)<σ(3): only 1,2,3 but then no descent at 3.
        // Wait: position 3 means σ(3)>σ(4), but n=3 has no position 4.
        // So for n=3, k=3: no descents required -> any increasing perm? No:
        // descent at pos 3 requires n>3. For n=3, pos 3 is the last pos, no comparison.
        // So: asc at 1, asc at 2, that's it. Only [1,2,3].
        assert_eq!(k_alternating_permutations(3, 3).len(), 1);
        // n=4, k=3: asc at 1, asc at 2, desc at 3 (σ(3)>σ(4)).
        // Count manually: σ(1)<σ(2)<σ(3)>σ(4). This is 132-shape.
        // C(3,1)=3 choices for σ(4), rest increasing. So 3.
        assert_eq!(k_alternating_permutations(4, 3).len(), 3);
    }

    #[test]
    fn test_k_alt_descent_set() {
        // Verify that k-alternating perms have correct descent set
        for k in 1..=4u8 {
            for n in 1..=7u8 {
                for w in k_alternating_permutations(n, k) {
                    for i in 0..(w.len() - 1) {
                        let pos = i + 1; // 1-indexed
                        if pos % k as usize == 0 {
                            assert!(
                                w[i] > w[i + 1],
                                "k={} n={} w={:?} expected descent at pos {}",
                                k,
                                n,
                                w,
                                pos
                            );
                        } else {
                            assert!(
                                w[i] < w[i + 1],
                                "k={} n={} w={:?} expected ascent at pos {}",
                                k,
                                n,
                                w,
                                pos
                            );
                        }
                    }
                }
            }
        }
    }
}
