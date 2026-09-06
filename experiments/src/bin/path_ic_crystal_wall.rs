//! Exact path-IC wall certificates for the flagged-UIG atom/crystal project.
//!
//! Project-agnostic inverse `P`-RS tracing lives in `sym-poly-core::p_rs`.
//! This binary supplies the path order, the intrinsic ternary insertion-word
//! grammar, the four-color duplicate carry, and the local excluded-wall test.

use clap::{Parser, ValueEnum};
use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use sym_poly_core::{
    inverse_column_step_with_trace, reverse_complement, InverseColumnEvent, InverseColumnEventKind,
    InverseColumnStep, PInsertionOrder, Partition, Tableau,
};

#[derive(Debug, Clone, Copy, ValueEnum)]
enum OutputFormat {
    Text,
    Json,
}

#[derive(Debug, Parser)]
#[command(about = "Certify the four-color path-IC inverse-P-RS wall")]
struct Arguments {
    /// Largest path order to check (the calibrated exact range ends at 17).
    #[arg(long, default_value_t = 17)]
    max_n: usize,

    /// Standard-output format.
    #[arg(long, value_enum, default_value_t = OutputFormat::Text)]
    format: OutputFormat,

    /// Optional detailed per-length-type TSV output.
    #[arg(long)]
    tsv_out: Option<PathBuf>,

    /// Also certify insertion-word independence of vertices 1, 2, and 3.
    #[arg(long)]
    verify_frozen_core: bool,
}

#[derive(Debug, Clone, Copy)]
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

type Column = Vec<u32>;
type Insertion = Vec<Column>;
type InverseTrace = (Vec<Option<u32>>, Vec<InverseColumnStep>);
type DynError = Box<dyn Error>;

#[derive(Debug, Clone)]
struct RejectedTableau {
    tableau: Tableau,
    raw_recording: Insertion,
    number_of_ones: usize,
    reason: RejectionReason,
}

#[derive(Debug, Clone, Copy)]
enum RejectionReason {
    RowTwoStartsThree,
    RowOneHasMultipleTwos,
}

#[derive(Debug, Default, Clone)]
struct TypeSummary {
    n: usize,
    s: usize,
    r: usize,
    d: usize,
    insertion_words: usize,
    excluded_tableaux: usize,
    excluded_incidences: usize,
    row_two_starts_three: usize,
    row_one_has_multiple_twos: usize,
    ladder_copy: usize,
    ladder_move: usize,
    starts_at_boundary: usize,
    starts_before_boundary: usize,
    boundary_signature_classes: usize,
    merge_patterns: BTreeSet<String>,
    double_source_merges: usize,
    q_coefficient: String,
}

#[derive(Debug, Default)]
struct TotalSummary {
    length_types: usize,
    insertion_words: usize,
    excluded_tableaux: usize,
    excluded_incidences: usize,
    ladder_copy: usize,
    ladder_move: usize,
    starts_at_boundary: usize,
    starts_before_boundary: usize,
    merge_pattern_classes: usize,
    double_source_merges: usize,
}

#[derive(Debug, Default)]
struct FrozenCoreSummary {
    length_types: usize,
    insertion_words: usize,
    source_tableaux: usize,
    repair_tableaux: usize,
    source_incidences: usize,
    repair_incidences: usize,
    positioned_core_checks: usize,
    core_emission_checks: usize,
    standard_max_n: usize,
    source_standard_tableaux: usize,
    repair_standard_tableaux: usize,
    source_standard_incidences: usize,
    repair_standard_incidences: usize,
    local_move_edges: usize,
    local_move_types: BTreeSet<LocalMove>,
    max_canonical_distance: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct FrozenCoreSignature {
    positioned_core: [usize; 3],
    emissions: Vec<(u32, usize, InverseColumnEventKind, usize, usize)>,
}

type LocalMove = (Vec<u8>, Vec<u8>);

#[derive(Debug)]
struct LocalMoveGraphSummary {
    edges: usize,
    move_types: BTreeSet<LocalMove>,
    max_canonical_distance: usize,
}

fn intrinsic_words(s: usize, r: usize, d: usize) -> Vec<Vec<u8>> {
    if d == 0 {
        return Vec::new();
    }
    let target = [s - 1, r - 1, d];
    let mut counts = [0, 0, 1];
    let mut word = vec![2];
    let mut output = Vec::new();
    intrinsic_words_recursive(target, &mut counts, &mut word, &mut output);
    output
}

fn intrinsic_words_recursive(
    target: [usize; 3],
    counts: &mut [usize; 3],
    word: &mut Vec<u8>,
    output: &mut Vec<Vec<u8>>,
) {
    if counts == &target {
        output.push(word.clone());
        return;
    }
    let choices: &[u8] = if counts[2] > counts[1] {
        &[1]
    } else if counts[1] > counts[0] {
        &[0]
    } else {
        &[0, 1, 2]
    };
    let last = *word.last().expect("an intrinsic word starts with 2");
    for &letter in choices {
        let index = letter as usize;
        if letter == last || counts[index] == target[index] {
            continue;
        }
        counts[index] += 1;
        word.push(letter);
        intrinsic_words_recursive(target, counts, word, output);
        word.pop();
        counts[index] -= 1;
    }
}

fn canonical_word(s: usize, r: usize, d: usize) -> Vec<u8> {
    let t = s - r;
    let p = r + t - d - 1;
    let u = d - t - 1;
    let mut word = vec![2, 1, 0];
    for _ in 0..p {
        word.extend([1, 0]);
    }
    for _ in 0..t {
        word.extend([2, 0]);
    }
    for _ in 0..u {
        word.extend([2, 1, 0]);
    }
    word
}

fn local_move_window(left: &[u8], right: &[u8]) -> Option<LocalMove> {
    let differing: Vec<_> = left
        .iter()
        .zip(right)
        .enumerate()
        .filter_map(|(index, (a, b))| (a != b).then_some(index))
        .collect();
    let (&first, &last) = (differing.first()?, differing.last()?);
    if first == 0 || last - first + 1 > 3 {
        return None;
    }
    let mut left_sorted = left[first..=last].to_vec();
    let mut right_sorted = right[first..=last].to_vec();
    left_sorted.sort_unstable();
    right_sorted.sort_unstable();
    if left_sorted != right_sorted {
        return None;
    }
    let blocks = (left[first..=last].to_vec(), right[first..=last].to_vec());
    Some(if blocks.0 < blocks.1 {
        blocks
    } else {
        (blocks.1, blocks.0)
    })
}

fn local_move_graph(
    words: &[Vec<u8>],
    canonical: &[u8],
) -> Result<LocalMoveGraphSummary, DynError> {
    let root = words
        .iter()
        .position(|word| word == canonical)
        .ok_or_else(|| format!("canonical word {canonical:?} is not intrinsic"))?;
    let mut adjacency = vec![Vec::new(); words.len()];
    let mut move_types = BTreeSet::new();
    let mut edges = 0;
    for left in 0..words.len() {
        for right in left + 1..words.len() {
            if let Some(move_type) = local_move_window(&words[left], &words[right]) {
                adjacency[left].push(right);
                adjacency[right].push(left);
                move_types.insert(move_type);
                edges += 1;
            }
        }
    }
    let mut distances = vec![None; words.len()];
    distances[root] = Some(0usize);
    let mut queue = VecDeque::from([root]);
    while let Some(current) = queue.pop_front() {
        let next_distance = distances[current].expect("a queued word has a distance") + 1;
        for &next in &adjacency[current] {
            if distances[next].is_none() {
                distances[next] = Some(next_distance);
                queue.push_back(next);
            }
        }
    }
    if distances.iter().any(Option::is_none) {
        return Err(
            format!("intrinsic words are disconnected from canonical word {canonical:?}").into(),
        );
    }
    Ok(LocalMoveGraphSummary {
        edges,
        move_types,
        max_canonical_distance: distances.into_iter().flatten().max().unwrap_or(0),
    })
}

fn repair_insertion(word: &[u8], n: usize) -> Insertion {
    let mut groups = vec![vec![1], Vec::new(), Vec::new()];
    for (label, &group) in (4..=n as u32).zip(word) {
        groups[group as usize].push(label);
    }
    for group in &mut groups {
        group.sort_unstable_by(|left, right| right.cmp(left));
    }
    vec![
        groups[0].clone(),
        groups[1].iter().copied().chain([3]).collect(),
        groups[2].iter().copied().chain([2]).collect(),
    ]
}

fn source_insertion(word: &[u8], n: usize) -> Insertion {
    let mut groups = vec![vec![1], Vec::new(), Vec::new()];
    for (label, &group) in (4..=n as u32).zip(word) {
        groups[group as usize].push(label);
    }
    for group in &mut groups {
        group.sort_unstable_by(|left, right| right.cmp(left));
    }
    vec![
        groups[0].clone(),
        groups[1].iter().copied().chain([2]).collect(),
        groups[2].clone(),
        vec![3],
    ]
}

fn duplicate_carry(tableau: &Tableau) -> (Tableau, Vec<u32>) {
    let mut rows = tableau.rows().to_vec();
    assert_eq!(rows.len(), 3);
    let mut moved = rows[2].pop().expect("the third row is nonempty");
    let mut trace = vec![moved];
    rows.push(vec![moved]);
    let mut lower = 2;
    while lower > 0 && rows[lower].contains(&moved) {
        let upper = lower - 1;
        let replacement = rows[upper]
            .iter()
            .copied()
            .filter(|&entry| entry < moved)
            .max()
            .expect("a semistandard carry has a smaller upper entry");
        let upper_index = rows[upper]
            .iter()
            .rposition(|&entry| entry == replacement)
            .expect("the chosen replacement occurs");
        let lower_index = rows[lower]
            .iter()
            .position(|&entry| entry == moved)
            .expect("the duplicated carried entry occurs");
        rows[upper][upper_index] = moved;
        rows[lower][lower_index] = replacement;
        rows[upper].sort_unstable();
        rows[lower].sort_unstable();
        moved = replacement;
        trace.push(moved);
        lower = upper;
    }
    (Tableau::new(rows), trace)
}

fn is_pre_base_candidate(tableau: &Tableau, d: usize) -> Option<Vec<u32>> {
    if tableau.right_key_weight_via_ssaf(4).parts()[0] != 0 || tableau.rows()[2].last() != Some(&4)
    {
        return None;
    }
    let (carried, trace) = duplicate_carry(tableau);
    if !matches!(trace.as_slice(), [4] | [4, 3] | [4, 3, 2])
        || !carried.is_semistandard()
        || carried.right_key_weight_via_ssaf(4).parts()[0] != d as u32
    {
        return None;
    }
    Some(trace)
}

fn rejected_tableaux(s: usize, r: usize, d: usize) -> Vec<RejectedTableau> {
    let shape = Partition::from_sorted(vec![s as u32, r as u32, d as u32 + 1]);
    Tableau::semistandard_tableaux(&shape, 4)
        .into_iter()
        .filter_map(|tableau| {
            let trace = is_pre_base_candidate(&tableau, d)?;
            if trace != [4, 3, 2] {
                return None;
            }
            let row_two_starts_two = tableau.rows()[1][0] == 2;
            let row_one_twos = tableau.rows()[0]
                .iter()
                .filter(|&&entry| entry == 2)
                .count();
            if row_two_starts_two && row_one_twos <= 1 {
                return None;
            }
            let reason = if !row_two_starts_two {
                RejectionReason::RowTwoStartsThree
            } else {
                RejectionReason::RowOneHasMultipleTwos
            };
            let number_of_ones = tableau
                .rows()
                .iter()
                .flatten()
                .filter(|&&entry| entry == 1)
                .count();
            Some(RejectedTableau {
                raw_recording: raw_recording(&tableau),
                tableau,
                number_of_ones,
                reason,
            })
        })
        .collect()
}

fn raw_recording(tableau: &Tableau) -> Insertion {
    let transposed = transpose_rows(tableau.rows());
    let mut counts = BTreeMap::new();
    for &entry in transposed.iter().flatten() {
        *counts.entry(entry).or_insert(0u32) += 1;
    }
    let mut next_label = 1u32;
    let mut starts = BTreeMap::new();
    for (&entry, &count) in &counts {
        starts.insert(entry, next_label);
        next_label += count;
    }
    let mut seen = BTreeMap::new();
    let standardized_rows: Vec<Vec<u32>> = transposed
        .iter()
        .map(|row| {
            row.iter()
                .map(|entry| {
                    let used = seen.entry(*entry).or_insert(0u32);
                    let label = starts[entry] + *used;
                    *used += 1;
                    label
                })
                .collect()
        })
        .collect();
    let evacuated = Tableau::new(standardized_rows).evacuation();
    rows_to_columns(evacuated.rows())
}

fn transpose_rows(rows: &[Vec<u32>]) -> Vec<Vec<u32>> {
    let width = rows.first().map_or(0, Vec::len);
    (0..width)
        .map(|column| {
            rows.iter()
                .filter_map(|row| row.get(column).copied())
                .collect()
        })
        .collect()
}

fn rows_to_columns(rows: &[Vec<u32>]) -> Insertion {
    let width = rows.first().map_or(0, Vec::len);
    (0..width)
        .map(|column| {
            let mut entries: Vec<u32> = rows
                .iter()
                .filter_map(|row| row.get(column).copied())
                .collect();
            entries.reverse();
            entries
        })
        .collect()
}

fn inverse_with_trace(
    insertion: &Insertion,
    recording: &Insertion,
    n: usize,
) -> Result<InverseTrace, DynError> {
    if insertion.len() != recording.len() {
        return Err(format!(
            "insertion and recording column counts differ: {} != {}, insertion={insertion:?}, recording={recording:?}",
            insertion.len(),
            recording.len()
        )
        .into());
    }
    let mut current = vec![None; n];
    let mut stages = Vec::new();
    for column in (0..insertion.len()).rev() {
        let drag_positions: BTreeSet<usize> = recording[column]
            .iter()
            .map(|&label| n + 1 - label as usize)
            .collect();
        let source = reverse_complement(&current, n as u32)?;
        let column_word: Vec<Option<u32>> = insertion[column].iter().copied().map(Some).collect();
        let chain = reverse_complement(&column_word, n as u32)?
            .into_iter()
            .collect::<Option<Vec<_>>>()
            .expect("an insertion column has no empty entries");
        let stage = inverse_column_step_with_trace(&source, &chain, &drag_positions, &PathOrder)?;
        if !stage.remaining_chain.is_empty() {
            return Err(format!(
                "inverse column {column} left active chain {:?}",
                stage.remaining_chain
            )
            .into());
        }
        current = reverse_complement(&stage.previous_word, n as u32)?;
        stages.push(stage);
    }
    Ok((current, stages))
}

fn frozen_core_signature(
    insertion: &Insertion,
    recording: &Insertion,
    n: usize,
    core_columns: &[(usize, u32); 3],
) -> Result<FrozenCoreSignature, DynError> {
    let (inverse_word, stages) = inverse_with_trace(insertion, recording, n)?;
    let mut positioned_core = [0; 3];
    for vertex in 1..=3u32 {
        positioned_core[vertex as usize - 1] = inverse_word
            .iter()
            .position(|&entry| entry == Some(vertex))
            .map(|position| position + 1)
            .ok_or_else(|| format!("inverse word does not contain core vertex {vertex}"))?;
    }

    let mut emissions = Vec::new();
    for &(column, vertex) in core_columns {
        let stage_index = insertion.len() - 1 - column;
        let distinguished = n as u32 + 1 - vertex;
        let event = stages[stage_index]
            .events
            .iter()
            .find(|event| event.emitted.contains(&distinguished))
            .ok_or_else(|| {
                format!("column {column} never emits complemented core vertex {distinguished}")
            })?;
        emissions.push((
            vertex,
            column,
            event.kind,
            event.first_position,
            event.last_position,
        ));
    }
    emissions.sort_unstable();
    Ok(FrozenCoreSignature {
        positioned_core,
        emissions,
    })
}

fn verify_frozen_core_tableaux(
    n: usize,
    shape: &Partition,
    tableaux: &[Tableau],
    insertions: &[Insertion],
    core_columns: &[(usize, u32); 3],
    side: &str,
) -> Result<usize, DynError> {
    for tableau in tableaux {
        let recording = raw_recording(tableau);
        let expected = frozen_core_signature(&insertions[0], &recording, n, core_columns)?;
        for insertion in insertions.iter().skip(1) {
            let actual = frozen_core_signature(insertion, &recording, n, core_columns)?;
            if actual != expected {
                return Err(format!(
                    "frozen-core failure on {side}: n={n}, shape={:?}, tableau={:?}, expected={expected:?}, actual={actual:?}",
                    shape.parts(),
                    tableau.rows(),
                )
                .into());
            }
        }
    }
    Ok(tableaux.len() * insertions.len())
}

fn verify_frozen_core(max_n: usize) -> Result<FrozenCoreSummary, DynError> {
    let mut summary = FrozenCoreSummary {
        standard_max_n: max_n.min(11),
        ..FrozenCoreSummary::default()
    };
    for n in 6..=max_n {
        for s in 2..n {
            for r in 2..=s {
                let Some(d) = n.checked_sub(1 + s + r) else {
                    continue;
                };
                if d == 0 || d >= r {
                    continue;
                }
                let words = intrinsic_words(s, r, d);
                if words.is_empty() {
                    continue;
                }
                let local_moves = local_move_graph(&words, &canonical_word(s, r, d))?;
                let sources: Vec<_> = words.iter().map(|word| source_insertion(word, n)).collect();
                let repairs: Vec<_> = words.iter().map(|word| repair_insertion(word, n)).collect();
                let source_shape = Partition::from_sorted(vec![s as u32, r as u32, d as u32, 1]);
                let repair_shape = Partition::from_sorted(vec![s as u32, r as u32, d as u32 + 1]);
                let source_tableaux = Tableau::semistandard_tableaux(&source_shape, 4);
                let repair_tableaux = Tableau::semistandard_tableaux(&repair_shape, 4);
                let source_incidences = verify_frozen_core_tableaux(
                    n,
                    &source_shape,
                    &source_tableaux,
                    &sources,
                    &[(0, 1), (1, 2), (3, 3)],
                    "source",
                )?;
                let repair_incidences = verify_frozen_core_tableaux(
                    n,
                    &repair_shape,
                    &repair_tableaux,
                    &repairs,
                    &[(0, 1), (1, 3), (2, 2)],
                    "repair",
                )?;
                summary.length_types += 1;
                summary.insertion_words += words.len();
                summary.local_move_edges += local_moves.edges;
                summary.local_move_types.extend(local_moves.move_types);
                summary.max_canonical_distance = summary
                    .max_canonical_distance
                    .max(local_moves.max_canonical_distance);
                summary.source_tableaux += source_tableaux.len();
                summary.repair_tableaux += repair_tableaux.len();
                summary.source_incidences += source_incidences;
                summary.repair_incidences += repair_incidences;

                if n <= summary.standard_max_n {
                    let source_standard = Tableau::standard_tableaux(&source_shape);
                    let repair_standard = Tableau::standard_tableaux(&repair_shape);
                    summary.source_standard_incidences += verify_frozen_core_tableaux(
                        n,
                        &source_shape,
                        &source_standard,
                        &sources,
                        &[(0, 1), (1, 2), (3, 3)],
                        "standard source",
                    )?;
                    summary.repair_standard_incidences += verify_frozen_core_tableaux(
                        n,
                        &repair_shape,
                        &repair_standard,
                        &repairs,
                        &[(0, 1), (1, 3), (2, 2)],
                        "standard repair",
                    )?;
                    summary.source_standard_tableaux += source_standard.len();
                    summary.repair_standard_tableaux += repair_standard.len();
                }
            }
        }
    }
    summary.positioned_core_checks = summary.source_incidences + summary.repair_incidences;
    summary.core_emission_checks = 3 * summary.positioned_core_checks;
    if max_n == 17 {
        assert_eq!(summary.length_types, 27);
        assert_eq!(summary.insertion_words, 282);
        assert_eq!(summary.source_tableaux, 6_915);
        assert_eq!(summary.repair_tableaux, 9_991);
        assert_eq!(summary.source_incidences, 91_115);
        assert_eq!(summary.repair_incidences, 107_697);
        assert_eq!(summary.positioned_core_checks, 198_812);
        assert_eq!(summary.core_emission_checks, 596_436);
        assert_eq!(summary.standard_max_n, 11);
        assert_eq!(summary.source_standard_tableaux, 2_621);
        assert_eq!(summary.repair_standard_tableaux, 1_013);
        assert_eq!(summary.source_standard_incidences, 6_749);
        assert_eq!(summary.repair_standard_incidences, 2_441);
        assert_eq!(summary.local_move_edges, 528);
        assert_eq!(summary.local_move_types.len(), 12);
        assert_eq!(summary.max_canonical_distance, 6);
    }
    Ok(summary)
}

fn q_coefficient(words: &[Vec<u8>]) -> String {
    let mut counts = BTreeMap::new();
    for word in words {
        let ascents = word.windows(2).filter(|pair| pair[0] < pair[1]).count();
        *counts.entry(2 + ascents).or_insert(0usize) += 1;
    }
    counts
        .into_iter()
        .map(|(degree, coefficient)| match (coefficient, degree) {
            (1, 0) => "1".to_string(),
            (1, 1) => "q".to_string(),
            (1, degree) => format!("q^{degree}"),
            (coefficient, 0) => coefficient.to_string(),
            (coefficient, 1) => format!("{coefficient}*q"),
            (coefficient, degree) => format!("{coefficient}*q^{degree}"),
        })
        .collect::<Vec<_>>()
        .join("+")
}

fn ladder_merge_pattern(
    event: &InverseColumnEvent,
    distinguished_maximum: u32,
) -> Result<String, DynError> {
    let first_source = *event
        .source
        .first()
        .ok_or_else(|| "a ladder event has no source letter".to_string())?;
    let row = event
        .active_before
        .iter()
        .rposition(|&active| active < first_source)
        .ok_or_else(|| "a ladder source has no active predecessor".to_string())?;
    if event.active_before[row] + 1 != first_source {
        return Err(format!(
            "ladder does not start with adjacent active/source letters: {event:?}"
        )
        .into());
    }

    let height = (0..event.active_before.len() - row)
        .filter(|&height| {
            let mut values = event.active_before[row..=row + height].to_vec();
            values.extend(&event.source);
            PathOrder.is_ladder(&values)
        })
        .max()
        .ok_or_else(|| "the traced ladder has no participating active interval".to_string())?;
    let participating_active = &event.active_before[row..=row + height];
    if participating_active.contains(&distinguished_maximum) {
        return Err("the distinguished maximum participates in the boundary ladder".into());
    }

    let mut tagged_values: Vec<(u32, char)> = participating_active
        .iter()
        .copied()
        .map(|value| (value, 'A'))
        .chain(event.source.iter().copied().map(|value| (value, 'S')))
        .collect();
    tagged_values.sort_unstable();
    if tagged_values
        .windows(2)
        .any(|pair| pair[0].0 + 1 != pair[1].0)
    {
        return Err("the normalized boundary merge is not consecutive".into());
    }
    Ok(tagged_values.into_iter().map(|(_, tag)| tag).collect())
}

fn check_merge_language(
    pattern: &str,
    reason: RejectionReason,
    kind: InverseColumnEventKind,
) -> Result<bool, DynError> {
    if !pattern.starts_with("AS") || pattern.as_bytes().windows(2).any(|pair| pair == b"AA") {
        return Err(format!("invalid active/source ladder language: {pattern}").into());
    }
    let double_sources = pattern
        .as_bytes()
        .windows(2)
        .filter(|pair| *pair == b"SS")
        .count();
    match reason {
        RejectionReason::RowTwoStartsThree if double_sources != 0 => {
            return Err(format!("row-two-starts-three merge is not alternating: {pattern}").into());
        }
        RejectionReason::RowOneHasMultipleTwos if double_sources > 1 => {
            return Err(format!("multiple-two merge has two source doublings: {pattern}").into());
        }
        _ => {}
    }
    let expected_last = match kind {
        InverseColumnEventKind::LadderCopy => 'A',
        InverseColumnEventKind::LadderMove => 'S',
        _ => return Err("merge-language check requires a ladder event".into()),
    };
    if !pattern.ends_with(expected_last) {
        return Err(format!(
            "ladder branch and final merge letter disagree: kind={kind:?}, pattern={pattern}"
        )
        .into());
    }
    Ok(double_sources == 1)
}

fn verify_type(
    n: usize,
    s: usize,
    r: usize,
    d: usize,
    words: &[Vec<u8>],
) -> Result<TypeSummary, DynError> {
    let rejected = rejected_tableaux(s, r, d);
    let mut summary = TypeSummary {
        n,
        s,
        r,
        d,
        insertion_words: words.len(),
        excluded_tableaux: rejected.len(),
        q_coefficient: q_coefficient(words),
        ..TypeSummary::default()
    };
    let mut signatures = BTreeSet::new();

    for word in words {
        let insertion = repair_insertion(word, n);
        for rejected_tableau in &rejected {
            let (inverse_word, stages) =
                inverse_with_trace(&insertion, &rejected_tableau.raw_recording, n)?;
            let final_stage = stages.last().expect("a repair insertion has columns");
            let boundary = rejected_tableau.number_of_ones;
            let certificate = final_stage
                .boundary_certificate(boundary, n as u32)
                .ok_or_else(|| {
                    format!(
                        "no boundary event: n={n}, word={word:?}, tableau={:?}, boundary={boundary}",
                        rejected_tableau.tableau.rows()
                    )
                })?;
            let position = inverse_word
                .iter()
                .position(|&entry| entry == Some(1))
                .map(|index| index + 1)
                .ok_or_else(|| "inverse word does not contain vertex 1".to_string())?;

            if !certificate.event_kind.is_ladder()
                || !certificate.active_before
                || !certificate.active_after
                || certificate.emitted
                || position <= boundary
            {
                return Err(format!(
                    "wall failure: n={n}, shape=({s},{r},{}), word={word:?}, tableau={:?}, position={position}, boundary={boundary}, certificate={certificate:?}",
                    d + 1,
                    rejected_tableau.tableau.rows()
                )
                .into());
            }
            let event = &final_stage.events[certificate.event_index];
            let merge_pattern = ladder_merge_pattern(event, n as u32)?;
            let has_double_source = check_merge_language(
                &merge_pattern,
                rejected_tableau.reason,
                certificate.event_kind,
            )?;

            summary.excluded_incidences += 1;
            summary.double_source_merges += usize::from(has_double_source);
            summary.merge_patterns.insert(merge_pattern);
            match rejected_tableau.reason {
                RejectionReason::RowTwoStartsThree => summary.row_two_starts_three += 1,
                RejectionReason::RowOneHasMultipleTwos => {
                    summary.row_one_has_multiple_twos += 1;
                }
            }
            match certificate.event_kind {
                InverseColumnEventKind::LadderCopy => summary.ladder_copy += 1,
                InverseColumnEventKind::LadderMove => summary.ladder_move += 1,
                _ => unreachable!("the wall check accepts only ladder events"),
            }
            if certificate.event_first_position == boundary {
                summary.starts_at_boundary += 1;
            } else {
                assert!(certificate.event_first_position < boundary);
                summary.starts_before_boundary += 1;
            }
            signatures.insert((
                certificate.event_kind,
                certificate.event_first_position as isize - boundary as isize,
                certificate.event_last_position as isize - boundary as isize,
            ));
        }
    }
    summary.boundary_signature_classes = signatures.len();
    if summary.excluded_incidences != words.len() * rejected.len() {
        return Err("excluded incidence count is not a Cartesian product".into());
    }
    Ok(summary)
}

fn write_tsv(path: &PathBuf, summaries: &[TypeSummary]) -> Result<(), DynError> {
    let mut output = BufWriter::new(File::create(path)?);
    writeln!(
        output,
        "n\ts\tr\td\tinsertion_words\texcluded_tableaux\texcluded_incidences\trow_two_starts_three\trow_one_has_multiple_twos\tladder_copy\tladder_move\tstarts_at_boundary\tstarts_before_boundary\tboundary_signature_classes\tmerge_pattern_classes\tdouble_source_merges\tq_coefficient\tmaximum_active_through_boundary\tposition_strictly_delayed\tfailures"
    )?;
    for row in summaries {
        writeln!(
            output,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\ttrue\ttrue\t0",
            row.n,
            row.s,
            row.r,
            row.d,
            row.insertion_words,
            row.excluded_tableaux,
            row.excluded_incidences,
            row.row_two_starts_three,
            row.row_one_has_multiple_twos,
            row.ladder_copy,
            row.ladder_move,
            row.starts_at_boundary,
            row.starts_before_boundary,
            row.boundary_signature_classes,
            row.merge_patterns.len(),
            row.double_source_merges,
            row.q_coefficient,
        )?;
    }
    Ok(())
}

fn print_json(
    total: &TotalSummary,
    summaries: &[TypeSummary],
    frozen_core: Option<&FrozenCoreSummary>,
) {
    println!("{{");
    println!("  \"length_types\": {},", total.length_types);
    println!("  \"insertion_words\": {},", total.insertion_words);
    println!("  \"excluded_tableaux\": {},", total.excluded_tableaux);
    println!("  \"excluded_incidences\": {},", total.excluded_incidences);
    println!("  \"ladder_copy\": {},", total.ladder_copy);
    println!("  \"ladder_move\": {},", total.ladder_move);
    println!(
        "  \"merge_pattern_classes\": {},",
        total.merge_pattern_classes
    );
    println!(
        "  \"double_source_merges\": {},",
        total.double_source_merges
    );
    println!("  \"maximum_active_through_boundary\": true,");
    println!("  \"position_strictly_delayed\": true,");
    println!("  \"types\": [");
    for (index, row) in summaries.iter().enumerate() {
        let comma = if index + 1 == summaries.len() {
            ""
        } else {
            ","
        };
        println!(
            "    {{\"n\":{},\"s\":{},\"r\":{},\"d\":{},\"insertion_words\":{},\"excluded_tableaux\":{},\"excluded_incidences\":{},\"ladder_copy\":{},\"ladder_move\":{},\"merge_pattern_classes\":{},\"double_source_merges\":{},\"q_coefficient\":\"{}\"}}{}",
            row.n,
            row.s,
            row.r,
            row.d,
            row.insertion_words,
            row.excluded_tableaux,
            row.excluded_incidences,
            row.ladder_copy,
            row.ladder_move,
            row.merge_patterns.len(),
            row.double_source_merges,
            row.q_coefficient,
            comma,
        );
    }
    println!("  ],");
    if let Some(core) = frozen_core {
        println!("  \"frozen_core\": {{");
        println!("    \"length_types\": {},", core.length_types);
        println!("    \"insertion_words\": {},", core.insertion_words);
        println!("    \"source_tableaux\": {},", core.source_tableaux);
        println!("    \"repair_tableaux\": {},", core.repair_tableaux);
        println!("    \"source_incidences\": {},", core.source_incidences);
        println!("    \"repair_incidences\": {},", core.repair_incidences);
        println!(
            "    \"positioned_core_checks\": {},",
            core.positioned_core_checks
        );
        println!(
            "    \"core_emission_checks\": {},",
            core.core_emission_checks
        );
        println!("    \"standard_max_n\": {},", core.standard_max_n);
        println!(
            "    \"source_standard_tableaux\": {},",
            core.source_standard_tableaux
        );
        println!(
            "    \"repair_standard_tableaux\": {},",
            core.repair_standard_tableaux
        );
        println!(
            "    \"source_standard_incidences\": {},",
            core.source_standard_incidences
        );
        println!(
            "    \"repair_standard_incidences\": {},",
            core.repair_standard_incidences
        );
        println!("    \"local_move_edges\": {},", core.local_move_edges);
        println!("    \"local_move_types\": {},", core.local_move_types.len());
        println!(
            "    \"max_canonical_distance\": {}",
            core.max_canonical_distance
        );
        println!("  }},");
    }
    println!("  \"failures\": 0");
    println!("}}");
}

fn main() -> Result<(), DynError> {
    let arguments = Arguments::parse();
    if !(6..=17).contains(&arguments.max_n) {
        return Err("--max-n must lie in 6..=17 for the calibrated scan".into());
    }

    let mut summaries = Vec::new();
    let mut total = TotalSummary::default();
    let mut merge_patterns = BTreeSet::new();
    for n in 6..=arguments.max_n {
        for s in 2..n {
            for r in 2..=s {
                let Some(d) = n.checked_sub(1 + s + r) else {
                    continue;
                };
                if d == 0 || d >= r {
                    continue;
                }
                let words = intrinsic_words(s, r, d);
                if words.is_empty() {
                    continue;
                }
                let summary = verify_type(n, s, r, d, &words)?;
                total.length_types += 1;
                total.insertion_words += summary.insertion_words;
                total.excluded_tableaux += summary.excluded_tableaux;
                total.excluded_incidences += summary.excluded_incidences;
                total.ladder_copy += summary.ladder_copy;
                total.ladder_move += summary.ladder_move;
                total.starts_at_boundary += summary.starts_at_boundary;
                total.starts_before_boundary += summary.starts_before_boundary;
                total.double_source_merges += summary.double_source_merges;
                merge_patterns.extend(summary.merge_patterns.iter().cloned());
                summaries.push(summary);
            }
        }
    }
    total.merge_pattern_classes = merge_patterns.len();

    if arguments.max_n == 17 {
        assert_eq!(total.length_types, 27);
        assert_eq!(total.insertion_words, 282);
        assert_eq!(total.excluded_tableaux, 989);
        assert_eq!(total.excluded_incidences, 7852);
        assert_eq!(total.ladder_copy, 5060);
        assert_eq!(total.ladder_move, 2792);
        assert_eq!(total.starts_at_boundary, 4002);
        assert_eq!(total.starts_before_boundary, 3850);
        assert_eq!(total.merge_pattern_classes, 27);
        assert_eq!(total.double_source_merges, 90);
    }

    let frozen_core = arguments
        .verify_frozen_core
        .then(|| verify_frozen_core(arguments.max_n))
        .transpose()?;

    if let Some(path) = &arguments.tsv_out {
        write_tsv(path, &summaries)?;
    }
    match arguments.format {
        OutputFormat::Text => {
            println!(
                "path_ic_crystal_wall\tlength_types={}\tinsertion_words={}\texcluded_tableaux={}\texcluded_incidences={}\tladder_copy_move={}/{}\tdouble_source_merges={}\tmaximum_active_through_boundary=all\tposition_strictly_delayed=all\tfailures=0",
                total.length_types,
                total.insertion_words,
                total.excluded_tableaux,
                total.excluded_incidences,
                total.ladder_copy,
                total.ladder_move,
                total.double_source_merges,
            );
            if let Some(core) = &frozen_core {
                println!(
                    "path_ic_frozen_core\tlength_types={}\tinsertion_words={}\tsource_tableaux={}\trepair_tableaux={}\tsource_incidences={}\trepair_incidences={}\tpositioned_core_checks={}\tcore_emission_checks={}\tstandard_max_n={}\tsource_standard_tableaux={}\trepair_standard_tableaux={}\tsource_standard_incidences={}\trepair_standard_incidences={}\tlocal_move_edges={}\tlocal_move_types={}\tmax_canonical_distance={}\tfailures=0",
                    core.length_types,
                    core.insertion_words,
                    core.source_tableaux,
                    core.repair_tableaux,
                    core.source_incidences,
                    core.repair_incidences,
                    core.positioned_core_checks,
                    core.core_emission_checks,
                    core.standard_max_n,
                    core.source_standard_tableaux,
                    core.repair_standard_tableaux,
                    core.source_standard_incidences,
                    core.repair_standard_incidences,
                    core.local_move_edges,
                    core.local_move_types.len(),
                    core.max_canonical_distance,
                );
            }
        }
        OutputFormat::Json => print_json(&total, &summaries, frozen_core.as_ref()),
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn smallest_intrinsic_word_builds_the_expected_repair() {
        let words = intrinsic_words(2, 2, 1);
        assert_eq!(words, vec![vec![2, 1, 0]]);
        assert_eq!(
            repair_insertion(&words[0], 6),
            vec![vec![6, 1], vec![5, 3], vec![4, 2]]
        );
    }

    #[test]
    fn raw_recording_matches_the_smallest_excluded_tableau() {
        let tableau = Tableau::new(vec![vec![1, 2], vec![3, 3], vec![4, 4]]);
        assert_eq!(
            raw_recording(&tableau),
            vec![vec![2, 1], vec![4, 3], vec![6, 5]]
        );
    }

    #[test]
    fn smallest_excluded_wall_is_a_ladder_move() {
        let words = intrinsic_words(2, 2, 1);
        let summary = verify_type(6, 2, 2, 1, &words).unwrap();
        assert_eq!(summary.excluded_tableaux, 1);
        assert_eq!(summary.excluded_incidences, 1);
        assert_eq!(summary.ladder_copy, 0);
        assert_eq!(summary.ladder_move, 1);
        assert_eq!(summary.starts_at_boundary, 1);
        assert_eq!(summary.q_coefficient, "q^2");
        assert_eq!(summary.merge_patterns.len(), 1);
        assert_eq!(summary.double_source_merges, 0);
    }

    #[test]
    fn diagonal_ic_words_have_the_same_positioned_core() {
        let words = intrinsic_words(3, 3, 2);
        assert_eq!(words.len(), 2);
        let tableau = Tableau::new(vec![vec![1, 1, 1], vec![2, 2, 2], vec![3, 3, 3]]);
        let recording = raw_recording(&tableau);

        let repair_signatures: Vec<_> = words
            .iter()
            .map(|word| {
                frozen_core_signature(
                    &repair_insertion(word, 9),
                    &recording,
                    9,
                    &[(0, 1), (1, 3), (2, 2)],
                )
                .unwrap()
            })
            .collect();
        assert_eq!(repair_signatures[0], repair_signatures[1]);
        assert_eq!(repair_signatures[0].positioned_core, [3, 9, 6]);

        let source_tableau = Tableau::new(vec![vec![1, 1, 1], vec![2, 2, 2], vec![3, 3], vec![4]]);
        let source_recording = raw_recording(&source_tableau);
        let source_signatures: Vec<_> = words
            .iter()
            .map(|word| {
                frozen_core_signature(
                    &source_insertion(word, 9),
                    &source_recording,
                    9,
                    &[(0, 1), (1, 2), (3, 3)],
                )
                .unwrap()
            })
            .collect();
        assert_eq!(source_signatures[0], source_signatures[1]);
    }
}
