//! Reversible packed-state prototype for the ordinary vertical-strip LLT/HRV identity.
//!
//! A source object is a packed strict coloring together with a subset of its
//! normal ascents.  A target object is an orientation together with a strictly
//! increasing color column on every highest-reachable-vertex fibre.  The
//! map follows the positive shifted recursion objectwise through column,
//! product, and unicellular steps.  At the first bounce seed it uses canonical
//! rank matching within the exact ordered color-content and q-degree bucket.
//! Every step has an explicit inverse and a replayable certificate; the bounce
//! rank rule is intentionally labelled as a fallback rather than presented as
//! the missing local Alexandersson--Sulzgruber branch map.

use std::collections::{BTreeMap, BTreeSet};

use clap::Parser;
use sym_poly_sym::llt::highest_reachable_vertices;
use sym_poly_sym::unit_interval_edges;

#[derive(Parser, Debug)]
#[command(about = "Prototype the ordinary marked-coloring to HRV-column bijection")]
struct Args {
    /// Exhaustively check every ordinary decorated Schroeder path through this rank.
    #[arg(long, default_value_t = 5)]
    max_n: usize,

    /// Check only this path, written with N, D, E.
    #[arg(long)]
    path: Option<String>,

    /// Print the complete canonical map for the undecorated K3 path NNNEEE.
    #[arg(long)]
    k3_map: bool,

    /// With --path, print the certificate for this packed coloring, e.g. 1,2,1.
    #[arg(long)]
    trace_coloring: Option<String>,

    /// Bit mask in the path's printed normal-edge order for the marked ascent set.
    #[arg(long, default_value_t = 0)]
    trace_mask: u64,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
enum Step {
    North,
    Diagonal,
    East,
}

impl Step {
    fn as_char(self) -> char {
        match self {
            Self::North => 'N',
            Self::Diagonal => 'D',
            Self::East => 'E',
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct PathData {
    steps: Vec<Step>,
    n: usize,
    normal_edges: Vec<(usize, usize)>,
    strict_edges: Vec<(usize, usize)>,
}

impl PathData {
    fn from_steps(steps: Vec<Step>) -> Result<Self, String> {
        let mut x = 0usize;
        let mut y = 0usize;
        let mut normal_edges = Vec::new();
        let mut strict_edges = Vec::new();
        for &step in &steps {
            match step {
                Step::North => y += 1,
                Step::Diagonal => {
                    if x >= y {
                        return Err("diagonal step is not strictly above the diagonal".into());
                    }
                    for v in x + 1..y {
                        normal_edges.push((x, v));
                    }
                    strict_edges.push((x, y));
                    x += 1;
                    y += 1;
                }
                Step::East => {
                    if x >= y {
                        return Err("east step is not strictly above the diagonal".into());
                    }
                    for v in x + 1..y {
                        normal_edges.push((x, v));
                    }
                    x += 1;
                }
            }
        }
        if x != y {
            return Err(format!("path ends at ({x},{y}), not on the diagonal"));
        }
        normal_edges.sort_unstable();
        normal_edges.dedup();
        strict_edges.sort_unstable();
        strict_edges.dedup();
        let path = Self {
            steps,
            n: x,
            normal_edges,
            strict_edges,
        };
        path.validate_library_edges()?;
        Ok(path)
    }

    fn parse(word: &str) -> Result<Self, String> {
        let steps = word
            .chars()
            .map(|ch| match ch.to_ascii_uppercase() {
                'N' => Ok(Step::North),
                'D' => Ok(Step::Diagonal),
                'E' => Ok(Step::East),
                other => Err(format!("invalid path step {other:?}")),
            })
            .collect::<Result<Vec<_>, _>>()?;
        Self::from_steps(steps)
    }

    fn word(&self) -> String {
        self.steps.iter().map(|step| step.as_char()).collect()
    }

    fn natural_edges(&self) -> Vec<(usize, usize)> {
        let mut edges = self.normal_edges.clone();
        edges.extend(self.strict_edges.iter().copied());
        edges.sort_unstable();
        edges.dedup();
        edges
    }

    fn area_sequence(&self) -> Vec<u8> {
        let natural = self.natural_edges();
        (0..self.n)
            .map(|v| {
                u8::try_from(natural.iter().filter(|&&(_, head)| head == v).count())
                    .expect("prototype ranks fit in u8")
            })
            .collect()
    }

    fn validate_library_edges(&self) -> Result<(), String> {
        let expected = unit_interval_edges(&self.area_sequence());
        let actual = self.natural_edges();
        if expected != actual {
            return Err(format!(
                "path/library edge mismatch for {}: path={actual:?}, library={expected:?}",
                self.word()
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
struct Signature {
    /// Multiplicity of each packed color, in increasing color order.
    content: Vec<usize>,
    /// Number of marked ascents / ascending normal orientation edges.
    q_degree: usize,
}

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
struct SourceState {
    coloring: Vec<u8>,
    marked_mask: u64,
}

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
struct TargetState {
    orientation_mask: u64,
    /// Canonical encoding of the fibre columns: values increase on each fibre
    /// when its vertices are read in increasing vertex order.
    column_word: Vec<u8>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct Certificate {
    path_word: String,
    signature: Signature,
    bucket_rank: usize,
    bucket_size: usize,
    source: SourceState,
    target: TargetState,
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct RecursiveCertificate {
    source: SourceState,
    target: TargetState,
    trace: Vec<String>,
}

#[derive(Clone, Debug)]
enum Reduction {
    Empty,
    Column,
    Product {
        left: PathData,
        right: PathData,
        split_vertex: usize,
    },
    Unicellular {
        strict_child: PathData,
        absent_child: PathData,
        distinguished_edge: (usize, usize),
    },
    BounceFallback,
}

#[derive(Clone, Debug)]
struct CanonicalBijection {
    source_buckets: BTreeMap<Signature, Vec<SourceState>>,
    target_buckets: BTreeMap<Signature, Vec<TargetState>>,
}

impl CanonicalBijection {
    fn build(path: &PathData) -> Result<Self, String> {
        let mut source_buckets: BTreeMap<Signature, Vec<SourceState>> = BTreeMap::new();
        for source in enumerate_sources(path) {
            let signature = source_signature(&source);
            source_buckets.entry(signature).or_default().push(source);
        }

        let mut target_buckets: BTreeMap<Signature, Vec<TargetState>> = BTreeMap::new();
        for target in enumerate_targets(path)? {
            let signature = target_signature(&target);
            target_buckets.entry(signature).or_default().push(target);
        }

        for bucket in source_buckets.values_mut() {
            bucket.sort();
        }
        for bucket in target_buckets.values_mut() {
            bucket.sort();
        }
        if source_buckets.keys().collect::<Vec<_>>() != target_buckets.keys().collect::<Vec<_>>() {
            return Err(format!(
                "signature support mismatch for path {}",
                path.word()
            ));
        }
        for (signature, sources) in &source_buckets {
            let targets = &target_buckets[signature];
            if sources.len() != targets.len() {
                return Err(format!(
                    "bucket mismatch for path {}, signature {signature:?}: {} sources, {} targets",
                    path.word(),
                    sources.len(),
                    targets.len()
                ));
            }
        }

        Ok(Self {
            source_buckets,
            target_buckets,
        })
    }

    fn map(&self, path: &PathData, source: &SourceState) -> Result<Certificate, String> {
        validate_source(path, source)?;
        let signature = source_signature(source);
        let sources = self
            .source_buckets
            .get(&signature)
            .ok_or_else(|| "missing source signature bucket".to_string())?;
        let rank = sources
            .binary_search(source)
            .map_err(|_| "source absent from canonical bucket".to_string())?;
        let targets = &self.target_buckets[&signature];
        let target = targets[rank].clone();
        validate_target(path, &target)?;
        Ok(Certificate {
            path_word: path.word(),
            signature,
            bucket_rank: rank,
            bucket_size: sources.len(),
            source: source.clone(),
            target,
        })
    }

    fn inverse(&self, path: &PathData, target: &TargetState) -> Result<SourceState, String> {
        validate_target(path, target)?;
        let signature = target_signature(target);
        let targets = self
            .target_buckets
            .get(&signature)
            .ok_or_else(|| "missing target signature bucket".to_string())?;
        let rank = targets
            .binary_search(target)
            .map_err(|_| "target absent from canonical bucket".to_string())?;
        Ok(self.source_buckets[&signature][rank].clone())
    }

    fn replay(&self, path: &PathData, certificate: &Certificate) -> Result<(), String> {
        if certificate.path_word != path.word() {
            return Err("certificate path mismatch".into());
        }
        let remapped = self.map(path, &certificate.source)?;
        if &remapped != certificate {
            return Err("certificate does not replay to the recorded target".into());
        }
        let recovered = self.inverse(path, &certificate.target)?;
        if recovered != certificate.source {
            return Err("inverse trace did not recover the source".into());
        }
        Ok(())
    }

    fn len(&self) -> usize {
        self.source_buckets.values().map(Vec::len).sum()
    }
}

#[derive(Default)]
struct RecursiveMatcher {
    fallback_cache: BTreeMap<String, CanonicalBijection>,
}

impl RecursiveMatcher {
    fn map(
        &mut self,
        path: &PathData,
        source: &SourceState,
    ) -> Result<RecursiveCertificate, String> {
        validate_source(path, source)?;
        let mut trace = Vec::new();
        let target = self.map_inner(path, source, &mut trace)?;
        validate_target(path, &target)?;
        if source_signature(source) != target_signature(&target) {
            return Err("recursive map did not preserve exact weight".into());
        }
        Ok(RecursiveCertificate {
            source: source.clone(),
            target,
            trace,
        })
    }

    fn map_inner(
        &mut self,
        path: &PathData,
        source: &SourceState,
        trace: &mut Vec<String>,
    ) -> Result<TargetState, String> {
        match reduction(path)? {
            Reduction::Empty => {
                trace.push("empty base: the empty coloring maps to the empty HRV filling".into());
                Ok(TargetState {
                    orientation_mask: 0,
                    column_word: Vec::new(),
                })
            }
            Reduction::Column => {
                trace.push(format!(
                    "column base {}: the strict coloring becomes its single increasing HRV column",
                    path.word()
                ));
                Ok(TargetState {
                    orientation_mask: 0,
                    column_word: source.coloring.clone(),
                })
            }
            Reduction::Product {
                left,
                right,
                split_vertex,
            } => {
                trace.push(format!(
                    "product {} = {} · {} at vertex split {}",
                    path.word(),
                    left.word(),
                    right.word(),
                    split_vertex
                ));
                let (left_word, left_alphabet) = standardize_word(&source.coloring[..split_vertex]);
                let (right_word, right_alphabet) =
                    standardize_word(&source.coloring[split_vertex..]);
                let left_source = SourceState {
                    coloring: left_word,
                    marked_mask: restrict_mask(path, source.marked_mask, &left, 0)?,
                };
                let right_source = SourceState {
                    coloring: right_word,
                    marked_mask: restrict_mask(path, source.marked_mask, &right, split_vertex)?,
                };
                let left_target = self.map_inner(&left, &left_source, trace)?;
                let right_target = self.map_inner(&right, &right_source, trace)?;
                let selected = shifted_selected_edges(&left, left_target.orientation_mask, 0)
                    .into_iter()
                    .chain(shifted_selected_edges(
                        &right,
                        right_target.orientation_mask,
                        split_vertex,
                    ))
                    .collect::<Vec<_>>();
                let mut column_word = lift_word(&left_target.column_word, &left_alphabet)?;
                column_word.extend(lift_word(&right_target.column_word, &right_alphabet)?);
                Ok(TargetState {
                    orientation_mask: mask_for_edges(path, &selected)?,
                    column_word,
                })
            }
            Reduction::Unicellular {
                strict_child,
                absent_child,
                distinguished_edge,
            } => {
                let selected = mask_contains_edge(path, source.marked_mask, distinguished_edge)?;
                let child = if selected {
                    trace.push(format!(
                        "shifted unicellular selected branch at {}: {} -> {}, consume one q mark",
                        format_edge(distinguished_edge),
                        path.word(),
                        strict_child.word()
                    ));
                    &strict_child
                } else {
                    trace.push(format!(
                        "shifted unicellular absent branch at {}: {} -> {}",
                        format_edge(distinguished_edge),
                        path.word(),
                        absent_child.word()
                    ));
                    &absent_child
                };
                let child_source = SourceState {
                    coloring: source.coloring.clone(),
                    marked_mask: restrict_mask(path, source.marked_mask, child, 0)?,
                };
                let child_target = self.map_inner(child, &child_source, trace)?;
                let mut selected_edges =
                    shifted_selected_edges(child, child_target.orientation_mask, 0);
                if selected {
                    selected_edges.push(distinguished_edge);
                }
                Ok(TargetState {
                    orientation_mask: mask_for_edges(path, &selected_edges)?,
                    column_word: child_target.column_word,
                })
            }
            Reduction::BounceFallback => {
                let key = path.word();
                if !self.fallback_cache.contains_key(&key) {
                    self.fallback_cache
                        .insert(key.clone(), CanonicalBijection::build(path)?);
                }
                let fallback = &self.fallback_cache[&key];
                let certificate = fallback.map(path, source)?;
                trace.push(format!(
                    "bounce seed {}: canonical exact-bucket fallback rank {}/{} (recursion lift not yet supplied)",
                    key,
                    certificate.bucket_rank + 1,
                    certificate.bucket_size
                ));
                Ok(certificate.target)
            }
        }
    }

    fn inverse(&mut self, path: &PathData, target: &TargetState) -> Result<SourceState, String> {
        validate_target(path, target)?;
        let source = self.inverse_inner(path, target)?;
        validate_source(path, &source)?;
        if source_signature(&source) != target_signature(target) {
            return Err("recursive inverse did not preserve exact weight".into());
        }
        Ok(source)
    }

    fn inverse_inner(
        &mut self,
        path: &PathData,
        target: &TargetState,
    ) -> Result<SourceState, String> {
        match reduction(path)? {
            Reduction::Empty | Reduction::Column => Ok(SourceState {
                coloring: target.column_word.clone(),
                marked_mask: 0,
            }),
            Reduction::Product {
                left,
                right,
                split_vertex,
            } => {
                let (left_word, left_alphabet) =
                    standardize_word(&target.column_word[..split_vertex]);
                let (right_word, right_alphabet) =
                    standardize_word(&target.column_word[split_vertex..]);
                let left_target = TargetState {
                    orientation_mask: restrict_mask(path, target.orientation_mask, &left, 0)?,
                    column_word: left_word,
                };
                let right_target = TargetState {
                    orientation_mask: restrict_mask(
                        path,
                        target.orientation_mask,
                        &right,
                        split_vertex,
                    )?,
                    column_word: right_word,
                };
                let left_source = self.inverse_inner(&left, &left_target)?;
                let right_source = self.inverse_inner(&right, &right_target)?;
                let mut coloring = lift_word(&left_source.coloring, &left_alphabet)?;
                coloring.extend(lift_word(&right_source.coloring, &right_alphabet)?);
                let selected = shifted_selected_edges(&left, left_source.marked_mask, 0)
                    .into_iter()
                    .chain(shifted_selected_edges(
                        &right,
                        right_source.marked_mask,
                        split_vertex,
                    ))
                    .collect::<Vec<_>>();
                Ok(SourceState {
                    coloring,
                    marked_mask: mask_for_edges(path, &selected)?,
                })
            }
            Reduction::Unicellular {
                strict_child,
                absent_child,
                distinguished_edge,
            } => {
                let selected =
                    mask_contains_edge(path, target.orientation_mask, distinguished_edge)?;
                let child = if selected {
                    &strict_child
                } else {
                    &absent_child
                };
                let child_target = TargetState {
                    orientation_mask: restrict_mask(path, target.orientation_mask, child, 0)?,
                    column_word: target.column_word.clone(),
                };
                let child_source = self.inverse_inner(child, &child_target)?;
                let mut selected_edges = shifted_selected_edges(child, child_source.marked_mask, 0);
                if selected {
                    selected_edges.push(distinguished_edge);
                }
                Ok(SourceState {
                    coloring: child_source.coloring,
                    marked_mask: mask_for_edges(path, &selected_edges)?,
                })
            }
            Reduction::BounceFallback => {
                let key = path.word();
                if !self.fallback_cache.contains_key(&key) {
                    self.fallback_cache
                        .insert(key.clone(), CanonicalBijection::build(path)?);
                }
                self.fallback_cache[&key].inverse(path, target)
            }
        }
    }

    fn replay(
        &mut self,
        path: &PathData,
        certificate: &RecursiveCertificate,
    ) -> Result<(), String> {
        let remapped = self.map(path, &certificate.source)?;
        if remapped != *certificate {
            return Err("recursive certificate trace did not replay".into());
        }
        if self.inverse(path, &certificate.target)? != certificate.source {
            return Err("recursive inverse did not recover the source".into());
        }
        Ok(())
    }
}

fn reduction(path: &PathData) -> Result<Reduction, String> {
    if path.n == 0 {
        return Ok(Reduction::Empty);
    }
    let east_positions = path
        .steps
        .iter()
        .enumerate()
        .filter_map(|(index, step)| (*step == Step::East).then_some(index))
        .collect::<Vec<_>>();
    if east_positions.len() == 1 {
        return Ok(Reduction::Column);
    }
    let first_east = *east_positions
        .first()
        .ok_or_else(|| "positive-rank path has no east step".to_string())?;
    let (x, y) = endpoint(&path.steps[..first_east]);
    if y == x + 1 {
        let left = PathData::from_steps(path.steps[..=first_east].to_vec())?;
        let right = PathData::from_steps(path.steps[first_east + 1..].to_vec())?;
        return Ok(Reduction::Product {
            split_vertex: left.n,
            left,
            right,
        });
    }
    if first_east > 0 && path.steps[first_east - 1] == Step::North {
        let prefix = &path.steps[..first_east - 1];
        let suffix = &path.steps[first_east + 1..];
        let (u, v) = endpoint(prefix);
        let mut strict_steps = prefix.to_vec();
        strict_steps.push(Step::Diagonal);
        strict_steps.extend_from_slice(suffix);
        let mut absent_steps = prefix.to_vec();
        absent_steps.extend([Step::East, Step::North]);
        absent_steps.extend_from_slice(suffix);
        return Ok(Reduction::Unicellular {
            strict_child: PathData::from_steps(strict_steps)?,
            absent_child: PathData::from_steps(absent_steps)?,
            distinguished_edge: (u, v),
        });
    }
    Ok(Reduction::BounceFallback)
}

fn endpoint(steps: &[Step]) -> (usize, usize) {
    steps
        .iter()
        .fold((0usize, 0usize), |(x, y), step| match step {
            Step::North => (x, y + 1),
            Step::Diagonal => (x + 1, y + 1),
            Step::East => (x + 1, y),
        })
}

fn standardize_word(word: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let alphabet = word
        .iter()
        .copied()
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>();
    let ranks = alphabet
        .iter()
        .enumerate()
        .map(|(rank, &color)| (color, rank as u8))
        .collect::<BTreeMap<_, _>>();
    (word.iter().map(|color| ranks[color]).collect(), alphabet)
}

fn lift_word(word: &[u8], alphabet: &[u8]) -> Result<Vec<u8>, String> {
    word.iter()
        .map(|&color| {
            alphabet
                .get(color as usize)
                .copied()
                .ok_or_else(|| "standardized color is outside its recorded alphabet".to_string())
        })
        .collect()
}

fn shifted_selected_edges(path: &PathData, mask: u64, shift: usize) -> Vec<(usize, usize)> {
    path.normal_edges
        .iter()
        .enumerate()
        .filter(|(index, _)| mask & (1u64 << index) != 0)
        .map(|(_, &(u, v))| (u + shift, v + shift))
        .collect()
}

fn mask_for_edges(path: &PathData, selected: &[(usize, usize)]) -> Result<u64, String> {
    let selected = selected.iter().copied().collect::<BTreeSet<_>>();
    if selected
        .iter()
        .any(|edge| path.normal_edges.binary_search(edge).is_err())
    {
        return Err(format!(
            "selected edge is not normal in path {}",
            path.word()
        ));
    }
    Ok(path
        .normal_edges
        .iter()
        .enumerate()
        .fold(0u64, |mask, (index, edge)| {
            if selected.contains(edge) {
                mask | (1u64 << index)
            } else {
                mask
            }
        }))
}

fn restrict_mask(
    parent: &PathData,
    parent_mask: u64,
    child: &PathData,
    child_shift: usize,
) -> Result<u64, String> {
    let selected = shifted_selected_edges(parent, parent_mask, 0)
        .into_iter()
        .filter_map(|(u, v)| {
            (u >= child_shift && v >= child_shift).then_some((u - child_shift, v - child_shift))
        })
        .filter(|edge| child.normal_edges.binary_search(edge).is_ok())
        .collect::<Vec<_>>();
    mask_for_edges(child, &selected)
}

fn mask_contains_edge(path: &PathData, mask: u64, edge: (usize, usize)) -> Result<bool, String> {
    let index = path.normal_edges.binary_search(&edge).map_err(|_| {
        format!(
            "distinguished edge {edge:?} is not normal in {}",
            path.word()
        )
    })?;
    Ok(mask & (1u64 << index) != 0)
}

fn format_edge((u, v): (usize, usize)) -> String {
    format!("{}{}", u + 1, v + 1)
}

fn enumerate_sources(path: &PathData) -> Vec<SourceState> {
    let mut result = Vec::new();
    for coloring in packed_words(path.n) {
        if !increases_on_edges(&coloring, &path.strict_edges) {
            continue;
        }
        let ascent_mask = edge_ascent_mask(&coloring, &path.normal_edges);
        let mut marked_mask = ascent_mask;
        loop {
            result.push(SourceState {
                coloring: coloring.clone(),
                marked_mask,
            });
            if marked_mask == 0 {
                break;
            }
            marked_mask = (marked_mask - 1) & ascent_mask;
        }
    }
    result
}

fn enumerate_targets(path: &PathData) -> Result<Vec<TargetState>, String> {
    if path.normal_edges.len() >= u64::BITS as usize {
        return Err("too many normal edges for the prototype mask".into());
    }
    let words = packed_words(path.n);
    let mut result = Vec::new();
    for orientation_mask in 0..(1u64 << path.normal_edges.len()) {
        let hrv = hrv_for_mask(path, orientation_mask)?;
        let fibres = fibres(&hrv);
        for column_word in &words {
            if fibres.iter().all(|fibre| {
                fibre
                    .windows(2)
                    .all(|vertices| column_word[vertices[0]] < column_word[vertices[1]])
            }) {
                result.push(TargetState {
                    orientation_mask,
                    column_word: column_word.clone(),
                });
            }
        }
    }
    Ok(result)
}

fn validate_source(path: &PathData, source: &SourceState) -> Result<(), String> {
    validate_packed_word(&source.coloring, path.n)?;
    if !increases_on_edges(&source.coloring, &path.strict_edges) {
        return Err("source coloring violates a strict edge".into());
    }
    let ascents = edge_ascent_mask(&source.coloring, &path.normal_edges);
    if source.marked_mask & !ascents != 0 {
        return Err("source marks a normal edge that is not a coloring ascent".into());
    }
    Ok(())
}

fn validate_target(path: &PathData, target: &TargetState) -> Result<(), String> {
    validate_packed_word(&target.column_word, path.n)?;
    if target.orientation_mask >> path.normal_edges.len() != 0 {
        return Err("target orientation mask has an out-of-range edge".into());
    }
    let hrv = hrv_for_mask(path, target.orientation_mask)?;
    for fibre in fibres(&hrv) {
        if !fibre
            .windows(2)
            .all(|vertices| target.column_word[vertices[0]] < target.column_word[vertices[1]])
        {
            return Err("target column is not strictly increasing".into());
        }
    }
    Ok(())
}

fn source_signature(source: &SourceState) -> Signature {
    Signature {
        content: color_content(&source.coloring),
        q_degree: source.marked_mask.count_ones() as usize,
    }
}

fn target_signature(target: &TargetState) -> Signature {
    Signature {
        content: color_content(&target.column_word),
        q_degree: target.orientation_mask.count_ones() as usize,
    }
}

fn hrv_for_mask(path: &PathData, mask: u64) -> Result<Vec<usize>, String> {
    let mut upward_edges = path.strict_edges.clone();
    upward_edges.extend(
        path.normal_edges
            .iter()
            .enumerate()
            .filter(|(index, _)| mask & (1u64 << index) != 0)
            .map(|(_, &edge)| edge),
    );
    highest_reachable_vertices(path.n, &upward_edges)
        .ok_or_else(|| "ordinary upward graph unexpectedly invalid or cyclic".to_string())
}

fn fibres(hrv: &[usize]) -> Vec<Vec<usize>> {
    let mut by_sink: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (vertex, &sink) in hrv.iter().enumerate() {
        by_sink.entry(sink).or_default().push(vertex);
    }
    by_sink.into_values().collect()
}

fn edge_ascent_mask(word: &[u8], edges: &[(usize, usize)]) -> u64 {
    edges
        .iter()
        .enumerate()
        .fold(0u64, |mask, (index, &(u, v))| {
            if word[u] < word[v] {
                mask | (1u64 << index)
            } else {
                mask
            }
        })
}

fn increases_on_edges(word: &[u8], edges: &[(usize, usize)]) -> bool {
    edges.iter().all(|&(u, v)| word[u] < word[v])
}

fn color_content(word: &[u8]) -> Vec<usize> {
    let max_color = word
        .iter()
        .copied()
        .max()
        .map_or(0, |color| color as usize + 1);
    let mut content = vec![0usize; max_color];
    for &color in word {
        content[color as usize] += 1;
    }
    content
}

fn validate_packed_word(word: &[u8], n: usize) -> Result<(), String> {
    if word.len() != n {
        return Err(format!("word length {} differs from rank {n}", word.len()));
    }
    if word.is_empty() {
        return Ok(());
    }
    let colors: BTreeSet<_> = word.iter().copied().collect();
    let expected: BTreeSet<_> = (0..colors.len() as u8).collect();
    if colors != expected {
        return Err("word is not order-packed".into());
    }
    Ok(())
}

fn packed_words(n: usize) -> Vec<Vec<u8>> {
    if n == 0 {
        return vec![Vec::new()];
    }
    let mut result = Vec::new();
    for colors in 1..=n {
        let mut word = vec![0u8; n];
        enumerate_words_with_alphabet(0, colors as u8, &mut word, &mut result);
    }
    result
}

fn enumerate_words_with_alphabet(
    index: usize,
    colors: u8,
    word: &mut [u8],
    output: &mut Vec<Vec<u8>>,
) {
    if index == word.len() {
        let mut seen = vec![false; colors as usize];
        for &color in word.iter() {
            seen[color as usize] = true;
        }
        if seen.into_iter().all(|value| value) {
            output.push(word.to_vec());
        }
        return;
    }
    for color in 0..colors {
        word[index] = color;
        enumerate_words_with_alphabet(index + 1, colors, word, output);
    }
}

fn all_paths(n: usize) -> Vec<PathData> {
    let mut paths = Vec::new();
    let mut steps = Vec::new();
    enumerate_paths(n, 0, 0, &mut steps, &mut paths);
    paths
}

fn enumerate_paths(
    n: usize,
    x: usize,
    y: usize,
    steps: &mut Vec<Step>,
    output: &mut Vec<PathData>,
) {
    if x == n && y == n {
        output.push(PathData::from_steps(steps.clone()).expect("generated path is valid"));
        return;
    }
    if y < n {
        steps.push(Step::North);
        enumerate_paths(n, x, y + 1, steps, output);
        steps.pop();
    }
    if x < n && y < n && x < y {
        steps.push(Step::Diagonal);
        enumerate_paths(n, x + 1, y + 1, steps, output);
        steps.pop();
    }
    if x < n && x < y {
        steps.push(Step::East);
        enumerate_paths(n, x + 1, y, steps, output);
        steps.pop();
    }
}

fn same_orientation_obstruction(
    path: &PathData,
) -> Result<Option<(Signature, u64, usize, usize)>, String> {
    let mut source_counts: BTreeMap<(Signature, u64), usize> = BTreeMap::new();
    for source in enumerate_sources(path) {
        *source_counts
            .entry((source_signature(&source), source.marked_mask))
            .or_default() += 1;
    }
    let mut target_counts: BTreeMap<(Signature, u64), usize> = BTreeMap::new();
    for target in enumerate_targets(path)? {
        *target_counts
            .entry((target_signature(&target), target.orientation_mask))
            .or_default() += 1;
    }
    let keys: BTreeSet<_> = source_counts
        .keys()
        .chain(target_counts.keys())
        .cloned()
        .collect();
    Ok(keys.into_iter().find_map(|(signature, mask)| {
        let sources = source_counts
            .get(&(signature.clone(), mask))
            .copied()
            .unwrap_or(0);
        let targets = target_counts
            .get(&(signature.clone(), mask))
            .copied()
            .unwrap_or(0);
        (sources != targets).then_some((signature, mask, sources, targets))
    }))
}

fn checked_path(path: &PathData) -> Result<(usize, usize), String> {
    let bijection = CanonicalBijection::build(path)?;
    let all_targets = enumerate_targets(path)?
        .into_iter()
        .collect::<BTreeSet<_>>();
    let mut recursive_targets = BTreeSet::new();
    let mut recursive_matcher = RecursiveMatcher::default();
    let mut certificates = 0usize;
    for sources in bijection.source_buckets.values() {
        for source in sources {
            let certificate = bijection.map(path, source)?;
            if source_signature(&certificate.source) != target_signature(&certificate.target) {
                return Err("certificate failed weight preservation".into());
            }
            bijection.replay(path, &certificate)?;
            let recursive_certificate = recursive_matcher.map(path, source)?;
            recursive_matcher.replay(path, &recursive_certificate)?;
            if !recursive_targets.insert(recursive_certificate.target) {
                return Err("hybrid recursive map is not injective".into());
            }
            certificates += 1;
        }
    }
    if certificates != bijection.len() {
        return Err("certificate traversal count mismatch".into());
    }
    if recursive_targets != all_targets {
        return Err("hybrid recursive map is not surjective onto the HRV-column targets".into());
    }
    Ok((bijection.source_buckets.len(), certificates))
}

fn format_word(word: &[u8]) -> String {
    if word.is_empty() {
        return "-".into();
    }
    word.iter()
        .map(|color| (usize::from(*color) + 1).to_string())
        .collect::<Vec<_>>()
        .join("")
}

fn format_edges(mask: u64, edges: &[(usize, usize)]) -> String {
    let selected = edges
        .iter()
        .enumerate()
        .filter(|(index, _)| mask & (1u64 << index) != 0)
        .map(|(_, &(u, v))| format!("{}{}", u + 1, v + 1))
        .collect::<Vec<_>>();
    if selected.is_empty() {
        "-".into()
    } else {
        selected.join(",")
    }
}

fn format_columns(path: &PathData, target: &TargetState) -> Result<String, String> {
    let hrv = hrv_for_mask(path, target.orientation_mask)?;
    let rendered = fibres(&hrv)
        .into_iter()
        .map(|fibre| {
            let sink = hrv[fibre[0]] + 1;
            let colors = fibre
                .iter()
                .map(|&vertex| (usize::from(target.column_word[vertex]) + 1).to_string())
                .collect::<Vec<_>>()
                .join("<");
            format!("{sink}:[{colors}]")
        })
        .collect::<Vec<_>>()
        .join(" ");
    Ok(rendered)
}

fn print_k3_map() -> Result<(), String> {
    let path = PathData::parse("NNNEEE")?;
    let bijection = CanonicalBijection::build(&path)?;
    let mut recursive_matcher = RecursiveMatcher::default();
    println!("K3 normal-edge order: {:?}", path.normal_edges);
    println!("source coloring\tmarked A\tqdeg\tcontent\ttarget eta ascents\tHRV columns\ttrace");
    for (signature, sources) in &bijection.source_buckets {
        for source in sources {
            let certificate = recursive_matcher.map(&path, source)?;
            println!(
                "{}\t{}\t{}\t{:?}\t{}\t{}\t{}",
                format_word(&source.coloring),
                format_edges(source.marked_mask, &path.normal_edges),
                signature.q_degree,
                signature.content,
                format_edges(certificate.target.orientation_mask, &path.normal_edges),
                format_columns(&path, &certificate.target)?,
                certificate.trace.join(" | "),
            );
        }
    }
    if let Some((signature, mask, sources, targets)) = same_orientation_obstruction(&path)? {
        println!(
            "K3 keep-orientation obstruction: signature={signature:?}, eta={}, source_count={sources}, target_count={targets}",
            format_edges(mask, &path.normal_edges)
        );
    }
    Ok(())
}

fn parse_coloring(text: &str) -> Result<Vec<u8>, String> {
    text.split(',')
        .map(|part| {
            let value = part
                .trim()
                .parse::<u8>()
                .map_err(|error| format!("invalid color {part:?}: {error}"))?;
            value
                .checked_sub(1)
                .ok_or_else(|| "colors must be positive".to_string())
        })
        .collect()
}

fn print_trace(path: &PathData, coloring_text: &str, marked_mask: u64) -> Result<(), String> {
    let source = SourceState {
        coloring: parse_coloring(coloring_text)?,
        marked_mask,
    };
    let mut matcher = RecursiveMatcher::default();
    let certificate = matcher.map(path, &source)?;
    matcher.replay(path, &certificate)?;
    println!("path: {}", path.word());
    println!("normal-edge order: {:?}", path.normal_edges);
    println!(
        "state expansion: packed coloring {} with marked ascents {}",
        format_word(&source.coloring),
        format_edges(source.marked_mask, &path.normal_edges)
    );
    println!(
        "exact invariant: content {:?}, q-degree {}",
        source_signature(&certificate.source).content,
        source_signature(&certificate.source).q_degree
    );
    for (index, step) in certificate.trace.iter().enumerate() {
        println!("trace {}: {step}", index + 1);
    }
    println!(
        "target: eta ascents {}; HRV columns {}",
        format_edges(certificate.target.orientation_mask, &path.normal_edges),
        format_columns(path, &certificate.target)?
    );
    println!("reverse certificate: target -> same bucket/rank -> original source [OK]");
    Ok(())
}

fn run(args: Args) -> Result<(), String> {
    if args.k3_map {
        print_k3_map()?;
    }
    if let Some(word) = args.path {
        let path = PathData::parse(&word)?;
        if let Some(coloring) = args.trace_coloring.as_deref() {
            print_trace(&path, coloring, args.trace_mask)?;
        } else {
            let (buckets, certificates) = checked_path(&path)?;
            println!(
                "checked path {}: rank {}, {} signatures, {} reversible certificates",
                path.word(),
                path.n,
                buckets,
                certificates
            );
        }
        return Ok(());
    }

    let mut total_paths = 0usize;
    let mut total_buckets = 0usize;
    let mut total_certificates = 0usize;
    for n in 0..=args.max_n {
        let paths = all_paths(n);
        let mut rank_certificates = 0usize;
        for path in &paths {
            let (buckets, certificates) = checked_path(path)?;
            total_buckets += buckets;
            total_certificates += certificates;
            rank_certificates += certificates;
        }
        total_paths += paths.len();
        println!(
            "rank {n}: {} paths, {} reversible certificates",
            paths.len(),
            rank_certificates
        );
    }
    println!(
        "total: {total_paths} ordinary decorated paths, {total_buckets} exact invariant buckets, {total_certificates} reversible certificates"
    );
    Ok(())
}

fn main() {
    if let Err(error) = run(Args::parse()) {
        eprintln!("error: {error}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn k3_map_is_reversible_and_weight_preserving() {
        let path = PathData::parse("NNNEEE").unwrap();
        let bijection = CanonicalBijection::build(&path).unwrap();
        assert!(!bijection.source_buckets.is_empty());
        for sources in bijection.source_buckets.values() {
            for source in sources {
                let certificate = bijection.map(&path, source).unwrap();
                assert_eq!(
                    source_signature(&certificate.source),
                    target_signature(&certificate.target)
                );
                bijection.replay(&path, &certificate).unwrap();
            }
        }
    }

    #[test]
    fn k3_hybrid_recursion_is_a_bijection() {
        let path = PathData::parse("NNNEEE").unwrap();
        checked_path(&path).unwrap();
    }

    #[test]
    fn k3_requires_mixing_orientations() {
        let path = PathData::parse("NNNEEE").unwrap();
        assert!(same_orientation_obstruction(&path).unwrap().is_some());
    }

    #[test]
    fn every_ordinary_decorated_path_through_rank_four_checks() {
        for n in 0..=4 {
            for path in all_paths(n) {
                checked_path(&path).unwrap_or_else(|error| panic!("{}: {error}", path.word()));
            }
        }
    }
}
