//! Reversible packed-state prototype for the ordinary vertical-strip LLT/HRV identity.
//!
//! A source object is a packed strict coloring together with a subset of its
//! normal ascents.  A target object is an orientation together with a strictly
//! increasing color column on every highest-reachable-vertex fibre.  The
//! map follows the positive shifted recursion objectwise.  Single NN and DN
//! bounces use the published coloring involutions and explicit orientation
//! correspondences.  Generalized NN/DN/ND bounces use canonical rank matching
//! only between the exact smaller recursion branches, because the published
//! many-bounce proof proceeds by cancellation rather than an object map.  Every
//! step has an explicit inverse and a replayable certificate.

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

    /// Print the complete recursive map for the undecorated K3 path NNNEEE.
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
    FirstBounce {
        successor: PathData,
        base: PathData,
        low: usize,
        mid: usize,
        top: usize,
    },
    SecondBounce {
        successor: PathData,
        base: PathData,
        corner: usize,
        low: usize,
        mid: usize,
        top: usize,
    },
    GeneralizedBounce(GeneralizedBounce),
    BounceFallback,
}

#[derive(Clone, Debug)]
enum GeneralizedBounce {
    First {
        successor: PathData,
    },
    Second {
        successor: PathData,
    },
    Third {
        selected: PathData,
        optional: PathData,
    },
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
enum GeneralizedBranch {
    First,
    Second,
    ThirdSelected,
    ThirdOptional,
}

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
struct BranchSource {
    branch: GeneralizedBranch,
    factor: bool,
    state: SourceState,
}

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
struct BranchTarget {
    branch: GeneralizedBranch,
    factor: bool,
    state: TargetState,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum FibreTransport {
    Identity,
    Swap(usize, usize),
}

#[derive(Clone, Copy, Debug)]
struct BounceVertices {
    low: usize,
    mid: usize,
    top: usize,
}

#[derive(Clone, Debug)]
struct CanonicalBijection {
    source_buckets: BTreeMap<Signature, Vec<SourceState>>,
    target_buckets: BTreeMap<Signature, Vec<TargetState>>,
}

#[derive(Clone, Debug)]
struct GeneralizedBijection {
    reduction: GeneralizedBounce,
    parent_sources: BTreeMap<Signature, Vec<SourceState>>,
    branch_sources: BTreeMap<Signature, Vec<BranchSource>>,
    parent_targets: BTreeMap<Signature, Vec<TargetState>>,
    branch_targets: BTreeMap<Signature, Vec<BranchTarget>>,
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

impl GeneralizedBounce {
    fn path_for(&self, branch: GeneralizedBranch) -> Result<&PathData, String> {
        match (self, branch) {
            (Self::First { successor }, GeneralizedBranch::First)
            | (Self::Second { successor }, GeneralizedBranch::Second) => Ok(successor),
            (Self::Third { selected, .. }, GeneralizedBranch::ThirdSelected) => Ok(selected),
            (Self::Third { optional, .. }, GeneralizedBranch::ThirdOptional) => Ok(optional),
            _ => Err("generalized-bounce branch does not belong to this reduction".into()),
        }
    }

    fn source_objects(&self) -> Vec<BranchSource> {
        match self {
            Self::First { successor } => enumerate_sources(successor)
                .into_iter()
                .flat_map(|state| {
                    [false, true].into_iter().map(move |factor| BranchSource {
                        branch: GeneralizedBranch::First,
                        factor,
                        state: state.clone(),
                    })
                })
                .collect(),
            Self::Second { successor } => enumerate_sources(successor)
                .into_iter()
                .map(|state| BranchSource {
                    branch: GeneralizedBranch::Second,
                    factor: false,
                    state,
                })
                .collect(),
            Self::Third { selected, optional } => {
                let mut objects = enumerate_sources(selected)
                    .into_iter()
                    .map(|state| BranchSource {
                        branch: GeneralizedBranch::ThirdSelected,
                        factor: false,
                        state,
                    })
                    .collect::<Vec<_>>();
                objects.extend(enumerate_sources(optional).into_iter().flat_map(|state| {
                    [false, true].into_iter().map(move |factor| BranchSource {
                        branch: GeneralizedBranch::ThirdOptional,
                        factor,
                        state: state.clone(),
                    })
                }));
                objects
            }
        }
    }

    fn target_objects(&self) -> Result<Vec<BranchTarget>, String> {
        Ok(match self {
            Self::First { successor } => enumerate_targets(successor)?
                .into_iter()
                .flat_map(|state| {
                    [false, true].into_iter().map(move |factor| BranchTarget {
                        branch: GeneralizedBranch::First,
                        factor,
                        state: state.clone(),
                    })
                })
                .collect(),
            Self::Second { successor } => enumerate_targets(successor)?
                .into_iter()
                .map(|state| BranchTarget {
                    branch: GeneralizedBranch::Second,
                    factor: false,
                    state,
                })
                .collect(),
            Self::Third { selected, optional } => {
                let mut objects = enumerate_targets(selected)?
                    .into_iter()
                    .map(|state| BranchTarget {
                        branch: GeneralizedBranch::ThirdSelected,
                        factor: false,
                        state,
                    })
                    .collect::<Vec<_>>();
                objects.extend(enumerate_targets(optional)?.into_iter().flat_map(|state| {
                    [false, true].into_iter().map(move |factor| BranchTarget {
                        branch: GeneralizedBranch::ThirdOptional,
                        factor,
                        state: state.clone(),
                    })
                }));
                objects
            }
        })
    }
}

fn generalized_offset(branch: GeneralizedBranch, factor: bool) -> usize {
    match branch {
        GeneralizedBranch::First | GeneralizedBranch::ThirdOptional => usize::from(factor),
        GeneralizedBranch::Second => 0,
        GeneralizedBranch::ThirdSelected => 1,
    }
}

fn branch_source_signature(branch: &BranchSource) -> Signature {
    let mut signature = source_signature(&branch.state);
    signature.q_degree += generalized_offset(branch.branch, branch.factor);
    signature
}

fn branch_target_signature(branch: &BranchTarget) -> Signature {
    let mut signature = target_signature(&branch.state);
    signature.q_degree += generalized_offset(branch.branch, branch.factor);
    signature
}

fn sorted_buckets<T: Ord>(
    objects: impl IntoIterator<Item = T>,
    signature: impl Fn(&T) -> Signature,
) -> BTreeMap<Signature, Vec<T>> {
    let mut buckets: BTreeMap<Signature, Vec<T>> = BTreeMap::new();
    for object in objects {
        buckets.entry(signature(&object)).or_default().push(object);
    }
    for bucket in buckets.values_mut() {
        bucket.sort();
    }
    buckets
}

fn validate_bucket_sizes<A, B>(
    path: &PathData,
    side: &str,
    left: &BTreeMap<Signature, Vec<A>>,
    right: &BTreeMap<Signature, Vec<B>>,
) -> Result<(), String> {
    let keys = left
        .keys()
        .chain(right.keys())
        .cloned()
        .collect::<BTreeSet<_>>();
    for signature in keys {
        let left_size = left.get(&signature).map_or(0, Vec::len);
        let right_size = right.get(&signature).map_or(0, Vec::len);
        if left_size != right_size {
            return Err(format!(
                "generalized-bounce {side} obstruction for {}, signature {signature:?}: parent {left_size}, branches {right_size}",
                path.word()
            ));
        }
    }
    Ok(())
}

impl GeneralizedBijection {
    fn build(path: &PathData, reduction: GeneralizedBounce) -> Result<Self, String> {
        let parent_sources = sorted_buckets(enumerate_sources(path), source_signature);
        let branch_sources = sorted_buckets(reduction.source_objects(), branch_source_signature);
        validate_bucket_sizes(path, "source", &parent_sources, &branch_sources)?;
        let parent_targets = sorted_buckets(enumerate_targets(path)?, target_signature);
        let branch_targets = sorted_buckets(reduction.target_objects()?, branch_target_signature);
        validate_bucket_sizes(path, "target", &parent_targets, &branch_targets)?;
        Ok(Self {
            reduction,
            parent_sources,
            branch_sources,
            parent_targets,
            branch_targets,
        })
    }

    fn source_to_branch(
        &self,
        source: &SourceState,
    ) -> Result<(BranchSource, usize, usize), String> {
        let signature = source_signature(source);
        let parents = self
            .parent_sources
            .get(&signature)
            .ok_or_else(|| "missing generalized parent-source bucket".to_string())?;
        let rank = parents
            .binary_search(source)
            .map_err(|_| "source absent from generalized parent bucket".to_string())?;
        Ok((
            self.branch_sources[&signature][rank].clone(),
            rank,
            parents.len(),
        ))
    }

    fn branch_to_source(&self, branch: &BranchSource) -> Result<SourceState, String> {
        let signature = branch_source_signature(branch);
        let branches = self
            .branch_sources
            .get(&signature)
            .ok_or_else(|| "missing generalized branch-source bucket".to_string())?;
        let rank = branches
            .binary_search(branch)
            .map_err(|_| "branch source absent from generalized bucket".to_string())?;
        Ok(self.parent_sources[&signature][rank].clone())
    }

    fn target_to_branch(&self, target: &TargetState) -> Result<BranchTarget, String> {
        let signature = target_signature(target);
        let parents = self
            .parent_targets
            .get(&signature)
            .ok_or_else(|| "missing generalized parent-target bucket".to_string())?;
        let rank = parents
            .binary_search(target)
            .map_err(|_| "target absent from generalized parent bucket".to_string())?;
        Ok(self.branch_targets[&signature][rank].clone())
    }

    fn branch_to_target(&self, branch: &BranchTarget) -> Result<TargetState, String> {
        let signature = branch_target_signature(branch);
        let branches = self
            .branch_targets
            .get(&signature)
            .ok_or_else(|| "missing generalized branch-target bucket".to_string())?;
        let rank = branches
            .binary_search(branch)
            .map_err(|_| "branch target absent from generalized bucket".to_string())?;
        Ok(self.parent_targets[&signature][rank].clone())
    }
}

#[derive(Default)]
struct RecursiveMatcher {
    generalized_cache: BTreeMap<String, GeneralizedBijection>,
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
            Reduction::FirstBounce {
                successor,
                base,
                low,
                mid,
                top,
            } => {
                let vertices = BounceVertices { low, mid, top };
                let (factor, successor_source, coloring_transport) =
                    first_bounce_coloring_forward(path, &successor, &base, vertices, source)?;
                trace.push(format!(
                    "first bounce {} -> {}: paper coloring {}, external q-bit {}; recurse on successor",
                    path.word(),
                    successor.word(),
                    coloring_transport,
                    usize::from(factor)
                ));
                let successor_target = self.map_inner(&successor, &successor_source, trace)?;
                let successor_background =
                    restrict_mask(&successor, successor_target.orientation_mask, &base, 0)?;
                let (source_bit, source_background) =
                    bounce_unmatch(&base, low, mid, top, factor, successor_background)?;
                let source_orientation =
                    insert_orientation_bit(path, &base, source_background, (mid, top), source_bit)?;
                let fibre_transport = detect_fibre_transport(
                    path,
                    source_orientation,
                    &successor,
                    successor_target.orientation_mask,
                    low,
                    mid,
                )?;
                let column_word = transport_columns(
                    &successor,
                    successor_target.orientation_mask,
                    &successor_target.column_word,
                    path,
                    source_orientation,
                    fibre_transport,
                )?;
                trace.push(format!(
                    "first-bounce orientation unmatch: {} fibre conjugacy",
                    fibre_transport.describe()
                ));
                Ok(TargetState {
                    orientation_mask: source_orientation,
                    column_word,
                })
            }
            Reduction::SecondBounce {
                successor,
                base,
                corner,
                low,
                mid,
                top,
            } => {
                let vertices = BounceVertices { low, mid, top };
                let (successor_source, coloring_transport) = second_bounce_coloring_forward(
                    path, &successor, &base, corner, vertices, source,
                )?;
                trace.push(format!(
                    "second bounce {} -> {}: paper coloring {}; recurse on successor",
                    path.word(),
                    successor.word(),
                    coloring_transport
                ));
                let successor_target = self.map_inner(&successor, &successor_source, trace)?;
                let successor_bit = mask_contains_edge(
                    &successor,
                    successor_target.orientation_mask,
                    (corner, low),
                )?;
                let successor_background =
                    restrict_mask(&successor, successor_target.orientation_mask, &base, 0)?;
                let (source_bit, source_background) =
                    second_unmatch(&base, low, mid, top, successor_bit, successor_background)?;
                let source_orientation =
                    insert_orientation_bit(path, &base, source_background, (mid, top), source_bit)?;
                let fibre_transport = detect_fibre_transport(
                    path,
                    source_orientation,
                    &successor,
                    successor_target.orientation_mask,
                    low,
                    mid,
                )?;
                let column_word = transport_columns(
                    &successor,
                    successor_target.orientation_mask,
                    &successor_target.column_word,
                    path,
                    source_orientation,
                    fibre_transport,
                )?;
                trace.push(format!(
                    "second-bounce orientation unmatch: {} fibre conjugacy",
                    fibre_transport.describe()
                ));
                Ok(TargetState {
                    orientation_mask: source_orientation,
                    column_word,
                })
            }
            Reduction::GeneralizedBounce(generalized) => {
                let key = path.word();
                if !self.generalized_cache.contains_key(&key) {
                    self.generalized_cache
                        .insert(key.clone(), GeneralizedBijection::build(path, generalized)?);
                }
                let (branch, rank, bucket_size, child) = {
                    let bijection = &self.generalized_cache[&key];
                    let (branch, rank, bucket_size) = bijection.source_to_branch(source)?;
                    let child = bijection.reduction.path_for(branch.branch)?.clone();
                    (branch, rank, bucket_size, child)
                };
                trace.push(format!(
                    "generalized bounce {}: canonical recursion-branch match {:?}, external q-offset {}, rank {}/{}",
                    key,
                    branch.branch,
                    generalized_offset(branch.branch, branch.factor),
                    rank + 1,
                    bucket_size
                ));
                let child_target = self.map_inner(&child, &branch.state, trace)?;
                let branch_target = BranchTarget {
                    branch: branch.branch,
                    factor: branch.factor,
                    state: child_target,
                };
                self.generalized_cache[&key].branch_to_target(&branch_target)
            }
            Reduction::BounceFallback => Err(format!(
                "unsupported bounce seed {}; no global source/target fallback is allowed",
                path.word()
            )),
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
            Reduction::FirstBounce {
                successor,
                base,
                low,
                mid,
                top,
            } => {
                let source_bit = mask_contains_edge(path, target.orientation_mask, (mid, top))?;
                let source_background = restrict_mask(path, target.orientation_mask, &base, 0)?;
                let (factor, successor_background) =
                    bounce_match(&base, low, mid, top, source_bit, source_background)?;
                let successor_orientation = successor_background;
                let fibre_transport = detect_fibre_transport(
                    path,
                    target.orientation_mask,
                    &successor,
                    successor_orientation,
                    low,
                    mid,
                )?;
                let successor_word = transport_columns(
                    path,
                    target.orientation_mask,
                    &target.column_word,
                    &successor,
                    successor_orientation,
                    fibre_transport,
                )?;
                let successor_target = TargetState {
                    orientation_mask: successor_orientation,
                    column_word: successor_word,
                };
                let successor_source = self.inverse_inner(&successor, &successor_target)?;
                first_bounce_coloring_inverse(
                    path,
                    &successor,
                    &base,
                    BounceVertices { low, mid, top },
                    factor,
                    &successor_source,
                )
            }
            Reduction::SecondBounce {
                successor,
                base,
                corner,
                low,
                mid,
                top,
            } => {
                let source_bit = mask_contains_edge(path, target.orientation_mask, (mid, top))?;
                let source_background = restrict_mask(path, target.orientation_mask, &base, 0)?;
                let (successor_bit, successor_background) =
                    second_match(&base, low, mid, top, source_bit, source_background)?;
                let successor_orientation = insert_orientation_bit(
                    &successor,
                    &base,
                    successor_background,
                    (corner, low),
                    successor_bit,
                )?;
                let fibre_transport = detect_fibre_transport(
                    path,
                    target.orientation_mask,
                    &successor,
                    successor_orientation,
                    low,
                    mid,
                )?;
                let successor_word = transport_columns(
                    path,
                    target.orientation_mask,
                    &target.column_word,
                    &successor,
                    successor_orientation,
                    fibre_transport,
                )?;
                let successor_target = TargetState {
                    orientation_mask: successor_orientation,
                    column_word: successor_word,
                };
                let successor_source = self.inverse_inner(&successor, &successor_target)?;
                second_bounce_coloring_inverse(
                    path,
                    &successor,
                    &base,
                    corner,
                    BounceVertices { low, mid, top },
                    &successor_source,
                )
            }
            Reduction::GeneralizedBounce(generalized) => {
                let key = path.word();
                if !self.generalized_cache.contains_key(&key) {
                    self.generalized_cache
                        .insert(key.clone(), GeneralizedBijection::build(path, generalized)?);
                }
                let (branch, child) = {
                    let bijection = &self.generalized_cache[&key];
                    let branch = bijection.target_to_branch(target)?;
                    let child = bijection.reduction.path_for(branch.branch)?.clone();
                    (branch, child)
                };
                let child_source = self.inverse_inner(&child, &branch.state)?;
                let branch_source = BranchSource {
                    branch: branch.branch,
                    factor: branch.factor,
                    state: child_source,
                };
                self.generalized_cache[&key].branch_to_source(&branch_source)
            }
            Reduction::BounceFallback => Err(format!(
                "unsupported bounce seed {}; no global source/target fallback is allowed",
                path.word()
            )),
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
    if first_east > 0 && path.steps[first_east - 1] == Step::Diagonal {
        let prefix = &path.steps[..first_east - 1];
        let suffix = &path.steps[first_east + 1..];
        let (low, top) = endpoint(prefix);
        let mid = low + 1;
        let mut candidates = Vec::new();
        for split in 0..prefix.len().saturating_sub(1) {
            let u = &prefix[..split];
            let v = &prefix[split + 2..];
            if endpoint(u).1 != low {
                continue;
            }
            match (prefix[split], prefix[split + 1]) {
                (Step::North, Step::North) => {
                    let mut successor_steps = prefix.to_vec();
                    successor_steps.push(Step::East);
                    successor_steps.push(Step::Diagonal);
                    successor_steps.extend_from_slice(suffix);
                    let mut base_steps = prefix.to_vec();
                    base_steps.extend([Step::East, Step::East, Step::North]);
                    base_steps.extend_from_slice(suffix);
                    candidates.push(Reduction::FirstBounce {
                        successor: PathData::from_steps(successor_steps)?,
                        base: PathData::from_steps(base_steps)?,
                        low,
                        mid,
                        top,
                    });
                }
                (Step::Diagonal, Step::North) => {
                    if endpoint(u).0 >= endpoint(u).1 {
                        continue;
                    }
                    let mut successor_steps = u.to_vec();
                    successor_steps.extend([Step::North, Step::Diagonal]);
                    successor_steps.extend_from_slice(v);
                    successor_steps.extend([Step::East, Step::Diagonal]);
                    successor_steps.extend_from_slice(suffix);
                    let mut base_steps = u.to_vec();
                    base_steps.extend([Step::East, Step::North, Step::North]);
                    base_steps.extend_from_slice(v);
                    base_steps.extend([Step::East, Step::East, Step::North]);
                    base_steps.extend_from_slice(suffix);
                    candidates.push(Reduction::SecondBounce {
                        successor: PathData::from_steps(successor_steps)?,
                        base: PathData::from_steps(base_steps)?,
                        corner: endpoint(u).0,
                        low,
                        mid,
                        top,
                    });
                }
                _ => {}
            }
        }
        if candidates.len() > 1 {
            return Err(format!(
                "ambiguous single-bounce decomposition for {}",
                path.word()
            ));
        }
        if let Some(candidate) = candidates.pop() {
            validate_bounce_reduction(path, &candidate)?;
            return Ok(candidate);
        }
        if let Some(generalized) = generalized_bounce_reduction(path, first_east)? {
            return Ok(Reduction::GeneralizedBounce(generalized));
        }
    }
    Ok(Reduction::BounceFallback)
}

fn generalized_bounce_reduction(
    path: &PathData,
    first_east: usize,
) -> Result<Option<GeneralizedBounce>, String> {
    let original = &path.steps[..first_east];
    if original.last() != Some(&Step::Diagonal) {
        return Ok(None);
    }
    let suffix = &path.steps[first_east + 1..];
    let mut current = original.to_vec();
    loop {
        let body = current
            .get(..current.len().saturating_sub(1))
            .ok_or_else(|| "bounce scanner reached an empty seed".to_string())?;
        let split = body.iter().filter(|&&step| step == Step::Diagonal).count();
        if split + 1 >= body.len() || split + 1 >= original.len() - 1 {
            return Err(format!("malformed bounce seed in {}", path.word()));
        }
        let first = body[split];
        let second = body[split + 1];
        if first == Step::Diagonal && second == Step::Diagonal {
            current = body[..split].to_vec();
            current.push(Step::Diagonal);
            continue;
        }
        if original[split] != first || original[split + 1] != second {
            return Err("bounce scanner lost its original-prefix coordinates".into());
        }
        let u = &original[..split];
        let v = &original[split + 2..original.len() - 1];
        let make_path = |middle: [Step; 2]| -> Result<PathData, String> {
            let mut steps = u.to_vec();
            steps.extend(middle);
            steps.extend_from_slice(v);
            steps.extend([Step::East, Step::Diagonal]);
            steps.extend_from_slice(suffix);
            PathData::from_steps(steps)
        };
        let reduction = match (first, second) {
            (Step::North, Step::North) => GeneralizedBounce::First {
                successor: make_path([Step::North, Step::North])?,
            },
            (Step::Diagonal, Step::North) => GeneralizedBounce::Second {
                successor: make_path([Step::North, Step::Diagonal])?,
            },
            (Step::North, Step::Diagonal) => GeneralizedBounce::Third {
                selected: make_path([Step::North, Step::Diagonal])?,
                optional: make_path([Step::Diagonal, Step::North])?,
            },
            _ => return Err("bounce scanner ended in a forbidden terminal pair".into()),
        };
        let parent_height = endpoint(original).1;
        for successor in match &reduction {
            GeneralizedBounce::First { successor } | GeneralizedBounce::Second { successor } => {
                vec![successor]
            }
            GeneralizedBounce::Third { selected, optional } => vec![selected, optional],
        } {
            let east = successor
                .steps
                .iter()
                .position(|step| *step == Step::East)
                .ok_or_else(|| "generalized successor has no east step".to_string())?;
            if endpoint(&successor.steps[..east]).1 >= parent_height {
                return Err("generalized bounce did not lower first-east height".into());
            }
        }
        return Ok(Some(reduction));
    }
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

impl FibreTransport {
    fn apply(self, vertex: usize) -> usize {
        match self {
            Self::Identity => vertex,
            Self::Swap(left, right) if vertex == left => right,
            Self::Swap(left, right) if vertex == right => left,
            Self::Swap(_, _) => vertex,
        }
    }

    fn describe(self) -> String {
        match self {
            Self::Identity => "identity".into(),
            Self::Swap(left, right) => format!("swap {}<->{}", left + 1, right + 1),
        }
    }
}

fn edge_set(edges: &[(usize, usize)]) -> BTreeSet<(usize, usize)> {
    edges.iter().copied().collect()
}

fn edge_set_with(
    edges: &[(usize, usize)],
    additions: &[(usize, usize)],
) -> BTreeSet<(usize, usize)> {
    edges.iter().chain(additions).copied().collect()
}

fn swap_vertex(vertex: usize, low: usize, mid: usize) -> usize {
    if vertex == low {
        mid
    } else if vertex == mid {
        low
    } else {
        vertex
    }
}

fn swap_edge(edge: (usize, usize), low: usize, mid: usize) -> (usize, usize) {
    (swap_vertex(edge.0, low, mid), swap_vertex(edge.1, low, mid))
}

fn validate_background_swap(base: &PathData, low: usize, mid: usize) -> Result<(), String> {
    let distinguished = (low, mid);
    if base.normal_edges.binary_search(&distinguished).is_err() {
        return Err("bounce background omits its distinguished normal edge".into());
    }
    let swapped_normal = base
        .normal_edges
        .iter()
        .copied()
        .filter(|&edge| edge != distinguished)
        .map(|edge| swap_edge(edge, low, mid))
        .collect::<BTreeSet<_>>();
    let expected_normal = base
        .normal_edges
        .iter()
        .copied()
        .filter(|&edge| edge != distinguished)
        .collect::<BTreeSet<_>>();
    let swapped_strict = base
        .strict_edges
        .iter()
        .copied()
        .map(|edge| swap_edge(edge, low, mid))
        .collect::<BTreeSet<_>>();
    if swapped_normal != expected_normal || swapped_strict != edge_set(&base.strict_edges) {
        return Err("bounce background is not invariant under the adjacent vertex swap".into());
    }
    Ok(())
}

fn validate_bounce_reduction(path: &PathData, reduction: &Reduction) -> Result<(), String> {
    match reduction {
        Reduction::FirstBounce {
            successor,
            base,
            low,
            mid,
            top,
        } => {
            if *mid != *low + 1
                || edge_set(&path.normal_edges)
                    != edge_set_with(&base.normal_edges, &[(*mid, *top)])
                || edge_set(&path.strict_edges)
                    != edge_set_with(&base.strict_edges, &[(*low, *top)])
                || edge_set(&successor.normal_edges) != edge_set(&base.normal_edges)
                || edge_set(&successor.strict_edges)
                    != edge_set_with(&base.strict_edges, &[(*mid, *top)])
            {
                return Err("first-bounce graph decomposition failed its edge certificate".into());
            }
            validate_background_swap(base, *low, *mid)
        }
        Reduction::SecondBounce {
            successor,
            base,
            corner,
            low,
            mid,
            top,
        } => {
            if *mid != *low + 1
                || *corner >= *low
                || edge_set(&path.normal_edges)
                    != edge_set_with(&base.normal_edges, &[(*mid, *top)])
                || edge_set(&path.strict_edges)
                    != edge_set_with(&base.strict_edges, &[(*corner, *low), (*low, *top)])
                || edge_set(&successor.normal_edges)
                    != edge_set_with(&base.normal_edges, &[(*corner, *low)])
                || edge_set(&successor.strict_edges)
                    != edge_set_with(&base.strict_edges, &[(*corner, *mid), (*mid, *top)])
            {
                return Err("second-bounce graph decomposition failed its edge certificate".into());
            }
            validate_background_swap(base, *low, *mid)
        }
        _ => Err("non-bounce reduction passed to bounce validator".into()),
    }
}

fn selected_edges(path: &PathData, mask: u64) -> Vec<(usize, usize)> {
    shifted_selected_edges(path, mask, 0)
}

fn swap_word_vertices(word: &[u8], low: usize, mid: usize) -> Vec<u8> {
    let mut swapped = word.to_vec();
    swapped.swap(low, mid);
    swapped
}

fn word_ascends(word: &[u8], edge: (usize, usize)) -> bool {
    word[edge.0] < word[edge.1]
}

fn first_bounce_coloring_forward(
    path: &PathData,
    successor: &PathData,
    base: &PathData,
    vertices: BounceVertices,
    source: &SourceState,
) -> Result<(bool, SourceState, &'static str), String> {
    let BounceVertices { low, mid, top } = vertices;
    let identity =
        source.coloring[low] < source.coloring[top] && source.coloring[mid] < source.coloring[top];
    let coloring = if identity {
        source.coloring.clone()
    } else {
        swap_word_vertices(&source.coloring, low, mid)
    };
    let selected = selected_edges(path, source.marked_mask);
    let mut factor = false;
    let mut successor_selected = Vec::new();
    if identity {
        for edge in selected {
            if edge == (mid, top) {
                factor = true;
            } else {
                successor_selected.push(edge);
            }
        }
    } else {
        for edge in selected {
            if edge == (low, mid) {
                factor = true;
            } else if edge == (mid, top) {
                return Err("nonascending first-bounce fresh edge was marked".into());
            } else {
                if base.normal_edges.binary_search(&edge).is_err() {
                    return Err("first-bounce mark lies outside the background".into());
                }
                successor_selected.push(swap_edge(edge, low, mid));
            }
        }
    }
    let result = SourceState {
        coloring,
        marked_mask: mask_for_edges(successor, &successor_selected)?,
    };
    validate_source(successor, &result)?;
    if source.marked_mask.count_ones() != result.marked_mask.count_ones() + u32::from(factor) {
        return Err("first-bounce coloring lift changed marked degree".into());
    }
    Ok((factor, result, if identity { "identity" } else { "swap" }))
}

fn first_bounce_coloring_inverse(
    path: &PathData,
    successor: &PathData,
    base: &PathData,
    vertices: BounceVertices,
    factor: bool,
    source: &SourceState,
) -> Result<SourceState, String> {
    let BounceVertices { low, mid, top } = vertices;
    let identity =
        source.coloring[low] < source.coloring[top] && source.coloring[mid] < source.coloring[top];
    let coloring = if identity {
        source.coloring.clone()
    } else {
        swap_word_vertices(&source.coloring, low, mid)
    };
    let mut selected = Vec::new();
    if identity {
        selected.extend(selected_edges(successor, source.marked_mask));
        if factor {
            selected.push((mid, top));
        }
    } else {
        for edge in selected_edges(successor, source.marked_mask) {
            if edge == (low, mid) {
                return Err("nonascending successor distinguished edge was marked".into());
            }
            if base.normal_edges.binary_search(&edge).is_err() {
                return Err("first-bounce inverse mark lies outside the background".into());
            }
            selected.push(swap_edge(edge, low, mid));
        }
        if factor {
            selected.push((low, mid));
        }
    }
    let result = SourceState {
        coloring,
        marked_mask: mask_for_edges(path, &selected)?,
    };
    validate_source(path, &result)?;
    if result.marked_mask.count_ones() != source.marked_mask.count_ones() + u32::from(factor) {
        return Err("first-bounce coloring inverse changed marked degree".into());
    }
    Ok(result)
}

fn map_special_marks(
    from_word: &[u8],
    from_selected: &BTreeSet<(usize, usize)>,
    from_special: &[(usize, usize)],
    to_word: &[u8],
    to_special: &[(usize, usize)],
) -> Result<Vec<(usize, usize)>, String> {
    let from_ascents = from_special
        .iter()
        .copied()
        .filter(|&edge| word_ascends(from_word, edge))
        .collect::<Vec<_>>();
    let to_ascents = to_special
        .iter()
        .copied()
        .filter(|&edge| word_ascends(to_word, edge))
        .collect::<Vec<_>>();
    if from_ascents.len() != to_ascents.len() {
        return Err("second-bounce special ascent counts do not match".into());
    }
    Ok(from_ascents
        .iter()
        .zip(to_ascents)
        .filter_map(|(from, to)| from_selected.contains(from).then_some(to))
        .collect())
}

fn second_bounce_coloring_forward(
    path: &PathData,
    successor: &PathData,
    base: &PathData,
    corner: usize,
    vertices: BounceVertices,
    source: &SourceState,
) -> Result<(SourceState, &'static str), String> {
    let BounceVertices { low, mid, top } = vertices;
    let inside_low = source.coloring[corner] < source.coloring[low]
        && source.coloring[low] < source.coloring[top];
    let inside_mid = source.coloring[corner] < source.coloring[mid]
        && source.coloring[mid] < source.coloring[top];
    let identity = inside_low && inside_mid;
    let coloring = if identity {
        source.coloring.clone()
    } else {
        swap_word_vertices(&source.coloring, low, mid)
    };
    let from_selected = selected_edges(path, source.marked_mask)
        .into_iter()
        .collect::<BTreeSet<_>>();
    let mut selected = Vec::new();
    if identity {
        for edge in &from_selected {
            selected.push(if *edge == (mid, top) {
                (corner, low)
            } else {
                *edge
            });
        }
    } else {
        for &edge in &from_selected {
            if edge == (low, mid) || edge == (mid, top) {
                continue;
            }
            if base.normal_edges.binary_search(&edge).is_err() {
                return Err("second-bounce mark lies outside the background".into());
            }
            selected.push(swap_edge(edge, low, mid));
        }
        selected.extend(map_special_marks(
            &source.coloring,
            &from_selected,
            &[(low, mid), (mid, top)],
            &coloring,
            &[(low, mid), (corner, low)],
        )?);
    }
    let result = SourceState {
        coloring,
        marked_mask: mask_for_edges(successor, &selected)?,
    };
    validate_source(successor, &result)?;
    if source_signature(source) != source_signature(&result) {
        return Err("second-bounce coloring lift changed exact weight".into());
    }
    Ok((result, if identity { "identity" } else { "swap" }))
}

fn second_bounce_coloring_inverse(
    path: &PathData,
    successor: &PathData,
    base: &PathData,
    corner: usize,
    vertices: BounceVertices,
    source: &SourceState,
) -> Result<SourceState, String> {
    let BounceVertices { low, mid, top } = vertices;
    let inside_low = source.coloring[corner] < source.coloring[low]
        && source.coloring[low] < source.coloring[top];
    let inside_mid = source.coloring[corner] < source.coloring[mid]
        && source.coloring[mid] < source.coloring[top];
    let identity = inside_low && inside_mid;
    let coloring = if identity {
        source.coloring.clone()
    } else {
        swap_word_vertices(&source.coloring, low, mid)
    };
    let from_selected = selected_edges(successor, source.marked_mask)
        .into_iter()
        .collect::<BTreeSet<_>>();
    let mut selected = Vec::new();
    if identity {
        for edge in &from_selected {
            selected.push(if *edge == (corner, low) {
                (mid, top)
            } else {
                *edge
            });
        }
    } else {
        for &edge in &from_selected {
            if edge == (low, mid) || edge == (corner, low) {
                continue;
            }
            if base.normal_edges.binary_search(&edge).is_err() {
                return Err("second-bounce inverse mark lies outside the background".into());
            }
            selected.push(swap_edge(edge, low, mid));
        }
        selected.extend(map_special_marks(
            &source.coloring,
            &from_selected,
            &[(low, mid), (corner, low)],
            &coloring,
            &[(low, mid), (mid, top)],
        )?);
    }
    let result = SourceState {
        coloring,
        marked_mask: mask_for_edges(path, &selected)?,
    };
    validate_source(path, &result)?;
    if source_signature(source) != source_signature(&result) {
        return Err("second-bounce coloring inverse changed exact weight".into());
    }
    Ok(result)
}

fn background_values(
    base: &PathData,
    mask: u64,
    low: usize,
    mid: usize,
    top: usize,
) -> Result<(usize, usize, usize), String> {
    let mut upward = base.strict_edges.clone();
    upward.extend(
        base.normal_edges
            .iter()
            .enumerate()
            .filter(|(index, edge)| **edge != (low, mid) && mask & (1u64 << index) != 0)
            .map(|(_, &edge)| edge),
    );
    let hrv = highest_reachable_vertices(base.n, &upward)
        .ok_or_else(|| "bounce background upward graph is invalid".to_string())?;
    Ok((hrv[low], hrv[mid], hrv[top]))
}

fn set_background_bit(
    base: &PathData,
    mask: u64,
    edge: (usize, usize),
    value: bool,
) -> Result<u64, String> {
    let index = base
        .normal_edges
        .binary_search(&edge)
        .map_err(|_| format!("background does not contain distinguished edge {edge:?}"))?;
    Ok(if value {
        mask | (1u64 << index)
    } else {
        mask & !(1u64 << index)
    })
}

fn swap_background_orientation(
    base: &PathData,
    mask: u64,
    low: usize,
    mid: usize,
) -> Result<u64, String> {
    let distinguished = (low, mid);
    let selected = selected_edges(base, mask)
        .into_iter()
        .map(|edge| {
            if edge == distinguished {
                edge
            } else {
                swap_edge(edge, low, mid)
            }
        })
        .collect::<Vec<_>>();
    mask_for_edges(base, &selected)
}

fn bounce_match(
    base: &PathData,
    low: usize,
    mid: usize,
    top: usize,
    factor: bool,
    orientation: u64,
) -> Result<(bool, u64), String> {
    let beta = mask_contains_edge(base, orientation, (low, mid))?;
    let (a, c, h) = background_values(base, orientation, low, mid, top)?;
    if h <= c {
        if beta || h <= a {
            Ok((factor, orientation))
        } else {
            Ok((
                factor,
                set_background_bit(
                    base,
                    swap_background_orientation(base, orientation, low, mid)?,
                    (low, mid),
                    false,
                )?,
            ))
        }
    } else if factor {
        Ok((
            beta,
            set_background_bit(base, orientation, (low, mid), true)?,
        ))
    } else {
        Ok((
            beta,
            set_background_bit(
                base,
                swap_background_orientation(base, orientation, low, mid)?,
                (low, mid),
                false,
            )?,
        ))
    }
}

fn bounce_unmatch(
    base: &PathData,
    low: usize,
    mid: usize,
    top: usize,
    factor: bool,
    orientation: u64,
) -> Result<(bool, u64), String> {
    let beta = mask_contains_edge(base, orientation, (low, mid))?;
    let (a, c, h) = background_values(base, orientation, low, mid, top)?;
    if beta {
        if h <= c {
            Ok((factor, orientation))
        } else {
            Ok((
                true,
                set_background_bit(base, orientation, (low, mid), factor)?,
            ))
        }
    } else if h <= a {
        if h <= c {
            Ok((factor, orientation))
        } else {
            Ok((
                factor,
                swap_background_orientation(base, orientation, low, mid)?,
            ))
        }
    } else {
        Ok((
            false,
            set_background_bit(
                base,
                swap_background_orientation(base, orientation, low, mid)?,
                (low, mid),
                factor,
            )?,
        ))
    }
}

fn second_match(
    base: &PathData,
    low: usize,
    mid: usize,
    top: usize,
    source_bit: bool,
    orientation: u64,
) -> Result<(bool, u64), String> {
    let beta = mask_contains_edge(base, orientation, (low, mid))?;
    let (a, c, h) = background_values(base, orientation, low, mid, top)?;
    if source_bit == beta {
        if source_bit {
            Ok((true, orientation))
        } else {
            Ok((
                false,
                set_background_bit(
                    base,
                    swap_background_orientation(base, orientation, low, mid)?,
                    (low, mid),
                    false,
                )?,
            ))
        }
    } else if a <= h && c <= h {
        if source_bit {
            Ok((
                false,
                set_background_bit(base, orientation, (low, mid), true)?,
            ))
        } else {
            Ok((
                true,
                set_background_bit(
                    base,
                    swap_background_orientation(base, orientation, low, mid)?,
                    (low, mid),
                    false,
                )?,
            ))
        }
    } else if source_bit {
        if c <= a {
            Ok((true, orientation))
        } else {
            Ok((
                false,
                set_background_bit(
                    base,
                    swap_background_orientation(base, orientation, low, mid)?,
                    (low, mid),
                    true,
                )?,
            ))
        }
    } else if a <= c {
        Ok((false, orientation))
    } else {
        Ok((
            true,
            set_background_bit(
                base,
                swap_background_orientation(base, orientation, low, mid)?,
                (low, mid),
                false,
            )?,
        ))
    }
}

fn second_unmatch(
    base: &PathData,
    low: usize,
    mid: usize,
    top: usize,
    successor_bit: bool,
    orientation: u64,
) -> Result<(bool, u64), String> {
    let beta = mask_contains_edge(base, orientation, (low, mid))?;
    let (a, c, h) = background_values(base, orientation, low, mid, top)?;
    if successor_bit == beta {
        if successor_bit {
            Ok((true, orientation))
        } else {
            Ok((
                false,
                set_background_bit(
                    base,
                    swap_background_orientation(base, orientation, low, mid)?,
                    (low, mid),
                    false,
                )?,
            ))
        }
    } else if a <= h && c <= h {
        if successor_bit {
            Ok((
                false,
                set_background_bit(
                    base,
                    swap_background_orientation(base, orientation, low, mid)?,
                    (low, mid),
                    true,
                )?,
            ))
        } else {
            Ok((
                true,
                set_background_bit(base, orientation, (low, mid), false)?,
            ))
        }
    } else if successor_bit {
        if c <= a {
            Ok((true, orientation))
        } else {
            Ok((
                false,
                set_background_bit(
                    base,
                    swap_background_orientation(base, orientation, low, mid)?,
                    (low, mid),
                    true,
                )?,
            ))
        }
    } else if a <= c {
        Ok((false, orientation))
    } else {
        Ok((
            true,
            set_background_bit(
                base,
                swap_background_orientation(base, orientation, low, mid)?,
                (low, mid),
                false,
            )?,
        ))
    }
}

fn insert_orientation_bit(
    path: &PathData,
    base: &PathData,
    base_mask: u64,
    edge: (usize, usize),
    bit: bool,
) -> Result<u64, String> {
    let mut selected = selected_edges(base, base_mask);
    if bit {
        selected.push(edge);
    }
    mask_for_edges(path, &selected)
}

fn detect_fibre_transport(
    source_path: &PathData,
    source_mask: u64,
    target_path: &PathData,
    target_mask: u64,
    low: usize,
    mid: usize,
) -> Result<FibreTransport, String> {
    let source_hrv = hrv_for_mask(source_path, source_mask)?;
    let target_hrv = hrv_for_mask(target_path, target_mask)?;
    if source_hrv == target_hrv {
        return Ok(FibreTransport::Identity);
    }
    let transport = FibreTransport::Swap(low, mid);
    if (0..source_path.n)
        .all(|vertex| target_hrv[transport.apply(vertex)] == transport.apply(source_hrv[vertex]))
    {
        Ok(transport)
    } else {
        Err("matched bounce orientations have no proved identity/swap HRV conjugacy".into())
    }
}

fn transport_columns(
    source_path: &PathData,
    source_mask: u64,
    source_word: &[u8],
    target_path: &PathData,
    target_mask: u64,
    transport: FibreTransport,
) -> Result<Vec<u8>, String> {
    let source = TargetState {
        orientation_mask: source_mask,
        column_word: source_word.to_vec(),
    };
    validate_target(source_path, &source)?;
    if source_path.n != target_path.n {
        return Err("cannot transport columns between different ranks".into());
    }
    let source_hrv = hrv_for_mask(source_path, source_mask)?;
    let target_hrv = hrv_for_mask(target_path, target_mask)?;
    if !(0..source_path.n)
        .all(|vertex| target_hrv[transport.apply(vertex)] == transport.apply(source_hrv[vertex]))
    {
        return Err("column transport does not conjugate the two HRV maps".into());
    }
    let mut target_word = vec![0u8; source_path.n];
    let mut assigned = vec![false; source_path.n];
    for fibre in fibres(&source_hrv) {
        let mut colors = fibre
            .iter()
            .map(|&vertex| source_word[vertex])
            .collect::<Vec<_>>();
        colors.sort_unstable();
        let mut target_vertices = fibre
            .iter()
            .map(|&vertex| transport.apply(vertex))
            .collect::<Vec<_>>();
        target_vertices.sort_unstable();
        for (vertex, color) in target_vertices.into_iter().zip(colors) {
            if std::mem::replace(&mut assigned[vertex], true) {
                return Err("column transport assigned a vertex twice".into());
            }
            target_word[vertex] = color;
        }
    }
    if assigned.iter().any(|assigned| !assigned) {
        return Err("column transport omitted a vertex".into());
    }
    let target = TargetState {
        orientation_mask: target_mask,
        column_word: target_word,
    };
    validate_target(target_path, &target)?;
    if color_content(source_word) != color_content(&target.column_word) {
        return Err("column transport changed color content".into());
    }
    Ok(target.column_word)
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

fn checked_path(path: &PathData) -> Result<(usize, usize, usize, usize, usize, usize), String> {
    let bijection = CanonicalBijection::build(path)?;
    let all_targets = enumerate_targets(path)?
        .into_iter()
        .collect::<BTreeSet<_>>();
    let mut recursive_targets = BTreeSet::new();
    let mut recursive_matcher = RecursiveMatcher::default();
    let mut certificates = 0usize;
    let mut fallback_certificates = 0usize;
    let mut first_bounce_certificates = 0usize;
    let mut second_bounce_certificates = 0usize;
    let mut generalized_bounce_certificates = 0usize;
    for sources in bijection.source_buckets.values() {
        for source in sources {
            let certificate = bijection.map(path, source)?;
            if source_signature(&certificate.source) != target_signature(&certificate.target) {
                return Err("certificate failed weight preservation".into());
            }
            bijection.replay(path, &certificate)?;
            let recursive_certificate = recursive_matcher.map(path, source)?;
            recursive_matcher.replay(path, &recursive_certificate)?;
            fallback_certificates += usize::from(
                recursive_certificate
                    .trace
                    .iter()
                    .any(|step| step.contains("canonical exact-bucket fallback")),
            );
            first_bounce_certificates += usize::from(
                recursive_certificate
                    .trace
                    .iter()
                    .any(|step| step.starts_with("first bounce")),
            );
            second_bounce_certificates += usize::from(
                recursive_certificate
                    .trace
                    .iter()
                    .any(|step| step.starts_with("second bounce")),
            );
            generalized_bounce_certificates += usize::from(
                recursive_certificate
                    .trace
                    .iter()
                    .any(|step| step.starts_with("generalized bounce")),
            );
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
    Ok((
        bijection.source_buckets.len(),
        certificates,
        fallback_certificates,
        first_bounce_certificates,
        second_bounce_certificates,
        generalized_bounce_certificates,
    ))
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
    println!("reverse certificate: replayed inverse recursion -> original source [OK]");
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
            let (buckets, certificates, fallback, first_bounce, second_bounce, generalized) =
                checked_path(&path)?;
            println!(
                "checked path {}: rank {}, {} signatures, {} reversible certificates; residual-fallback {}, first-bounce {}, second-bounce {}, generalized-branch {}",
                path.word(), path.n, buckets, certificates, fallback, first_bounce, second_bounce,
                generalized
            );
        }
        return Ok(());
    }

    let mut total_paths = 0usize;
    let mut total_buckets = 0usize;
    let mut total_certificates = 0usize;
    let mut total_fallback = 0usize;
    let mut total_first_bounce = 0usize;
    let mut total_second_bounce = 0usize;
    let mut total_generalized_bounce = 0usize;
    for n in 0..=args.max_n {
        let paths = all_paths(n);
        let mut rank_certificates = 0usize;
        for path in &paths {
            let (buckets, certificates, fallback, first_bounce, second_bounce, generalized) =
                checked_path(path)?;
            total_buckets += buckets;
            total_certificates += certificates;
            rank_certificates += certificates;
            total_fallback += fallback;
            total_first_bounce += first_bounce;
            total_second_bounce += second_bounce;
            total_generalized_bounce += generalized;
        }
        total_paths += paths.len();
        println!(
            "rank {n}: {} paths, {} reversible certificates",
            paths.len(),
            rank_certificates
        );
    }
    println!(
        "total: {total_paths} ordinary decorated paths, {total_buckets} exact invariant buckets, {total_certificates} reversible certificates; residual-fallback {total_fallback}, first-bounce {total_first_bounce}, second-bounce {total_second_bounce}, generalized-branch {total_generalized_bounce}"
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
        let (_, _, fallback, first_bounce, _, _) = checked_path(&path).unwrap();
        assert_eq!(fallback, 0);
        assert!(first_bounce > 0);
    }

    #[test]
    fn single_second_bounce_is_explicit_and_reversible() {
        let path = PathData::parse("NDNDEE").unwrap();
        let (_, _, fallback, _, second_bounce, _) = checked_path(&path).unwrap();
        assert_eq!(fallback, 0);
        assert!(second_bounce > 0);
    }

    #[test]
    fn third_bounce_uses_reversible_recursion_branch_matching() {
        let path = PathData::parse("NNDDEE").unwrap();
        let (_, _, fallback, _, _, generalized) = checked_path(&path).unwrap();
        assert_eq!(fallback, 0);
        assert!(generalized > 0);
    }

    #[test]
    fn k3_requires_mixing_orientations() {
        let path = PathData::parse("NNNEEE").unwrap();
        assert!(same_orientation_obstruction(&path).unwrap().is_some());
    }

    #[test]
    fn every_ordinary_decorated_path_through_rank_four_checks() {
        let mut generalized = 0usize;
        for n in 0..=4 {
            for path in all_paths(n) {
                let (_, _, fallback, _, _, path_generalized) =
                    checked_path(&path).unwrap_or_else(|error| panic!("{}: {error}", path.word()));
                assert_eq!(fallback, 0, "{} used the forbidden fallback", path.word());
                generalized += path_generalized;
            }
        }
        assert!(generalized > 0);
    }
}
