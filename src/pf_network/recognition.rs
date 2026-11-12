//! Bounded support search and PF-network recognition.

use super::*;

/// All primitive step types within the bounds, in increasing order.
pub fn primitive_step_types(max_step_length: usize, max_rise: usize) -> Vec<PfStepType> {
    let mut types = Vec::new();
    for length in 1..=max_step_length {
        for rise in 0..=max_rise {
            let step = PfStepType { length, rise };
            if is_primitive(step) {
                types.push(step);
            }
        }
    }
    types
}

/// Search for a height-invariant planar network whose rows are the given
/// polynomials.  See the module documentation for the setting, the method and
/// what a certificate proves.
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::pf_network::{recognize_pf_network, PfNetworkOptions, PfNetworkRecognition};
/// use polytool::Polynomial;
///
/// let rows: Vec<Polynomial<BigInt>> = [vec![1, 1], vec![1, 2, 1], vec![1, 3, 3, 1]]
///     .iter()
///     .map(|r| Polynomial::from_i64_coeffs(r))
///     .collect();
/// let result = recognize_pf_network(&rows, &PfNetworkOptions::default()).unwrap();
/// assert!(result.is_certified());
///
/// // 1 + t + t^2 is not real-rooted, so no network exists at any bound.
/// let bad = vec![Polynomial::from_i64_coeffs(&[1, 1, 1])];
/// let result = recognize_pf_network(&bad, &PfNetworkOptions::default()).unwrap();
/// assert!(matches!(result, PfNetworkRecognition::ImpossibleWithinBounds { .. }));
/// ```
pub fn recognize_pf_network(
    rows: &[Polynomial<BigInt>],
    options: &PfNetworkOptions,
) -> Result<PfNetworkRecognition, PfNetworkRecognitionError> {
    if rows.is_empty() {
        return Err(PfNetworkRecognitionError::EmptyInput);
    }
    if options.first_row_index > 1 {
        return Err(PfNetworkRecognitionError::UnsupportedFirstRowIndex {
            found: options.first_row_index,
        });
    }
    if options.max_step_length == 0 {
        return Err(PfNetworkRecognitionError::ZeroStepLength);
    }
    let mut vectors: Vec<Vec<BigInt>> = Vec::with_capacity(rows.len() + 1);
    if options.first_row_index == 1 {
        vectors.push(vec![BigInt::one()]);
    } else if rows[0].coeffs() != [BigInt::one()] {
        return Err(PfNetworkRecognitionError::FirstRowNotOne {
            found: rows[0].coeffs().to_vec(),
        });
    }
    vectors.extend(rows.iter().map(|p| p.coeffs().to_vec()));

    for (row, v) in vectors.iter().enumerate() {
        if let Some((k, c)) = v.iter().enumerate().find(|(_, c)| **c < BigInt::zero()) {
            return Ok(PfNetworkRecognition::ImpossibleWithinBounds {
                supports_checked: 0,
                obstruction: format!(
                    "row {} has negative coefficient {} at t^{}; path weights are nonnegative",
                    row, c, k
                ),
            });
        }
    }
    for (row, v) in vectors.iter().enumerate() {
        if !v.is_empty() && !is_real_rooted_bigint_coeffs(v) {
            return Ok(PfNetworkRecognition::ImpossibleWithinBounds {
                supports_checked: 0,
                obstruction: format!(
                    "row {} is not real-rooted, but every row of a height-invariant planar \
                     network is a PF sequence (no network exists at any bound)",
                    row
                ),
            });
        }
    }

    let row_indices: Vec<usize> = (0..rows.len())
        .map(|i| i + options.first_row_index)
        .collect();
    let mut search = SupportSearch::new(&vectors, options);
    match search.run() {
        SearchOutcome::Found(columns) => {
            let network = HeightInvariantNetwork {
                steps: columns
                    .into_iter()
                    .map(|column| {
                        column
                            .into_iter()
                            .filter(|(_, w)| !w.is_zero())
                            .map(|(step, weight)| PfStep {
                                length: step.length,
                                rise: step.rise,
                                weight,
                            })
                            .collect()
                    })
                    .collect(),
            };
            let certificate = build_pf_rows_certificate(network, &row_indices, rows)
                .map_err(PfNetworkRecognitionError::CertificateRejected)?;
            Ok(PfNetworkRecognition::Certified { certificate })
        }
        SearchOutcome::Impossible(obstruction) => {
            Ok(PfNetworkRecognition::ImpossibleWithinBounds {
                supports_checked: search.nodes,
                obstruction,
            })
        }
        SearchOutcome::Inconclusive(reason) => Ok(PfNetworkRecognition::Inconclusive { reason }),
    }
}

type ColumnWeights = Vec<(PfStepType, Q)>;

pub(super) enum SearchOutcome {
    Found(Vec<ColumnWeights>),
    Impossible(String),
    Inconclusive(String),
}

/// A minimal feasible support of one row problem, with its weights.
#[derive(Clone)]
struct Support {
    types: Vec<usize>,
    weights: Vec<Q>,
}

pub(super) struct SupportSearch<'a> {
    vectors: &'a [Vec<BigInt>],
    types: Vec<PfStepType>,
    max_length: usize,
    integer: bool,
    max_nodes: usize,
    max_integer_nodes: usize,
    lp_cache: HashMap<(usize, Vec<usize>), ColumnOutcome>,
    minimal_cache: HashMap<(usize, Vec<usize>), Vec<Support>>,
    failed: std::collections::HashSet<(usize, Vec<Vec<usize>>)>,
    nodes: usize,
    budget_exhausted: bool,
    inconclusive_rows: usize,
}

impl<'a> SupportSearch<'a> {
    pub(super) fn new(vectors: &'a [Vec<BigInt>], options: &PfNetworkOptions) -> Self {
        Self {
            vectors,
            types: primitive_step_types(options.max_step_length, options.max_rise),
            max_length: options.max_step_length,
            integer: options.weights == BrentiWeightDomain::Integer,
            max_nodes: options.max_search_nodes,
            max_integer_nodes: options.max_integer_search_nodes,
            lp_cache: HashMap::new(),
            minimal_cache: HashMap::new(),
            failed: std::collections::HashSet::new(),
            nodes: 0,
            budget_exhausted: false,
            inconclusive_rows: 0,
        }
    }

    fn last_row(&self) -> usize {
        self.vectors.len() - 1
    }

    pub(super) fn run(&mut self) -> SearchOutcome {
        for n in 1..=self.last_row() {
            let usable: Vec<usize> = (0..self.types.len())
                .filter(|&i| self.types[i].length <= n)
                .collect();
            if self.row(n, &usable) == ColumnOutcome::Infeasible {
                return SearchOutcome::Impossible(format!(
                    "row {} is not a nonnegative combination of t^rise * P_(n - length) over \
                     the step types within the bounds, even ignoring planarity",
                    n
                ));
            }
        }
        let mut chosen: Vec<Support> = Vec::with_capacity(self.last_row());
        if self.extend(&mut chosen) {
            let columns = chosen
                .iter()
                .map(|support| {
                    support
                        .types
                        .iter()
                        .map(|&i| self.types[i])
                        .zip(support.weights.iter().cloned())
                        .collect()
                })
                .collect();
            return SearchOutcome::Found(columns);
        }
        if self.budget_exhausted {
            SearchOutcome::Inconclusive(format!(
                "support search budget of {} nodes exhausted",
                self.max_nodes
            ))
        } else if self.inconclusive_rows > 0 {
            SearchOutcome::Inconclusive(format!(
                "integer weight search budget of {} nodes exhausted for some row",
                self.max_integer_nodes
            ))
        } else {
            let domain = if self.integer { "integer" } else { "rational" };
            SearchOutcome::Impossible(format!(
                "no planar choice of steps with length <= {} and rise <= {} admits \
                 nonnegative {} weights",
                self.max_length,
                self.types.iter().map(|t| t.rise).max().unwrap_or(0),
                domain
            ))
        }
    }

    /// Depth-first search over the columns; `chosen[m - 1]` is the support
    /// into column `m`.  Only minimal feasible supports are tried: replacing a
    /// support by a feasible subset keeps its row and removes obstacles for
    /// later columns, so this loses no solutions.
    fn extend(&mut self, chosen: &mut Vec<Support>) -> bool {
        let n = chosen.len() + 1;
        if n > self.last_row() {
            return true;
        }
        // Only edges into columns m > n - max_length can meet later edges.
        let first_relevant = (n + 1).saturating_sub(self.max_length).max(1);
        let key = (
            n,
            chosen[first_relevant - 1..]
                .iter()
                .map(|s| s.types.clone())
                .collect::<Vec<_>>(),
        );
        if self.failed.contains(&key) {
            return false;
        }
        let candidates: Vec<usize> = (0..self.types.len())
            .filter(|&i| {
                let step = self.types[i];
                step.length <= n
                    && chosen[first_relevant - 1..]
                        .iter()
                        .enumerate()
                        .all(|(offset, support)| {
                            let m = first_relevant + offset;
                            support.types.iter().all(|&j| {
                                let other = self.types[j];
                                find_conflict(
                                    (m - other.length) as i64,
                                    other,
                                    (n - step.length) as i64,
                                    step,
                                    false,
                                )
                                .is_none()
                            })
                        })
            })
            .collect();
        let inconclusive_before = self.inconclusive_rows;
        let supports = self.minimal_supports(n, candidates);
        for support in supports {
            if self.nodes >= self.max_nodes {
                self.budget_exhausted = true;
                return false;
            }
            self.nodes += 1;
            chosen.push(support);
            if self.extend(chosen) {
                return true;
            }
            chosen.pop();
            if self.budget_exhausted {
                return false;
            }
        }
        if self.inconclusive_rows == inconclusive_before {
            self.failed.insert(key);
        }
        false
    }

    /// The inclusion-minimal subsets of `candidates` that are planar within
    /// column `n` and admit nonnegative weights for row `n`.
    fn minimal_supports(&mut self, n: usize, candidates: Vec<usize>) -> Vec<Support> {
        let key = (n, candidates);
        if let Some(found) = self.minimal_cache.get(&key) {
            return found.clone();
        }
        let candidates = &key.1;
        let mut cliques: Vec<Vec<usize>> = vec![Vec::new()];
        for &i in candidates {
            let extended: Vec<Vec<usize>> = cliques
                .iter()
                .filter(|clique| clique.iter().all(|&j| self.same_column_compatible(n, i, j)))
                .map(|clique| {
                    let mut next = clique.clone();
                    next.push(i);
                    next
                })
                .collect();
            cliques.extend(extended);
            if cliques.len() > self.max_nodes {
                self.budget_exhausted = true;
                break;
            }
        }
        cliques.sort_by_key(Vec::len);
        let mut minimal: Vec<Support> = Vec::new();
        for clique in cliques {
            if minimal
                .iter()
                .any(|m| m.types.iter().all(|t| clique.contains(t)))
            {
                continue;
            }
            match self.row(n, &clique) {
                ColumnOutcome::Feasible(weights) => minimal.push(Support {
                    types: clique,
                    weights,
                }),
                ColumnOutcome::Infeasible => {}
                ColumnOutcome::Inconclusive => self.inconclusive_rows += 1,
            }
        }
        self.minimal_cache.insert(key, minimal.clone());
        minimal
    }

    fn same_column_compatible(&self, n: usize, i: usize, j: usize) -> bool {
        let (a, b) = (self.types[i], self.types[j]);
        find_conflict((n - a.length) as i64, a, (n - b.length) as i64, b, i == j).is_none()
    }

    /// Decide `P_n = sum_i w_i t^{rise_i} P_{n - length_i}`, `w >= 0`, over
    /// the step types `usable` (all with `length <= n`).
    fn row(&mut self, n: usize, usable: &[usize]) -> ColumnOutcome {
        let key = (n, usable.to_vec());
        if let Some(outcome) = self.lp_cache.get(&key) {
            return outcome.clone();
        }
        let target = &self.vectors[n];
        let mut height = target.len();
        for &i in usable {
            let step = self.types[i];
            height = height.max(self.vectors[n - step.length].len() + step.rise);
        }
        let matrix: Vec<Vec<Q>> = (0..height)
            .map(|k| {
                usable
                    .iter()
                    .map(|&i| {
                        let step = self.types[i];
                        k.checked_sub(step.rise)
                            .and_then(|j| self.vectors[n - step.length].get(j))
                            .map_or_else(Q::zero, |c| Q::from_integer(c.clone()))
                    })
                    .collect()
            })
            .collect();
        let rhs: Vec<Q> = (0..height)
            .map(|k| {
                target
                    .get(k)
                    .map_or_else(Q::zero, |c| Q::from_integer(c.clone()))
            })
            .collect();
        let outcome = solve_nonnegative_column(
            &matrix,
            &rhs,
            usable.len(),
            self.integer,
            self.max_integer_nodes,
        );
        self.lp_cache.insert(key, outcome.clone());
        outcome
    }
}
