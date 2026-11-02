//! Lazy enumeration of candidate recurrence shapes in increasing complexity.

use super::*;
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::collections::HashSet;

type CandidateSortKey = (
    usize,
    usize,
    bool,
    bool,
    bool,
    usize,
    usize,
    usize,
    usize,
    usize,
    usize,
    (LeaderChart, usize, u128),
);

/// Extra weight of a left-hand chart per degree of its fixed monomial
/// (`candidate_complexity` documents the weights).
const LEADER_CHART_WEIGHT: usize = 1;

#[derive(Clone, Copy, Debug, Hash, Eq, PartialEq, Ord, PartialOrd)]
enum CandidateKind {
    Ordinary,
    Inhomogeneous,
    DenominatorIndexOnly,
    DenominatorPositiveVariableDegree,
}

#[derive(Clone, Copy, Debug, Hash, Eq, PartialEq, Ord, PartialOrd)]
struct CandidateNode {
    kind: CandidateKind,
    rec_len: usize,
    diff_deg: usize,
    idx_deg: usize,
    var_deg: usize,
    alternating_sign: bool,
    extra_idx_deg: usize,
    extra_var_deg: usize,
    /// The normalization chart of the left-hand factor (always the constant
    /// chart outside the denominator kinds, and without charts).
    chart: LeaderChart,
    /// The power `e` of the left-hand monomial `x^e` (0 without
    /// `try_laurent`).
    laurent: usize,
}

pub(super) struct CandidateIterator<'a> {
    search: &'a AdaptiveSearchOptions,
    /// The number of rows the enumeration was built for.
    rows: usize,
    /// Also enumerate the left-hand variants: monomials `x^e`
    /// (`e <= max_laurent_x_power`) and the non-constant charts.
    with_charts: bool,
    min_rec_len: usize,
    max_rec_len: usize,
    heap: BinaryHeap<Reverse<(CandidateSortKey, CandidateNode)>>,
    visited: HashSet<CandidateNode>,
    total_candidates: usize,
}

fn inclusive_count(minimum: usize, maximum: usize) -> usize {
    maximum
        .checked_sub(minimum)
        .and_then(|difference| difference.checked_add(1))
        .unwrap_or(0)
}

impl CandidateNode {
    fn to_options(self, modular_prefilter: bool) -> RecurrenceOptions {
        let (homogeneous, inhomo_idx_deg, inhomo_var_deg, denom_idx_deg, denom_var_deg) = match self
            .kind
        {
            CandidateKind::Ordinary => (true, 0, 0, 0, 0),
            CandidateKind::Inhomogeneous => (false, self.extra_idx_deg, self.extra_var_deg, 0, 0),
            CandidateKind::DenominatorIndexOnly
            | CandidateKind::DenominatorPositiveVariableDegree => {
                (true, 0, 0, self.extra_idx_deg, self.extra_var_deg)
            }
        };
        RecurrenceOptions {
            rec_len: self.rec_len,
            var_deg: self.var_deg,
            idx_deg: self.idx_deg,
            diff_deg: self.diff_deg,
            homogeneous,
            inhomo_var_deg,
            inhomo_idx_deg,
            denom_var_deg,
            denom_idx_deg,
            alternating_sign: self.alternating_sign,
            modular_prefilter,
            leader_x_power: self.laurent,
        }
    }
}

/// Number of left-hand charts of all denominator shapes `(a, v) != (0, 0)`
/// with `a <= max_idx`, `v <= max_var`: the constant chart plus the charts
/// `x^e` (`1 <= e <= v`) and `n^i` (`1 <= i <= a`) of each shape.
fn chart_count(max_idx: usize, max_var: usize) -> usize {
    let mut count = 0usize;
    for a in 0..=max_idx {
        for v in 0..=max_var {
            if a == 0 && v == 0 {
                continue;
            }
            count = count.saturating_add(1 + a + v);
        }
    }
    count
}

impl<'a> CandidateIterator<'a> {
    pub(super) fn new(m: usize, search: &'a AdaptiveSearchOptions) -> Self {
        Self::build(m, search, false)
    }

    /// The enumeration including the non-constant left-hand charts of every
    /// denominator shape (module documentation of `adaptive`).
    pub(super) fn with_charts(m: usize, search: &'a AdaptiveSearchOptions) -> Self {
        Self::build(m, search, true)
    }

    /// The number of left-hand monomials `x^e` tried per shape.
    fn laurent_count(&self) -> usize {
        if self.with_charts {
            self.search.max_laurent_x_power.saturating_add(1)
        } else {
            1
        }
    }

    fn denominator_shape_count(&self) -> usize {
        if self.with_charts {
            chart_count(self.search.max_denom_idx_deg, self.search.max_denom_var_deg)
                .saturating_mul(self.laurent_count())
        } else {
            self.search
                .max_denom_idx_deg
                .saturating_add(1)
                .saturating_mul(self.search.max_denom_var_deg.saturating_add(1))
                .saturating_sub(1)
        }
    }

    fn build(m: usize, search: &'a AdaptiveSearchOptions, with_charts: bool) -> Self {
        let min_rec_len = search.min_rec_len.max(1).min(m.saturating_sub(1));
        let max_rec_len = search.max_rec_len.min(m.saturating_sub(1));
        let base_count = inclusive_count(min_rec_len, max_rec_len)
            .saturating_mul(inclusive_count(search.min_diff_deg, search.max_diff_deg))
            .saturating_mul(inclusive_count(search.min_idx_deg, search.max_idx_deg))
            .saturating_mul(inclusive_count(search.min_var_deg, search.max_var_deg))
            .saturating_mul(if search.try_alternating_sign { 2 } else { 1 });
        let inhomogeneous_count = if search.try_inhomogeneous {
            inclusive_count(search.min_inhomo_idx_deg, search.max_inhomo_idx_deg).saturating_mul(
                inclusive_count(search.min_inhomo_var_deg, search.max_inhomo_var_deg),
            )
        } else {
            0
        };
        let laurent_count = if with_charts {
            search.max_laurent_x_power.saturating_add(1)
        } else {
            1
        };
        let denominator_count = if search.try_denominator {
            if with_charts {
                chart_count(search.max_denom_idx_deg, search.max_denom_var_deg)
                    .saturating_mul(laurent_count)
            } else {
                search
                    .max_denom_idx_deg
                    .saturating_add(1)
                    .saturating_mul(search.max_denom_var_deg.saturating_add(1))
                    .saturating_sub(1)
            }
        } else {
            0
        };
        let total_candidates = base_count.saturating_mul(
            laurent_count
                .saturating_add(inhomogeneous_count)
                .saturating_add(denominator_count),
        );
        let mut iterator = Self {
            search,
            rows: m,
            with_charts,
            min_rec_len,
            max_rec_len,
            heap: BinaryHeap::new(),
            visited: HashSet::new(),
            total_candidates,
        };
        if total_candidates == 0 {
            return iterator;
        }

        iterator.push(CandidateNode {
            kind: CandidateKind::Ordinary,
            rec_len: min_rec_len,
            diff_deg: search.min_diff_deg,
            idx_deg: search.min_idx_deg,
            var_deg: search.min_var_deg,
            alternating_sign: false,
            extra_idx_deg: 0,
            extra_var_deg: 0,
            chart: LeaderChart::CONSTANT,
            laurent: 0,
        });
        if inhomogeneous_count > 0 {
            iterator.push(CandidateNode {
                kind: CandidateKind::Inhomogeneous,
                rec_len: min_rec_len,
                diff_deg: search.min_diff_deg,
                idx_deg: search.min_idx_deg,
                var_deg: search.min_var_deg,
                alternating_sign: false,
                extra_idx_deg: search.min_inhomo_idx_deg,
                extra_var_deg: search.min_inhomo_var_deg,
                chart: LeaderChart::CONSTANT,
                laurent: 0,
            });
        }
        if denominator_count > 0 && search.max_denom_idx_deg > 0 {
            iterator.push(CandidateNode {
                kind: CandidateKind::DenominatorIndexOnly,
                rec_len: min_rec_len,
                diff_deg: search.min_diff_deg,
                idx_deg: search.min_idx_deg,
                var_deg: search.min_var_deg,
                alternating_sign: false,
                extra_idx_deg: 1,
                extra_var_deg: 0,
                chart: LeaderChart::CONSTANT,
                laurent: 0,
            });
        }
        if denominator_count > 0 && search.max_denom_var_deg > 0 {
            iterator.push(CandidateNode {
                kind: CandidateKind::DenominatorPositiveVariableDegree,
                rec_len: min_rec_len,
                diff_deg: search.min_diff_deg,
                idx_deg: search.min_idx_deg,
                var_deg: search.min_var_deg,
                alternating_sign: false,
                extra_idx_deg: 0,
                extra_var_deg: 1,
                chart: LeaderChart::CONSTANT,
                laurent: 0,
            });
        }
        iterator
    }

    pub(super) fn total_candidates(&self) -> usize {
        self.total_candidates
    }

    /// The number of left-hand variants (candidates other than the
    /// historical constant-chart, `e = 0` ones) of this enumeration.
    pub(super) fn variant_candidates(&self) -> usize {
        let historical = Self::build(self.rows, self.search, false).total_candidates;
        self.total_candidates.saturating_sub(historical)
    }

    pub(super) fn denominator_candidates(&self) -> usize {
        if !self.search.try_denominator {
            return 0;
        }
        let base_count = inclusive_count(self.min_rec_len, self.max_rec_len)
            .saturating_mul(inclusive_count(
                self.search.min_diff_deg,
                self.search.max_diff_deg,
            ))
            .saturating_mul(inclusive_count(
                self.search.min_idx_deg,
                self.search.max_idx_deg,
            ))
            .saturating_mul(inclusive_count(
                self.search.min_var_deg,
                self.search.max_var_deg,
            ))
            .saturating_mul(if self.search.try_alternating_sign {
                2
            } else {
                1
            });
        base_count.saturating_mul(self.denominator_shape_count())
    }

    fn raw_ordinal(&self, node: CandidateNode) -> u128 {
        let diff_count =
            inclusive_count(self.search.min_diff_deg, self.search.max_diff_deg) as u128;
        let idx_count = inclusive_count(self.search.min_idx_deg, self.search.max_idx_deg) as u128;
        let var_count = inclusive_count(self.search.min_var_deg, self.search.max_var_deg) as u128;
        let alternating_count = if self.search.try_alternating_sign {
            2
        } else {
            1
        };
        let inhomogeneous_count = if self.search.try_inhomogeneous {
            inclusive_count(
                self.search.min_inhomo_idx_deg,
                self.search.max_inhomo_idx_deg,
            )
            .saturating_mul(inclusive_count(
                self.search.min_inhomo_var_deg,
                self.search.max_inhomo_var_deg,
            ))
        } else {
            0
        };
        let denominator_count = if self.search.try_denominator {
            self.search
                .max_denom_idx_deg
                .saturating_add(1)
                .saturating_mul(self.search.max_denom_var_deg.saturating_add(1))
                .saturating_sub(1)
        } else {
            0
        };
        let block_count = 1usize
            .saturating_add(inhomogeneous_count)
            .saturating_add(denominator_count) as u128;
        let mut base = node.rec_len.saturating_sub(self.min_rec_len) as u128;
        base = base
            .saturating_mul(diff_count)
            .saturating_add(node.diff_deg.saturating_sub(self.search.min_diff_deg) as u128);
        base = base
            .saturating_mul(idx_count)
            .saturating_add(node.idx_deg.saturating_sub(self.search.min_idx_deg) as u128);
        base = base
            .saturating_mul(var_count)
            .saturating_add(node.var_deg.saturating_sub(self.search.min_var_deg) as u128);
        base = base
            .saturating_mul(alternating_count)
            .saturating_add(u128::from(node.alternating_sign));
        let within = match node.kind {
            CandidateKind::Ordinary => 0,
            CandidateKind::Inhomogeneous => {
                let variable_count = inclusive_count(
                    self.search.min_inhomo_var_deg,
                    self.search.max_inhomo_var_deg,
                ) as u128;
                1u128
                    .saturating_add(
                        ((node.extra_idx_deg - self.search.min_inhomo_idx_deg) as u128)
                            .saturating_mul(variable_count),
                    )
                    .saturating_add((node.extra_var_deg - self.search.min_inhomo_var_deg) as u128)
            }
            CandidateKind::DenominatorIndexOnly
            | CandidateKind::DenominatorPositiveVariableDegree => {
                let index_count = self.search.max_denom_idx_deg.saturating_add(1) as u128;
                1u128
                    .saturating_add(inhomogeneous_count as u128)
                    .saturating_add(
                        (node.extra_var_deg as u128)
                            .saturating_mul(index_count)
                            .saturating_add(node.extra_idx_deg as u128)
                            .saturating_sub(1),
                    )
            }
        };
        base.saturating_mul(block_count).saturating_add(within)
    }

    fn key(&self, node: CandidateNode) -> CandidateSortKey {
        let opts = node.to_options(self.search.modular_prefilter);
        let complexity = candidate_complexity(&opts);
        let has_denominator = opts.denom_var_deg > 0 || opts.denom_idx_deg > 0;
        let chart_weight = LEADER_CHART_WEIGHT
            .saturating_mul(node.chart.index_power.saturating_add(node.chart.x_power));
        (
            complexity.weighted_unknowns.saturating_add(chart_weight),
            complexity.raw_unknowns,
            opts.alternating_sign,
            has_denominator,
            !opts.homogeneous,
            opts.diff_deg,
            opts.rec_len,
            opts.idx_deg,
            opts.var_deg,
            opts.denom_idx_deg.saturating_add(opts.denom_var_deg),
            opts.inhomo_idx_deg.saturating_add(opts.inhomo_var_deg),
            (node.chart, node.laurent, self.raw_ordinal(node)),
        )
    }

    fn push(&mut self, node: CandidateNode) {
        if self.visited.insert(node) {
            self.heap.push(Reverse((self.key(node), node)));
        }
    }

    fn push_neighbors(&mut self, node: CandidateNode) {
        let mut neighbor = node;
        if node.rec_len < self.max_rec_len {
            neighbor.rec_len += 1;
            self.push(neighbor);
        }
        neighbor = node;
        if node.diff_deg < self.search.max_diff_deg {
            neighbor.diff_deg += 1;
            self.push(neighbor);
        }
        neighbor = node;
        if node.idx_deg < self.search.max_idx_deg {
            neighbor.idx_deg += 1;
            self.push(neighbor);
        }
        neighbor = node;
        if node.var_deg < self.search.max_var_deg {
            neighbor.var_deg += 1;
            self.push(neighbor);
        }
        neighbor = node;
        if !node.alternating_sign && self.search.try_alternating_sign {
            neighbor.alternating_sign = true;
            self.push(neighbor);
        }
        neighbor = node;
        match node.kind {
            CandidateKind::Ordinary => self.push_laurent_neighbor(node),
            CandidateKind::Inhomogeneous => {
                if node.extra_idx_deg < self.search.max_inhomo_idx_deg {
                    neighbor.extra_idx_deg += 1;
                    self.push(neighbor);
                }
                neighbor = node;
                if node.extra_var_deg < self.search.max_inhomo_var_deg {
                    neighbor.extra_var_deg += 1;
                    self.push(neighbor);
                }
            }
            CandidateKind::DenominatorIndexOnly => {
                if node.extra_idx_deg < self.search.max_denom_idx_deg {
                    neighbor.extra_idx_deg += 1;
                    self.push(neighbor);
                }
                self.push_chart_neighbors(node);
            }
            CandidateKind::DenominatorPositiveVariableDegree => {
                if node.extra_idx_deg < self.search.max_denom_idx_deg {
                    neighbor.extra_idx_deg += 1;
                    self.push(neighbor);
                }
                neighbor = node;
                if node.extra_var_deg < self.search.max_denom_var_deg {
                    neighbor.extra_var_deg += 1;
                    self.push(neighbor);
                }
                self.push_chart_neighbors(node);
            }
        }
    }

    /// The next chart along the `x` axis (`x^e -> x^(e+1)`, from the
    /// constant chart) or the `n` axis, inside the node's denominator shape.
    /// The next left-hand monomial `x^e -> x^(e+1)`.
    fn push_laurent_neighbor(&mut self, node: CandidateNode) {
        if self.with_charts && node.laurent < self.search.max_laurent_x_power {
            let mut neighbor = node;
            neighbor.laurent += 1;
            self.push(neighbor);
        }
    }

    fn push_chart_neighbors(&mut self, node: CandidateNode) {
        if !self.with_charts {
            return;
        }
        self.push_laurent_neighbor(node);
        let chart = node.chart;
        if chart.index_power == 0 && chart.x_power < node.extra_var_deg {
            let mut neighbor = node;
            neighbor.chart.x_power += 1;
            self.push(neighbor);
        }
        if chart.x_power == 0 && chart.index_power < node.extra_idx_deg {
            let mut neighbor = node;
            neighbor.chart.index_power += 1;
            self.push(neighbor);
        }
    }

    /// The next candidate with its left-hand chart.
    pub(super) fn next_candidate(&mut self) -> Option<(RecurrenceOptions, LeaderChart)> {
        let Reverse((_key, node)) = self.heap.pop()?;
        self.push_neighbors(node);
        Some((node.to_options(self.search.modular_prefilter), node.chart))
    }
}

impl Iterator for CandidateIterator<'_> {
    type Item = RecurrenceOptions;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_candidate().map(|(opts, _)| opts)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lazy_candidate_iterator_matches_eager_stable_order() {
        fn eager(m: usize, search: &AdaptiveSearchOptions) -> Vec<RecurrenceOptions> {
            let min_rec_len = search.min_rec_len.max(1).min(m.saturating_sub(1));
            let max_rec_len = search.max_rec_len.min(m.saturating_sub(1));
            let mut candidates = Vec::new();
            if min_rec_len > max_rec_len {
                return candidates;
            }
            for rec_len in min_rec_len..=max_rec_len {
                for diff_deg in search.min_diff_deg..=search.max_diff_deg {
                    for idx_deg in search.min_idx_deg..=search.max_idx_deg {
                        for var_deg in search.min_var_deg..=search.max_var_deg {
                            for alternating_sign in [false, true]
                                .into_iter()
                                .take(if search.try_alternating_sign { 2 } else { 1 })
                            {
                                candidates.push(RecurrenceOptions {
                                    rec_len,
                                    var_deg,
                                    idx_deg,
                                    diff_deg,
                                    homogeneous: true,
                                    inhomo_var_deg: 0,
                                    inhomo_idx_deg: 0,
                                    denom_var_deg: 0,
                                    denom_idx_deg: 0,
                                    alternating_sign,
                                    modular_prefilter: search.modular_prefilter,
                                    leader_x_power: 0,
                                });
                                if search.try_inhomogeneous {
                                    for inhomo_idx_deg in
                                        search.min_inhomo_idx_deg..=search.max_inhomo_idx_deg
                                    {
                                        for inhomo_var_deg in
                                            search.min_inhomo_var_deg..=search.max_inhomo_var_deg
                                        {
                                            candidates.push(RecurrenceOptions {
                                                rec_len,
                                                var_deg,
                                                idx_deg,
                                                diff_deg,
                                                homogeneous: false,
                                                inhomo_var_deg,
                                                inhomo_idx_deg,
                                                denom_var_deg: 0,
                                                denom_idx_deg: 0,
                                                alternating_sign,
                                                modular_prefilter: search.modular_prefilter,
                                                leader_x_power: 0,
                                            });
                                        }
                                    }
                                }
                                if search.try_denominator {
                                    for denom_var_deg in 0..=search.max_denom_var_deg {
                                        for denom_idx_deg in 0..=search.max_denom_idx_deg {
                                            if denom_var_deg == 0 && denom_idx_deg == 0 {
                                                continue;
                                            }
                                            candidates.push(RecurrenceOptions {
                                                rec_len,
                                                var_deg,
                                                idx_deg,
                                                diff_deg,
                                                homogeneous: true,
                                                inhomo_var_deg: 0,
                                                inhomo_idx_deg: 0,
                                                denom_var_deg,
                                                denom_idx_deg,
                                                alternating_sign,
                                                modular_prefilter: search.modular_prefilter,
                                                leader_x_power: 0,
                                            });
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            candidates.sort_by_key(|opts| {
                let complexity = candidate_complexity(opts);
                (
                    complexity.weighted_unknowns,
                    complexity.raw_unknowns,
                    opts.alternating_sign,
                    opts.denom_var_deg > 0 || opts.denom_idx_deg > 0,
                    !opts.homogeneous,
                    opts.diff_deg,
                    opts.rec_len,
                    opts.idx_deg,
                    opts.var_deg,
                    opts.denom_idx_deg.saturating_add(opts.denom_var_deg),
                    opts.inhomo_idx_deg.saturating_add(opts.inhomo_var_deg),
                )
            });
            candidates
        }

        fn signatures(candidates: impl IntoIterator<Item = RecurrenceOptions>) -> Vec<String> {
            candidates
                .into_iter()
                .map(|opts| {
                    format!(
                        "{}:{}:{}:{}:{}:{}:{}:{}:{}:{}",
                        opts.rec_len,
                        opts.diff_deg,
                        opts.idx_deg,
                        opts.var_deg,
                        opts.alternating_sign,
                        opts.homogeneous,
                        opts.inhomo_idx_deg,
                        opts.inhomo_var_deg,
                        opts.denom_idx_deg,
                        opts.denom_var_deg,
                    )
                })
                .collect()
        }

        let search = AdaptiveSearchOptions {
            min_rec_len: 1,
            max_rec_len: 3,
            min_var_deg: 0,
            max_var_deg: 2,
            min_idx_deg: 0,
            max_idx_deg: 2,
            min_diff_deg: 0,
            max_diff_deg: 1,
            try_inhomogeneous: true,
            min_inhomo_var_deg: 0,
            max_inhomo_var_deg: 2,
            min_inhomo_idx_deg: 0,
            max_inhomo_idx_deg: 2,
            try_denominator: true,
            try_alternating_sign: true,
            max_denom_var_deg: 2,
            max_denom_idx_deg: 2,
            ..Default::default()
        };
        let lazy = CandidateIterator::new(8, &search);
        assert_eq!(lazy.total_candidates(), eager(8, &search).len());
        assert_eq!(signatures(lazy), signatures(eager(8, &search)));
    }

    #[test]
    fn charted_enumeration_extends_the_historical_order() {
        let search = AdaptiveSearchOptions {
            min_rec_len: 1,
            max_rec_len: 2,
            max_var_deg: 2,
            max_idx_deg: 1,
            max_diff_deg: 1,
            try_denominator: true,
            max_denom_var_deg: 2,
            max_denom_idx_deg: 2,
            try_alternating_sign: true,
            ..Default::default()
        };
        let mut charted = CandidateIterator::with_charts(8, &search);
        let total = charted.total_candidates();
        let mut all = Vec::new();
        while let Some(candidate) = charted.next_candidate() {
            all.push(candidate);
        }
        assert_eq!(all.len(), total);
        for (opts, chart) in &all {
            assert!(chart.index_power <= opts.denom_idx_deg);
            assert!(chart.x_power <= opts.denom_var_deg);
            assert!(chart.index_power == 0 || chart.x_power == 0);
            assert!(opts.leader_x_power <= search.max_laurent_x_power);
        }
        // The historical candidates keep their relative order.
        let historical: Vec<String> = all
            .iter()
            .filter(|(opts, chart)| chart.is_constant() && opts.leader_x_power == 0)
            .map(|(opts, _)| format!("{opts:?}"))
            .collect();
        let expected: Vec<String> = CandidateIterator::new(8, &search)
            .map(|opts| format!("{opts:?}"))
            .collect();
        assert_eq!(historical, expected);
        let variants = CandidateIterator::with_charts(8, &search).variant_candidates();
        assert_eq!(variants, total - expected.len());
    }
}
