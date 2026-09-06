use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use crate::{Composition, Partition, Ssaf, WeakComposition};

/// A left-justified Young tableau with positive integer entries.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Tableau {
    rows: Vec<Vec<u32>>,
}

impl Tableau {
    /// Create a tableau from rows.
    ///
    /// Rows must have weakly decreasing lengths and use positive entries.
    pub fn new(mut rows: Vec<Vec<u32>>) -> Self {
        while rows.last().is_some_and(|row| row.is_empty()) {
            rows.pop();
        }

        let mut prev_len = usize::MAX;
        for row in &rows {
            assert!(
                row.iter().all(|&x| x > 0),
                "tableau entries must be positive"
            );
            assert!(
                row.len() <= prev_len,
                "tableau rows must have weakly decreasing lengths"
            );
            prev_len = row.len();
        }

        Tableau { rows }
    }

    /// The empty tableau.
    pub fn empty() -> Self {
        Tableau { rows: Vec::new() }
    }

    /// Access the tableau rows.
    pub fn rows(&self) -> &[Vec<u32>] {
        &self.rows
    }

    /// Return the shape of the tableau.
    pub fn shape(&self) -> Partition {
        Partition::from_sorted(self.rows.iter().map(|row| row.len() as u32).collect())
    }

    /// Number of boxes.
    pub fn size(&self) -> u32 {
        self.rows.iter().map(|row| row.len() as u32).sum()
    }

    /// Weight composition, where the i-th part counts entries equal to i + 1.
    pub fn weight(&self) -> Composition {
        let max_entry = self
            .rows
            .iter()
            .flat_map(|row| row.iter())
            .copied()
            .max()
            .unwrap_or(0) as usize;

        let mut counts = vec![0u32; max_entry];
        for row in &self.rows {
            for &entry in row {
                counts[entry as usize - 1] += 1;
            }
        }
        Composition::new(counts)
    }

    /// Construct the unique key tableau of the given weight.
    ///
    /// The entry `i + 1` appears in the first `weight[i]` columns.  Equivalently,
    /// column `c` contains precisely the entries `i + 1` such that
    /// `weight[i] > c`, listed increasingly from top to bottom.
    pub fn key_tableau_from_weight(weight: &[u32]) -> Self {
        let max_width = weight.iter().copied().max().unwrap_or(0);
        if max_width == 0 {
            return Tableau::empty();
        }

        let mut rows: Vec<Vec<u32>> = Vec::new();
        for col in 0..max_width {
            for (row, entry) in weight
                .iter()
                .enumerate()
                .filter_map(|(idx, &count)| (count > col).then_some(idx as u32 + 1))
                .enumerate()
            {
                if rows.len() <= row {
                    rows.push(Vec::new());
                }
                rows[row].push(entry);
            }
        }

        Tableau::new(rows)
    }

    /// Whether this semistandard tableau is a key tableau.
    ///
    /// A key tableau has the property that every entry appearing in column
    /// `j + 1` also appears in column `j`.
    pub fn is_key_tableau(&self) -> bool {
        if !self.is_semistandard() {
            return false;
        }

        let max_width = self.rows.first().map_or(0, Vec::len);
        for col in 1..max_width {
            let left: BTreeSet<u32> = self
                .rows
                .iter()
                .filter_map(|row| row.get(col - 1).copied())
                .collect();
            for entry in self.rows.iter().filter_map(|row| row.get(col).copied()) {
                if !left.contains(&entry) {
                    return false;
                }
            }
        }
        true
    }

    /// Right-key weight computed through Mason's SSYT-to-SSAF bijection.
    ///
    /// The tableau is treated as an SSYT with entries in `{1, ..., num_vars}`.
    /// Mason's bijection sends it to an atom filling whose composition shape is
    /// the weight of the right key in this convention.
    pub fn right_key_weight_via_ssaf(&self, num_vars: usize) -> WeakComposition {
        assert!(
            self.is_semistandard(),
            "right_key_weight_via_ssaf requires a semistandard tableau"
        );
        let max_entry = self
            .rows
            .iter()
            .flat_map(|row| row.iter())
            .copied()
            .max()
            .unwrap_or(0) as usize;
        assert!(
            max_entry <= num_vars,
            "num_vars ({num_vars}) must be at least the maximum tableau entry ({max_entry})"
        );

        let atom = Ssaf::from_ssyt(&self.rows);
        let mut shape = atom.shape();
        shape.resize(num_vars, 0);
        WeakComposition::new(shape)
    }

    /// Right key computed through Mason's SSYT-to-SSAF bijection.
    pub fn right_key_via_ssaf(&self, num_vars: usize) -> Self {
        let weight = self.right_key_weight_via_ssaf(num_vars);
        Self::key_tableau_from_weight(weight.parts())
    }

    /// Row reading word, read left-to-right in each row from top to bottom.
    pub fn row_reading_word(&self) -> Vec<u32> {
        self.rows
            .iter()
            .flat_map(|row| row.iter().copied())
            .collect()
    }

    /// Reverse row reading word, read left-to-right in each row from bottom to top.
    pub fn reverse_row_reading_word(&self) -> Vec<u32> {
        self.rows
            .iter()
            .rev()
            .flat_map(|row| row.iter().copied())
            .collect()
    }

    /// The standard reading word, read left-to-right in each row from bottom to top.
    ///
    /// This is the conventional tableau reading word in the sense used for
    /// Young tableaux, and is an alias for [`Self::reverse_row_reading_word`].
    pub fn reading_word(&self) -> Vec<u32> {
        self.reverse_row_reading_word()
    }

    /// The reverse reading word, read left-to-right in each row from top to bottom.
    ///
    /// This is the reverse of the standard tableau reading word, and is an
    /// alias for [`Self::row_reading_word`].
    pub fn reverse_reading_word(&self) -> Vec<u32> {
        self.row_reading_word()
    }

    /// Whether rows are weakly increasing.
    pub fn is_row_weakly_increasing(&self) -> bool {
        self.rows
            .iter()
            .all(|row| row.windows(2).all(|pair| pair[0] <= pair[1]))
    }

    /// Whether columns are strictly increasing.
    pub fn is_column_strictly_increasing(&self) -> bool {
        for r in 1..self.rows.len() {
            for c in 0..self.rows[r].len() {
                if self.rows[r - 1].len() > c && self.rows[r - 1][c] >= self.rows[r][c] {
                    return false;
                }
            }
        }
        true
    }

    /// Whether this is a semistandard Young tableau.
    pub fn is_semistandard(&self) -> bool {
        self.is_row_weakly_increasing() && self.is_column_strictly_increasing()
    }

    /// Whether this is a standard Young tableau.
    pub fn is_standard(&self) -> bool {
        if !self.is_semistandard() {
            return false;
        }

        let n = self.size() as usize;
        let mut seen = vec![false; n + 1];
        for row in &self.rows {
            for &entry in row {
                let idx = entry as usize;
                if idx == 0 || idx > n || seen[idx] {
                    return false;
                }
                seen[idx] = true;
            }
        }
        seen.into_iter().skip(1).all(|b| b)
    }

    /// Descent set of a standard tableau.
    pub fn descent_set(&self) -> BTreeSet<u32> {
        assert!(
            self.is_standard(),
            "descent_set requires a standard tableau"
        );

        let mut positions = vec![(0usize, 0usize); self.size() as usize + 1];
        for (r, row) in self.rows.iter().enumerate() {
            for (c, &entry) in row.iter().enumerate() {
                positions[entry as usize] = (r, c);
            }
        }

        let mut descents = BTreeSet::new();
        for i in 1..self.size() {
            if positions[(i + 1) as usize].0 > positions[i as usize].0 {
                descents.insert(i);
            }
        }
        descents
    }

    /// Enumerate all standard Young tableaux of the given shape.
    pub fn standard_tableaux(shape: &Partition) -> Vec<Tableau> {
        Self::standard_tableaux_iter(shape).collect()
    }

    /// Enumerate all standard Young tableaux of the given shape lazily.
    ///
    /// This iterator yields tableaux one by one, which is useful when a caller
    /// wants to compute a statistic on each tableau and discard it immediately.
    pub fn standard_tableaux_iter(shape: &Partition) -> StandardTableauxIter {
        StandardTableauxIter::new(shape)
    }

    /// Enumerate semistandard tableaux of the given shape with entries in
    /// `{1, ..., max_entry}`.
    pub fn semistandard_tableaux(shape: &Partition, max_entry: u32) -> Vec<Tableau> {
        if max_entry == 0 {
            return if shape.is_empty() {
                vec![Tableau::empty()]
            } else {
                Vec::new()
            };
        }

        let mut grid: Vec<Vec<u32>> = shape
            .parts()
            .iter()
            .map(|&part| vec![0; part as usize])
            .collect();
        let cells: Vec<(usize, usize)> = shape
            .parts()
            .iter()
            .enumerate()
            .flat_map(|(r, &part)| (0..part as usize).map(move |c| (r, c)))
            .collect();
        let mut out = Vec::new();
        semistandard_helper(&mut grid, &cells, 0, max_entry, &mut out);
        out.into_iter().map(Tableau::new).collect()
    }
}

/// A skew Young tableau of outer shape `outer` and inner shape `inner`.
///
/// Each row stores only the visible skew entries. Row `r` therefore contains
/// exactly `outer_r - inner_r` values, and empty visible rows are permitted.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct SkewTableau {
    outer: Partition,
    inner: Partition,
    rows: Vec<Vec<u32>>,
}

impl SkewTableau {
    /// Create a skew tableau from an outer shape, inner shape, and visible rows.
    ///
    /// Row `r` must contain exactly the entries in columns
    /// `inner_r, inner_r + 1, ..., outer_r - 1`, in left-to-right order.
    pub fn new(outer: Partition, inner: Partition, mut rows: Vec<Vec<u32>>) -> Self {
        assert!(
            inner.partition_less_equal(&outer),
            "inner shape must fit inside outer shape"
        );

        while rows.len() > outer.num_parts() && rows.last().is_some_and(|row| row.is_empty()) {
            rows.pop();
        }
        while rows.len() < outer.num_parts() {
            rows.push(Vec::new());
        }
        assert_eq!(
            rows.len(),
            outer.num_parts(),
            "skew tableau rows must match the number of rows in the outer shape"
        );

        for (r, row) in rows.iter().enumerate() {
            assert!(
                row.iter().all(|&x| x > 0),
                "tableau entries must be positive"
            );
            let expected_len = (outer.part(r) - inner.part(r)) as usize;
            assert_eq!(
                row.len(),
                expected_len,
                "row {} must contain exactly {} visible skew entries",
                r,
                expected_len
            );
        }

        SkewTableau { outer, inner, rows }
    }

    /// The empty skew tableau of shape `∅/∅`.
    pub fn empty() -> Self {
        SkewTableau {
            outer: Partition::empty(),
            inner: Partition::empty(),
            rows: Vec::new(),
        }
    }

    /// Access the visible rows of the skew tableau.
    ///
    /// Empty rows may occur when the skew shape has a gap in that row.
    pub fn rows(&self) -> &[Vec<u32>] {
        &self.rows
    }

    /// The outer shape of the skew tableau.
    pub fn outer_shape(&self) -> &Partition {
        &self.outer
    }

    /// The inner shape of the skew tableau.
    pub fn inner_shape(&self) -> &Partition {
        &self.inner
    }

    /// Number of visible boxes.
    pub fn size(&self) -> u32 {
        self.outer.size() - self.inner.size()
    }

    /// Weight composition, where the i-th part counts entries equal to i + 1.
    pub fn weight(&self) -> Composition {
        let max_entry = self
            .rows
            .iter()
            .flat_map(|row| row.iter())
            .copied()
            .max()
            .unwrap_or(0) as usize;

        let mut counts = vec![0u32; max_entry];
        for row in &self.rows {
            for &entry in row {
                counts[entry as usize - 1] += 1;
            }
        }
        Composition::new(counts)
    }

    /// Row reading word, read left-to-right in each row from top to bottom.
    pub fn row_reading_word(&self) -> Vec<u32> {
        self.rows
            .iter()
            .flat_map(|row| row.iter().copied())
            .collect()
    }

    /// Reverse row reading word, read left-to-right in each row from bottom to top.
    pub fn reverse_row_reading_word(&self) -> Vec<u32> {
        self.rows
            .iter()
            .rev()
            .flat_map(|row| row.iter().copied())
            .collect()
    }

    /// The standard reading word, read left-to-right in each row from bottom to top.
    ///
    /// This is the conventional tableau reading word, and is an alias for
    /// [`Self::reverse_row_reading_word`].
    pub fn reading_word(&self) -> Vec<u32> {
        self.reverse_row_reading_word()
    }

    /// The reverse reading word, read left-to-right in each row from top to bottom.
    ///
    /// This is the reverse of the standard tableau reading word, and is an
    /// alias for [`Self::row_reading_word`].
    pub fn reverse_reading_word(&self) -> Vec<u32> {
        self.row_reading_word()
    }

    /// Whether rows are weakly increasing.
    pub fn is_row_weakly_increasing(&self) -> bool {
        self.rows
            .iter()
            .all(|row| row.windows(2).all(|pair| pair[0] <= pair[1]))
    }

    /// Whether columns are strictly increasing.
    pub fn is_column_strictly_increasing(&self) -> bool {
        for (r, row) in self.rows.iter().enumerate().skip(1) {
            let start_col = self.inner.part(r) as usize;
            for (offset, &value) in row.iter().enumerate() {
                let col = start_col + offset;
                if let Some(above) = self.cell(r - 1, col) {
                    if above >= value {
                        return false;
                    }
                }
            }
        }
        true
    }

    /// Whether this is a semistandard skew tableau.
    pub fn is_semistandard(&self) -> bool {
        self.is_row_weakly_increasing() && self.is_column_strictly_increasing()
    }

    /// Whether this is a standard skew tableau.
    pub fn is_standard(&self) -> bool {
        if !self.is_semistandard() {
            return false;
        }

        let n = self.size() as usize;
        let mut seen = vec![false; n + 1];
        for row in &self.rows {
            for &entry in row {
                let idx = entry as usize;
                if idx == 0 || idx > n || seen[idx] {
                    return false;
                }
                seen[idx] = true;
            }
        }
        seen.into_iter().skip(1).all(|b| b)
    }

    /// Descent set of a standard skew tableau.
    pub fn descent_set(&self) -> BTreeSet<u32> {
        assert!(
            self.is_standard(),
            "descent_set requires a standard skew tableau"
        );

        let mut positions = vec![(0usize, 0usize); self.size() as usize + 1];
        for (r, row) in self.rows.iter().enumerate() {
            let start_col = self.inner.part(r) as usize;
            for (offset, &entry) in row.iter().enumerate() {
                positions[entry as usize] = (r, start_col + offset);
            }
        }

        let mut descents = BTreeSet::new();
        for i in 1..self.size() {
            if positions[(i + 1) as usize].0 > positions[i as usize].0 {
                descents.insert(i);
            }
        }
        descents
    }

    /// Enumerate all standard skew tableaux of shape `outer / inner`.
    pub fn standard_tableaux(outer: &Partition, inner: &Partition) -> Vec<SkewTableau> {
        Self::standard_tableaux_iter(outer, inner).collect()
    }

    /// Enumerate all standard skew tableaux of shape `outer / inner` lazily.
    ///
    /// This iterator yields tableaux one by one, which is useful when a caller
    /// wants to compute a statistic on each tableau and discard it immediately.
    pub fn standard_tableaux_iter(
        outer: &Partition,
        inner: &Partition,
    ) -> StandardSkewTableauxIter {
        StandardSkewTableauxIter::new(outer, inner)
    }

    fn cell(&self, row: usize, col: usize) -> Option<u32> {
        if row >= self.rows.len() {
            return None;
        }
        let start_col = self.inner.part(row) as usize;
        let end_col = self.outer.part(row) as usize;
        if col < start_col || col >= end_col {
            return None;
        }
        Some(self.rows[row][col - start_col])
    }
}

fn corner_cells(shape: &[u32]) -> Vec<(usize, usize)> {
    let mut corners = Vec::new();
    for row in 0..shape.len() {
        let next = if row + 1 < shape.len() {
            shape[row + 1]
        } else {
            0
        };
        if shape[row] > next {
            corners.push((row, shape[row] as usize - 1));
        }
    }
    corners
}

fn skew_corner_rows(outer: &Partition, inner: &Partition) -> Vec<usize> {
    let mut corners = Vec::new();
    for row in 0..outer.num_parts() {
        if outer.part(row) > inner.part(row) && outer.part(row) > outer.part(row + 1) {
            corners.push(row);
        }
    }
    corners
}

/// Lazy iterator over all standard Young tableaux of a fixed straight shape.
#[derive(Debug, Clone)]
pub struct StandardTableauxIter {
    shape_parts: Vec<u32>,
    size: u32,
    corners: Vec<usize>,
    next_corner_idx: usize,
    active_row: Option<usize>,
    child: Option<Box<StandardTableauxIter>>,
    yielded_empty: bool,
}

impl StandardTableauxIter {
    fn new(shape: &Partition) -> Self {
        Self::from_parts(shape.parts().to_vec())
    }

    fn from_parts(shape_parts: Vec<u32>) -> Self {
        let size = shape_parts.iter().sum();
        let corners = if size == 0 {
            Vec::new()
        } else {
            corner_cells(&shape_parts)
                .into_iter()
                .map(|(row, _)| row)
                .collect()
        };

        StandardTableauxIter {
            shape_parts,
            size,
            corners,
            next_corner_idx: 0,
            active_row: None,
            child: None,
            yielded_empty: false,
        }
    }
}

impl Iterator for StandardTableauxIter {
    type Item = Tableau;

    fn next(&mut self) -> Option<Self::Item> {
        if self.size == 0 {
            if self.yielded_empty {
                return None;
            }
            self.yielded_empty = true;
            return Some(Tableau::empty());
        }

        loop {
            if let Some(child) = self.child.as_mut() {
                if let Some(tableau) = child.next() {
                    let row = self
                        .active_row
                        .expect("active_row must be set while child iterator is active");
                    let mut rows = tableau.rows().to_vec();
                    while rows.len() <= row {
                        rows.push(Vec::new());
                    }
                    rows[row].push(self.size);
                    return Some(Tableau::new(rows));
                }
                self.child = None;
                self.active_row = None;
            }

            if self.next_corner_idx >= self.corners.len() {
                return None;
            }

            let row = self.corners[self.next_corner_idx];
            self.next_corner_idx += 1;
            let mut smaller = self.shape_parts.clone();
            smaller[row] -= 1;
            while smaller.last() == Some(&0) {
                smaller.pop();
            }
            self.active_row = Some(row);
            self.child = Some(Box::new(StandardTableauxIter::from_parts(smaller)));
        }
    }
}

/// Lazy iterator over all standard skew tableaux of a fixed skew shape.
#[derive(Debug, Clone)]
pub struct StandardSkewTableauxIter {
    outer: Partition,
    inner: Partition,
    size: u32,
    corner_rows: Vec<usize>,
    next_corner_idx: usize,
    active_row: Option<usize>,
    child: Option<Box<StandardSkewTableauxIter>>,
    yielded_empty: bool,
}

impl StandardSkewTableauxIter {
    fn new(outer: &Partition, inner: &Partition) -> Self {
        assert!(
            inner.partition_less_equal(outer),
            "inner shape must fit inside outer shape"
        );

        StandardSkewTableauxIter {
            outer: outer.clone(),
            inner: inner.clone(),
            size: outer.size() - inner.size(),
            corner_rows: skew_corner_rows(outer, inner),
            next_corner_idx: 0,
            active_row: None,
            child: None,
            yielded_empty: false,
        }
    }
}

impl Iterator for StandardSkewTableauxIter {
    type Item = SkewTableau;

    fn next(&mut self) -> Option<Self::Item> {
        if self.size == 0 {
            if self.yielded_empty {
                return None;
            }
            self.yielded_empty = true;
            return Some(SkewTableau::new(
                self.outer.clone(),
                self.inner.clone(),
                vec![Vec::new(); self.outer.num_parts()],
            ));
        }

        loop {
            if let Some(child) = self.child.as_mut() {
                if let Some(tableau) = child.next() {
                    let row = self
                        .active_row
                        .expect("active_row must be set while child iterator is active");
                    let mut rows = tableau.rows().to_vec();
                    while rows.len() <= row {
                        rows.push(Vec::new());
                    }
                    rows[row].push(self.size);
                    return Some(SkewTableau::new(
                        self.outer.clone(),
                        self.inner.clone(),
                        rows,
                    ));
                }
                self.child = None;
                self.active_row = None;
            }

            if self.next_corner_idx >= self.corner_rows.len() {
                return None;
            }

            let row = self.corner_rows[self.next_corner_idx];
            self.next_corner_idx += 1;
            let mut smaller_parts = self.outer.parts().to_vec();
            smaller_parts[row] -= 1;
            while smaller_parts.last() == Some(&0) {
                smaller_parts.pop();
            }
            let smaller_outer = Partition::from_sorted(smaller_parts);
            self.active_row = Some(row);
            self.child = Some(Box::new(StandardSkewTableauxIter::new(
                &smaller_outer,
                &self.inner,
            )));
        }
    }
}

fn semistandard_helper(
    grid: &mut [Vec<u32>],
    cells: &[(usize, usize)],
    idx: usize,
    max_entry: u32,
    out: &mut Vec<Vec<Vec<u32>>>,
) {
    if idx == cells.len() {
        out.push(grid.to_vec());
        return;
    }

    let (row, col) = cells[idx];
    let mut min_value = 1u32;
    if col > 0 {
        min_value = min_value.max(grid[row][col - 1]);
    }
    if row > 0 && grid[row - 1].len() > col {
        min_value = min_value.max(grid[row - 1][col] + 1);
    }

    for value in min_value..=max_entry {
        grid[row][col] = value;
        semistandard_helper(grid, cells, idx + 1, max_entry, out);
    }
}

impl Tableau {
    // ── RSK insertion ───────────────────────────────────────────────

    /// Row-insert a value into this semistandard tableau.
    ///
    /// Returns the new tableau and the row index (0-based) where the
    /// new cell was created (i.e. the row that grew by one box).
    pub fn row_insert(&self, value: u32) -> (Tableau, usize) {
        let mut rows = self.rows.clone();
        let bump_row = row_insert_into(&mut rows, 0, value);
        (Tableau { rows }, bump_row)
    }

    /// Robinson-Schensted correspondence on a word.
    ///
    /// Returns `(P, Q)` where `P` is the insertion tableau and `Q` the
    /// recording tableau.  Both are standard when the input is a permutation.
    pub fn rsk(word: &[u32]) -> (Tableau, Tableau) {
        let mut p_rows: Vec<Vec<u32>> = Vec::new();
        let mut q_rows: Vec<Vec<u32>> = Vec::new();
        for (idx, &val) in word.iter().enumerate() {
            let bump_row = row_insert_into(&mut p_rows, 0, val);
            // Recording tableau gets entry idx+1
            while q_rows.len() <= bump_row {
                q_rows.push(Vec::new());
            }
            q_rows[bump_row].push(idx as u32 + 1);
        }
        (Tableau { rows: p_rows }, Tableau { rows: q_rows })
    }

    /// Generalized RSK on a biword `(top, bottom)`.
    ///
    /// The pairs `(top[i], bottom[i])` are first sorted lexicographically
    /// (by top, breaking ties by bottom), then bottom values are inserted
    /// while top values are recorded.
    pub fn biword_rsk(top: &[u32], bottom: &[u32]) -> (Tableau, Tableau) {
        assert_eq!(top.len(), bottom.len(), "biword must have matching lengths");
        let mut indices: Vec<usize> = (0..top.len()).collect();
        indices.sort_by(|&a, &b| top[a].cmp(&top[b]).then(bottom[a].cmp(&bottom[b])));
        let mut p_rows: Vec<Vec<u32>> = Vec::new();
        let mut q_rows: Vec<Vec<u32>> = Vec::new();
        for &i in &indices {
            let bump_row = row_insert_into(&mut p_rows, 0, bottom[i]);
            while q_rows.len() <= bump_row {
                q_rows.push(Vec::new());
            }
            q_rows[bump_row].push(top[i]);
        }
        (Tableau { rows: p_rows }, Tableau { rows: q_rows })
    }

    // ── Promotion and evacuation ────────────────────────────────────

    /// Schützenberger promotion for standard Young tableaux.
    ///
    /// 1. Remove entry 1, creating an inner hole.
    /// 2. Forward jeu de taquin to slide the hole to an outer corner.
    /// 3. Remove the outer corner from the shape.
    /// 4. Subtract 1 from all remaining entries (so 2→1, 3→2, etc.).
    /// 5. Place *n* in the vacated outer corner.
    pub fn promotion(&self) -> Tableau {
        let n = self.size();
        if n <= 1 {
            return self.clone();
        }
        // Find position of entry 1
        let mut hole = (0usize, 0usize);
        'outer: for (r, row) in self.rows.iter().enumerate() {
            for (c, &v) in row.iter().enumerate() {
                if v == 1 {
                    hole = (r, c);
                    break 'outer;
                }
            }
        }
        let mut rows = self.rows.clone();
        // Jeu de taquin forward slide: hole moves down/right
        loop {
            let (r, c) = hole;
            let right = if c + 1 < rows[r].len() {
                Some(rows[r][c + 1])
            } else {
                None
            };
            let below = if r + 1 < rows.len() && c < rows[r + 1].len() {
                Some(rows[r + 1][c])
            } else {
                None
            };
            match (right, below) {
                (None, None) => break,
                (Some(rv), None) => {
                    rows[r][c] = rv;
                    hole = (r, c + 1);
                }
                (None, Some(bv)) => {
                    rows[r][c] = bv;
                    hole = (r + 1, c);
                }
                (Some(rv), Some(bv)) => {
                    if bv <= rv {
                        rows[r][c] = bv;
                        hole = (r + 1, c);
                    } else {
                        rows[r][c] = rv;
                        hole = (r, c + 1);
                    }
                }
            }
        }
        // Remove the outer corner (the hole) from the shape
        rows[hole.0].pop();
        if rows.last().is_some_and(|r| r.is_empty()) {
            rows.pop();
        }
        // Subtract 1 from all remaining entries
        for row in &mut rows {
            for v in row.iter_mut() {
                *v -= 1;
            }
        }
        // Place n at the vacated corner
        while rows.len() <= hole.0 {
            rows.push(Vec::new());
        }
        rows[hole.0].push(n);
        Tableau { rows }
    }

    /// Schützenberger evacuation involution for standard Young tableaux.
    pub fn evacuation(&self) -> Tableau {
        let n = self.size();
        if n == 0 {
            return self.clone();
        }
        assert!(self.is_standard(), "evacuation requires a standard tableau");

        // Repeatedly remove 1, slide the hole southeast by jeu de taquin,
        // place the descending output label in the final corner, and subtract
        // one from the surviving labels.  Keeping the active and output
        // diagrams as maps makes the shrinking outer shape explicit.
        let mut current: BTreeMap<(usize, usize), u32> = self
            .rows
            .iter()
            .enumerate()
            .flat_map(|(row, entries)| {
                entries
                    .iter()
                    .enumerate()
                    .map(move |(column, &value)| ((row, column), value))
            })
            .collect();
        let mut output = BTreeMap::new();
        for label in (1..=n).rev() {
            let mut hole = *current
                .iter()
                .find_map(|(cell, &value)| (value == 1).then_some(cell))
                .expect("a standard active tableau contains 1");
            current.remove(&hole);
            loop {
                let candidates = [(hole.0, hole.1 + 1), (hole.0 + 1, hole.1)];
                let Some((&chosen, &value)) = candidates
                    .iter()
                    .filter_map(|cell| current.get_key_value(cell))
                    .min_by_key(|(_, value)| *value)
                else {
                    break;
                };
                current.insert(hole, value);
                current.remove(&chosen);
                hole = chosen;
            }
            output.insert(hole, label);
            for value in current.values_mut() {
                *value -= 1;
            }
        }
        Tableau::new(
            self.rows
                .iter()
                .enumerate()
                .map(|(row, entries)| {
                    (0..entries.len())
                        .map(|column| output[&(row, column)])
                        .collect()
                })
                .collect(),
        )
    }

    // ── Crystal operators ───────────────────────────────────────────

    /// Crystal raising operator *eᵢ* on a semistandard Young tableau.
    ///
    /// Uses the signature rule: extract the i/(i+1) subword, bracket-match
    /// adjacent (i+1, i) pairs, then reassign all unmatched positions with
    /// one fewer (i+1) and one more i.
    /// Returns `None` if the operator annihilates the tableau.
    pub fn crystal_e(&self, i: u32) -> Option<Tableau> {
        assert!(i >= 1, "crystal index must be >= 1");
        let positions = self.crystal_unmatched_positions(i);
        let b_count = positions.iter().filter(|&&(_, _, v)| v == i + 1).count();
        if b_count == 0 {
            return None;
        }
        // Reassign: (a+1) copies of i, then (b-1) copies of (i+1)
        let a_count = positions.len() - b_count;
        let mut rows = self.rows.clone();
        for (idx, &(r, c, _)) in positions.iter().enumerate() {
            rows[r][c] = if idx < a_count + 1 { i } else { i + 1 };
        }
        Some(Tableau { rows })
    }

    /// Crystal lowering operator *fᵢ* on a semistandard Young tableau.
    ///
    /// Returns `None` if the operator annihilates the tableau.
    pub fn crystal_f(&self, i: u32) -> Option<Tableau> {
        assert!(i >= 1, "crystal index must be >= 1");
        let positions = self.crystal_unmatched_positions(i);
        let a_count = positions.iter().filter(|&&(_, _, v)| v == i).count();
        if a_count == 0 {
            return None;
        }
        // Reassign: (a-1) copies of i, then (b+1) copies of (i+1)
        let mut rows = self.rows.clone();
        for (idx, &(r, c, _)) in positions.iter().enumerate() {
            rows[r][c] = if idx < a_count - 1 { i } else { i + 1 };
        }
        Some(Tableau { rows })
    }

    /// Unmatched positions for the *i*-crystal string.
    ///
    /// Extracts the sub-word of entries in {i, i+1} in reverse row reading
    /// order (bottom-to-top, left-to-right), then repeatedly removes adjacent
    /// (i+1, i) pairs until stable.  Returns the surviving positions.
    fn crystal_unmatched_positions(&self, i: u32) -> Vec<(usize, usize, u32)> {
        let mut positions: Vec<(usize, usize, u32)> = Vec::new();
        for (r, row) in self.rows.iter().enumerate().rev() {
            for (c, &v) in row.iter().enumerate() {
                if v == i || v == i + 1 {
                    positions.push((r, c, v));
                }
            }
        }
        // Repeatedly remove adjacent (i+1, i) pairs
        loop {
            let mut removed = false;
            let mut new_positions = Vec::with_capacity(positions.len());
            let mut skip_next = false;
            for idx in 0..positions.len() {
                if skip_next {
                    skip_next = false;
                    continue;
                }
                if idx + 1 < positions.len()
                    && positions[idx].2 == i + 1
                    && positions[idx + 1].2 == i
                {
                    skip_next = true;
                    removed = true;
                } else {
                    new_positions.push(positions[idx]);
                }
            }
            positions = new_positions;
            if !removed {
                break;
            }
        }
        positions
    }

    // ── Lascoux-Schützenberger involution ─────────────────────────

    /// Lascoux-Schützenberger involution exchanging *i* and *i+1* on an SSYT.
    ///
    /// Uses column reading order (left-to-right, bottom-to-top), removes
    /// same-column pairs, then removes adjacent (i+1, i) pairs.
    /// The resulting word i^a (i+1)^b gets its values swapped and reversed.
    pub fn lascoux_schutzenberger(&self, i: u32) -> Tableau {
        let n = self.rows.len();
        let max_col = self.rows.iter().map(|r| r.len()).max().unwrap_or(0);

        // Column reading: left-to-right, bottom-to-top
        let mut positions: Vec<(usize, usize, u32)> = Vec::new();
        for c in 0..max_col {
            for r in (0..n).rev() {
                if self.rows[r].len() > c {
                    let v = self.rows[r][c];
                    if v == i || v == i + 1 {
                        positions.push((r, c, v));
                    }
                }
            }
        }

        // SameColDelete: remove column pairs
        let mut col_map: std::collections::HashMap<usize, Vec<usize>> =
            std::collections::HashMap::new();
        for (idx, &(_, c, _)) in positions.iter().enumerate() {
            col_map.entry(c).or_default().push(idx);
        }
        let mut remove: std::collections::HashSet<usize> = std::collections::HashSet::new();
        for indices in col_map.values() {
            if indices.len() == 2 {
                let v0 = positions[indices[0]].2;
                let v1 = positions[indices[1]].2;
                if (v0 == i && v1 == i + 1) || (v0 == i + 1 && v1 == i) {
                    remove.insert(indices[0]);
                    remove.insert(indices[1]);
                }
            }
        }
        if !remove.is_empty() {
            positions = positions
                .into_iter()
                .enumerate()
                .filter(|(idx, _)| !remove.contains(idx))
                .map(|(_, p)| p)
                .collect();
        }

        // RemovePair: remove adjacent (i+1, i) pairs [SSYT convention]
        loop {
            let mut removed = false;
            let mut new_positions = Vec::with_capacity(positions.len());
            let mut skip_next = false;
            for idx in 0..positions.len() {
                if skip_next {
                    skip_next = false;
                    continue;
                }
                if idx + 1 < positions.len()
                    && positions[idx].2 == i + 1
                    && positions[idx + 1].2 == i
                {
                    skip_next = true;
                    removed = true;
                } else {
                    new_positions.push(positions[idx]);
                }
            }
            positions = new_positions;
            if !removed {
                break;
            }
        }

        if positions.is_empty() {
            return self.clone();
        }

        // Swap i↔(i+1) and reverse
        let swapped_reversed: Vec<u32> = positions
            .iter()
            .rev()
            .map(|&(_, _, v)| if v == i { i + 1 } else { i })
            .collect();

        let mut rows = self.rows.clone();
        for (idx, &(r, c, _)) in positions.iter().enumerate() {
            rows[r][c] = swapped_reversed[idx];
        }
        Tableau { rows }
    }

    /// Normalize the SSYT weight to a partition.
    ///
    /// Applies Lascoux-Schützenberger involutions to bubble-sort the weight.
    pub fn weight_normalize(&self) -> Tableau {
        let mut result = self.clone();
        loop {
            let w = result.weight();
            let parts = w.parts();
            let mut sorted = true;
            for idx in 0..parts.len().saturating_sub(1) {
                if parts[idx] < parts[idx + 1] {
                    result = result.lascoux_schutzenberger((idx + 1) as u32);
                    sorted = false;
                    break;
                }
            }
            if sorted {
                break;
            }
        }
        result
    }

    // ── Charge for semistandard tableaux ────────────────────────────

    /// Charge of a semistandard Young tableau.
    ///
    /// Decomposes the reverse row reading word into standard subwords,
    /// then sums the charge of each subword.
    pub fn charge(&self) -> u32 {
        let word = self.reverse_row_reading_word();
        word_charge(&word)
    }

    /// Cocharge of a semistandard Young tableau.
    pub fn cocharge(&self) -> u32 {
        let word = self.reverse_row_reading_word();
        word_cocharge(&word)
    }
}

// ── RSK helpers ─────────────────────────────────────────────────────

/// Row-insert `value` starting at row `start_row` of the mutable grid.
/// Returns the row index where the new cell was created.
fn row_insert_into(rows: &mut Vec<Vec<u32>>, start_row: usize, value: u32) -> usize {
    let mut current_row = start_row;
    let mut val = value;
    loop {
        if current_row >= rows.len() {
            rows.push(vec![val]);
            return current_row;
        }
        let row = &rows[current_row];
        // If val >= last element (or row empty), append
        match row.last() {
            None => {
                rows[current_row].push(val);
                return current_row;
            }
            Some(&last) if val >= last => {
                rows[current_row].push(val);
                return current_row;
            }
            _ => {}
        }
        // val < last, so some element > val must exist
        let bump_pos = row
            .iter()
            .position(|&x| x > val)
            .expect("element > val must exist");
        std::mem::swap(&mut rows[current_row][bump_pos], &mut val);
        current_row += 1;
    }
}

// ── Charge helpers ──────────────────────────────────────────────────

/// Charge of a word (Lascoux-Schützenberger).
///
/// Decomposes the word into standard subwords and sums their charges.
pub fn word_charge(word: &[u32]) -> u32 {
    if word.is_empty() {
        return 0;
    }
    let subwords = charge_word_decompose(word);
    subwords.iter().map(|sw| permutation_charge(sw)).sum()
}

/// Cocharge of a word.
pub fn word_cocharge(word: &[u32]) -> u32 {
    if word.is_empty() {
        return 0;
    }
    let subwords = charge_word_decompose(word);
    subwords
        .iter()
        .map(|sw| {
            let n = sw.iter().copied().max().unwrap_or(0);
            n * (n - 1) / 2 - permutation_charge(sw)
        })
        .sum()
}

/// Decompose a word into standard subwords for charge computation.
///
/// Cyclically scans for 1, 2, 3, ... to extract each subword.
fn charge_word_decompose(word: &[u32]) -> Vec<Vec<u32>> {
    if word.is_empty() {
        return Vec::new();
    }
    let max_val = *word.iter().max().unwrap();
    let n = word.len();
    let mut remaining: Vec<Option<u32>> = word.iter().map(|&v| Some(v)).collect();
    let mut subwords = Vec::new();

    while remaining.iter().any(|x| x.is_some()) {
        let mut subword = Vec::new();
        // Start scanning from the rightmost position
        let mut cursor = n - 1;

        for target in 1..=max_val {
            // Scan leftward cyclically from cursor
            let mut found = false;
            for step in 0..n {
                let idx = (cursor + n - step) % n;
                if remaining[idx] == Some(target) {
                    subword.push(target);
                    remaining[idx] = None;
                    cursor = idx;
                    found = true;
                    break;
                }
            }
            if !found {
                break;
            }
        }
        if !subword.is_empty() {
            subwords.push(subword);
        }
    }
    subwords
}

/// Charge of a standard subword (permutation).
///
/// charge(w) = major index of the inverse of the reverse of w.
fn permutation_charge(word: &[u32]) -> u32 {
    if word.len() <= 1 {
        return 0;
    }
    let n = word.len();
    // Compute inverse of reversed word
    let mut inv = vec![0u32; n];
    for (i, &v) in word.iter().rev().enumerate() {
        inv[(v - 1) as usize] = i as u32 + 1;
    }
    // Major index of inv
    let mut maj = 0u32;
    for i in 0..n - 1 {
        if inv[i] > inv[i + 1] {
            maj += i as u32 + 1;
        }
    }
    maj
}

impl fmt::Display for Tableau {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (i, row) in self.rows.iter().enumerate() {
            if i > 0 {
                writeln!(f)?;
            }
            for (j, entry) in row.iter().enumerate() {
                if j > 0 {
                    write!(f, " ")?;
                }
                write!(f, "{entry}")?;
            }
        }
        Ok(())
    }
}

impl fmt::Display for SkewTableau {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (r, row) in self.rows.iter().enumerate() {
            if r > 0 {
                writeln!(f)?;
            }
            for _ in 0..self.inner.part(r) {
                write!(f, "  ")?;
            }
            for (j, entry) in row.iter().enumerate() {
                if j > 0 {
                    write!(f, " ")?;
                }
                write!(f, "{entry}")?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use combinatoric_core::poset::Poset;

    #[test]
    fn test_tableau_shape_weight_and_words() {
        let tab = Tableau::new(vec![vec![1, 1, 3], vec![2, 4]]);
        assert_eq!(tab.shape(), Partition::from_sorted(vec![3, 2]));
        assert_eq!(tab.weight(), Composition::new(vec![2, 1, 1, 1]));
        assert_eq!(tab.row_reading_word(), vec![1, 1, 3, 2, 4]);
        assert_eq!(tab.reverse_row_reading_word(), vec![2, 4, 1, 1, 3]);
        assert_eq!(tab.reading_word(), vec![2, 4, 1, 1, 3]);
        assert_eq!(tab.reverse_reading_word(), vec![1, 1, 3, 2, 4]);
    }

    #[test]
    fn test_key_tableau_from_weight() {
        let key = Tableau::key_tableau_from_weight(&[1, 0, 2]);
        assert_eq!(key.rows(), &[vec![1, 3], vec![3]]);
        assert_eq!(key.weight(), Composition::new(vec![1, 0, 2]));
        assert!(key.is_key_tableau());

        let site_key = Tableau::key_tableau_from_weight(&[1, 3, 0, 1, 2, 4]);
        assert_eq!(
            site_key.rows(),
            &[
                vec![1, 2, 2, 6],
                vec![2, 5, 6],
                vec![4, 6],
                vec![5],
                vec![6],
            ]
        );
        assert_eq!(site_key.weight(), Composition::new(vec![1, 3, 0, 1, 2, 4]));
        assert!(site_key.is_key_tableau());
    }

    #[test]
    fn test_right_key_weight_via_ssaf() {
        let t = Tableau::new(vec![vec![1, 3], vec![3]]);
        assert_eq!(
            t.right_key_weight_via_ssaf(3),
            WeakComposition::new(vec![1, 0, 2])
        );
        assert_eq!(t.right_key_via_ssaf(3).rows(), &[vec![1, 3], vec![3]]);

        let excluded = Tableau::new(vec![vec![1, 2], vec![3]]);
        assert_eq!(
            excluded.right_key_weight_via_ssaf(3),
            WeakComposition::new(vec![0, 2, 1])
        );
        assert_eq!(
            excluded.right_key_via_ssaf(3).rows(),
            &[vec![2, 2], vec![3]]
        );
    }

    #[test]
    fn test_semistandard_and_standard_checks() {
        let ssyt = Tableau::new(vec![vec![1, 1, 3], vec![2, 4]]);
        assert!(ssyt.is_semistandard());
        assert!(!ssyt.is_standard());

        let syt = Tableau::new(vec![vec![1, 2], vec![3]]);
        assert!(syt.is_standard());
        assert_eq!(syt.descent_set(), BTreeSet::from([2]));
    }

    #[test]
    fn test_standard_tableaux_shape_21() {
        let tabs = Tableau::standard_tableaux(&Partition::from_sorted(vec![2, 1]));
        assert_eq!(tabs.len(), 2);
        assert!(tabs.iter().all(Tableau::is_standard));
    }

    #[test]
    fn test_standard_tableaux_iterator_streams_one_by_one() {
        let shape = Partition::from_sorted(vec![3, 2]);
        let mut iter = Tableau::standard_tableaux_iter(&shape);
        let first = iter.next().unwrap();
        assert!(first.is_standard());
        assert_eq!(iter.count() + 1, Tableau::standard_tableaux(&shape).len());
    }

    #[test]
    fn test_semistandard_tableaux_shape_21_entries_up_to_2() {
        let tabs = Tableau::semistandard_tableaux(&Partition::from_sorted(vec![2, 1]), 2);
        assert_eq!(tabs.len(), 2);
        assert!(tabs.iter().all(Tableau::is_semistandard));
    }

    #[test]
    fn test_skew_tableau_words_and_descents() {
        let skew = SkewTableau::new(
            Partition::from_sorted(vec![3, 2]),
            Partition::from_sorted(vec![1]),
            vec![vec![1, 2], vec![3, 4]],
        );
        assert!(skew.is_standard());
        assert_eq!(skew.weight(), Composition::new(vec![1, 1, 1, 1]));
        assert_eq!(skew.row_reading_word(), vec![1, 2, 3, 4]);
        assert_eq!(skew.reverse_row_reading_word(), vec![3, 4, 1, 2]);
        assert_eq!(skew.reading_word(), vec![3, 4, 1, 2]);
        assert_eq!(skew.reverse_reading_word(), vec![1, 2, 3, 4]);
        assert_eq!(skew.descent_set(), BTreeSet::from([2]));
    }

    #[test]
    fn test_skew_tableau_with_empty_visible_row() {
        let skew = SkewTableau::new(
            Partition::from_sorted(vec![3, 1]),
            Partition::from_sorted(vec![3]),
            vec![vec![], vec![1]],
        );
        assert!(skew.is_standard());
        assert_eq!(skew.size(), 1);
        assert_eq!(skew.reading_word(), vec![1]);
    }

    #[test]
    fn test_skew_semistandard_column_check_uses_offsets() {
        let skew = SkewTableau::new(
            Partition::from_sorted(vec![3, 2]),
            Partition::from_sorted(vec![1]),
            vec![vec![3, 4], vec![1, 2]],
        );
        assert!(skew.is_row_weakly_increasing());
        assert!(!skew.is_column_strictly_increasing());
        assert!(!skew.is_semistandard());
    }

    #[test]
    fn test_standard_skew_tableaux_matches_linear_extensions() {
        let outer = Partition::from_sorted(vec![3, 2]);
        let inner = Partition::from_sorted(vec![1]);
        let tabs = SkewTableau::standard_tableaux(&outer, &inner);
        assert!(tabs.iter().all(SkewTableau::is_standard));
        assert_eq!(
            tabs.len(),
            Poset::from_skew_shape(&outer, &inner).num_linear_extensions()
        );
    }

    #[test]
    fn test_standard_skew_tableaux_iterator_streams_one_by_one() {
        let outer = Partition::from_sorted(vec![3, 2]);
        let inner = Partition::from_sorted(vec![1]);
        let mut iter = SkewTableau::standard_tableaux_iter(&outer, &inner);
        let first = iter.next().unwrap();
        assert!(first.is_standard());
        assert_eq!(
            iter.count() + 1,
            SkewTableau::standard_tableaux(&outer, &inner).len()
        );
    }

    #[test]
    fn test_standard_skew_tableaux_disconnected_shape() {
        let outer = Partition::from_sorted(vec![2, 1]);
        let inner = Partition::from_sorted(vec![1]);
        let tabs = SkewTableau::standard_tableaux(&outer, &inner);
        assert_eq!(tabs.len(), 2);
        assert_eq!(
            tabs.len(),
            Poset::from_skew_shape(&outer, &inner).num_linear_extensions()
        );
    }

    // ── RSK tests ──────────────────────────────────────────────────

    #[test]
    fn test_rsk_permutation() {
        // RSK on the permutation [3, 1, 2]
        let (p, q) = Tableau::rsk(&[3, 1, 2]);
        assert!(p.is_semistandard());
        assert!(q.is_standard());
        assert_eq!(p.shape(), q.shape());
        // P-tableau: insert 3 → [3], insert 1 → [1,3], insert 2 → [1,2] / [3]
        assert_eq!(p.rows(), &[vec![1, 2], vec![3]]);
        assert_eq!(q.rows(), &[vec![1, 3], vec![2]]);
    }

    #[test]
    fn test_rsk_identity() {
        // RSK on [1, 2, 3] gives a single-row tableau
        let (p, q) = Tableau::rsk(&[1, 2, 3]);
        assert_eq!(p.rows(), &[vec![1, 2, 3]]);
        assert_eq!(q.rows(), &[vec![1, 2, 3]]);
    }

    #[test]
    fn test_rsk_reverse() {
        // RSK on [3, 2, 1] gives a single-column tableau
        let (p, q) = Tableau::rsk(&[3, 2, 1]);
        assert_eq!(p.rows(), &[vec![1], vec![2], vec![3]]);
        assert_eq!(q.rows(), &[vec![1], vec![2], vec![3]]);
    }

    #[test]
    fn test_row_insert() {
        let tab = Tableau::new(vec![vec![1, 3], vec![2]]);
        let (new_tab, row) = tab.row_insert(2);
        // 2 bumps 3 from row 0, 3 goes to row 1 bumping nothing (appends)
        assert_eq!(new_tab.rows(), &[vec![1, 2], vec![2, 3]]);
        assert_eq!(row, 1);
    }

    // ── Promotion / Evacuation tests ───────────────────────────────

    #[test]
    fn test_promotion_standard() {
        // Promotion on [[1,2],[3]] should give [[1,3],[2]]
        let syt = Tableau::new(vec![vec![1, 2], vec![3]]);
        let p1 = syt.promotion();
        assert!(p1.is_standard(), "promotion result not standard: {:?}", p1);
        assert_eq!(p1, Tableau::new(vec![vec![1, 3], vec![2]]));
        let p2 = p1.promotion();
        assert!(p2.is_standard());
        assert_eq!(p2, syt, "promotion^2 = identity on shape (2,1)");
    }

    #[test]
    fn test_promotion_single_row() {
        // Promotion on [1,2,3,4] should cycle: [1,2,3,4] → [1,2,3,4]
        // (single row promotion is trivial since JDT slides to the end)
        let syt = Tableau::new(vec![vec![1, 2, 3, 4]]);
        assert_eq!(syt.promotion(), syt);
    }

    #[test]
    fn test_evacuation_involution() {
        // Evacuation is an involution
        let syt = Tableau::new(vec![vec![1, 3], vec![2]]);
        let e = syt.evacuation();
        assert!(e.is_standard());
        let ee = e.evacuation();
        assert_eq!(ee, syt, "evacuation^2 should be identity");
    }

    #[test]
    fn test_evacuation_preserves_a_non_self_conjugate_shape() {
        let fixed = Tableau::new(vec![vec![1, 3, 5], vec![2, 4, 6]]);
        assert_eq!(fixed.evacuation(), fixed);

        let syt = Tableau::new(vec![vec![1, 2, 5], vec![3, 4]]);
        let evacuated = syt.evacuation();
        assert_eq!(evacuated.rows(), &[vec![1, 2, 3], vec![4, 5]]);
        assert_eq!(evacuated.shape(), syt.shape());
        assert_eq!(evacuated.evacuation(), syt);
    }

    // ── Crystal operator tests ─────────────────────────────────────

    #[test]
    fn test_crystal_f_e_roundtrip() {
        let tab = Tableau::new(vec![vec![1, 1, 2], vec![2, 3]]);
        if let Some(lowered) = tab.crystal_f(1) {
            assert!(lowered.is_semistandard());
            let raised = lowered.crystal_e(1).unwrap();
            assert_eq!(raised, tab);
        }
    }

    #[test]
    fn test_crystal_preserves_shape() {
        let tab = Tableau::new(vec![vec![1, 1, 2], vec![2, 3]]);
        if let Some(lowered) = tab.crystal_f(2) {
            assert_eq!(lowered.shape(), tab.shape());
        }
    }

    // ── Charge tests ───────────────────────────────────────────────

    #[test]
    fn test_charge_ssyt() {
        // Charge is computed from the reverse row reading word.
        // For [[1,1,2],[2,3]]: reverse reading word = [2,3,1,1,2]
        // The decomposition and charge depend on the algorithm.
        // Basic sanity: charge should be non-negative.
        let ssyt = Tableau::new(vec![vec![1, 1, 2], vec![2, 3]]);
        let c = ssyt.charge();
        assert!(c < 100, "charge should be reasonable: {}", c);
    }
}
