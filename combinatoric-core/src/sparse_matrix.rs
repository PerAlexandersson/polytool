//! Sparse matrices with explicit shape, deterministic CSR storage, and a
//! mutable BigInt form for exact chain-complex algorithms.

use num_bigint::BigInt;
use num_traits::Zero;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::ops::{Add, Mul};

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum SparseMatrixError {
    RowOutOfRange { row: usize, rows: usize },
    ColumnOutOfRange { column: usize, columns: usize },
    RaggedRow { row: usize, expected: usize, actual: usize },
    DimensionMismatch { left: (usize, usize), right: (usize, usize) },
    DenseBudgetExceeded { entries: usize, budget: usize },
    AllocationOverflow,
}

impl fmt::Display for SparseMatrixError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::RowOutOfRange { row, rows } => write!(f, "row {row} is outside 0..{rows}"),
            Self::ColumnOutOfRange { column, columns } => {
                write!(f, "column {column} is outside 0..{columns}")
            }
            Self::RaggedRow { row, expected, actual } => write!(
                f,
                "dense row {row} has length {actual}, expected {expected}"
            ),
            Self::DimensionMismatch { left, right } => write!(
                f,
                "incompatible matrix dimensions {left:?} and {right:?}"
            ),
            Self::DenseBudgetExceeded { entries, budget } => write!(
                f,
                "dense conversion needs {entries} entries, above budget {budget}"
            ),
            Self::AllocationOverflow => write!(f, "matrix dimensions overflow usize"),
        }
    }
}

impl std::error::Error for SparseMatrixError {}

/// Immutable compressed-sparse-row matrix.  Empty dimensions are preserved:
/// a 0-by-n matrix is distinct from an n-by-0 matrix.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SparseMatrix<C> {
    rows: usize,
    columns: usize,
    row_offsets: Vec<usize>,
    column_indices: Vec<usize>,
    values: Vec<C>,
}

impl<C> SparseMatrix<C> {
    pub fn rows(&self) -> usize { self.rows }
    pub fn columns(&self) -> usize { self.columns }
    pub fn shape(&self) -> (usize, usize) { (self.rows, self.columns) }
    pub fn nnz(&self) -> usize { self.values.len() }

    pub fn row(&self, row: usize) -> Result<impl Iterator<Item = (usize, &C)>, SparseMatrixError> {
        if row >= self.rows {
            return Err(SparseMatrixError::RowOutOfRange { row, rows: self.rows });
        }
        let start = self.row_offsets[row];
        let end = self.row_offsets[row + 1];
        Ok(self.column_indices[start..end]
            .iter()
            .copied()
            .zip(self.values[start..end].iter()))
    }

    pub fn rows_iter(&self) -> impl Iterator<Item = impl Iterator<Item = (usize, &C)>> {
        (0..self.rows).map(|row| {
            let start = self.row_offsets[row];
            let end = self.row_offsets[row + 1];
            self.column_indices[start..end]
                .iter()
                .copied()
                .zip(self.values[start..end].iter())
        })
    }
}

impl<C> SparseMatrix<C>
where
    C: Clone + Zero + Add<Output = C>,
{
    pub fn builder(rows: usize, columns: usize) -> SparseMatrixBuilder<C> {
        SparseMatrixBuilder::new(rows, columns)
    }

    pub fn zero(rows: usize, columns: usize) -> Self {
        Self {
            rows,
            columns,
            row_offsets: vec![0; rows.saturating_add(1)],
            column_indices: Vec::new(),
            values: Vec::new(),
        }
    }

    pub fn from_dense(dense: &[Vec<C>]) -> Result<Self, SparseMatrixError> {
        let rows = dense.len();
        let columns = dense.first().map_or(0, Vec::len);
        let mut builder = Self::builder(rows, columns);
        for (row, values) in dense.iter().enumerate() {
            if values.len() != columns {
                return Err(SparseMatrixError::RaggedRow {
                    row,
                    expected: columns,
                    actual: values.len(),
                });
            }
            for (column, value) in values.iter().enumerate() {
                builder.add(row, column, value.clone())?;
            }
        }
        Ok(builder.finish())
    }

    pub fn transpose(&self) -> Self {
        let mut builder = Self::builder(self.columns, self.rows);
        for (row, entries) in self.rows_iter().enumerate() {
            for (column, value) in entries {
                // Original indices are valid by the CSR invariant.
                builder.add(column, row, value.clone()).expect("valid transpose index");
            }
        }
        builder.finish()
    }

    pub fn apply(&self, vector: &[C]) -> Result<Vec<C>, SparseMatrixError>
    where
        C: Mul<Output = C>,
    {
        if vector.len() != self.columns {
            return Err(SparseMatrixError::DimensionMismatch {
                left: (self.rows, self.columns),
                right: (vector.len(), 1),
            });
        }
        Ok(self.rows_iter()
            .map(|row| row.fold(C::zero(), |sum, (column, value)| {
                sum + value.clone() * vector[column].clone()
            }))
            .collect())
    }

    pub fn to_dense_with_budget(&self, budget: usize) -> Result<Vec<Vec<C>>, SparseMatrixError> {
        let entries = self.rows.checked_mul(self.columns).ok_or(SparseMatrixError::AllocationOverflow)?;
        if entries > budget {
            return Err(SparseMatrixError::DenseBudgetExceeded { entries, budget });
        }
        let mut dense = vec![vec![C::zero(); self.columns]; self.rows];
        for (row, entries) in self.rows_iter().enumerate() {
            for (column, value) in entries {
                dense[row][column] = value.clone();
            }
        }
        Ok(dense)
    }

    /// Test `left * right = 0` without allocating a product.  The witness is
    /// the first deterministic nonzero entry found.
    pub fn compose_is_zero(
        left: &Self,
        right: &Self,
    ) -> Result<Result<(), (usize, usize, C)>, SparseMatrixError>
    where
        C: Mul<Output = C>,
    {
        if left.columns != right.rows {
            return Err(SparseMatrixError::DimensionMismatch {
                left: left.shape(),
                right: right.shape(),
            });
        }
        let left_transpose = left.transpose();
        let right_transpose = right.transpose();
        for column in 0..right.columns {
            let mut accum: BTreeMap<usize, C> = BTreeMap::new();
            for (middle, coefficient) in right_transpose.row(column)? {
                for (row, left_coefficient) in left_transpose.row(middle)? {
                    let old = accum.remove(&row).unwrap_or_else(C::zero);
                    let value = old + left_coefficient.clone() * coefficient.clone();
                    if !value.is_zero() {
                        accum.insert(row, value);
                    }
                }
            }
            if let Some((row, value)) = accum.into_iter().next() {
                return Ok(Err((row, column, value)));
            }
        }
        Ok(Ok(()))
    }

}

/// Deterministic triplet builder.  Duplicate coefficients are summed and
/// exact zeroes disappear on `finish`.
#[derive(Clone, Debug)]
pub struct SparseMatrixBuilder<C> {
    rows: usize,
    columns: usize,
    entries: BTreeMap<(usize, usize), C>,
}

impl<C> SparseMatrixBuilder<C>
where
    C: Clone + Zero + Add<Output = C>,
{
    pub fn new(rows: usize, columns: usize) -> Self {
        Self { rows, columns, entries: BTreeMap::new() }
    }

    pub fn add(&mut self, row: usize, column: usize, value: C) -> Result<(), SparseMatrixError> {
        if row >= self.rows {
            return Err(SparseMatrixError::RowOutOfRange { row, rows: self.rows });
        }
        if column >= self.columns {
            return Err(SparseMatrixError::ColumnOutOfRange { column, columns: self.columns });
        }
        let old = self.entries.remove(&(row, column)).unwrap_or_else(C::zero);
        let sum = old + value;
        if !sum.is_zero() {
            self.entries.insert((row, column), sum);
        }
        Ok(())
    }

    pub fn finish(self) -> SparseMatrix<C> {
        let mut row_offsets = vec![0; self.rows + 1];
        for &(row, _) in self.entries.keys() {
            row_offsets[row + 1] += 1;
        }
        for row in 0..self.rows {
            row_offsets[row + 1] += row_offsets[row];
        }
        let mut column_indices = Vec::with_capacity(self.entries.len());
        let mut values = Vec::with_capacity(self.entries.len());
        for ((_, column), value) in self.entries {
            column_indices.push(column);
            values.push(value);
        }
        SparseMatrix { rows: self.rows, columns: self.columns, row_offsets, column_indices, values }
    }
}

/// Mutable sparse BigInt matrix with row maps and reverse column incidence.
/// It is designed for bounded exact cancellation, not as a universal sparse
/// numerical backend.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MutableSparseMatrix {
    rows: usize,
    columns: usize,
    row_maps: Vec<BTreeMap<usize, BigInt>>,
    column_rows: Vec<BTreeSet<usize>>,
    nnz: usize,
}

impl MutableSparseMatrix {
    pub fn new(rows: usize, columns: usize) -> Self {
        Self {
            rows,
            columns,
            row_maps: vec![BTreeMap::new(); rows],
            column_rows: vec![BTreeSet::new(); columns],
            nnz: 0,
        }
    }

    pub fn from_sparse(matrix: &SparseMatrix<BigInt>) -> Self {
        let mut result = Self::new(matrix.rows, matrix.columns);
        for (row, entries) in matrix.rows_iter().enumerate() {
            for (column, value) in entries {
                result.set(row, column, value.clone()).expect("valid CSR index");
            }
        }
        result
    }

    pub fn rows(&self) -> usize { self.rows }
    pub fn columns(&self) -> usize { self.columns }
    pub fn nnz(&self) -> usize { self.nnz }

    pub fn get(&self, row: usize, column: usize) -> Result<BigInt, SparseMatrixError> {
        self.check_index(row, column)?;
        Ok(self.row_maps[row].get(&column).cloned().unwrap_or_else(BigInt::zero))
    }

    pub fn row_entries(&self, row: usize) -> Result<impl Iterator<Item = (usize, &BigInt)>, SparseMatrixError> {
        if row >= self.rows {
            return Err(SparseMatrixError::RowOutOfRange { row, rows: self.rows });
        }
        Ok(self.row_maps[row].iter().map(|(&column, value)| (column, value)))
    }

    pub fn column_rows(&self, column: usize) -> Result<impl Iterator<Item = usize> + '_, SparseMatrixError> {
        if column >= self.columns {
            return Err(SparseMatrixError::ColumnOutOfRange { column, columns: self.columns });
        }
        Ok(self.column_rows[column].iter().copied())
    }

    pub fn set(&mut self, row: usize, column: usize, value: BigInt) -> Result<(), SparseMatrixError> {
        self.check_index(row, column)?;
        let existed = self.row_maps[row].contains_key(&column);
        if value.is_zero() {
            if existed {
                self.row_maps[row].remove(&column);
                self.column_rows[column].remove(&row);
                self.nnz -= 1;
            }
        } else if existed {
            self.row_maps[row].insert(column, value);
        } else {
            self.row_maps[row].insert(column, value);
            self.column_rows[column].insert(row);
            self.nnz += 1;
        }
        Ok(())
    }

    pub fn add_to(&mut self, row: usize, column: usize, value: &BigInt) -> Result<(), SparseMatrixError> {
        let next = self.get(row, column)? + value;
        self.set(row, column, next)
    }

    pub fn clear_row(&mut self, row: usize) -> Result<(), SparseMatrixError> {
        let columns = self.row_entries(row)?.map(|(column, _)| column).collect::<Vec<_>>();
        for column in columns { self.set(row, column, BigInt::zero())?; }
        Ok(())
    }

    pub fn clear_column(&mut self, column: usize) -> Result<(), SparseMatrixError> {
        let rows = self.column_rows(column)?.collect::<Vec<_>>();
        for row in rows { self.set(row, column, BigInt::zero())?; }
        Ok(())
    }

    pub fn to_sparse(&self) -> SparseMatrix<BigInt> {
        let mut builder = SparseMatrix::builder(self.rows, self.columns);
        for (row, entries) in self.row_maps.iter().enumerate() {
            for (&column, value) in entries {
                builder.add(row, column, value.clone()).expect("mutable indices are valid");
            }
        }
        builder.finish()
    }

    pub fn validate(&self) -> Result<(), String> {
        let mut counted = 0;
        for (row, entries) in self.row_maps.iter().enumerate() {
            for (&column, value) in entries {
                if value.is_zero() || column >= self.columns || !self.column_rows[column].contains(&row) {
                    return Err(format!("invalid sparse edge ({row}, {column})"));
                }
                counted += 1;
            }
        }
        if counted != self.nnz { return Err("nnz counter is inconsistent".into()); }
        for (column, rows) in self.column_rows.iter().enumerate() {
            for &row in rows {
                if row >= self.rows || !self.row_maps[row].contains_key(&column) {
                    return Err(format!("invalid reverse sparse edge ({row}, {column})"));
                }
            }
        }
        Ok(())
    }

    fn check_index(&self, row: usize, column: usize) -> Result<(), SparseMatrixError> {
        if row >= self.rows { return Err(SparseMatrixError::RowOutOfRange { row, rows: self.rows }); }
        if column >= self.columns { return Err(SparseMatrixError::ColumnOutOfRange { column, columns: self.columns }); }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn csr_aggregates_duplicates_and_preserves_empty_shape() {
        let builder = SparseMatrix::<BigInt>::builder(0, 3);
        assert_eq!(builder.finish().shape(), (0, 3));
        let mut builder = SparseMatrix::<BigInt>::builder(2, 3);
        builder.add(1, 2, 4.into()).unwrap();
        builder.add(1, 2, (-4).into()).unwrap();
        assert_eq!(builder.finish().nnz(), 0);
    }

    #[test]
    fn transpose_apply_and_composition_are_exact() {
        let matrix = SparseMatrix::<BigInt>::from_dense(&vec![vec![1.into(), 2.into()], vec![0.into(), 3.into()]]).unwrap();
        assert_eq!(matrix.apply(&[4.into(), 5.into()]).unwrap(), vec![14.into(), 15.into()]);
        assert_eq!(matrix.transpose().transpose(), matrix);
        let zero = SparseMatrix::zero(2, 2);
        assert_eq!(SparseMatrix::compose_is_zero(&matrix, &zero).unwrap(), Ok(()));
    }

    #[test]
    fn mutable_form_keeps_reverse_incidence_and_zeroes() {
        let mut matrix = MutableSparseMatrix::new(3, 2);
        matrix.set(2, 1, 7.into()).unwrap();
        matrix.add_to(2, 1, &(-7).into()).unwrap();
        assert_eq!(matrix.nnz(), 0);
        matrix.set(0, 0, 1.into()).unwrap();
        matrix.set(2, 0, 2.into()).unwrap();
        matrix.clear_column(0).unwrap();
        assert_eq!(matrix.nnz(), 0);
        matrix.validate().unwrap();
    }
}
