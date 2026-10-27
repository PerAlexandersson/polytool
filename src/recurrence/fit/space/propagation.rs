//! Polynomial jet-state certificates for monic polynomial-index sources.
//!
//! We retain exact row-operation multipliers. Replay checks polynomial
//! identities, without searching for a dependence or dividing at integer ranks.
use super::*;
use crate::recurrence::RecurrenceError;
use num_traits::ToPrimitive;
type PolynomialMatrix = Vec<Vec<BivarPoly>>;

/// Monic source `P[n+k] = sum_i coefficients[i](n,X) P[n+i]`.
/// Each coefficient is indexed by [n power][X power], with exact strings.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct IndexedLagSource {
    /// `coefficients[i]` is `a_i(n, X)` as an array indexed `[n power][X power]` of
    /// exact rational strings; the source order `k` is `coefficients.len()`.
    pub coefficients: Vec<ExactRows>,
    /// Initial rows: `seeds[j]` is `P[origin + j]` as ascending-power rational strings.
    /// Exactly `k` seeds are required.
    pub seeds: ExactRows,
    /// Original index of `seeds[0]`; the rule is used for base indices `n >= origin`.
    /// Defaults to zero.
    #[serde(default)]
    pub origin: usize,
    /// Free-form note on where the source came from. It is informational and is never
    /// used in a computation.
    #[serde(default)]
    pub provenance: Option<String>,
}
impl IndexedLagSource {
    /// Canonical exact generator, in the original index convention.
    pub fn generate_rows(
        &self,
        count: usize,
        limits: &RecurrenceSpaceLimits,
    ) -> Result<ExactRows, RecurrenceError> {
        let request = RecurrenceSpaceRequest {
            source: None,
            source_recurrence: None,
            indexed_source: Some(self.clone()),
            coefficient_rescaling: None,
            include_leader_coverage: false,
            rows: None,
            row_origin: self.origin,
            domain_start: self.origin,
            shapes: vec![],
            limits: limits.clone(),
        };
        let input = prepare(&request)?;
        let mut rows = source_rows(&input, count, limits)?;
        rows.truncate(count);
        Ok(encode_polys(&rows))
    }
}

/// Polynomial coefficients of a common residual annihilator.
/// `sum_j coefficients[j](n,X) E[n+j] = 0` for every shaped residual.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RecurrencePropagationCertificate {
    /// `coefficients[j]` is the polynomial multiplying `E[n + j]`, as an array indexed
    /// `[n power][X power]`. The last entry is the leader.
    pub coefficients: Vec<ExactRows>,
    /// Dimension of the jet state: the source order times one more than the largest
    /// derivative order in the support.
    pub state_dimension: usize,
    /// Number of distinct `(lag, derivative_order, index_degree)` triples in the support.
    pub basic_columns: usize,
    /// Nonnegative integer ranks, as decimal strings, at which the leader vanishes
    /// identically in `X`. `None` means the bounded exact analysis was not performed.
    pub exceptional_ranks: Option<Vec<String>>,
    /// How `regular_from` was derived: `exact_integer_zeros` or `coefficient_root_bound`.
    pub cutoff_method: String,
    /// Base index from which the leader does not vanish identically in `X` at any rank:
    /// at least the domain start and beyond every exceptional rank (or, when those are
    /// unknown, beyond a Cauchy bound on the roots).
    pub regular_from: usize,
}

pub(super) struct Arithmetic<'a> {
    pub(super) limits: &'a RecurrenceSpaceLimits,
    remaining: usize,
}
impl<'a> Arithmetic<'a> {
    pub(super) fn new(limits: &'a RecurrenceSpaceLimits) -> Self {
        Self {
            limits,
            remaining: limits.max_operations,
        }
    }
    fn spend(&mut self, n: usize) -> Result<(), RecurrenceError> {
        self.remaining = self
            .remaining
            .checked_sub(n)
            .ok_or(RecurrenceError::BudgetExhausted(
                "polynomial propagation arithmetic budget".into(),
            ))?;
        Ok(())
    }
    fn check(&self, p: &BivarPoly) -> Result<(), RecurrenceError> {
        if p.coeffs.len().saturating_sub(1) > self.limits.max_degree
            || p.coeffs
                .iter()
                .any(|r| r.len().saturating_sub(1) > self.limits.max_degree)
        {
            return Err(RecurrenceError::BudgetExhausted(
                "polynomial propagation degree budget".into(),
            ));
        }
        if p.coeffs.iter().map(Vec::len).sum::<usize>() > self.limits.max_matrix_entries {
            return Err(RecurrenceError::BudgetExhausted(
                "polynomial propagation entry budget".into(),
            ));
        }
        if p.coeffs
            .iter()
            .flatten()
            .any(|q| q.numer().bits().max(q.denom().bits()) > self.limits.max_coefficient_bits)
        {
            return Err(RecurrenceError::BudgetExhausted(
                "polynomial propagation coefficient budget".into(),
            ));
        }
        Ok(())
    }
    pub(super) fn add(
        &mut self,
        a: &BivarPoly,
        b: &BivarPoly,
        minus: bool,
    ) -> Result<BivarPoly, RecurrenceError> {
        let n = a.coeffs.len().max(b.coeffs.len());
        let x = a
            .coeffs
            .iter()
            .chain(&b.coeffs)
            .map(Vec::len)
            .max()
            .unwrap_or(0);
        if n.checked_mul(x)
            .is_none_or(|v| v > self.limits.max_matrix_entries)
        {
            return Err(RecurrenceError::BudgetExhausted(
                "polynomial propagation entry budget".into(),
            ));
        }
        self.spend(n.checked_mul(x).ok_or(RecurrenceError::BudgetExhausted(
            "polynomial propagation entry budget".into(),
        ))?)?;
        let mut out = vec![vec![Q::zero(); x]; n];
        for (i, row) in a.coeffs.iter().enumerate() {
            for (j, c) in row.iter().enumerate() {
                out[i][j] += c;
            }
        }
        for (i, row) in b.coeffs.iter().enumerate() {
            for (j, c) in row.iter().enumerate() {
                if minus {
                    out[i][j] -= c;
                } else {
                    out[i][j] += c;
                }
            }
        }
        let p = trim(out);
        self.check(&p)?;
        Ok(p)
    }
    pub(super) fn mul(
        &mut self,
        a: &BivarPoly,
        b: &BivarPoly,
    ) -> Result<BivarPoly, RecurrenceError> {
        if a.is_zero() || b.is_zero() {
            return Ok(zero());
        }
        let na = a.coeffs.len();
        let nb = b.coeffs.len();
        let xa = a.coeffs.iter().map(Vec::len).max().unwrap_or(0);
        let xb = b.coeffs.iter().map(Vec::len).max().unwrap_or(0);
        let n = na.checked_add(nb).and_then(|v| v.checked_sub(1)).ok_or(
            RecurrenceError::BudgetExhausted("polynomial propagation degree budget".into()),
        )?;
        let x = xa.checked_add(xb).and_then(|v| v.checked_sub(1)).ok_or(
            RecurrenceError::BudgetExhausted("polynomial propagation degree budget".into()),
        )?;
        if n.saturating_sub(1) > self.limits.max_degree
            || x.saturating_sub(1) > self.limits.max_degree
            || n.checked_mul(x)
                .is_none_or(|v| v > self.limits.max_matrix_entries)
        {
            return Err(RecurrenceError::BudgetExhausted(
                "polynomial propagation degree/entry budget".into(),
            ));
        }
        let work = na
            .checked_mul(xa)
            .and_then(|v| v.checked_mul(nb))
            .and_then(|v| v.checked_mul(xb))
            .ok_or(RecurrenceError::BudgetExhausted(
                "polynomial propagation arithmetic budget".into(),
            ))?;
        self.spend(work)?;
        let mut out = vec![vec![Q::zero(); x]; n];
        for (i, ar) in a.coeffs.iter().enumerate() {
            for (j, ac) in ar.iter().enumerate().filter(|(_, c)| !c.is_zero()) {
                for (l, br) in b.coeffs.iter().enumerate() {
                    for (m, bc) in br.iter().enumerate().filter(|(_, c)| !c.is_zero()) {
                        out[i + l][j + m] += ac * bc;
                    }
                }
            }
        }
        let p = trim(out);
        self.check(&p)?;
        Ok(p)
    }
    fn shift(&mut self, p: &BivarPoly, amount: usize) -> Result<BivarPoly, RecurrenceError> {
        self.spend(
            p.coeffs
                .len()
                .checked_mul(p.coeffs.iter().map(Vec::len).sum())
                .and_then(|v| v.checked_mul(2))
                .ok_or(RecurrenceError::BudgetExhausted(
                    "polynomial propagation arithmetic budget".into(),
                ))?,
        )?;
        let q = trim(
            p.shift_index(
                i64::try_from(amount)
                    .map_err(|_| RecurrenceError::LimitExceeded("shift overflow".into()))?,
            )
            .coeffs,
        );
        self.check(&q)?;
        Ok(q)
    }
    fn derivative(&mut self, p: &BivarPoly, d: usize) -> Result<BivarPoly, RecurrenceError> {
        self.spend(
            p.coeffs
                .iter()
                .map(Vec::len)
                .sum::<usize>()
                .checked_mul(d.saturating_add(1))
                .ok_or(RecurrenceError::BudgetExhausted(
                    "polynomial propagation arithmetic budget".into(),
                ))?,
        )?;
        let q = trim(
            p.coeffs
                .iter()
                .map(|r| Poly::new(r.clone()).nth_derivative(d).coeffs().to_vec())
                .collect(),
        );
        self.check(&q)?;
        Ok(q)
    }
}
fn trim(mut rows: Vec<Vec<Q>>) -> BivarPoly {
    for row in &mut rows {
        while row.last().is_some_and(Zero::is_zero) {
            row.pop();
        }
    }
    while rows.last().is_some_and(Vec::is_empty) {
        rows.pop();
    }
    BivarPoly { coeffs: rows }
}
fn zero() -> BivarPoly {
    BivarPoly { coeffs: vec![] }
}
fn one() -> BivarPoly {
    BivarPoly {
        coeffs: vec![vec![Q::one()]],
    }
}
fn encode(p: &BivarPoly) -> ExactRows {
    encode_matrix(&p.coeffs)
}
pub(super) fn parse(
    p: &ExactRows,
    lim: &RecurrenceSpaceLimits,
) -> Result<BivarPoly, RecurrenceError> {
    if p.len().saturating_sub(1) > lim.max_degree
        || p.iter().any(|r| r.len().saturating_sub(1) > lim.max_degree)
    {
        return Err(RecurrenceError::BudgetExhausted(
            "source coefficient degree budget".into(),
        ));
    }
    Ok(trim(parse_matrix(p, lim)?))
}
pub(super) fn evaluate(
    p: &BivarPoly,
    n: usize,
    lim: &RecurrenceSpaceLimits,
) -> Result<Poly, RecurrenceError> {
    let rank = Q::from_integer(BigInt::from(n));
    let mut out = Poly::zero();
    for row in p.coeffs.iter().rev() {
        out = out.scale(&rank) + Poly::new(row.clone());
        poly_budget(&out, lim)?;
    }
    Ok(out)
}
fn matrix_budget(
    rows: &[Vec<BivarPoly>],
    lim: &RecurrenceSpaceLimits,
) -> Result<(), RecurrenceError> {
    let mut size = 0usize;
    for p in rows.iter().flatten() {
        size = size
            .checked_add(1 + p.coeffs.iter().map(Vec::len).sum::<usize>())
            .ok_or(RecurrenceError::BudgetExhausted(
                "propagation matrix entry budget".into(),
            ))?;
    }
    if size > lim.max_matrix_entries {
        return Err(RecurrenceError::BudgetExhausted(
            "propagation matrix entry budget".into(),
        ));
    }
    Ok(())
}
// Check aggregate retained coefficients before storing each new polynomial.
// A single arithmetic temporary is separately bounded by Arithmetic::check.
fn replace_entry(
    entry: &mut BivarPoly,
    value: BivarPoly,
    retained: &mut usize,
    lim: &RecurrenceSpaceLimits,
) -> Result<(), RecurrenceError> {
    let old = entry.coeffs.iter().map(Vec::len).sum::<usize>();
    let new = value.coeffs.iter().map(Vec::len).sum::<usize>();
    let size = retained
        .checked_sub(old)
        .and_then(|n| n.checked_add(new))
        .ok_or(RecurrenceError::BudgetExhausted(
            "propagation matrix entry budget".into(),
        ))?;
    if size > lim.max_matrix_entries {
        return Err(RecurrenceError::BudgetExhausted(
            "propagation matrix entry budget".into(),
        ));
    }
    *retained = size;
    *entry = value;
    Ok(())
}
fn transition(
    a: &[BivarPoly],
    s: usize,
    ar: &mut Arithmetic<'_>,
) -> Result<Vec<Vec<BivarPoly>>, RecurrenceError> {
    let k = a.len();
    let q = k
        .checked_mul(s.checked_add(1).ok_or(RecurrenceError::BudgetExhausted(
            "jet dimension budget".into(),
        ))?)
        .ok_or(RecurrenceError::BudgetExhausted(
            "jet dimension budget".into(),
        ))?;
    if q > ar.limits.max_columns
        || q.checked_mul(q)
            .is_none_or(|v| v > ar.limits.max_matrix_entries)
    {
        return Err(RecurrenceError::BudgetExhausted(
            "jet dimension budget".into(),
        ));
    }
    let mut t = vec![vec![zero(); q]; q];
    let mut retained = q * q;
    for d in 0..=s {
        for i in 0..k - 1 {
            replace_entry(
                &mut t[d * k + i][d * k + i + 1],
                one(),
                &mut retained,
                ar.limits,
            )?;
        }
        let mut binomial = BigInt::one();
        for e in 0..=d {
            let factor = BivarPoly {
                coeffs: vec![vec![Q::from_integer(binomial.clone())]],
            };
            for (i, c) in a.iter().enumerate() {
                let derivative = ar.derivative(c, d - e)?;
                let value = ar.mul(&factor, &derivative)?;
                replace_entry(
                    &mut t[d * k + k - 1][e * k + i],
                    value,
                    &mut retained,
                    ar.limits,
                )?;
            }
            if e < d {
                binomial = binomial * BigInt::from(d - e) / BigInt::from(e + 1);
            }
        }
    }
    matrix_budget(&t, ar.limits)?;
    Ok(t)
}
fn advance(
    rows: &[Vec<BivarPoly>],
    t: &[Vec<BivarPoly>],
    ar: &mut Arithmetic<'_>,
) -> Result<Vec<Vec<BivarPoly>>, RecurrenceError> {
    let q = t.len();
    let mut out = vec![vec![zero(); q]; rows.len()];
    let mut retained = q * rows.len();
    for (i, row) in rows.iter().enumerate() {
        for (l, p) in row.iter().enumerate().filter(|(_, p)| !p.is_zero()) {
            let shifted = ar.shift(p, 1)?;
            for (j, c) in t[l].iter().enumerate().filter(|(_, c)| !c.is_zero()) {
                let product = ar.mul(&shifted, c)?;
                let value = ar.add(&out[i][j], &product, false)?;
                replace_entry(&mut out[i][j], value, &mut retained, ar.limits)?;
            }
        }
    }
    matrix_budget(&out, ar.limits)?;
    Ok(out)
}
fn setup(
    a: &[BivarPoly],
    cols: &[RecurrenceMonomial],
    ar: &mut Arithmetic<'_>,
) -> Result<(PolynomialMatrix, PolynomialMatrix), RecurrenceError> {
    let basic: BTreeSet<_> = cols
        .iter()
        .map(|c| (c.lag, c.derivative_order, c.index_degree))
        .collect();
    let s = basic.iter().map(|c| c.1).max().unwrap_or(0);
    let t = transition(a, s, ar)?;
    let k = a.len();
    let q = t.len();
    if basic
        .len()
        .checked_mul(q)
        .and_then(|v| v.checked_mul(v))
        .is_none_or(|v| v > ar.limits.max_matrix_entries / 2)
    {
        return Err(RecurrenceError::BudgetExhausted(
            "propagation dependence dimension budget".into(),
        ));
    }
    let mut v = Vec::new();
    for (lag, d, u) in basic {
        if lag > ar.limits.max_rows {
            return Err(RecurrenceError::BudgetExhausted(
                "propagation lag budget".into(),
            ));
        }
        let mut row = vec![zero(); q];
        row[d * k] = one();
        for _ in 0..lag {
            row = advance(&[row], &t, ar)?.remove(0);
        }
        let mut power = vec![vec![]; u + 1];
        power[u] = vec![Q::one()];
        let factor = BivarPoly { coeffs: power };
        for p in &mut row {
            *p = ar.mul(p, &factor)?;
        }
        v.push(row);
        matrix_budget(&v, ar.limits)?;
    }
    matrix_budget(&v, ar.limits)?;
    Ok((t, v))
}
fn finish(
    coefficients: Vec<BivarPoly>,
    q: usize,
    m: usize,
    start: usize,
) -> Result<RecurrencePropagationCertificate, RecurrenceError> {
    let leading = coefficients.last().ok_or(RecurrenceError::InvalidInput(
        "empty propagation relation".into(),
    ))?;
    if leading.is_zero() {
        return Err(RecurrenceError::InvalidInput(
            "zero propagation leader".into(),
        ));
    }
    let width = leading.coeffs.iter().map(Vec::len).max().unwrap_or(0);
    let polynomials: Vec<_> = (0..width)
        .map(|x| {
            Poly::new(
                leading
                    .coeffs
                    .iter()
                    .map(|row| row.get(x).cloned().unwrap_or_else(Q::zero))
                    .collect(),
            )
        })
        .collect();
    let roots = integer_zeros_of_gcd(&polynomials);
    let mut regular = start;
    let cutoff_method = if let Some(roots) = &roots {
        for root in roots {
            let n = root
                .parse::<BigInt>()
                .map_err(|_| RecurrenceError::InvalidInput("invalid exceptional rank".into()))?
                .to_usize()
                .ok_or(RecurrenceError::BudgetExhausted(
                    "exceptional rank row budget".into(),
                ))?;
            regular = regular.max(n.checked_add(1).ok_or(RecurrenceError::BudgetExhausted(
                "exceptional rank row budget".into(),
            ))?);
        }
        "exact_integer_zeros"
    } else {
        // A zero of the whole leader is a zero of every nonzero X-coefficient.
        // Cauchy's absolute root bound on any one coefficient is sufficient.
        let mut best: Option<BigInt> = None;
        for p in polynomials.iter().filter(|p| !p.is_zero()) {
            let lc = p.leading_coefficient().unwrap().abs();
            let ratio = p
                .coeffs()
                .iter()
                .take(p.coeffs().len() - 1)
                .map(|a| a.abs() / &lc)
                .max()
                .unwrap_or_else(Q::zero);
            let bound = ratio.numer().div_floor(ratio.denom()) + BigInt::from(2);
            best = Some(best.map_or(bound.clone(), |b| b.min(bound)));
        }
        regular = regular.max(
            best.ok_or(RecurrenceError::InvalidInput(
                "zero propagation leader".into(),
            ))?
            .to_usize()
            .ok_or(RecurrenceError::BudgetExhausted(
                "exceptional rank row budget".into(),
            ))?,
        );
        "coefficient_root_bound"
    };
    Ok(RecurrencePropagationCertificate {
        coefficients: coefficients.iter().map(encode).collect(),
        state_dimension: q,
        basic_columns: m,
        exceptional_ranks: roots,
        cutoff_method: cutoff_method.into(),
        regular_from: regular,
    })
}

pub(super) fn build(
    a: &[BivarPoly],
    cols: &[RecurrenceMonomial],
    start: usize,
    lim: &RecurrenceSpaceLimits,
) -> Result<RecurrencePropagationCertificate, RecurrenceError> {
    let mut ar = Arithmetic::new(lim);
    let (t, mut w) = setup(a, cols, &mut ar)?;
    let q = t.len();
    let m = w.len();
    let dimension = q * m;
    let mut basis: Vec<(usize, Vec<BivarPoly>, Vec<BivarPoly>)> = Vec::new();
    for order in 0..=dimension {
        let mut row: Vec<_> = w.iter().flatten().cloned().collect();
        let mut relation = vec![zero(); order + 1];
        relation[order] = one();
        let mut retained = row.len()
            + relation.len()
            + row
                .iter()
                .chain(&relation)
                .map(|p| p.coeffs.iter().map(Vec::len).sum::<usize>())
                .sum::<usize>();
        for (pivot, previous, combination) in &basis {
            if row[*pivot].is_zero() {
                continue;
            }
            let factor = row[*pivot].clone();
            let divisor = previous[*pivot].clone();
            for (a, b) in row.iter_mut().zip(previous) {
                let first = ar.mul(&divisor, a)?;
                let second = ar.mul(&factor, b)?;
                let value = ar.add(&first, &second, true)?;
                replace_entry(a, value, &mut retained, lim)?;
            }
            for (i, a) in relation.iter_mut().enumerate() {
                let first = ar.mul(&divisor, a)?;
                let value = if let Some(b) = combination.get(i) {
                    let second = ar.mul(&factor, b)?;
                    ar.add(&first, &second, true)?
                } else {
                    first
                };
                replace_entry(a, value, &mut retained, lim)?;
            }
        }
        let Some(pivot) = row.iter().position(|p| !p.is_zero()) else {
            return finish(relation, q, m, start);
        };
        basis.push((pivot, row, relation));
        let mut entries = 0usize;
        for (_, r, c) in &basis {
            for p in r.iter().chain(c) {
                entries = entries
                    .checked_add(1 + p.coeffs.iter().map(Vec::len).sum::<usize>())
                    .ok_or(RecurrenceError::BudgetExhausted(
                        "propagation basis entry budget".into(),
                    ))?;
            }
        }
        if entries > lim.max_matrix_entries {
            return Err(RecurrenceError::BudgetExhausted(
                "propagation basis entry budget".into(),
            ));
        }
        w = advance(&w, &t, &mut ar)?;
    }
    Err(RecurrenceError::InvalidInput(
        "polynomial dependence dimension invariant failed".into(),
    ))
}

pub(super) fn verify(
    a: &[BivarPoly],
    cols: &[RecurrenceMonomial],
    start: usize,
    cert: &RecurrencePropagationCertificate,
    lim: &RecurrenceSpaceLimits,
) -> Result<(), RecurrenceError> {
    let mut ar = Arithmetic::new(lim);
    let (t, mut w) = setup(a, cols, &mut ar)?;
    let q = t.len();
    let m = w.len();
    if cert.coefficients.is_empty() || cert.coefficients.len() - 1 > q * m {
        return Err(RecurrenceError::VerificationFailed(
            "invalid propagation order".into(),
        ));
    }
    let coeffs = cert
        .coefficients
        .iter()
        .map(|p| parse(p, lim))
        .collect::<Result<Vec<_>, _>>()?;
    let mut residual = vec![vec![zero(); q]; m];
    let mut retained = q * m;
    for (j, c) in coeffs.iter().enumerate() {
        for (r, wr) in residual.iter_mut().zip(&w) {
            for (entry, p) in r.iter_mut().zip(wr) {
                let term = ar.mul(c, p)?;
                let value = ar.add(entry, &term, false)?;
                replace_entry(entry, value, &mut retained, lim)?;
            }
        }
        if j + 1 < coeffs.len() {
            w = advance(&w, &t, &mut ar)?;
        }
    }
    if residual.iter().flatten().any(|p| !p.is_zero()) {
        return Err(RecurrenceError::VerificationFailed(
            "invalid polynomial propagation identity".into(),
        ));
    }
    if &finish(coeffs, q, m, start)? != cert {
        return Err(RecurrenceError::VerificationFailed(
            "propagation cutoff or metadata mismatch".into(),
        ));
    }
    Ok(())
}
