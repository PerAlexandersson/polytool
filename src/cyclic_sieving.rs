//! Exact root-of-unity checks for cyclic sieving.
//!
//! The implementation avoids floating point arithmetic.  For a primitive
//! `r`-th root of unity, `P(zeta) = m` if and only if `P(q) - m` is divisible by
//! the cyclotomic polynomial `Phi_r(q)`.
//!
//! The module also decomposes an integral polynomial into cyclic orbits
//! ([`cyclic_sieving_orbit_decomposition_bigint`]), applies the
//! Alexandersson--Amini existence criterion for a cyclic action, converts
//! between orbit sizes and fixed-point counts, and checks bicyclic sieving
//! ([`bicyclic_sieving_report_bigint`]).  Reference: P. Alexandersson and
//! N. Amini, *The cone of cyclic sieving phenomena*, Discrete Math. 342 (2019)
//! 1581-1601, arXiv:1804.01447.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Zero};
use std::collections::BTreeMap;

/// Exact evaluation data for one root of unity.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RootOfUnityEvaluation {
    pub group_order: usize,
    pub power: usize,
    pub root_order: usize,
    pub integer_value: Option<BigInt>,
    pub remainder: Vec<BigInt>,
}

/// One CSP equality check against a supplied fixed-point count.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CyclicSievingPowerCheck {
    pub power: usize,
    pub root_order: usize,
    pub expected_fixed_points: BigInt,
    pub evaluation: RootOfUnityEvaluation,
    pub holds: bool,
}

/// CSP report for one polynomial and one cyclic group order.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CyclicSievingReport {
    pub order: usize,
    pub coefficients: Vec<BigInt>,
    pub fixed_counts: Option<Vec<BigInt>>,
    pub evaluations: Vec<RootOfUnityEvaluation>,
    pub checks: Vec<CyclicSievingPowerCheck>,
    pub holds: Option<bool>,
}

/// Sequence-level CSP report for one row and several candidate group orders.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CyclicSievingSequenceItem {
    pub row: usize,
    pub index: isize,
    pub coefficients: Vec<BigInt>,
    pub candidate_orders: Vec<CyclicSievingReport>,
}

fn trim(mut p: Vec<BigInt>) -> Vec<BigInt> {
    while p.last().is_some_and(|c| c.is_zero()) {
        p.pop();
    }
    p
}

fn degree(p: &[BigInt]) -> Option<usize> {
    p.iter().rposition(|c| !c.is_zero())
}

fn divisors(n: usize) -> Vec<usize> {
    let mut out = Vec::new();
    for d in 1..=n {
        if n % d == 0 {
            out.push(d);
        }
    }
    out
}

#[cfg(test)]
fn poly_mul(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![BigInt::zero(); a.len() + b.len() - 1];
    for (i, ca) in a.iter().enumerate() {
        if ca.is_zero() {
            continue;
        }
        for (j, cb) in b.iter().enumerate() {
            if !cb.is_zero() {
                out[i + j] += ca * cb;
            }
        }
    }
    trim(out)
}

fn exact_div_monic(mut numerator: Vec<BigInt>, divisor: &[BigInt]) -> Option<Vec<BigInt>> {
    numerator = trim(numerator);
    let divisor = trim(divisor.to_vec());
    let dd = degree(&divisor)?;
    if divisor[dd] != BigInt::one() {
        return None;
    }
    let nd = match degree(&numerator) {
        Some(nd) => nd,
        None => return Some(Vec::new()),
    };
    if nd < dd {
        return Some(Vec::new());
    }

    let mut quotient = vec![BigInt::zero(); nd - dd + 1];
    for shift in (0..=nd - dd).rev() {
        let coeff = numerator[dd + shift].clone();
        if coeff.is_zero() {
            continue;
        }
        quotient[shift] = coeff.clone();
        for (j, dcoeff) in divisor.iter().enumerate().take(dd + 1) {
            numerator[j + shift] -= &coeff * dcoeff;
        }
    }
    if trim(numerator).is_empty() {
        Some(trim(quotient))
    } else {
        None
    }
}

fn remainder_monic(mut numerator: Vec<BigInt>, divisor: &[BigInt]) -> Vec<BigInt> {
    numerator = trim(numerator);
    let divisor = trim(divisor.to_vec());
    let Some(dd) = degree(&divisor) else {
        panic!("division by zero polynomial");
    };
    assert_eq!(divisor[dd], BigInt::one(), "divisor must be monic");
    while let Some(nd) = degree(&numerator) {
        if nd < dd {
            break;
        }
        let shift = nd - dd;
        let coeff = numerator[nd].clone();
        if coeff.is_zero() {
            numerator.pop();
            continue;
        }
        for (j, dcoeff) in divisor.iter().enumerate().take(dd + 1) {
            numerator[j + shift] -= &coeff * dcoeff;
        }
        numerator = trim(numerator);
    }
    trim(numerator)
}

fn cyclotomic_cached(n: usize, cache: &mut BTreeMap<usize, Vec<BigInt>>) -> Vec<BigInt> {
    if let Some(poly) = cache.get(&n) {
        return poly.clone();
    }
    let poly = if n == 1 {
        vec![BigInt::from(-1), BigInt::one()]
    } else {
        let mut p = vec![BigInt::zero(); n + 1];
        p[0] = BigInt::from(-1);
        p[n] = BigInt::one();
        for d in divisors(n).into_iter().filter(|&d| d < n) {
            let phi_d = cyclotomic_cached(d, cache);
            p = exact_div_monic(p, &phi_d).expect("cyclotomic division is exact");
        }
        p
    };
    cache.insert(n, poly.clone());
    poly
}

/// Return the `n`-th cyclotomic polynomial in ascending coefficient order.
pub fn cyclotomic_polynomial_bigint(n: usize) -> Vec<BigInt> {
    assert!(n >= 1, "cyclotomic index must be positive");
    cyclotomic_cached(n, &mut BTreeMap::new())
}

fn evaluate_at_one(coeffs: &[BigInt]) -> BigInt {
    coeffs.iter().fold(BigInt::zero(), |acc, coeff| acc + coeff)
}

/// Exactly reduce `P(q)` at `q = zeta_order^power`.
pub fn root_of_unity_evaluation_bigint(
    coeffs: &[BigInt],
    group_order: usize,
    power: usize,
) -> RootOfUnityEvaluation {
    assert!(group_order >= 1, "group order must be positive");
    let power = power % group_order;
    let root_order = if power == 0 {
        1
    } else {
        group_order / group_order.gcd(&power)
    };
    if root_order == 1 {
        let value = evaluate_at_one(coeffs);
        return RootOfUnityEvaluation {
            group_order,
            power,
            root_order,
            integer_value: Some(value.clone()),
            remainder: vec![value],
        };
    }

    let phi = cyclotomic_polynomial_bigint(root_order);
    let remainder = remainder_monic(coeffs.to_vec(), &phi);
    let integer_value = match remainder.as_slice() {
        [] => Some(BigInt::zero()),
        [constant] => Some(constant.clone()),
        _ => None,
    };
    RootOfUnityEvaluation {
        group_order,
        power,
        root_order,
        integer_value,
        remainder,
    }
}

/// Check one CSP candidate.  If `fixed_counts` is `None`, this only profiles
/// root-of-unity evaluations.
pub fn cyclic_sieving_report_bigint(
    coeffs: &[BigInt],
    order: usize,
    fixed_counts: Option<&[BigInt]>,
) -> Result<CyclicSievingReport, String> {
    if order == 0 {
        return Err("cyclic group order must be positive".to_string());
    }
    if let Some(counts) = fixed_counts {
        if counts.len() != order {
            return Err(format!(
                "expected {order} fixed-point counts for order {order}, got {}",
                counts.len()
            ));
        }
    }

    let coefficients = trim(coeffs.to_vec());
    let mut evaluations = Vec::with_capacity(order);
    let mut checks = Vec::new();
    for power in 0..order {
        let evaluation = root_of_unity_evaluation_bigint(&coefficients, order, power);
        if let Some(counts) = fixed_counts {
            let expected = counts[power].clone();
            let holds = evaluation.integer_value.as_ref() == Some(&expected);
            checks.push(CyclicSievingPowerCheck {
                power,
                root_order: evaluation.root_order,
                expected_fixed_points: expected,
                evaluation: evaluation.clone(),
                holds,
            });
        }
        evaluations.push(evaluation);
    }
    let holds = fixed_counts.map(|_| checks.iter().all(|check| check.holds));
    Ok(CyclicSievingReport {
        order,
        coefficients,
        fixed_counts: fixed_counts.map(|counts| counts.to_vec()),
        evaluations,
        checks,
        holds,
    })
}

/// Build sequence CSP profiles for candidate orders `index + offset`.
pub fn cyclic_sieving_sequence_reports_bigint(
    polynomials: &[Vec<BigInt>],
    first_index: isize,
    offsets: &[isize],
    fixed_counts: &BTreeMap<(isize, usize), Vec<BigInt>>,
) -> Vec<CyclicSievingSequenceItem> {
    polynomials
        .iter()
        .enumerate()
        .map(|(row, coeffs)| {
            let index = first_index + row as isize;
            let mut candidate_orders = Vec::new();
            for &offset in offsets {
                let order_index = index + offset;
                if order_index <= 0 {
                    continue;
                }
                let order = order_index as usize;
                let counts = fixed_counts.get(&(index, order)).map(Vec::as_slice);
                if let Ok(report) = cyclic_sieving_report_bigint(coeffs, order, counts) {
                    candidate_orders.push(report);
                }
            }
            CyclicSievingSequenceItem {
                row,
                index,
                coefficients: trim(coeffs.clone()),
                candidate_orders,
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Orbit decomposition, the Alexandersson-Amini criterion and bicyclic sieving
// ---------------------------------------------------------------------------
//
// Reference: P. Alexandersson and N. Amini, *The cone of cyclic sieving
// phenomena*, Discrete Math. 342 (2019) 1581-1601 (arXiv:1804.01447).  Their
// Desarmenien proposition (integrality), the lemma giving the bases `B_1`,
// `B_2` of `M(n)`, and the theorem on the existence of a cyclic action
// (`thm:cspCond`) are used below.  All arithmetic is exact integer arithmetic.

/// Mobius function `mu(n)` for `n >= 1`.
fn mobius(mut n: usize) -> i32 {
    let mut result = 1;
    let mut p = 2;
    while p * p <= n {
        if n % p == 0 {
            n /= p;
            if n % p == 0 {
                return 0;
            }
            result = -result;
        }
        p += 1;
    }
    if n > 1 {
        result = -result;
    }
    result
}

/// `gcd(j, n)` with the convention `gcd(0, n) = n`.
fn gcd_class(j: usize, n: usize) -> usize {
    j.gcd(&n)
}

fn big_mobius_term(mu: i32, value: &BigInt) -> BigInt {
    match mu {
        0 => BigInt::zero(),
        1 => value.clone(),
        _ => -value,
    }
}

/// Reduce `f(q)` modulo `q^n - 1`; the result has exactly `n` entries
/// `c_0, ..., c_{n-1}` (ascending order).
pub fn reduce_modulo_cyclic_bigint(coeffs: &[BigInt], n: usize) -> Vec<BigInt> {
    assert!(n >= 1, "group order must be positive");
    let mut out = vec![BigInt::zero(); n];
    for (i, c) in coeffs.iter().enumerate() {
        out[i % n] += c;
    }
    out
}

/// The residues `c_j` (`j = 0..n-1`) of `f mod q^n - 1` that share one value
/// of `gcd(j, n)`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ResidueClass {
    /// The divisor `e = gcd(j, n)` of `n` (with `gcd(0, n) = n`).
    pub gcd_class: usize,
    /// The indices `j` with `gcd(j, n) = e`, ascending.
    pub indices: Vec<usize>,
    /// The residues `c_j` for these indices.
    pub residues: Vec<BigInt>,
    /// The common value `c_e`, or `None` if the residues differ.
    pub common_value: Option<BigInt>,
}

/// One row of an orbit decomposition: orbits of size `orbit_size` (a divisor of
/// `n`).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct OrbitSizeEntry {
    /// The orbit size `k`, a divisor of `n`.
    pub orbit_size: usize,
    /// The index `r = n / k` of the orbit-basis polynomial `h_r`.
    pub basis_index: usize,
    /// The number `a_r` of orbits of size `k`.  For a polynomial that is not
    /// realizable this can be negative.
    pub orbit_count: BigInt,
    /// The number `S_k = k * a_r` of elements lying in orbits of size `k`.
    pub element_count: BigInt,
}

/// Orbit decomposition of a polynomial for a cyclic group `C_n`.
///
/// Let `c_0, ..., c_{n-1}` be the coefficients of `f mod q^n - 1`.
///
/// * **Integrality** (Desarmenien; Alexandersson--Amini, Prop. on `M(n)`):
///   `f(zeta_n^d)` is an integer for every `d` if and only if `c_j = c_{gcd(j,n)}`
///   for all `j` (with `gcd(0, n) = n`).
/// * **Orbit basis**: with `h_r(q) = sum_{i=0}^{n/r-1} q^{r i}` (the polynomial of
///   one orbit of size `n/r`), an integral `f` satisfies
///   `f = sum_{r | n} a_r h_r (mod q^n - 1)` with `a_r = sum_{e | r} mu(r/e) c_e`.
///   A `C_n`-action realizing the CSP must have exactly `a_r` orbits of size
///   `n/r`, and `S_k = k a_{n/k} = sum_{j | k} mu(k/j) f(zeta_n^j)` elements in
///   orbits of size `k`.
/// * **Fixed points**: `f(zeta_n^d) = sum_{r | n, (n/r) | d} (n/r) a_r`.
/// * **Cone condition** (Alexandersson--Amini, `thm:cspCond`): for
///   `f` in `N[q]` that is integral, a `C_n`-action with `(X, C_n, f)` a CSP
///   exists if and only if all `a_r >= 0`.  The theorem is stated for natural
///   coefficients; `natural_coefficients` and `cone_condition` are reported
///   separately.  (The construction in the proof only uses integrality and
///   `a_r >= 0`, which is what `realizable` reports.)
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct OrbitDecomposition {
    pub order: usize,
    /// The input polynomial with trailing zeros removed.
    pub coefficients: Vec<BigInt>,
    /// `f mod q^n - 1`: exactly `n` residues `c_0, ..., c_{n-1}`.
    pub reduced_residues: Vec<BigInt>,
    /// The residue classes `gcd(j, n) = e` for `e | n`, ascending in `e`.
    pub classes: Vec<ResidueClass>,
    /// `f(zeta_n^d) in Z` for all `d` (equivalently `c_j = c_{gcd(j,n)}`).
    pub integral: bool,
    /// The divisors `e` whose residue class is not constant; empty iff integral.
    pub violating_classes: Vec<usize>,
    /// Orbit data for each orbit size `k | n` (ascending in `k`); `None` unless
    /// the polynomial is integral.
    pub orbits: Option<Vec<OrbitSizeEntry>>,
    /// `f(zeta_n^d)` for `d = 0..n-1`, computed from the orbit basis; `None`
    /// unless integral.
    pub fixed_points: Option<Vec<BigInt>>,
    /// All input coefficients are nonnegative (the hypothesis `f in N[q]` of
    /// `thm:cspCond`).
    pub natural_coefficients: bool,
    /// All `a_r >= 0`, equivalently all `S_k >= 0`; `None` unless integral.
    pub cone_condition: Option<bool>,
    /// Some `C_n`-action has fixed-point counts `f(zeta_n^d)`: the polynomial is
    /// integral and satisfies the cone condition.
    pub realizable: bool,
}

/// Compute the orbit decomposition of `coeffs` for the cyclic group `C_order`.
pub fn cyclic_sieving_orbit_decomposition_bigint(
    coeffs: &[BigInt],
    order: usize,
) -> Result<OrbitDecomposition, String> {
    if order == 0 {
        return Err("cyclic group order must be positive".to_string());
    }
    let n = order;
    let coefficients = trim(coeffs.to_vec());
    let reduced_residues = reduce_modulo_cyclic_bigint(&coefficients, n);
    let divs = divisors(n);

    let classes: Vec<ResidueClass> = divs
        .iter()
        .map(|&e| {
            let indices: Vec<usize> = (0..n).filter(|&j| gcd_class(j, n) == e).collect();
            let residues: Vec<BigInt> = indices
                .iter()
                .map(|&j| reduced_residues[j].clone())
                .collect();
            let common_value = match residues.first() {
                Some(first) if residues.iter().all(|r| r == first) => Some(first.clone()),
                _ => None,
            };
            ResidueClass {
                gcd_class: e,
                indices,
                residues,
                common_value,
            }
        })
        .collect();
    let violating_classes: Vec<usize> = classes
        .iter()
        .filter(|class| class.common_value.is_none())
        .map(|class| class.gcd_class)
        .collect();
    let integral = violating_classes.is_empty();
    let natural_coefficients = coefficients
        .iter()
        .all(|c| c.sign() != num_bigint::Sign::Minus);

    let (orbits, fixed_points, cone_condition) = if integral {
        let values: BTreeMap<usize, BigInt> = classes
            .iter()
            .map(|class| (class.gcd_class, class.common_value.clone().unwrap()))
            .collect();
        // a_r = sum_{e | r} mu(r/e) c_e
        let mut orbits = Vec::with_capacity(divs.len());
        for &k in &divs {
            let r = n / k;
            let mut a_r = BigInt::zero();
            for &e in divs.iter().filter(|&&e| r % e == 0) {
                a_r += big_mobius_term(mobius(r / e), &values[&e]);
            }
            orbits.push(OrbitSizeEntry {
                orbit_size: k,
                basis_index: r,
                element_count: &a_r * BigInt::from(k),
                orbit_count: a_r,
            });
        }
        orbits.sort_by_key(|entry| entry.orbit_size);
        let fixed = fixed_points_from_entries(&orbits, n);
        let cone = orbits
            .iter()
            .all(|entry| entry.orbit_count.sign() != num_bigint::Sign::Minus);
        (Some(orbits), Some(fixed), Some(cone))
    } else {
        (None, None, None)
    };

    Ok(OrbitDecomposition {
        order,
        coefficients,
        reduced_residues,
        classes,
        integral,
        violating_classes,
        orbits,
        fixed_points,
        natural_coefficients,
        realizable: cone_condition == Some(true),
        cone_condition,
    })
}

/// `fix(g^d) = sum_{k | gcd(d, n)} S_k` for `d = 0..n-1`.
fn fixed_points_from_entries(entries: &[OrbitSizeEntry], n: usize) -> Vec<BigInt> {
    (0..n)
        .map(|d| {
            entries
                .iter()
                .filter(|entry| d % entry.orbit_size == 0)
                .fold(BigInt::zero(), |acc, entry| acc + &entry.element_count)
        })
        .collect()
}

/// Fixed-point counts `fix(g^d)`, `d = 0..n-1`, of a `C_n`-action with the
/// given orbit counts `(orbit_size, number_of_orbits)`.
///
/// Every orbit size must divide `n` and counts must be nonnegative; repeated
/// sizes are summed.  `fix(g^d) = sum_{k | d, k | n} k * (number of orbits of size k)`.
pub fn fixed_points_from_orbit_counts(
    orbit_counts: &[(usize, BigInt)],
    order: usize,
) -> Result<Vec<BigInt>, String> {
    if order == 0 {
        return Err("cyclic group order must be positive".to_string());
    }
    let mut merged: BTreeMap<usize, BigInt> = BTreeMap::new();
    for (size, count) in orbit_counts {
        if *size == 0 || order % size != 0 {
            return Err(format!(
                "orbit size {size} does not divide the order {order}"
            ));
        }
        if count.sign() == num_bigint::Sign::Minus {
            return Err(format!("negative number of orbits of size {size}"));
        }
        *merged.entry(*size).or_insert_with(BigInt::zero) += count;
    }
    let entries: Vec<OrbitSizeEntry> = merged
        .into_iter()
        .map(|(k, count)| OrbitSizeEntry {
            orbit_size: k,
            basis_index: order / k,
            element_count: &count * BigInt::from(k),
            orbit_count: count,
        })
        .collect();
    Ok(fixed_points_from_entries(&entries, order))
}

/// Orbit counts of an action given as a multiset of orbit sizes (each a divisor
/// of `order`).  Returns one entry for every divisor `k` of `order`, ascending.
pub fn orbit_counts_from_orbit_sizes(
    orbit_sizes: &[usize],
    order: usize,
) -> Result<Vec<OrbitSizeEntry>, String> {
    if order == 0 {
        return Err("cyclic group order must be positive".to_string());
    }
    let mut counts: BTreeMap<usize, BigInt> = divisors(order)
        .into_iter()
        .map(|k| (k, BigInt::zero()))
        .collect();
    for &size in orbit_sizes {
        match counts.get_mut(&size) {
            Some(count) => *count += 1,
            None => {
                return Err(format!(
                    "orbit size {size} does not divide the order {order}"
                ))
            }
        }
    }
    Ok(counts
        .into_iter()
        .map(|(k, count)| OrbitSizeEntry {
            orbit_size: k,
            basis_index: order / k,
            element_count: &count * BigInt::from(k),
            orbit_count: count,
        })
        .collect())
}

/// Fixed-point counts `fix(g^d)`, `d = 0..n-1`, of an action with the given
/// orbit sizes.
pub fn fixed_points_from_orbit_sizes(
    orbit_sizes: &[usize],
    order: usize,
) -> Result<Vec<BigInt>, String> {
    let entries = orbit_counts_from_orbit_sizes(orbit_sizes, order)?;
    Ok(fixed_points_from_entries(&entries, order))
}

/// Recover the orbit counts of a `C_n`-action from its fixed-point counts
/// `fix(g^d)`, `d = 0..n-1`, by Mobius inversion.
///
/// With `F_j = fix(g^j)` for `j | n` (`F_n = fix(g^0)` is the size of the set),
/// `S_k = sum_{j | k} mu(k/j) F_j` is the number of elements in orbits of size
/// `k`.  An error is returned if the table has the wrong length, if
/// `fix(g^d)` does not depend only on `gcd(d, n)`, or if some `S_k` is negative
/// or not divisible by `k` (so the table does not come from an action).
/// Returns one entry for every divisor of `n`, ascending.
pub fn orbit_counts_from_fixed_points(
    fixed_points: &[BigInt],
    order: usize,
) -> Result<Vec<OrbitSizeEntry>, String> {
    if order == 0 {
        return Err("cyclic group order must be positive".to_string());
    }
    let n = order;
    if fixed_points.len() != n {
        return Err(format!(
            "expected {n} fixed-point counts for order {n}, got {}",
            fixed_points.len()
        ));
    }
    let divs = divisors(n);
    // F_j for j | n, using fix(g^n) = fix(g^0).
    let f_at = |j: usize| &fixed_points[j % n];
    for (d, value) in fixed_points.iter().enumerate() {
        let g = gcd_class(d, n);
        if value != f_at(g) {
            return Err(format!(
                "fixed-point counts are not those of a C_{n}-action: fix(g^{d}) = {} but \
                 fix(g^{g}) = {} although gcd({d}, {n}) = {g}",
                value,
                f_at(g)
            ));
        }
    }
    let mut entries = Vec::with_capacity(divs.len());
    for &k in &divs {
        let mut s_k = BigInt::zero();
        for &j in divs.iter().filter(|&&j| k % j == 0) {
            s_k += big_mobius_term(mobius(k / j), f_at(j));
        }
        let kk = BigInt::from(k);
        if s_k.sign() == num_bigint::Sign::Minus {
            return Err(format!(
                "fixed-point counts are not those of a C_{n}-action: S_{k} = {s_k} is negative"
            ));
        }
        if !(&s_k % &kk).is_zero() {
            return Err(format!(
                "fixed-point counts are not those of a C_{n}-action: S_{k} = {s_k} is not \
                 divisible by {k}"
            ));
        }
        entries.push(OrbitSizeEntry {
            orbit_size: k,
            basis_index: n / k,
            orbit_count: &s_k / &kk,
            element_count: s_k,
        });
    }
    Ok(entries)
}

/// Result of checking a polynomial against a concrete cyclic action given by
/// its orbit sizes.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct OrbitSizesCheck {
    pub order: usize,
    /// The supplied multiset of orbit sizes.
    pub orbit_sizes: Vec<usize>,
    /// Number of elements of the set (the sum of the orbit sizes).
    pub num_elements: BigInt,
    /// Actual orbit counts of the action (every divisor of `n`).
    pub actual_orbits: Vec<OrbitSizeEntry>,
    /// Actual fixed-point counts `fix(g^d)`, `d = 0..n-1`.
    pub actual_fixed_points: Vec<BigInt>,
    /// The predicted decomposition of the polynomial.
    pub decomposition: OrbitDecomposition,
    /// Orbit sizes `k` where the predicted and actual orbit counts differ
    /// (empty if the polynomial is not integral).
    pub mismatched_orbit_sizes: Vec<usize>,
    /// CSP holds: the polynomial is integral and the predicted orbit counts
    /// equal the actual ones.  Equivalent to `f(zeta^d) = fix(g^d)` for all `d`.
    pub holds: bool,
    /// The ordinary fixed-point CSP report for `actual_fixed_points`.
    pub sieving_report: CyclicSievingReport,
}

/// Check the CSP for a polynomial and an action given by its orbit sizes.
pub fn cyclic_sieving_check_orbit_sizes(
    coeffs: &[BigInt],
    order: usize,
    orbit_sizes: &[usize],
) -> Result<OrbitSizesCheck, String> {
    let actual_orbits = orbit_counts_from_orbit_sizes(orbit_sizes, order)?;
    let actual_fixed_points = fixed_points_from_entries(&actual_orbits, order);
    let decomposition = cyclic_sieving_orbit_decomposition_bigint(coeffs, order)?;
    let sieving_report = cyclic_sieving_report_bigint(coeffs, order, Some(&actual_fixed_points))?;
    let mismatched_orbit_sizes: Vec<usize> = match &decomposition.orbits {
        Some(predicted) => predicted
            .iter()
            .zip(&actual_orbits)
            .filter(|(p, a)| p.orbit_count != a.orbit_count)
            .map(|(p, _)| p.orbit_size)
            .collect(),
        None => Vec::new(),
    };
    let holds = decomposition.integral && mismatched_orbit_sizes.is_empty();
    Ok(OrbitSizesCheck {
        order,
        orbit_sizes: orbit_sizes.to_vec(),
        num_elements: orbit_sizes.iter().map(|&k| BigInt::from(k)).sum(),
        actual_orbits,
        actual_fixed_points,
        decomposition,
        mismatched_orbit_sizes,
        holds,
        sieving_report,
    })
}

/// Largest `n1 * n2` accepted by [`bicyclic_sieving_report_bigint`].
pub const BICYCLIC_MAX_GROUP_SIZE: usize = 2500;

/// Bicyclic sieving report for `G = C_{n1} x C_{n2}` and `f(q, t)`.
///
/// `f[i][j]` is the coefficient of `q^i t^j`.  The polynomial is reduced modulo
/// `(q^{n1} - 1, t^{n2} - 1)`.  The values `f(zeta_{n1}^{d1}, zeta_{n2}^{d2})` are
/// computed exactly (as one reduction modulo the cyclotomic polynomial
/// `Phi_L`, `L = lcm(n1, n2)`).
///
/// If `c_{j1,j2} = c_{gcd(j1,n1), gcd(j2,n2)}` for all `j1, j2` (the polynomial
/// lies in the span of the products `h_{r1}(q) h_{r2}(t)`, the product basis
/// `B_2` of the Alexandersson--Amini multivariate basis theorem), the
/// coordinates are
/// `a_{r1,r2} = sum_{e1 | r1, e2 | r2} mu(r1/e1) mu(r2/e2) c_{e1,e2}`, and
/// `f(zeta^{d1}, omega^{d2}) = sum_{(n1/r1) | d1, (n2/r2) | d2} (n1/r1)(n2/r2) a_{r1,r2}`.
///
/// These product coordinates are **not** orbit counts and need not be
/// nonnegative: for the diagonal action of `C_2 x C_2` on a 2-element set,
/// `f = 1 + q t` and the coordinates `a_{1,1}, a_{1,2}, a_{2,1}, a_{2,2}` are
/// `1, -1, -1, 2`.  Products `h_{r1}(q) h_{r2}(t)` are the CSP polynomials of
/// product orbits `C_{n1}/H1 x C_{n2}/H2`; a general `G`-orbit `G/H` need not be
/// of this form.
///
/// Integrality at all points is **not** equivalent to product-basis form when
/// `gcd(n1, n2) > 1`: for `n1 = n2 = 3`,
/// `f = sum_{i,j} (zeta_3^{i+j} + zeta_3^{-(i+j)}) q^i t^j` has coefficients
/// `2` if `i + j = 0 (mod 3)` and `-1` otherwise, takes only the values `0` and
/// `9` on `mu_3 x mu_3`, yet is not in the product span.  (For coprime `n1, n2`
/// the two notions agree, because the Galois groups act independently.)  We
/// therefore report `integral` (from exact evaluation) and `product_basis`
/// separately; the CSP check against a fixed-point table uses the exact
/// evaluations only.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BicyclicSievingReport {
    pub order1: usize,
    pub order2: usize,
    pub coefficient_matrix: Vec<Vec<BigInt>>,
    /// `f mod (q^{n1} - 1, t^{n2} - 1)`, an `n1 x n2` matrix.
    pub reduced_matrix: Vec<Vec<BigInt>>,
    /// `f(zeta_{n1}^{d1}, zeta_{n2}^{d2})` for `d1 < n1`, `d2 < n2`; `None`
    /// where the value is not an integer.
    pub evaluations: Vec<Vec<Option<BigInt>>>,
    /// Every evaluation is an integer.
    pub integral: bool,
    /// All input coefficients are nonnegative.
    pub natural_coefficients: bool,
    /// Divisors `r1` of `n1`, ascending: the row index of `product_coordinates`
    /// (the row for `r1` belongs to `h_{r1}(q)`, a polynomial of one orbit of
    /// size `n1 / r1`).
    pub divisors1: Vec<usize>,
    /// Divisors `r2` of `n2`, ascending: the column index of `product_coordinates`.
    pub divisors2: Vec<usize>,
    /// Whether `c_{j1,j2} = c_{gcd(j1,n1), gcd(j2,n2)}` for all `j1, j2`.
    pub product_basis_form: bool,
    /// The coordinates `a_{r1,r2}` in the basis `h_{r1}(q) h_{r2}(t)`, rows
    /// indexed by `r1 in divisors1` and columns by `r2 in divisors2`; `Some` iff
    /// `product_basis_form`.
    pub product_coordinates: Option<Vec<Vec<BigInt>>>,
    /// `f(zeta^{d1}, omega^{d2})` recomputed from `product_coordinates`.
    pub product_formula_evaluations: Option<Vec<Vec<BigInt>>>,
    /// The supplied table `fix(g1^{d1} g2^{d2})`.
    pub fixed_table: Option<Vec<Vec<BigInt>>>,
    /// Positions `(d1, d2)` where the evaluation is not the supplied count.
    pub mismatches: Vec<(usize, usize)>,
    /// `Some(true)` iff a table was supplied and equals the evaluations.
    pub holds: Option<bool>,
}

/// Compute a bicyclic sieving report; see [`BicyclicSievingReport`].
pub fn bicyclic_sieving_report_bigint(
    coeff_matrix: &[Vec<BigInt>],
    order1: usize,
    order2: usize,
    fixed_table: Option<&[Vec<BigInt>]>,
) -> Result<BicyclicSievingReport, String> {
    if order1 == 0 || order2 == 0 {
        return Err("cyclic group orders must be positive".to_string());
    }
    if order1
        .checked_mul(order2)
        .is_none_or(|g| g > BICYCLIC_MAX_GROUP_SIZE)
    {
        return Err(format!(
            "group order {order1} x {order2} exceeds the limit {BICYCLIC_MAX_GROUP_SIZE}"
        ));
    }
    if let Some(table) = fixed_table {
        if table.len() != order1 || table.iter().any(|row| row.len() != order2) {
            return Err(format!(
                "fixed-point table must be a {order1} x {order2} matrix"
            ));
        }
    }
    let (n1, n2) = (order1, order2);
    let mut reduced = vec![vec![BigInt::zero(); n2]; n1];
    for (i, row) in coeff_matrix.iter().enumerate() {
        for (j, c) in row.iter().enumerate() {
            reduced[i % n1][j % n2] += c;
        }
    }
    let natural_coefficients = coeff_matrix
        .iter()
        .flatten()
        .all(|c| c.sign() != num_bigint::Sign::Minus);

    // Exact evaluations: f(zeta_{n1}^{d1}, zeta_{n2}^{d2}) = R(xi), xi = zeta_L,
    // R = sum c_{ij} x^{(i d1 L/n1 + j d2 L/n2) mod L}.
    let l = n1.lcm(&n2);
    let phi_l = cyclotomic_polynomial_bigint(l);
    let mut evaluations = vec![vec![None; n2]; n1];
    for (d1, eval_row) in evaluations.iter_mut().enumerate() {
        for (d2, slot) in eval_row.iter_mut().enumerate() {
            let (u1, u2) = (d1 * (l / n1), d2 * (l / n2));
            let mut poly = vec![BigInt::zero(); l];
            for (i, row) in reduced.iter().enumerate() {
                for (j, c) in row.iter().enumerate() {
                    if !c.is_zero() {
                        poly[(i * u1 + j * u2) % l] += c;
                    }
                }
            }
            *slot = match remainder_monic(poly, &phi_l).as_slice() {
                [] => Some(BigInt::zero()),
                [constant] => Some(constant.clone()),
                _ => None,
            };
        }
    }
    let integral = evaluations.iter().flatten().all(Option::is_some);

    let divisors1 = divisors(n1);
    let divisors2 = divisors(n2);
    let product_basis_form = (0..n1).all(|j1| {
        (0..n2).all(|j2| reduced[j1][j2] == reduced[gcd_class(j1, n1) % n1][gcd_class(j2, n2) % n2])
    });
    let (product_coordinates, product_formula_evaluations) = if product_basis_form {
        let c_at = |e1: usize, e2: usize| &reduced[e1 % n1][e2 % n2];
        let coords: Vec<Vec<BigInt>> = divisors1
            .iter()
            .map(|&r1| {
                divisors2
                    .iter()
                    .map(|&r2| {
                        let mut a = BigInt::zero();
                        for &e1 in divisors1.iter().filter(|&&e| r1 % e == 0) {
                            for &e2 in divisors2.iter().filter(|&&e| r2 % e == 0) {
                                let mu = mobius(r1 / e1) * mobius(r2 / e2);
                                a += big_mobius_term(mu, c_at(e1, e2));
                            }
                        }
                        a
                    })
                    .collect()
            })
            .collect();
        let formula: Vec<Vec<BigInt>> = (0..n1)
            .map(|d1| {
                (0..n2)
                    .map(|d2| {
                        let mut value = BigInt::zero();
                        for (x1, &r1) in divisors1.iter().enumerate() {
                            for (x2, &r2) in divisors2.iter().enumerate() {
                                // (n1/r1) must divide d1 and (n2/r2) must divide d2
                                let (k1, k2) = (n1 / r1, n2 / r2);
                                if d1 % k1 == 0 && d2 % k2 == 0 {
                                    value += &coords[x1][x2] * BigInt::from(k1 * k2);
                                }
                            }
                        }
                        value
                    })
                    .collect()
            })
            .collect();
        (Some(coords), Some(formula))
    } else {
        (None, None)
    };

    let mut mismatches = Vec::new();
    if let Some(table) = fixed_table {
        for (d1, row) in table.iter().enumerate() {
            for (d2, expected) in row.iter().enumerate() {
                if evaluations[d1][d2].as_ref() != Some(expected) {
                    mismatches.push((d1, d2));
                }
            }
        }
    }
    let holds = fixed_table.map(|_| mismatches.is_empty());
    Ok(BicyclicSievingReport {
        order1,
        order2,
        coefficient_matrix: coeff_matrix.to_vec(),
        reduced_matrix: reduced,
        evaluations,
        integral,
        natural_coefficients,
        divisors1,
        divisors2,
        product_basis_form,
        product_coordinates,
        product_formula_evaluations,
        fixed_table: fixed_table.map(|t| t.to_vec()),
        mismatches,
        holds,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn b(values: &[i64]) -> Vec<BigInt> {
        values.iter().map(|&v| BigInt::from(v)).collect()
    }

    #[test]
    fn cyclotomic_polynomials_are_exact() {
        assert_eq!(cyclotomic_polynomial_bigint(1), b(&[-1, 1]));
        assert_eq!(cyclotomic_polynomial_bigint(2), b(&[1, 1]));
        assert_eq!(cyclotomic_polynomial_bigint(3), b(&[1, 1, 1]));
        assert_eq!(cyclotomic_polynomial_bigint(4), b(&[1, 0, 1]));
        assert_eq!(
            poly_mul(&cyclotomic_polynomial_bigint(3), &b(&[-1, 1])),
            b(&[-1, 0, 0, 1])
        );
    }

    #[test]
    fn root_of_unity_integer_profiles() {
        let p = b(&[1, 1, 1, 1]);
        let at_one = root_of_unity_evaluation_bigint(&p, 4, 0);
        assert_eq!(at_one.integer_value, Some(BigInt::from(4)));

        let at_minus_one = root_of_unity_evaluation_bigint(&p, 4, 2);
        assert_eq!(at_minus_one.integer_value, Some(BigInt::from(0)));

        let at_i = root_of_unity_evaluation_bigint(&p, 4, 1);
        assert_eq!(at_i.integer_value, Some(BigInt::from(0)));
    }

    #[test]
    fn checks_basic_csp() {
        let p = b(&[1, 1]);
        let counts = b(&[2, 0]);
        let report = cyclic_sieving_report_bigint(&p, 2, Some(&counts)).unwrap();
        assert_eq!(report.holds, Some(true));
    }

    // ----- orbit decomposition and bicyclic sieving -----

    fn big(v: i64) -> BigInt {
        BigInt::from(v)
    }

    /// Gaussian binomial `[n choose k]_q`, ascending coefficients.
    fn gaussian_binomial(n: usize, k: usize) -> Vec<BigInt> {
        // [n,k] = [n-1,k-1] + q^k [n-1,k]
        let mut table = vec![vec![vec![BigInt::zero()]; k + 1]; n + 1];
        for row in table.iter_mut() {
            row[0] = vec![BigInt::one()];
        }
        for i in 1..=n {
            for j in 1..=k.min(i) {
                let mut sum = table[i - 1][j - 1].clone();
                let shifted = &table[i - 1][j];
                if !shifted.iter().all(Zero::is_zero) {
                    let mut padded = vec![BigInt::zero(); j];
                    padded.extend(shifted.iter().cloned());
                    if padded.len() > sum.len() {
                        sum.resize(padded.len(), BigInt::zero());
                    }
                    for (x, y) in padded.iter().enumerate() {
                        sum[x] += y;
                    }
                }
                table[i][j] = trim(sum);
            }
        }
        table[n][k].clone()
    }

    /// Orbit sizes of the rotation action of `C_n` on `k`-subsets of `[n]`.
    fn rotation_orbit_sizes(n: usize, k: usize) -> Vec<usize> {
        let mask = (1u32 << n) - 1;
        let rotate = |x: u32| ((x << 1) | (x >> (n - 1))) & mask;
        let mut seen = std::collections::BTreeSet::new();
        let mut sizes = Vec::new();
        for x in 0..(1u32 << n) {
            if x.count_ones() as usize != k || seen.contains(&x) {
                continue;
            }
            let mut size = 0;
            let mut y = x;
            loop {
                seen.insert(y);
                size += 1;
                y = rotate(y);
                if y == x {
                    break;
                }
            }
            sizes.push(size);
        }
        sizes
    }

    /// Fixed points `fix(g^d)` by direct enumeration.
    fn rotation_fixed_points(n: usize, k: usize) -> Vec<BigInt> {
        let mask = (1u32 << n) - 1;
        (0..n)
            .map(|d| {
                let rot = |x: u32| ((x << d) | (x >> (n - d))) & mask;
                let count = (0..(1u32 << n))
                    .filter(|&x| x.count_ones() as usize == k)
                    .filter(|&x| d == 0 || rot(x) == x)
                    .count();
                big(count as i64)
            })
            .collect()
    }

    #[test]
    fn gaussian_binomial_helper_is_correct() {
        assert_eq!(gaussian_binomial(4, 2), b(&[1, 1, 2, 1, 1]));
        assert_eq!(gaussian_binomial(5, 1), b(&[1, 1, 1, 1, 1]));
    }

    #[test]
    fn mobius_values() {
        let expected = [1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0];
        for (i, &e) in expected.iter().enumerate() {
            assert_eq!(mobius(i + 1), e, "mu({})", i + 1);
        }
    }

    #[test]
    fn q_binomial_rotation_orbit_decomposition() {
        for (n, k) in [
            (4, 2),
            (5, 2),
            (6, 2),
            (6, 3),
            (7, 3),
            (8, 4),
            (9, 3),
            (10, 5),
            (10, 4),
            (12, 6),
        ] {
            let f = gaussian_binomial(n, k);
            let sizes = rotation_orbit_sizes(n, k);
            let fixed = rotation_fixed_points(n, k);
            // orbit sizes -> fixed points agrees with enumeration
            assert_eq!(fixed_points_from_orbit_sizes(&sizes, n).unwrap(), fixed);

            let dec = cyclic_sieving_orbit_decomposition_bigint(&f, n).unwrap();
            assert!(dec.integral, "[{n},{k}] integral");
            assert!(dec.natural_coefficients);
            assert_eq!(dec.cone_condition, Some(true));
            assert!(dec.realizable);
            assert_eq!(dec.fixed_points.as_ref().unwrap(), &fixed);
            // predicted orbit counts equal the actual ones
            let actual = orbit_counts_from_orbit_sizes(&sizes, n).unwrap();
            assert_eq!(dec.orbits.as_ref().unwrap(), &actual);
            // fixed points -> orbit counts round trip
            assert_eq!(orbit_counts_from_fixed_points(&fixed, n).unwrap(), actual);

            let check = cyclic_sieving_check_orbit_sizes(&f, n, &sizes).unwrap();
            assert!(check.holds);
            assert_eq!(check.sieving_report.holds, Some(true));
            assert!(check.mismatched_orbit_sizes.is_empty());
            assert_eq!(check.num_elements, f.iter().sum::<BigInt>());
        }
    }

    #[test]
    fn q_binomial_wrong_action_is_rejected() {
        // [4 choose 2]_q = 1 + q + 2q^2 + q^3 + q^4, n = 4: orbits 4,2 (6 elements).
        let f = gaussian_binomial(4, 2);
        let good = cyclic_sieving_check_orbit_sizes(&f, 4, &[4, 2]).unwrap();
        assert!(good.holds);
        let bad = cyclic_sieving_check_orbit_sizes(&f, 4, &[4, 1, 1]).unwrap();
        assert!(!bad.holds);
        assert_eq!(bad.sieving_report.holds, Some(false));
        assert_eq!(bad.mismatched_orbit_sizes, vec![1, 2]);
        assert!(cyclic_sieving_check_orbit_sizes(&f, 4, &[3]).is_err());
    }

    #[test]
    fn one_plus_q_order_two() {
        let dec = cyclic_sieving_orbit_decomposition_bigint(&b(&[1, 1]), 2).unwrap();
        assert!(dec.integral);
        let orbits = dec.orbits.unwrap();
        // one orbit of size 2, none of size 1
        assert_eq!(orbits[0].orbit_size, 1);
        assert_eq!(orbits[0].orbit_count, big(0));
        assert_eq!(orbits[1].orbit_size, 2);
        assert_eq!(orbits[1].orbit_count, big(1));
        assert_eq!(orbits[1].element_count, big(2));
        assert_eq!(dec.fixed_points.unwrap(), b(&[2, 0]));
        assert!(dec.realizable);
    }

    #[test]
    fn non_integral_example() {
        // 1 + q mod q^3 - 1: c_1 = 1 but c_2 = 0 (both in class gcd = 1).
        let dec = cyclic_sieving_orbit_decomposition_bigint(&b(&[1, 1]), 3).unwrap();
        assert!(!dec.integral);
        assert_eq!(dec.violating_classes, vec![1]);
        assert!(dec.orbits.is_none() && dec.fixed_points.is_none());
        assert_eq!(dec.cone_condition, None);
        assert!(!dec.realizable);
        // agrees with exact evaluation: 1 + zeta_3 is not an integer
        assert!(root_of_unity_evaluation_bigint(&b(&[1, 1]), 3, 1)
            .integer_value
            .is_none());
        let check = cyclic_sieving_check_orbit_sizes(&b(&[1, 1]), 3, &[1, 1]).unwrap();
        assert!(!check.holds);
        assert_eq!(check.sieving_report.holds, Some(false));
    }

    #[test]
    fn paper_example_is_integral_natural_but_not_realizable() {
        // Alexandersson-Amini, ex:counterex.
        let f = b(&[10, 1, 0, 3, 0, 1]);
        let dec = cyclic_sieving_orbit_decomposition_bigint(&f, 6).unwrap();
        assert!(dec.integral);
        assert!(dec.natural_coefficients);
        // f(zeta_6^j), j = 1..6 is 8,12,5,12,8,15 (j = 6 is d = 0)
        let fp = dec.fixed_points.clone().unwrap();
        let by_j: Vec<BigInt> = (1..=6).map(|j| fp[j % 6].clone()).collect();
        assert_eq!(by_j, b(&[8, 12, 5, 12, 8, 15]));
        // S_k, k = 1..6
        let orbits = dec.orbits.unwrap();
        let s: Vec<BigInt> = (1..=6)
            .map(|k| {
                orbits
                    .iter()
                    .find(|e| e.orbit_size == k)
                    .map(|e| e.element_count.clone())
                    .unwrap_or_else(BigInt::zero)
            })
            .collect();
        // S_4 and S_5 are not defined (k does not divide 6); the paper's formula gives 0.
        assert_eq!(s, b(&[8, 4, -3, 0, 0, 6]));
        assert_eq!(dec.cone_condition, Some(false));
        assert!(!dec.realizable);
        let neg = orbits.iter().find(|e| e.orbit_size == 3).unwrap();
        assert_eq!(neg.orbit_count, big(-1));
        // the fixed-point table is not that of an action
        assert!(orbit_counts_from_fixed_points(&fp, 6).is_err());
    }

    #[test]
    fn orbit_counts_from_fixed_points_errors() {
        // wrong length
        assert!(orbit_counts_from_fixed_points(&b(&[2, 0, 0]), 2).is_err());
        // fix(g^1) depends only on gcd(1,4) = 1, so fix(g) = fix(g^3) must agree
        assert!(orbit_counts_from_fixed_points(&b(&[4, 1, 0, 0]), 4).is_err());
        // S_2 = fix(g^2) - fix(g) = 3 is odd
        assert!(orbit_counts_from_fixed_points(&b(&[5, 2, 5, 2]), 4)
            .unwrap_err()
            .contains("not divisible"));
        // negative S_2
        assert!(orbit_counts_from_fixed_points(&b(&[3, 2, 1, 2]), 4)
            .unwrap_err()
            .contains("negative"));
    }

    #[test]
    fn orbit_size_conversion_roundtrip() {
        let sizes = [1, 1, 2, 3, 3, 6, 6, 6];
        let fixed = fixed_points_from_orbit_sizes(&sizes, 6).unwrap();
        let counts = orbit_counts_from_fixed_points(&fixed, 6).unwrap();
        let expected = orbit_counts_from_orbit_sizes(&sizes, 6).unwrap();
        assert_eq!(counts, expected);
        let pairs: Vec<(usize, BigInt)> = expected
            .iter()
            .map(|e| (e.orbit_size, e.orbit_count.clone()))
            .collect();
        assert_eq!(fixed_points_from_orbit_counts(&pairs, 6).unwrap(), fixed);
        assert!(fixed_points_from_orbit_sizes(&[4], 6).is_err());
        assert!(fixed_points_from_orbit_sizes(&[0], 6).is_err());
        assert!(
            fixed_points_from_orbit_counts(&[(2, big(-1))], 6).is_err(),
            "negative counts rejected"
        );
    }

    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self, bound: u64) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (self.0 >> 33) % bound
        }
    }

    #[test]
    fn orbit_formula_matches_cyclotomic_evaluation_on_random_polynomials() {
        let mut rng = Lcg(12345);
        for trial in 0..200 {
            let n = 1 + rng.next(12) as usize;
            let len = 1 + rng.next(20) as usize;
            // Half of the trials are forced integral by building f from orbit polynomials.
            let coeffs: Vec<BigInt> = if trial % 2 == 0 {
                (0..len).map(|_| big(rng.next(7) as i64 - 2)).collect()
            } else {
                let mut f = vec![BigInt::zero(); n];
                for k in divisors(n) {
                    let m = rng.next(4) as i64 - 1;
                    // orbit of size k has polynomial h_{n/k}
                    let r = n / k;
                    for i in 0..k {
                        f[(r * i) % n] += m;
                    }
                }
                f
            };
            let dec = cyclic_sieving_orbit_decomposition_bigint(&coeffs, n).unwrap();
            let all_integral = (0..n).all(|d| {
                root_of_unity_evaluation_bigint(&coeffs, n, d)
                    .integer_value
                    .is_some()
            });
            assert_eq!(dec.integral, all_integral, "n={n} f={coeffs:?}");
            if trial % 2 == 1 {
                assert!(dec.integral);
            }
            if let Some(fixed) = &dec.fixed_points {
                for (d, value) in fixed.iter().enumerate() {
                    assert_eq!(
                        Some(value),
                        root_of_unity_evaluation_bigint(&coeffs, n, d)
                            .integer_value
                            .as_ref(),
                        "n={n} d={d} f={coeffs:?}"
                    );
                }
                // S_k = sum_{j | k} mu(k/j) f(zeta^j) and a_r reconstruct f mod q^n - 1
                let orbits = dec.orbits.as_ref().unwrap();
                let mut rebuilt = vec![BigInt::zero(); n];
                for entry in orbits {
                    for i in 0..entry.orbit_size {
                        rebuilt[(entry.basis_index * i) % n] += &entry.orbit_count;
                    }
                }
                assert_eq!(rebuilt, dec.reduced_residues);
            }
        }
    }

    #[test]
    fn integral_natural_polynomial_with_negative_orbit_count_via_remainders() {
        // Same example through the general API at order 6 with fixed counts of an
        // existing action never matching: CSP fails for every candidate action.
        let f = b(&[10, 1, 0, 3, 0, 1]);
        for sizes in [vec![1; 15], vec![6, 6, 3], vec![6, 6, 1, 1, 1]] {
            let check = cyclic_sieving_check_orbit_sizes(&f, 6, &sizes).unwrap();
            assert!(!check.holds, "{sizes:?}");
        }
    }

    #[test]
    fn bicyclic_one_plus_qt() {
        let f = vec![b(&[1, 0]), b(&[0, 1])];
        // diagonal action of C_2 x C_2 on two points, swapped by both generators
        let table = vec![b(&[2, 0]), b(&[0, 2])];
        let report = bicyclic_sieving_report_bigint(&f, 2, 2, Some(&table)).unwrap();
        assert!(report.integral);
        assert!(report.product_basis_form);
        assert_eq!(report.holds, Some(true));
        assert_eq!(report.divisors1, vec![1, 2]);
        assert_eq!(
            report.product_coordinates.as_ref().unwrap(),
            &vec![b(&[1, -1]), b(&[-1, 2])]
        );
        assert_eq!(report.product_formula_evaluations.as_ref().unwrap(), &table);
        let evaluations: Vec<Vec<Option<BigInt>>> = table
            .iter()
            .map(|row| row.iter().cloned().map(Some).collect())
            .collect();
        assert_eq!(report.evaluations, evaluations);

        let wrong = vec![b(&[2, 0]), b(&[0, 0])];
        let report = bicyclic_sieving_report_bigint(&f, 2, 2, Some(&wrong)).unwrap();
        assert_eq!(report.holds, Some(false));
        assert_eq!(report.mismatches, vec![(1, 1)]);
    }

    fn matrix_product(p: &[BigInt], q: &[BigInt]) -> Vec<Vec<BigInt>> {
        p.iter()
            .map(|a| q.iter().map(|c| a * c).collect())
            .collect()
    }

    #[test]
    fn genuine_bicyclic_csp_from_two_q_binomials() {
        for ((n1, k1), (n2, k2)) in [((4, 2), (6, 3)), ((5, 2), (4, 2)), ((6, 2), (6, 3))] {
            let f = matrix_product(&gaussian_binomial(n1, k1), &gaussian_binomial(n2, k2));
            let fix1 = rotation_fixed_points(n1, k1);
            let fix2 = rotation_fixed_points(n2, k2);
            let table: Vec<Vec<BigInt>> = fix1
                .iter()
                .map(|a| fix2.iter().map(|c| a * c).collect())
                .collect();
            let report = bicyclic_sieving_report_bigint(&f, n1, n2, Some(&table)).unwrap();
            assert!(report.integral);
            assert_eq!(report.holds, Some(true), "({n1},{k1}) x ({n2},{k2})");
            // evaluations factor
            let u1 =
                cyclic_sieving_orbit_decomposition_bigint(&gaussian_binomial(n1, k1), n1).unwrap();
            let u2 =
                cyclic_sieving_orbit_decomposition_bigint(&gaussian_binomial(n2, k2), n2).unwrap();
            assert_eq!(u1.fixed_points.unwrap(), fix1);
            assert_eq!(u2.fixed_points.unwrap(), fix2);
            // a product of univariate polynomials in product form has product coordinates
            if report.product_basis_form {
                assert_eq!(report.product_formula_evaluations.as_ref().unwrap(), &table);
            }
        }
    }

    #[test]
    fn bicyclic_product_coordinates_are_products_of_univariate_ones() {
        let f1 = gaussian_binomial(6, 3);
        let f2 = gaussian_binomial(4, 2);
        let f = matrix_product(&f1, &f2);
        let report = bicyclic_sieving_report_bigint(&f, 6, 4, None).unwrap();
        assert!(report.product_basis_form);
        let a1 = cyclic_sieving_orbit_decomposition_bigint(&f1, 6).unwrap();
        let a2 = cyclic_sieving_orbit_decomposition_bigint(&f2, 4).unwrap();
        let coords = report.product_coordinates.unwrap();
        for (x1, &r1) in report.divisors1.iter().enumerate() {
            for (x2, &r2) in report.divisors2.iter().enumerate() {
                let e1 = a1
                    .orbits
                    .as_ref()
                    .unwrap()
                    .iter()
                    .find(|e| e.basis_index == r1);
                let e2 = a2
                    .orbits
                    .as_ref()
                    .unwrap()
                    .iter()
                    .find(|e| e.basis_index == r2);
                let (e1, e2) = (e1.unwrap(), e2.unwrap());
                assert_eq!(coords[x1][x2], &e1.orbit_count * &e2.orbit_count);
            }
        }
    }

    #[test]
    fn bicyclic_integral_but_not_product_form_when_orders_share_a_factor() {
        // f = sum_{i,j} (zeta_3^{i+j} + zeta_3^{-(i+j)}) q^i t^j: 2 if i + j = 0 mod 3, else -1.
        let f: Vec<Vec<BigInt>> = (0..3)
            .map(|i| {
                (0..3)
                    .map(|j| if (i + j) % 3 == 0 { big(2) } else { big(-1) })
                    .collect()
            })
            .collect();
        let report = bicyclic_sieving_report_bigint(&f, 3, 3, None).unwrap();
        assert!(report.integral);
        assert!(!report.product_basis_form);
        assert!(report.product_coordinates.is_none());
        let values: std::collections::BTreeSet<BigInt> = report
            .evaluations
            .iter()
            .flatten()
            .map(|v| v.clone().unwrap())
            .collect();
        assert_eq!(values, [big(0), big(9)].into_iter().collect());
    }

    #[test]
    fn bicyclic_non_integral_and_bad_input() {
        let f = vec![b(&[1]), b(&[1])];
        let report = bicyclic_sieving_report_bigint(&f, 3, 1, None).unwrap();
        assert!(!report.integral);
        assert!(bicyclic_sieving_report_bigint(&f, 0, 1, None).is_err());
        assert!(bicyclic_sieving_report_bigint(&f, 2, 2, Some(&[b(&[1, 1])])).is_err());
        assert!(bicyclic_sieving_report_bigint(&f, 3, 1, Some(&[b(&[1])])).is_err());
        assert!(bicyclic_sieving_report_bigint(&f, 100, 100, None).is_err());
    }

    #[test]
    fn bicyclic_univariate_case_agrees_with_cyclic() {
        let f = gaussian_binomial(6, 3);
        let matrix: Vec<Vec<BigInt>> = f.iter().map(|c| vec![c.clone()]).collect();
        let report = bicyclic_sieving_report_bigint(&matrix, 6, 1, None).unwrap();
        let dec = cyclic_sieving_orbit_decomposition_bigint(&f, 6).unwrap();
        let col: Vec<BigInt> = report
            .evaluations
            .iter()
            .map(|row| row[0].clone().unwrap())
            .collect();
        assert_eq!(col, dec.fixed_points.unwrap());
    }
}
