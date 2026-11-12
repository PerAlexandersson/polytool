//! Cyclic orbit decompositions and fixed-point conversions.

use super::*;

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
pub(super) fn mobius(mut n: usize) -> i32 {
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
pub(super) fn gcd_class(j: usize, n: usize) -> usize {
    j.gcd(&n)
}

pub(super) fn big_mobius_term(mu: i32, value: &BigInt) -> BigInt {
    match mu {
        0 => BigInt::zero(),
        1 => value.clone(),
        _ => -value,
    }
}

/// Reduce `f(q)` modulo `q^n - 1`; the result has exactly `n` entries
/// `c_0, ..., c_{n-1}` (ascending order).
///
/// # Panics
///
/// Panics if `n == 0`.
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
    /// The group order `n`.
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
///
/// Returns `Err` only for `order == 0`.  A non-integral polynomial is not an
/// error: it is reported with `integral == false` and `orbits`,
/// `fixed_points`, `cone_condition` set to `None`.
pub fn cyclic_sieving_orbit_decomposition_bigint(
    coeffs: &[BigInt],
    order: usize,
) -> Result<OrbitDecomposition, CyclicSievingError> {
    if order == 0 {
        return Err(CyclicSievingError::NonpositiveOrder);
    }
    let n = order;
    let mut coefficients = coeffs.to_vec();
    strip_trailing_zeros(&mut coefficients);
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
) -> Result<Vec<BigInt>, CyclicSievingError> {
    if order == 0 {
        return Err(CyclicSievingError::NonpositiveOrder);
    }
    let mut merged: BTreeMap<usize, BigInt> = BTreeMap::new();
    for (size, count) in orbit_counts {
        if *size == 0 || order % size != 0 {
            return Err(CyclicSievingError::OrbitSizeDoesNotDivideOrder { size: *size, order });
        }
        if count.sign() == num_bigint::Sign::Minus {
            return Err(CyclicSievingError::NegativeOrbitCount { size: *size });
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
) -> Result<Vec<OrbitSizeEntry>, CyclicSievingError> {
    if order == 0 {
        return Err(CyclicSievingError::NonpositiveOrder);
    }
    let mut counts: BTreeMap<usize, BigInt> = divisors(order)
        .into_iter()
        .map(|k| (k, BigInt::zero()))
        .collect();
    for &size in orbit_sizes {
        match counts.get_mut(&size) {
            Some(count) => *count += 1,
            None => return Err(CyclicSievingError::OrbitSizeDoesNotDivideOrder { size, order }),
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
) -> Result<Vec<BigInt>, CyclicSievingError> {
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
) -> Result<Vec<OrbitSizeEntry>, CyclicSievingError> {
    if order == 0 {
        return Err(CyclicSievingError::NonpositiveOrder);
    }
    let n = order;
    if fixed_points.len() != n {
        return Err(CyclicSievingError::FixedPointCountLength {
            order: n,
            found: fixed_points.len(),
        });
    }
    let divs = divisors(n);
    // F_j for j | n, using fix(g^n) = fix(g^0).
    let f_at = |j: usize| &fixed_points[j % n];
    for (d, value) in fixed_points.iter().enumerate() {
        let g = gcd_class(d, n);
        if value != f_at(g) {
            return Err(CyclicSievingError::FixedPointsNotConstantOnGcdClasses {
                order: n,
                power: d,
                gcd: g,
                value: value.clone(),
                class_value: f_at(g).clone(),
            });
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
            return Err(CyclicSievingError::NegativeOrbitElementCount {
                order: n,
                orbit_size: k,
                element_count: s_k,
            });
        }
        if !(&s_k % &kk).is_zero() {
            return Err(CyclicSievingError::OrbitElementCountNotDivisible {
                order: n,
                orbit_size: k,
                element_count: s_k,
            });
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
    /// The group order `n`.
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
///
/// `orbit_sizes` lists the size of every orbit (a multiset of divisors of
/// `order`); an `Err` is returned for `order == 0` or an orbit size that does not
/// divide `order`.
pub fn check_cyclic_sieving_orbit_sizes(
    coeffs: &[BigInt],
    order: usize,
    orbit_sizes: &[usize],
) -> Result<OrbitSizesCheck, CyclicSievingError> {
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

/// Deprecated name of [`check_cyclic_sieving_orbit_sizes`]; it keeps the former
/// `String` error (the `Display` text of [`CyclicSievingError`]).
#[deprecated(
    since = "0.6.0",
    note = "renamed to `check_cyclic_sieving_orbit_sizes`"
)]
pub fn cyclic_sieving_check_orbit_sizes(
    coeffs: &[BigInt],
    order: usize,
    orbit_sizes: &[usize],
) -> Result<OrbitSizesCheck, String> {
    check_cyclic_sieving_orbit_sizes(coeffs, order, orbit_sizes).map_err(|err| err.to_string())
}
