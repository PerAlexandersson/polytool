use std::collections::BTreeMap;
use std::fmt;
use std::ops::{Add, Neg, Sub};

use sym_poly_core::{Partition, Ring};

use crate::basis::Basis;
use crate::transition;
use crate::z_coefficient_i64;

/// A symmetric function expressed in a fixed basis with coefficients in a ring C.
///
/// Stores a formal sum Σ c_λ * B_λ where B is one of the six classical bases
/// and λ ranges over partitions.
#[derive(Debug, Clone)]
pub struct SymmetricFunction<C: Ring> {
    basis: Basis,
    pub(crate) terms: BTreeMap<Partition, C>,
}

impl<C: Ring> SymmetricFunction<C> {
    // -----------------------------------------------------------------------
    // Constructors
    // -----------------------------------------------------------------------

    /// Create from a basis and a map of partition -> coefficient. Strips zeros.
    pub fn from_terms(basis: Basis, terms: BTreeMap<Partition, C>) -> Self {
        let mut sf = SymmetricFunction { basis, terms };
        sf.strip_zeros();
        sf
    }

    /// A single basis element B_λ with coefficient 1.
    pub fn basis_element(basis: Basis, partition: Partition) -> Self {
        let mut terms = BTreeMap::new();
        terms.insert(partition, C::one());
        SymmetricFunction { basis, terms }
    }

    /// A single basis element B_λ with given coefficient.
    pub fn scaled_basis_element(basis: Basis, partition: Partition, coeff: C) -> Self {
        if coeff.is_zero() {
            return Self::zero(basis);
        }
        let mut terms = BTreeMap::new();
        terms.insert(partition, coeff);
        SymmetricFunction { basis, terms }
    }

    /// The zero symmetric function.
    pub fn zero(basis: Basis) -> Self {
        SymmetricFunction {
            basis,
            terms: BTreeMap::new(),
        }
    }

    // -----------------------------------------------------------------------
    // Named constructors
    // -----------------------------------------------------------------------

    pub fn monomial_symmetric(partition: Partition) -> Self {
        Self::basis_element(Basis::Monomial, partition)
    }

    pub fn elementary_symmetric(partition: Partition) -> Self {
        Self::basis_element(Basis::Elementary, partition)
    }

    pub fn complete_h_symmetric(partition: Partition) -> Self {
        Self::basis_element(Basis::CompleteH, partition)
    }

    pub fn power_sum_symmetric(partition: Partition) -> Self {
        Self::basis_element(Basis::PowerSum, partition)
    }

    pub fn schur_symmetric(partition: Partition) -> Self {
        Self::basis_element(Basis::Schur, partition)
    }

    pub fn forgotten_symmetric(partition: Partition) -> Self {
        Self::basis_element(Basis::Forgotten, partition)
    }

    // -----------------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------------

    pub fn basis(&self) -> Basis {
        self.basis
    }

    pub fn terms(&self) -> &BTreeMap<Partition, C> {
        &self.terms
    }

    pub fn into_terms(self) -> BTreeMap<Partition, C> {
        self.terms
    }

    pub fn coefficient(&self, partition: &Partition) -> C {
        self.terms.get(partition).cloned().unwrap_or_else(C::zero)
    }

    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    pub fn degree(&self) -> Option<u32> {
        let mut deg = None;
        for p in self.terms.keys() {
            match deg {
                None => deg = Some(p.size()),
                Some(d) if d != p.size() => return None,
                _ => {}
            }
        }
        deg
    }

    pub fn positive_coefficients(&self) -> bool
    where
        C: Ord,
    {
        self.terms.values().all(|c| *c >= C::zero())
    }

    pub(crate) fn strip_zeros(&mut self) {
        self.terms.retain(|_, c| !c.is_zero());
    }

    // -----------------------------------------------------------------------
    // Scalar operations
    // -----------------------------------------------------------------------

    pub fn scale(&self, scalar: &C) -> Self {
        if scalar.is_zero() {
            return Self::zero(self.basis);
        }
        let terms = self
            .terms
            .iter()
            .map(|(p, c)| (p.clone(), c.clone() * scalar.clone()))
            .collect();
        Self::from_terms(self.basis, terms)
    }

    // -----------------------------------------------------------------------
    // Basis conversion
    // -----------------------------------------------------------------------

    pub fn to_basis(&self, target: Basis) -> Self {
        if self.basis == target {
            return self.clone();
        }
        transition::convert(self, target)
    }

    pub fn to_monomial_basis(&self) -> Self {
        self.to_basis(Basis::Monomial)
    }
    pub fn to_elementary_basis(&self) -> Self {
        self.to_basis(Basis::Elementary)
    }
    pub fn to_complete_h_basis(&self) -> Self {
        self.to_basis(Basis::CompleteH)
    }
    pub fn to_power_sum_basis(&self) -> Self {
        self.to_basis(Basis::PowerSum)
    }
    pub fn to_schur_basis(&self) -> Self {
        self.to_basis(Basis::Schur)
    }
    pub fn to_forgotten_basis(&self) -> Self {
        self.to_basis(Basis::Forgotten)
    }

    // -----------------------------------------------------------------------
    // Omega involution
    // -----------------------------------------------------------------------

    pub fn omega_involution(&self) -> Self {
        match self.basis {
            Basis::Elementary => SymmetricFunction {
                basis: Basis::CompleteH,
                terms: self.terms.clone(),
            },
            Basis::CompleteH => SymmetricFunction {
                basis: Basis::Elementary,
                terms: self.terms.clone(),
            },
            Basis::Monomial => SymmetricFunction {
                basis: Basis::Forgotten,
                terms: self.terms.clone(),
            },
            Basis::Forgotten => SymmetricFunction {
                basis: Basis::Monomial,
                terms: self.terms.clone(),
            },
            Basis::Schur => {
                let terms = self
                    .terms
                    .iter()
                    .map(|(p, c)| (p.conjugate_partition(), c.clone()))
                    .collect();
                SymmetricFunction {
                    basis: Basis::Schur,
                    terms,
                }
            }
            Basis::PowerSum => {
                let terms = self
                    .terms
                    .iter()
                    .map(|(p, c)| {
                        let sign_exp = p.size() - p.num_parts() as u32;
                        let sign = if sign_exp % 2 == 0 {
                            C::one()
                        } else {
                            C::minus_one()
                        };
                        (p.clone(), c.clone() * sign)
                    })
                    .collect();
                SymmetricFunction {
                    basis: Basis::PowerSum,
                    terms,
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // Multiplication
    // -----------------------------------------------------------------------

    pub fn multiply(&self, other: &Self) -> Self {
        assert_eq!(
            self.basis, other.basis,
            "multiply requires same basis; convert first"
        );

        if self.is_zero() || other.is_zero() {
            return Self::zero(self.basis);
        }

        if self.basis.is_multiplicative() {
            let mut result_terms = BTreeMap::new();
            for (p1, c1) in &self.terms {
                for (p2, c2) in &other.terms {
                    let joined = p1.partition_join(p2);
                    let coeff = c1.clone() * c2.clone();
                    let entry = result_terms.entry(joined).or_insert_with(C::zero);
                    *entry = entry.clone() + coeff;
                }
            }
            Self::from_terms(self.basis, result_terms)
        } else {
            let bridge = match self.basis {
                Basis::Schur | Basis::Monomial => Basis::CompleteH,
                Basis::Forgotten => Basis::Elementary,
                _ => unreachable!(),
            };
            let a = self.to_basis(bridge);
            let b = other.to_basis(bridge);
            let product = a.multiply(&b);
            product.to_basis(self.basis)
        }
    }

    // -----------------------------------------------------------------------
    // Hall inner product
    // -----------------------------------------------------------------------

    pub fn hall_inner_product(&self, other: &Self) -> C {
        let f = self.to_power_sum_basis();
        let g = other.to_power_sum_basis();
        let mut result = C::zero();
        for (p, cf) in &f.terms {
            if let Some(cg) = g.terms.get(p) {
                let z = C::from_i64(z_coefficient_i64(p));
                result = result + cf.clone() * cg.clone() * z;
            }
        }
        result
    }

    // -----------------------------------------------------------------------
    // Skew Schur functions
    // -----------------------------------------------------------------------

    pub fn skew_schur_function(lambda: &Partition, mu: &Partition) -> Self {
        if !mu.partition_less_equal(lambda) {
            return Self::zero(Basis::CompleteH);
        }
        if mu.is_empty() && lambda.is_empty() {
            return Self::basis_element(Basis::CompleteH, Partition::empty());
        }

        let h_size = lambda.num_parts().max(mu.num_parts());
        let e_size = lambda.part(0).max(mu.part(0)) as usize;

        if e_size < h_size {
            Self::skew_schur_jacobi_trudi_e(lambda, mu)
        } else {
            Self::skew_schur_jacobi_trudi_h(lambda, mu)
        }
    }

    pub(crate) fn skew_schur_jacobi_trudi_h(lambda: &Partition, mu: &Partition) -> Self {
        if mu.is_empty() {
            return Self::schur_symmetric(lambda.clone()).to_complete_h_basis();
        }

        let n = lambda.num_parts().max(mu.num_parts());
        if n == 0 {
            return Self::basis_element(Basis::CompleteH, Partition::empty());
        }

        let perms = all_permutations(n);
        let mut result = Self::zero(Basis::CompleteH);

        for (sigma, sign) in &perms {
            let mut indices = Vec::with_capacity(n);
            let mut valid = true;
            for i in 0..n {
                let Some(val) = jacobi_trudi_index(lambda.part(i), mu.part(sigma[i]), i, sigma[i])
                else {
                    valid = false;
                    break;
                };
                if val > 0 {
                    indices.push(val);
                }
            }
            if !valid {
                continue;
            }
            indices.sort_unstable_by(|a, b| b.cmp(a));
            let partition = Partition::from_sorted(indices);
            let coeff = C::from_i64(*sign as i64);
            let entry = result.terms.entry(partition).or_insert_with(C::zero);
            *entry = entry.clone() + coeff;
        }

        result.strip_zeros();
        result
    }

    pub(crate) fn skew_schur_jacobi_trudi_e(lambda: &Partition, mu: &Partition) -> Self {
        if !mu.partition_less_equal(lambda) {
            return Self::zero(Basis::Elementary);
        }
        if mu.is_empty() {
            return Self::schur_symmetric(lambda.clone()).to_elementary_basis();
        }

        let lam_conj = lambda.conjugate_partition();
        let mu_conj = mu.conjugate_partition();
        let n = lam_conj.num_parts().max(mu_conj.num_parts());
        if n == 0 {
            return Self::basis_element(Basis::Elementary, Partition::empty());
        }

        let perms = all_permutations(n);
        let mut result = Self::zero(Basis::Elementary);

        for (sigma, sign) in &perms {
            let mut indices = Vec::with_capacity(n);
            let mut valid = true;
            for i in 0..n {
                let Some(val) =
                    jacobi_trudi_index(lam_conj.part(i), mu_conj.part(sigma[i]), i, sigma[i])
                else {
                    valid = false;
                    break;
                };
                if val > 0 {
                    indices.push(val);
                }
            }
            if !valid {
                continue;
            }
            indices.sort_unstable_by(|a, b| b.cmp(a));
            let partition = Partition::from_sorted(indices);
            let coeff = C::from_i64(*sign as i64);
            let entry = result.terms.entry(partition).or_insert_with(C::zero);
            *entry = entry.clone() + coeff;
        }

        result.strip_zeros();
        result
    }

    // -----------------------------------------------------------------------
    // Plethysm
    // -----------------------------------------------------------------------

    /// Compute the plethysm `self[inner]`.
    ///
    /// The implementation uses the defining power-sum rules
    /// `p_k[p_mu] = p_{k mu}` and multiplicativity in the outer
    /// power-sum parts. Coefficients are treated as scalars; plethysms with
    /// parameter substitution such as `p_k[q X] = q^k p_k[X]` need a separate
    /// parameter-aware coefficient layer.
    ///
    /// The result is returned in the power-sum basis.
    pub fn plethysm(&self, inner: &Self) -> Self {
        let outer_p = self.to_power_sum_basis();
        let inner_p = inner.to_power_sum_basis();
        let mut result = Self::zero(Basis::PowerSum);

        for (partition, coeff) in &outer_p.terms {
            let mut term =
                Self::scaled_basis_element(Basis::PowerSum, Partition::empty(), coeff.clone());
            for &part in partition.parts() {
                term = term.multiply(&inner_p.plethysm_power_sum(part));
            }
            result = result + term;
        }

        result
    }

    /// Compute `self[p_k]`, the Adams-operator special case of plethysm.
    ///
    /// The result is returned in the power-sum basis.
    pub fn plethysm_power_sum(&self, k: u32) -> Self {
        self.checked_plethysm_power_sum(k)
            .expect("plethysm_power_sum requires k >= 1 and scaled parts to fit in u32")
    }

    /// Compute `self[p_k]`, returning `None` when `k` is zero or scaling a
    /// power-sum part would exceed the `u32` partition-index representation.
    ///
    /// The result is returned in the power-sum basis.
    pub fn checked_plethysm_power_sum(&self, k: u32) -> Option<Self> {
        if k == 0 {
            return None;
        }
        let f_p = self.to_power_sum_basis();
        let mut result_terms = BTreeMap::new();
        for (partition, coeff) in &f_p.terms {
            let new_parts: Vec<u32> = partition
                .parts()
                .iter()
                .map(|&p| p.checked_mul(k))
                .collect::<Option<_>>()?;
            let new_partition = Partition::from_sorted(new_parts);
            let entry = result_terms.entry(new_partition).or_insert_with(C::zero);
            *entry = entry.clone() + coeff.clone();
        }
        Some(Self::from_terms(Basis::PowerSum, result_terms))
    }

    // -----------------------------------------------------------------------
    // Specializations
    // -----------------------------------------------------------------------

    /// Trivial specialization: f(1, 1, ..., 1) with n copies.
    pub fn trivial_specialization(&self, n: u32) -> C {
        let in_schur = self.to_schur_basis();
        let mut result = C::zero();
        for (partition, coeff) in in_schur.terms() {
            let val = schur_trivial_spec::<C>(partition, n);
            result = result + coeff.clone() * val;
        }
        result
    }

    /// Principal specialization: f(1, q, q^2, ..., q^{n-1}).
    pub fn principal_specialization(&self, n: u32) -> Vec<C> {
        let in_schur = self.to_schur_basis();
        let mut result: Vec<C> = vec![];
        for (partition, coeff) in in_schur.terms() {
            let poly = schur_principal_spec::<C>(partition, n);
            let scaled: Vec<C> = poly.iter().map(|c| c.clone() * coeff.clone()).collect();
            result = poly_add::<C>(&result, &scaled);
        }
        while result.last().map_or(false, |c| c.is_zero()) {
            result.pop();
        }
        result
    }

    /// Stable principal specialization: f(1, q, q^2, ...).
    /// Returns (numerator, denominator) polynomials in q.
    pub fn stable_principal_specialization(&self) -> (Vec<C>, Vec<C>) {
        let in_schur = self.to_schur_basis();

        if in_schur.is_zero() {
            return (vec![], vec![C::one()]);
        }

        let mut lcd_mults: std::collections::BTreeMap<u32, usize> =
            std::collections::BTreeMap::new();
        for (partition, _) in in_schur.terms() {
            let mut term_mults: std::collections::BTreeMap<u32, usize> =
                std::collections::BTreeMap::new();
            for (r, c) in partition.diagram_boxes() {
                let h = partition.hook_length(r, c).unwrap();
                *term_mults.entry(h).or_insert(0) += 1;
            }
            for (h, m) in &term_mults {
                let entry = lcd_mults.entry(*h).or_insert(0);
                *entry = (*entry).max(*m);
            }
        }

        let mut denom = vec![C::one()];
        for (&h, &mult) in &lcd_mults {
            let mut factor = vec![C::zero(); h as usize + 1];
            factor[0] = C::one();
            factor[h as usize] = C::minus_one();
            for _ in 0..mult {
                denom = poly_mul::<C>(&denom, &factor);
            }
        }

        let mut numer = vec![];
        for (partition, coeff) in in_schur.terms() {
            let mut term_denom = vec![C::one()];
            for (r, c) in partition.diagram_boxes() {
                let h = partition.hook_length(r, c).unwrap();
                let mut factor = vec![C::zero(); h as usize + 1];
                factor[0] = C::one();
                factor[h as usize] = C::minus_one();
                term_denom = poly_mul::<C>(&term_denom, &factor);
            }

            let quotient = poly_exact_div::<C>(&denom, &term_denom);

            let b = partition_b(partition);
            let mut contrib = vec![C::zero(); b as usize + quotient.len()];
            for (i, c) in quotient.iter().enumerate() {
                contrib[b as usize + i] = c.clone() * coeff.clone();
            }
            numer = poly_add::<C>(&numer, &contrib);
        }

        while numer.last().map_or(false, |c| c.is_zero()) {
            numer.pop();
        }
        while denom.last().map_or(false, |c| c.is_zero()) {
            denom.pop();
        }
        (numer, denom)
    }
}

// ---------------------------------------------------------------------------
// Specialization helpers
// ---------------------------------------------------------------------------

fn partition_b(p: &Partition) -> u32 {
    p.parts()
        .iter()
        .enumerate()
        .map(|(i, &part)| i as u32 * part)
        .sum()
}

fn schur_trivial_spec<C: Ring>(partition: &Partition, n: u32) -> C {
    if partition.is_empty() {
        return C::one();
    }
    if partition.num_parts() > n as usize {
        return C::zero();
    }

    // Cancel the hook-content factors before embedding them in C.  This keeps
    // intermediate integer arithmetic bounded by individual u32-sized factors
    // and lets arbitrary-precision coefficient rings retain their precision.
    let mut numerators = Vec::new();
    let mut denominators = Vec::new();
    for (r, c) in partition.diagram_boxes() {
        let content = c as i64 - r as i64;
        let factor = i64::from(n) + content;
        debug_assert!(factor > 0);
        numerators.push(factor);
        denominators.push(i64::from(partition.hook_length(r, c).unwrap()));
    }

    for denominator in &mut denominators {
        for numerator in &mut numerators {
            let divisor = gcd_i64(*numerator, *denominator);
            *numerator /= divisor;
            *denominator /= divisor;
            if *denominator == 1 {
                break;
            }
        }
    }

    assert!(
        denominators.iter().all(|&factor| factor == 1),
        "hook-content formula did not simplify to an integer"
    );
    numerators
        .into_iter()
        .fold(C::one(), |value, factor| value * C::from_i64(factor))
}

fn gcd_i64(mut a: i64, mut b: i64) -> i64 {
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

fn schur_principal_spec<C: Ring>(partition: &Partition, n: u32) -> Vec<C> {
    if partition.is_empty() {
        return vec![C::one()];
    }
    if partition.num_parts() > n as usize {
        return vec![];
    }

    let b = partition_b(partition);
    let mut numer = vec![C::one()];
    let mut denom = vec![C::one()];

    for (r, c) in partition.diagram_boxes() {
        let content = c as i32 - r as i32;
        let npc = n as i32 + content;
        let hook = partition.hook_length(r, c).unwrap() as i32;
        if npc <= 0 {
            return vec![];
        }
        numer = poly_mul::<C>(&numer, &q_analog::<C>(npc as u32));
        denom = poly_mul::<C>(&denom, &q_analog::<C>(hook as u32));
    }

    let mut poly = poly_exact_div::<C>(&numer, &denom);

    if b > 0 {
        let mut shifted = vec![C::zero(); b as usize];
        shifted.extend(poly);
        poly = shifted;
    }
    poly
}

fn q_analog<C: Ring>(k: u32) -> Vec<C> {
    (0..k).map(|_| C::one()).collect()
}

// ---------------------------------------------------------------------------
// Polynomial arithmetic helpers
// ---------------------------------------------------------------------------

fn poly_add<C: Ring>(a: &[C], b: &[C]) -> Vec<C> {
    let len = a.len().max(b.len());
    (0..len)
        .map(|i| {
            let ai = if i < a.len() { a[i].clone() } else { C::zero() };
            let bi = if i < b.len() { b[i].clone() } else { C::zero() };
            ai + bi
        })
        .collect()
}

fn poly_mul<C: Ring>(a: &[C], b: &[C]) -> Vec<C> {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
    let n = a.len() + b.len() - 1;
    let mut result = vec![C::zero(); n];
    for (i, ai) in a.iter().enumerate() {
        if ai.is_zero() {
            continue;
        }
        for (j, bj) in b.iter().enumerate() {
            result[i + j] = result[i + j].clone() + ai.clone() * bj.clone();
        }
    }
    result
}

fn poly_exact_div<C: Ring>(numer: &[C], denom: &[C]) -> Vec<C> {
    if denom.is_empty() {
        panic!("polynomial division by zero");
    }
    let dn = {
        let mut d = numer.len();
        while d > 0 && numer[d - 1].is_zero() {
            d -= 1;
        }
        d
    };
    let dd = {
        let mut d = denom.len();
        while d > 0 && denom[d - 1].is_zero() {
            d -= 1;
        }
        d
    };
    if dn == 0 {
        return vec![];
    }
    if dd == 0 {
        panic!("polynomial division by zero");
    }
    if dn < dd {
        return vec![];
    }
    let mut rem: Vec<C> = numer.to_vec();
    let dq = dn - dd;
    let mut quot = vec![C::zero(); dq + 1];

    for i in (0..=dq).rev() {
        let lc = rem[i + dd - 1].clone();
        if lc.is_zero() {
            continue;
        }
        let lc_denom = denom[dd - 1].clone();
        let q = if lc_denom == C::one() {
            lc
        } else if lc_denom == C::minus_one() {
            -lc
        } else {
            panic!("poly_exact_div: denominator leading coefficient must be +/-1");
        };
        quot[i] = q.clone();
        for (j, c) in denom[..dd].iter().enumerate() {
            rem[i + j] = rem[i + j].clone() - q.clone() * c.clone();
        }
    }
    while quot.last().map_or(false, |c| c.is_zero()) {
        quot.pop();
    }
    quot
}

/// Generate all permutations of {0, ..., n-1} with their signs.
pub(crate) fn all_permutations(n: usize) -> Vec<(Vec<usize>, i8)> {
    if n == 0 {
        return vec![(vec![], 1)];
    }
    let mut result = Vec::new();
    let mut perm: Vec<usize> = (0..n).collect();
    let mut sign: i8 = 1;
    let mut c = vec![0usize; n];
    result.push((perm.clone(), sign));
    let mut i = 0;
    while i < n {
        if c[i] < i {
            if i % 2 == 0 {
                perm.swap(0, i);
            } else {
                perm.swap(c[i], i);
            }
            sign = -sign;
            result.push((perm.clone(), sign));
            c[i] += 1;
            i = 0;
        } else {
            c[i] = 0;
            i += 1;
        }
    }
    result
}

/// Return the nonnegative Jacobi--Trudi matrix index, if it is valid.
///
/// Partition parts use the full `u32` range, so this calculation must not use
/// a narrower signed type merely to account for the row and column shifts.
fn jacobi_trudi_index(lambda_part: u32, mu_part: u32, row: usize, column: usize) -> Option<u32> {
    let row = i128::try_from(row).expect("Jacobi--Trudi row index exceeds i128");
    let column = i128::try_from(column).expect("Jacobi--Trudi column index exceeds i128");
    let index = i128::from(lambda_part) - i128::from(mu_part) - row + column;
    if index < 0 {
        None
    } else {
        Some(u32::try_from(index).expect("Jacobi--Trudi index exceeds u32"))
    }
}

// ---------------------------------------------------------------------------
// Arithmetic trait impls
// ---------------------------------------------------------------------------

impl<C: Ring> Add for SymmetricFunction<C> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        assert_eq!(self.basis, rhs.basis, "cannot add different bases");
        let mut terms = self.terms;
        for (p, c) in rhs.terms {
            let entry = terms.entry(p).or_insert_with(C::zero);
            *entry = entry.clone() + c;
        }
        Self::from_terms(self.basis, terms)
    }
}

impl<C: Ring> Sub for SymmetricFunction<C> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        self + (-rhs)
    }
}

impl<C: Ring> Neg for SymmetricFunction<C> {
    type Output = Self;
    fn neg(self) -> Self {
        let terms = self.terms.into_iter().map(|(p, c)| (p, -c)).collect();
        SymmetricFunction {
            basis: self.basis,
            terms,
        }
    }
}

impl<C: Ring> PartialEq for SymmetricFunction<C> {
    fn eq(&self, other: &Self) -> bool {
        self.basis == other.basis && self.terms == other.terms
    }
}

impl<C: Ring> Eq for SymmetricFunction<C> {}

// ---------------------------------------------------------------------------
// Display
// ---------------------------------------------------------------------------

impl<C: Ring> fmt::Display for SymmetricFunction<C> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        }
        let sym = self.basis.symbol();
        let mut first = true;
        for (partition, coeff) in &self.terms {
            if !first {
                write!(f, " + ")?;
            }
            first = false;
            let one = C::one();
            let minus_one = C::minus_one();
            if *coeff == one {
                write!(f, "{}[{}]", sym, partition)?;
            } else if *coeff == minus_one {
                write!(f, "-{}[{}]", sym, partition)?;
            } else {
                write!(f, "{}*{}[{}]", coeff, sym, partition)?;
            }
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basis_element() {
        let s: SymmetricFunction<i64> =
            SymmetricFunction::schur_symmetric(Partition::new(vec![3, 1]));
        assert_eq!(s.basis(), Basis::Schur);
        assert_eq!(s.coefficient(&Partition::new(vec![3, 1])), 1);
        assert_eq!(s.coefficient(&Partition::new(vec![2, 2])), 0);
    }

    #[test]
    fn test_add() {
        let s1: SymmetricFunction<i64> =
            SymmetricFunction::schur_symmetric(Partition::new(vec![2, 1]));
        let s2 = SymmetricFunction::scaled_basis_element(Basis::Schur, Partition::new(vec![3]), 2);
        let sum = s1 + s2;
        assert_eq!(sum.coefficient(&Partition::new(vec![2, 1])), 1);
        assert_eq!(sum.coefficient(&Partition::new(vec![3])), 2);
    }

    #[test]
    fn test_multiplicative_product() {
        let e1: SymmetricFunction<i64> =
            SymmetricFunction::elementary_symmetric(Partition::new(vec![2]));
        let e2 = SymmetricFunction::elementary_symmetric(Partition::new(vec![3]));
        let prod = e1.multiply(&e2);
        assert_eq!(prod.coefficient(&Partition::new(vec![3, 2])), 1);
        assert_eq!(prod.terms().len(), 1);
    }

    #[test]
    fn test_omega_schur() {
        let s: SymmetricFunction<i64> =
            SymmetricFunction::schur_symmetric(Partition::new(vec![3, 1]));
        let omega_s = s.omega_involution();
        assert_eq!(omega_s.basis(), Basis::Schur);
        assert_eq!(omega_s.coefficient(&Partition::new(vec![2, 1, 1])), 1);
    }

    #[test]
    fn test_skew_schur_21_1() {
        let skew: SymmetricFunction<i64> = SymmetricFunction::skew_schur_function(
            &Partition::new(vec![2, 1]),
            &Partition::new(vec![1]),
        );
        let in_schur = skew.to_schur_basis();
        assert_eq!(in_schur.coefficient(&Partition::new(vec![2])), 1);
        assert_eq!(in_schur.coefficient(&Partition::new(vec![1, 1])), 1);
        assert_eq!(in_schur.terms().len(), 2);
    }

    #[test]
    fn test_skew_schur_jacobi_trudi_accepts_full_u32_parts() {
        let skew: SymmetricFunction<i64> = SymmetricFunction::skew_schur_function(
            &Partition::new(vec![u32::MAX]),
            &Partition::new(vec![1]),
        );

        assert_eq!(skew.basis(), Basis::CompleteH);
        assert_eq!(skew.coefficient(&Partition::new(vec![u32::MAX - 1])), 1);
    }

    #[test]
    fn test_trivial_specialization() {
        let s: SymmetricFunction<i64> =
            SymmetricFunction::schur_symmetric(Partition::new(vec![2, 1]));
        assert_eq!(s.trivial_specialization(3), 8);
    }

    #[test]
    fn test_trivial_specialization_bigint_does_not_narrow_through_i64() {
        use num_bigint::BigInt;

        let s: SymmetricFunction<BigInt> =
            SymmetricFunction::schur_symmetric(Partition::new(vec![50]));
        let expected = (1u32..=50).fold(BigInt::from(1), |value, k| {
            value * BigInt::from(49 + k) / BigInt::from(k)
        });
        assert_eq!(s.trivial_specialization(50), expected);
    }

    #[test]
    fn test_principal_specialization_s1() {
        let s1: SymmetricFunction<i64> =
            SymmetricFunction::schur_symmetric(Partition::new(vec![1]));
        let poly = s1.principal_specialization(4);
        assert_eq!(poly, vec![1, 1, 1, 1]);
    }

    #[test]
    fn test_stable_principal_specialization_s1() {
        let s1: SymmetricFunction<i64> =
            SymmetricFunction::schur_symmetric(Partition::new(vec![1]));
        let (numer, denom) = s1.stable_principal_specialization();
        assert_eq!(numer, vec![1]);
        assert_eq!(denom, vec![1, -1]);
    }

    #[test]
    fn test_plethysm_p2_on_s2() {
        use num_rational::Ratio;
        type Q = Ratio<i64>;
        let s2: SymmetricFunction<Q> = SymmetricFunction::schur_symmetric(Partition::new(vec![2]));
        let result = s2.plethysm_power_sum(2);
        let half = Q::new(1, 2);
        assert_eq!(result.coefficient(&Partition::new(vec![2, 2])), half);
        assert_eq!(result.coefficient(&Partition::new(vec![4])), half);
    }

    #[test]
    fn test_checked_plethysm_power_sum_rejects_part_overflow() {
        let p_max: SymmetricFunction<i64> =
            SymmetricFunction::power_sum_symmetric(Partition::new(vec![u32::MAX]));

        assert!(p_max.checked_plethysm_power_sum(2).is_none());
    }

    #[test]
    #[should_panic(expected = "scaled parts to fit in u32")]
    fn test_plethysm_power_sum_panics_on_part_overflow() {
        let p_max: SymmetricFunction<i64> =
            SymmetricFunction::power_sum_symmetric(Partition::new(vec![u32::MAX]));

        let _ = p_max.plethysm_power_sum(2);
    }

    #[test]
    fn test_plethysm_power_sum_defining_rule() {
        let p22: SymmetricFunction<i64> =
            SymmetricFunction::power_sum_symmetric(Partition::new(vec![2, 2]));
        let p43: SymmetricFunction<i64> =
            SymmetricFunction::power_sum_symmetric(Partition::new(vec![4, 3]));

        let result = p22.plethysm(&p43);
        assert_eq!(result.basis(), Basis::PowerSum);
        assert_eq!(result.coefficient(&Partition::new(vec![8, 8, 6, 6])), 1);
        assert_eq!(result.terms().len(), 1);
    }

    #[test]
    fn test_plethysm_generalizes_power_sum_special_case() {
        use num_rational::Ratio;
        type Q = Ratio<i64>;

        let s21: SymmetricFunction<Q> =
            SymmetricFunction::schur_symmetric(Partition::new(vec![2, 1]));
        let p3: SymmetricFunction<Q> =
            SymmetricFunction::power_sum_symmetric(Partition::new(vec![3]));

        assert_eq!(s21.plethysm(&p3), s21.plethysm_power_sum(3));
    }

    #[test]
    fn test_plethysm_complete_h_orbit_harmonics_small_cases() {
        use num_rational::Ratio;
        type Q = Ratio<i64>;

        let h2: SymmetricFunction<Q> =
            SymmetricFunction::complete_h_symmetric(Partition::new(vec![2]));
        let h3: SymmetricFunction<Q> =
            SymmetricFunction::complete_h_symmetric(Partition::new(vec![3]));

        let h2_h2 = h2.plethysm(&h2).to_schur_basis();
        assert_eq!(h2_h2.coefficient(&Partition::new(vec![4])), Q::from_i64(1));
        assert_eq!(
            h2_h2.coefficient(&Partition::new(vec![2, 2])),
            Q::from_i64(1)
        );
        assert_eq!(h2_h2.terms().len(), 2);

        let h3_h2 = h3.plethysm(&h2).to_schur_basis();
        assert_eq!(h3_h2.coefficient(&Partition::new(vec![6])), Q::from_i64(1));
        assert_eq!(
            h3_h2.coefficient(&Partition::new(vec![4, 2])),
            Q::from_i64(1)
        );
        assert_eq!(
            h3_h2.coefficient(&Partition::new(vec![2, 2, 2])),
            Q::from_i64(1)
        );
        assert_eq!(h3_h2.terms().len(), 3);

        let h2_h3 = h2.plethysm(&h3).to_schur_basis();
        assert_eq!(h2_h3.coefficient(&Partition::new(vec![6])), Q::from_i64(1));
        assert_eq!(
            h2_h3.coefficient(&Partition::new(vec![4, 2])),
            Q::from_i64(1)
        );
        assert_eq!(h2_h3.terms().len(), 2);

        let foulkes_difference = h3.plethysm(&h2) - h2.plethysm(&h3);
        let schur_difference = foulkes_difference.to_schur_basis();
        assert_eq!(
            schur_difference.coefficient(&Partition::new(vec![2, 2, 2])),
            Q::from_i64(1)
        );
        assert_eq!(schur_difference.terms().len(), 1);
    }
}
