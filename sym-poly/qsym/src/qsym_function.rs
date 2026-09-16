use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::ops::{Add, Neg, Sub};

use sym_poly_core::{Composition, Ring};

use crate::basis::QSymBasis;

/// A quasisymmetric function expressed in a fixed basis with coefficients in a ring C.
///
/// Stores a formal sum Σ c_α * B_α where B is a QSym basis and α ranges
/// over compositions.
#[derive(Debug, Clone)]
pub struct QSymFunction<C: Ring> {
    basis: QSymBasis,
    terms: BTreeMap<Composition, C>,
}

impl<C: Ring> QSymFunction<C> {
    // -----------------------------------------------------------------------
    // Constructors
    // -----------------------------------------------------------------------

    pub fn from_terms(basis: QSymBasis, terms: BTreeMap<Composition, C>) -> Self {
        let mut f = QSymFunction { basis, terms };
        f.strip_zeros();
        f
    }

    pub fn basis_element(basis: QSymBasis, comp: Composition) -> Self {
        let mut terms = BTreeMap::new();
        terms.insert(comp, C::one());
        QSymFunction { basis, terms }
    }

    pub fn scaled_basis_element(basis: QSymBasis, comp: Composition, coeff: C) -> Self {
        if coeff.is_zero() {
            return Self::zero(basis);
        }
        let mut terms = BTreeMap::new();
        terms.insert(comp, coeff);
        QSymFunction { basis, terms }
    }

    pub fn zero(basis: QSymBasis) -> Self {
        QSymFunction {
            basis,
            terms: BTreeMap::new(),
        }
    }

    // Named constructors
    pub fn monomial_qsym(comp: Composition) -> Self {
        Self::basis_element(QSymBasis::Monomial, comp)
    }

    pub fn fundamental_qsym(comp: Composition) -> Self {
        Self::basis_element(QSymBasis::Fundamental, comp)
    }

    pub fn quasisymmetric_schur(comp: Composition) -> Self {
        Self::basis_element(QSymBasis::QuasisymmetricSchur, comp)
    }

    pub fn dual_immaculate(comp: Composition) -> Self {
        Self::basis_element(QSymBasis::DualImmaculate, comp)
    }

    pub fn extended_schur(comp: Composition) -> Self {
        Self::basis_element(QSymBasis::ExtendedSchur, comp)
    }

    pub fn row_strict_extended_schur(comp: Composition) -> Self {
        Self::basis_element(QSymBasis::RowStrictExtendedSchur, comp)
    }

    pub fn flipped_extended_schur(comp: Composition) -> Self {
        Self::basis_element(QSymBasis::FlippedExtendedSchur, comp)
    }

    pub fn backward_extended_schur(comp: Composition) -> Self {
        Self::basis_element(QSymBasis::BackwardExtendedSchur, comp)
    }

    pub fn combinatorial_power_sum(comp: Composition) -> Self {
        Self::basis_element(QSymBasis::CombinatorialPowerSum, comp)
    }

    pub fn reverse_combinatorial_power_sum(comp: Composition) -> Self {
        Self::basis_element(QSymBasis::ReverseCombinatorialPowerSum, comp)
    }

    // -----------------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------------

    pub fn basis(&self) -> QSymBasis {
        self.basis
    }
    pub fn terms(&self) -> &BTreeMap<Composition, C> {
        &self.terms
    }
    pub fn into_terms(self) -> BTreeMap<Composition, C> {
        self.terms
    }

    pub fn coefficient(&self, comp: &Composition) -> C {
        self.terms.get(comp).cloned().unwrap_or_else(C::zero)
    }

    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    pub fn degree(&self) -> Option<u32> {
        let mut deg = None;
        for c in self.terms.keys() {
            match deg {
                None => deg = Some(c.size()),
                Some(d) if d != c.size() => return None,
                _ => {}
            }
        }
        deg
    }

    fn strip_zeros(&mut self) {
        self.terms.retain(|_, c| !c.is_zero());
    }

    pub fn scale(&self, scalar: &C) -> Self {
        if scalar.is_zero() {
            return Self::zero(self.basis);
        }
        let terms = self
            .terms
            .iter()
            .map(|(c, v)| (c.clone(), v.clone() * scalar.clone()))
            .collect();
        Self::from_terms(self.basis, terms)
    }

    // -----------------------------------------------------------------------
    // Basis conversion
    // -----------------------------------------------------------------------

    pub fn to_basis(&self, target: QSymBasis) -> Self {
        if self.basis == target {
            return self.clone();
        }
        crate::transition::convert(self, target)
    }

    pub fn to_monomial_basis(&self) -> Self {
        self.to_basis(QSymBasis::Monomial)
    }
    pub fn to_fundamental_basis(&self) -> Self {
        self.to_basis(QSymBasis::Fundamental)
    }
    pub fn to_quasisymmetric_schur_basis(&self) -> Self {
        self.to_basis(QSymBasis::QuasisymmetricSchur)
    }
    pub fn to_dual_immaculate_basis(&self) -> Self {
        self.to_basis(QSymBasis::DualImmaculate)
    }
    pub fn to_extended_schur_basis(&self) -> Self {
        self.to_basis(QSymBasis::ExtendedSchur)
    }
    pub fn to_row_strict_extended_schur_basis(&self) -> Self {
        self.to_basis(QSymBasis::RowStrictExtendedSchur)
    }
    pub fn to_flipped_extended_schur_basis(&self) -> Self {
        self.to_basis(QSymBasis::FlippedExtendedSchur)
    }
    pub fn to_backward_extended_schur_basis(&self) -> Self {
        self.to_basis(QSymBasis::BackwardExtendedSchur)
    }
    pub fn to_psi_basis(&self) -> Self {
        self.to_basis(QSymBasis::PowerSumPsi)
    }
    pub fn to_phi_basis(&self) -> Self {
        self.to_basis(QSymBasis::PowerSumPhi)
    }
    pub fn to_combinatorial_power_sum_basis(&self) -> Self {
        self.to_basis(QSymBasis::CombinatorialPowerSum)
    }
    pub fn to_reverse_combinatorial_power_sum_basis(&self) -> Self {
        self.to_basis(QSymBasis::ReverseCombinatorialPowerSum)
    }

    // -----------------------------------------------------------------------
    // Involutions
    // -----------------------------------------------------------------------

    /// The omega involution on QSym.
    ///
    /// On the fundamental basis, ω sends F_α to F_β where β is determined by:
    ///   Des(β) = [n-1] \ { n - s : s ∈ Des(α) }
    ///
    /// Key properties:
    /// - ω(Ψ_α) = (-1)^{n-ℓ(α)} Ψ_{α^r}  (reverse composition, sign)
    /// - ω restricts to the classical omega involution on Sym ⊂ QSym
    /// - ω is an involution: ω² = id
    pub fn omega_involution(&self) -> Self {
        let in_f = self.to_fundamental_basis();
        let mut result_terms: BTreeMap<Composition, C> = BTreeMap::new();

        for (alpha, coeff) in &in_f.terms {
            let beta = omega_on_composition(alpha);
            let entry = result_terms.entry(beta).or_insert_with(C::zero);
            *entry = entry.clone() + coeff.clone();
        }

        let result = Self::from_terms(QSymBasis::Fundamental, result_terms);
        if self.basis == QSymBasis::Fundamental {
            result
        } else {
            result.to_basis(self.basis)
        }
    }

    /// The `rho` involution on QSym.
    ///
    /// On the fundamental basis, `rho` sends `F_alpha` to `F_{rev(alpha)}`.
    pub fn rho_involution(&self) -> Self {
        let in_f = self.to_fundamental_basis();
        let mut result_terms: BTreeMap<Composition, C> = BTreeMap::new();

        for (alpha, coeff) in &in_f.terms {
            let beta = reverse_on_composition(alpha);
            let entry = result_terms.entry(beta).or_insert_with(C::zero);
            *entry = entry.clone() + coeff.clone();
        }

        let result = Self::from_terms(QSymBasis::Fundamental, result_terms);
        if self.basis == QSymBasis::Fundamental {
            result
        } else {
            result.to_basis(self.basis)
        }
    }

    /// The `psi` involution on QSym.
    ///
    /// On the fundamental basis, `psi` sends `F_alpha` to `F_beta`, where
    /// `Des(beta) = [n - 1] \ Des(alpha)`.
    ///
    /// This is unrelated to the type-1 quasisymmetric power-sum basis
    /// represented by [`QSymBasis::PowerSumPsi`].
    pub fn psi_involution(&self) -> Self {
        let in_f = self.to_fundamental_basis();
        let mut result_terms: BTreeMap<Composition, C> = BTreeMap::new();

        for (alpha, coeff) in &in_f.terms {
            let beta = psi_on_composition(alpha);
            let entry = result_terms.entry(beta).or_insert_with(C::zero);
            *entry = entry.clone() + coeff.clone();
        }

        let result = Self::from_terms(QSymBasis::Fundamental, result_terms);
        if self.basis == QSymBasis::Fundamental {
            result
        } else {
            result.to_basis(self.basis)
        }
    }

    // -----------------------------------------------------------------------
    // Multiplication (shuffle product in monomial basis)
    // -----------------------------------------------------------------------

    /// Multiply two QSym functions in the same basis.
    /// For the monomial basis, uses the quasi-shuffle (overlapping shuffle) product.
    pub fn multiply(&self, other: &Self) -> Self {
        assert_eq!(self.basis, other.basis, "multiply requires same basis");

        if self.is_zero() || other.is_zero() {
            return Self::zero(self.basis);
        }

        // Convert to monomial, multiply there, convert back
        let a = self.to_monomial_basis();
        let b = other.to_monomial_basis();

        let mut result_terms: BTreeMap<Composition, C> = BTreeMap::new();
        for (alpha, ca) in &a.terms {
            for (beta, cb) in &b.terms {
                let coeff = ca.clone() * cb.clone();
                // M_α * M_β = Σ_{γ} c_{α,β}^γ M_γ, where the coefficient
                // counts quasi-shuffles.  Accumulate paths directly: distinct
                // paths can yield the same composition, so building a vector of
                // all paths first needlessly allocates duplicate compositions.
                accumulate_quasi_shuffles(
                    alpha.parts(),
                    0,
                    beta.parts(),
                    0,
                    &mut Vec::with_capacity(alpha.num_parts() + beta.num_parts()),
                    &coeff,
                    &mut result_terms,
                );
            }
        }

        let result = Self::from_terms(QSymBasis::Monomial, result_terms);
        if self.basis == QSymBasis::Monomial {
            result
        } else {
            result.to_basis(self.basis)
        }
    }
}

/// Reverse a composition, corresponding to the rho involution on fundamental
/// quasisymmetric functions.
fn reverse_on_composition(alpha: &Composition) -> Composition {
    Composition::new(alpha.parts().iter().rev().copied().collect())
}

/// Apply the omega involution to a composition (via descent sets).
///
/// ω(F_α) = F_β where Des(β) = [n-1] \ { n - s : s ∈ Des(α) }.
fn omega_on_composition(alpha: &Composition) -> Composition {
    let n = alpha.size();
    if n <= 1 {
        return alpha.clone();
    }
    let des = alpha.composition_to_descent_set();
    // n - S
    let flipped: BTreeSet<u32> = des.iter().map(|&s| n - s).collect();
    // [n-1] \ flipped
    let complemented: BTreeSet<u32> = (1..n).filter(|s| !flipped.contains(s)).collect();
    Composition::from_descent_set(&complemented, n)
}

/// Apply the psi involution to a composition (via descent sets).
///
/// `psi(F_alpha) = F_beta` where `Des(beta) = [n - 1] \ Des(alpha)`.
fn psi_on_composition(alpha: &Composition) -> Composition {
    let n = alpha.size();
    if n <= 1 {
        return alpha.clone();
    }
    let des = alpha.composition_to_descent_set();
    let complemented: BTreeSet<u32> = (1..n).filter(|s| !des.contains(s)).collect();
    Composition::from_descent_set(&complemented, n)
}

/// Add all quasi-shuffle paths of the two remaining composition suffixes.
///
/// A quasi-shuffle interleaves parts, with a third option to combine the next
/// part from each side.  The output map records multiplicities, which matters
/// whenever different paths have equal resulting parts.
fn accumulate_quasi_shuffles<C: Ring>(
    a: &[u32],
    ai: usize,
    b: &[u32],
    bi: usize,
    current: &mut Vec<u32>,
    coefficient: &C,
    output: &mut BTreeMap<Composition, C>,
) {
    if ai >= a.len() && bi >= b.len() {
        let entry = output
            .entry(Composition::new(current.clone()))
            .or_insert_with(C::zero);
        *entry = entry.clone() + coefficient.clone();
        return;
    }

    // Option 1: take next from a
    if ai < a.len() {
        current.push(a[ai]);
        accumulate_quasi_shuffles(a, ai + 1, b, bi, current, coefficient, output);
        current.pop();
    }

    // Option 2: take next from b
    if bi < b.len() {
        current.push(b[bi]);
        accumulate_quasi_shuffles(a, ai, b, bi + 1, current, coefficient, output);
        current.pop();
    }

    // Option 3: combine a[ai] + b[bi]
    if ai < a.len() && bi < b.len() {
        current.push(
            a[ai]
                .checked_add(b[bi])
                .expect("quasi-shuffle part sum exceeds u32"),
        );
        accumulate_quasi_shuffles(a, ai + 1, b, bi + 1, current, coefficient, output);
        current.pop();
    }
}

// ---------------------------------------------------------------------------
// Arithmetic trait impls
// ---------------------------------------------------------------------------

impl<C: Ring> Add for QSymFunction<C> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        assert_eq!(self.basis, rhs.basis, "cannot add different bases");
        let mut terms = self.terms;
        for (c, v) in rhs.terms {
            let entry = terms.entry(c).or_insert_with(C::zero);
            *entry = entry.clone() + v;
        }
        Self::from_terms(self.basis, terms)
    }
}

impl<C: Ring> Sub for QSymFunction<C> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        self + (-rhs)
    }
}

impl<C: Ring> Neg for QSymFunction<C> {
    type Output = Self;
    fn neg(self) -> Self {
        let terms = self.terms.into_iter().map(|(c, v)| (c, -v)).collect();
        QSymFunction {
            basis: self.basis,
            terms,
        }
    }
}

impl<C: Ring> PartialEq for QSymFunction<C> {
    fn eq(&self, other: &Self) -> bool {
        self.basis == other.basis && self.terms == other.terms
    }
}
impl<C: Ring> Eq for QSymFunction<C> {}

impl<C: Ring> fmt::Display for QSymFunction<C> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        }
        let sym = self.basis.symbol();
        let mut first = true;
        for (comp, coeff) in &self.terms {
            if !first {
                write!(f, " + ")?;
            }
            first = false;
            if *coeff == C::one() {
                write!(f, "{}[{}]", sym, comp)?;
            } else if *coeff == C::minus_one() {
                write!(f, "-{}[{}]", sym, comp)?;
            } else {
                write!(f, "{}*{}[{}]", coeff, sym, comp)?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigInt;

    /// Deliberately naive reference implementation: it materializes every
    /// path, including duplicate compositions, before counting them.
    fn naive_quasi_shuffles(alpha: &Composition, beta: &Composition) -> Vec<Composition> {
        fn visit(
            a: &[u32],
            ai: usize,
            b: &[u32],
            bi: usize,
            current: &mut Vec<u32>,
            output: &mut Vec<Composition>,
        ) {
            if ai == a.len() && bi == b.len() {
                output.push(Composition::new(current.clone()));
                return;
            }
            if ai < a.len() {
                current.push(a[ai]);
                visit(a, ai + 1, b, bi, current, output);
                current.pop();
            }
            if bi < b.len() {
                current.push(b[bi]);
                visit(a, ai, b, bi + 1, current, output);
                current.pop();
            }
            if ai < a.len() && bi < b.len() {
                current.push(a[ai] + b[bi]);
                visit(a, ai + 1, b, bi + 1, current, output);
                current.pop();
            }
        }

        let mut output = Vec::new();
        visit(
            alpha.parts(),
            0,
            beta.parts(),
            0,
            &mut Vec::new(),
            &mut output,
        );
        output
    }

    fn naive_monomial_product(
        alpha: &Composition,
        beta: &Composition,
    ) -> BTreeMap<Composition, i64> {
        let mut output = BTreeMap::new();
        for gamma in naive_quasi_shuffles(alpha, beta) {
            *output.entry(gamma).or_insert(0) += 1;
        }
        output
    }

    #[test]
    fn test_basic() {
        let m: QSymFunction<i64> = QSymFunction::monomial_qsym(Composition::new(vec![2, 1]));
        assert_eq!(m.basis(), QSymBasis::Monomial);
        assert_eq!(m.coefficient(&Composition::new(vec![2, 1])), 1);
        assert!(!m.is_zero());
    }

    #[test]
    fn test_add() {
        let m1: QSymFunction<i64> = QSymFunction::monomial_qsym(Composition::new(vec![2, 1]));
        let m2: QSymFunction<i64> = QSymFunction::scaled_basis_element(
            QSymBasis::Monomial,
            Composition::new(vec![1, 2]),
            3,
        );
        let sum = m1 + m2;
        assert_eq!(sum.coefficient(&Composition::new(vec![2, 1])), 1);
        assert_eq!(sum.coefficient(&Composition::new(vec![1, 2])), 3);
    }

    #[test]
    fn test_quasi_shuffles_1_1() {
        // QShuffle of (1) and (1): (1,1), (1,1), (2)
        let a = Composition::new(vec![1]);
        let b = Composition::new(vec![1]);
        let shuffles = naive_quasi_shuffles(&a, &b);
        // (1,1) appears twice (a first, b first) and (2) once
        assert_eq!(shuffles.len(), 3);
    }

    #[test]
    fn test_monomial_product_accumulates_quasi_shuffle_multiplicities() {
        // Equal adjacent parts produce several paths with the same output.
        // Empty compositions are the multiplicative identity boundary case.
        let cases = [
            (vec![], vec![]),
            (vec![], vec![2, 1]),
            (vec![2, 1], vec![]),
            (vec![1], vec![1]),
            (vec![1, 1], vec![1]),
            (vec![2, 2], vec![2, 2]),
        ];

        for (left, right) in cases {
            let alpha = Composition::new(left);
            let beta = Composition::new(right);
            let product = QSymFunction::<i64>::monomial_qsym(alpha.clone())
                .multiply(&QSymFunction::monomial_qsym(beta.clone()));
            assert_eq!(product.terms(), &naive_monomial_product(&alpha, &beta));
        }
    }

    #[test]
    fn test_monomial_product_keeps_bigint_coefficients_exact() {
        let alpha = Composition::new(vec![1]);
        let coefficient: BigInt = BigInt::from(1_u8) << 100_u32;
        let product = QSymFunction::scaled_basis_element(
            QSymBasis::Monomial,
            alpha.clone(),
            coefficient.clone(),
        )
        .multiply(&QSymFunction::scaled_basis_element(
            QSymBasis::Monomial,
            alpha,
            coefficient.clone(),
        ));
        let square = coefficient.clone() * coefficient;

        assert_eq!(
            product.coefficient(&Composition::new(vec![1, 1])),
            BigInt::from(2_u8) * square.clone()
        );
        assert_eq!(product.coefficient(&Composition::new(vec![2])), square);
    }

    #[test]
    fn test_multiply_m1_m1() {
        // M_(1) * M_(1) = 2*M_(1,1) + M_(2)
        let m1: QSymFunction<i64> = QSymFunction::monomial_qsym(Composition::new(vec![1]));
        let prod = m1.multiply(&m1.clone());
        assert_eq!(prod.coefficient(&Composition::new(vec![1, 1])), 2);
        assert_eq!(prod.coefficient(&Composition::new(vec![2])), 1);
    }

    #[test]
    fn test_psi_involution_on_fundamental() {
        let f21: QSymFunction<i64> = QSymFunction::fundamental_qsym(Composition::new(vec![2, 1]));
        let psi_f21 = f21.psi_involution();
        assert_eq!(psi_f21.basis(), QSymBasis::Fundamental);
        assert_eq!(psi_f21.coefficient(&Composition::new(vec![1, 2])), 1);
        assert_eq!(psi_f21.terms().len(), 1);

        let f3: QSymFunction<i64> = QSymFunction::fundamental_qsym(Composition::new(vec![3]));
        let psi_f3 = f3.psi_involution();
        assert_eq!(psi_f3.coefficient(&Composition::new(vec![1, 1, 1])), 1);
        assert_eq!(psi_f3.terms().len(), 1);
    }

    #[test]
    fn test_psi_involution_roundtrip() {
        let f = QSymFunction::scaled_basis_element(
            QSymBasis::Fundamental,
            Composition::new(vec![1, 2]),
            2i64,
        ) + QSymFunction::fundamental_qsym(Composition::new(vec![3]));

        assert_eq!(f.clone().psi_involution().psi_involution(), f);
    }

    #[test]
    fn test_rho_involution_on_fundamental() {
        let f21: QSymFunction<i64> = QSymFunction::fundamental_qsym(Composition::new(vec![2, 1]));
        let rho_f21 = f21.rho_involution();
        assert_eq!(rho_f21.basis(), QSymBasis::Fundamental);
        assert_eq!(rho_f21.coefficient(&Composition::new(vec![1, 2])), 1);
        assert_eq!(rho_f21.terms().len(), 1);

        let f111: QSymFunction<i64> =
            QSymFunction::fundamental_qsym(Composition::new(vec![1, 1, 1]));
        let rho_f111 = f111.rho_involution();
        assert_eq!(rho_f111.coefficient(&Composition::new(vec![1, 1, 1])), 1);
        assert_eq!(rho_f111.terms().len(), 1);
    }

    #[test]
    fn test_rho_involution_roundtrip() {
        let f = QSymFunction::scaled_basis_element(
            QSymBasis::Fundamental,
            Composition::new(vec![1, 2]),
            2i64,
        ) + QSymFunction::fundamental_qsym(Composition::new(vec![3]));

        assert_eq!(f.clone().rho_involution().rho_involution(), f);
    }

    #[test]
    fn test_involution_relations() {
        for degree in 1..=5 {
            for alpha in Composition::integer_compositions(degree) {
                let f: QSymFunction<i64> = QSymFunction::fundamental_qsym(alpha);
                assert_eq!(
                    f.clone().psi_involution(),
                    f.clone().omega_involution().rho_involution()
                );
                assert_eq!(
                    f.clone().psi_involution(),
                    f.clone().rho_involution().omega_involution()
                );
            }
        }
    }
}
