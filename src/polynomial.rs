//! Dense univariate polynomial with generic coefficients.
//!
//! `Polynomial<C>` stores coefficients in ascending degree order: `coeffs[i]` is
//! the coefficient of t^i.  Trailing zeros are stripped automatically.
//!
//! Ring arithmetic, derivatives, powers, evaluation, and reversal delegate to
//! [`crate::vec_poly`], so fixed-width coefficient types (`i64`, `i128`,
//! `Ratio<i64>`) panic on overflow instead of wrapping; see [`CoeffRing`].

use crate::vec_poly;
use std::fmt;
use std::ops::{Add, Mul, Neg, Sub};

/// Trait for polynomial coefficients. Any commutative ring with identity.
///
/// Besides the required constants, the trait has borrowing arithmetic hooks
/// ([`add_ref`](Self::add_ref), [`mul_ref`](Self::mul_ref),
/// [`add_mul_assign`](Self::add_mul_assign), ...) with clone-based defaults.
/// [`crate::vec_poly`] and the arithmetic of [`Polynomial`] use them.
///
/// # Overflow
///
/// The implementations for the fixed-width types `i64`, `i128` and
/// `Ratio<i64>` use checked arithmetic in these hooks and **panic on
/// overflow** (in both debug and release builds) instead of wrapping.  Use
/// `BigInt` or `Ratio<BigInt>` coefficients when values may grow without
/// bound.
pub trait CoeffRing:
    Clone
    + Eq
    + fmt::Debug
    + fmt::Display
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Neg<Output = Self>
    + Sized
{
    /// The additive identity (0).
    fn zero() -> Self;
    /// The multiplicative identity (1).
    fn one() -> Self;
    /// Whether this element is 0.
    fn is_zero(&self) -> bool;
    /// Embed an integer.
    fn from_i64(n: i64) -> Self;

    // -- Borrowing arithmetic hooks --
    //
    // The coefficient-vector routines in `crate::vec_poly` and the arithmetic
    // of `Polynomial<C>` go through these methods.  The defaults clone and use
    // the by-value operators.  Implementations override them to avoid clones
    // (`BigInt`, `BigRational`) or to turn silent wrap-around into a panic
    // (`i64`, `i128`, `Ratio<i64>`).

    /// `self + other` without consuming either operand.
    fn add_ref(&self, other: &Self) -> Self {
        self.clone() + other.clone()
    }
    /// `self - other` without consuming either operand.
    fn sub_ref(&self, other: &Self) -> Self {
        self.clone() - other.clone()
    }
    /// `self * other` without consuming either operand.
    fn mul_ref(&self, other: &Self) -> Self {
        self.clone() * other.clone()
    }
    /// `-self` without consuming the operand.
    fn neg_ref(&self) -> Self {
        -self.clone()
    }
    /// `self += other`.
    fn add_assign_ref(&mut self, other: &Self) {
        *self = self.add_ref(other);
    }
    /// `self -= other`.
    fn sub_assign_ref(&mut self, other: &Self) {
        *self = self.sub_ref(other);
    }
    /// `self += a * b`.
    fn add_mul_assign(&mut self, a: &Self, b: &Self) {
        let product = a.mul_ref(b);
        self.add_assign_ref(&product);
    }
}

// -- CoeffRing impls for standard types --

/// Implements `CoeffRing` for a primitive signed integer type.  All arithmetic
/// hooks use checked operations and panic on overflow instead of wrapping.
macro_rules! impl_checked_int_coeff_ring {
    ($t:ty) => {
        impl CoeffRing for $t {
            fn zero() -> Self {
                0
            }
            fn one() -> Self {
                1
            }
            fn is_zero(&self) -> bool {
                *self == 0
            }
            fn from_i64(n: i64) -> Self {
                n.into()
            }
            fn add_ref(&self, other: &Self) -> Self {
                self.checked_add(*other)
                    .unwrap_or_else(|| overflow_panic(stringify!($t), "addition"))
            }
            fn sub_ref(&self, other: &Self) -> Self {
                self.checked_sub(*other)
                    .unwrap_or_else(|| overflow_panic(stringify!($t), "subtraction"))
            }
            fn mul_ref(&self, other: &Self) -> Self {
                self.checked_mul(*other)
                    .unwrap_or_else(|| overflow_panic(stringify!($t), "multiplication"))
            }
            fn neg_ref(&self) -> Self {
                self.checked_neg()
                    .unwrap_or_else(|| overflow_panic(stringify!($t), "negation"))
            }
        }
    };
}

#[cold]
#[track_caller]
fn overflow_panic(type_name: &str, operation: &str) -> ! {
    panic!("{type_name} overflow in polynomial coefficient {operation}; use BigInt coefficients")
}

impl_checked_int_coeff_ring!(i64);
impl_checked_int_coeff_ring!(i128);

impl CoeffRing for num_bigint::BigInt {
    fn zero() -> Self {
        num_bigint::BigInt::from(0)
    }
    fn one() -> Self {
        num_bigint::BigInt::from(1)
    }
    fn is_zero(&self) -> bool {
        num_traits::Zero::is_zero(self)
    }
    fn from_i64(n: i64) -> Self {
        num_bigint::BigInt::from(n)
    }
    fn add_ref(&self, other: &Self) -> Self {
        self + other
    }
    fn sub_ref(&self, other: &Self) -> Self {
        self - other
    }
    fn mul_ref(&self, other: &Self) -> Self {
        self * other
    }
    fn neg_ref(&self) -> Self {
        -self
    }
    fn add_assign_ref(&mut self, other: &Self) {
        *self += other;
    }
    fn sub_assign_ref(&mut self, other: &Self) {
        *self -= other;
    }
    fn add_mul_assign(&mut self, a: &Self, b: &Self) {
        *self += a * b;
    }
}

impl CoeffRing for num_rational::Ratio<num_bigint::BigInt> {
    fn zero() -> Self {
        num_rational::Ratio::from_integer(num_bigint::BigInt::from(0))
    }
    fn one() -> Self {
        num_rational::Ratio::from_integer(num_bigint::BigInt::from(1))
    }
    fn is_zero(&self) -> bool {
        num_traits::Zero::is_zero(self.numer())
    }
    fn from_i64(n: i64) -> Self {
        num_rational::Ratio::from_integer(num_bigint::BigInt::from(n))
    }
    fn add_ref(&self, other: &Self) -> Self {
        self + other
    }
    fn sub_ref(&self, other: &Self) -> Self {
        self - other
    }
    fn mul_ref(&self, other: &Self) -> Self {
        self * other
    }
    fn neg_ref(&self) -> Self {
        -self
    }
    fn add_assign_ref(&mut self, other: &Self) {
        *self += other;
    }
    fn sub_assign_ref(&mut self, other: &Self) {
        *self -= other;
    }
    fn add_mul_assign(&mut self, a: &Self, b: &Self) {
        *self += a * b;
    }
}

impl CoeffRing for num_rational::Ratio<i64> {
    fn zero() -> Self {
        num_rational::Ratio::new(0, 1)
    }
    fn one() -> Self {
        num_rational::Ratio::new(1, 1)
    }
    fn is_zero(&self) -> bool {
        *self.numer() == 0
    }
    fn from_i64(n: i64) -> Self {
        num_rational::Ratio::from_integer(n)
    }
    fn add_ref(&self, other: &Self) -> Self {
        num_traits::CheckedAdd::checked_add(self, other)
            .unwrap_or_else(|| overflow_panic("Ratio<i64>", "addition"))
    }
    fn sub_ref(&self, other: &Self) -> Self {
        num_traits::CheckedSub::checked_sub(self, other)
            .unwrap_or_else(|| overflow_panic("Ratio<i64>", "subtraction"))
    }
    fn mul_ref(&self, other: &Self) -> Self {
        num_traits::CheckedMul::checked_mul(self, other)
            .unwrap_or_else(|| overflow_panic("Ratio<i64>", "multiplication"))
    }
    fn neg_ref(&self) -> Self {
        let numer = self
            .numer()
            .checked_neg()
            .unwrap_or_else(|| overflow_panic("Ratio<i64>", "negation"));
        num_rational::Ratio::new_raw(numer, *self.denom())
    }
}

// -- FieldRing trait and impls --

/// Trait for coefficient types that form a field (supports exact division).
///
/// Extends [`CoeffRing`] with division, enabling polynomial GCD, exact division,
/// and related algorithms that require inverting leading coefficients.
pub trait FieldRing: CoeffRing {
    /// Exact division. Panics if `other` is zero.
    fn field_div(self, other: Self) -> Self;
}

impl FieldRing for num_rational::Ratio<num_bigint::BigInt> {
    fn field_div(self, other: Self) -> Self {
        self / other
    }
}

impl FieldRing for num_rational::Ratio<i64> {
    fn field_div(self, other: Self) -> Self {
        self / other
    }
}

// ---------------------------------------------------------------------------
// Polynomial type
// ---------------------------------------------------------------------------

/// A dense univariate polynomial over a coefficient ring C.
///
/// Coefficients are stored in ascending order: `coeffs[i]` = coefficient of t^i.
/// The zero polynomial has an empty coefficient vector.
#[derive(Debug, Clone, Eq)]
pub struct Polynomial<C: CoeffRing> {
    coeffs: Vec<C>,
}

impl<C: CoeffRing> Polynomial<C> {
    /// Create a polynomial from coefficients in ascending degree order.
    /// Strips trailing zeros.
    pub fn new(coeffs: Vec<C>) -> Self {
        let mut p = Polynomial { coeffs };
        p.strip_trailing_zeros();
        p
    }

    /// The zero polynomial.
    pub fn zero() -> Self {
        Polynomial { coeffs: vec![] }
    }

    /// The constant polynomial 1.
    pub fn one() -> Self {
        Polynomial {
            coeffs: vec![C::one()],
        }
    }

    /// The variable t (= 0 + 1*t).
    pub fn variable() -> Self {
        Polynomial {
            coeffs: vec![C::zero(), C::one()],
        }
    }

    /// A constant polynomial.
    pub fn constant(c: C) -> Self {
        if c.is_zero() {
            Self::zero()
        } else {
            Polynomial { coeffs: vec![c] }
        }
    }

    /// A monomial c * t^k.
    pub fn monomial(c: C, k: usize) -> Self {
        if c.is_zero() {
            return Self::zero();
        }
        let mut coeffs = vec![C::zero(); k + 1];
        coeffs[k] = c;
        Polynomial { coeffs }
    }

    /// Construct from i64 coefficients.
    pub fn from_i64_coeffs(coeffs: &[i64]) -> Self {
        Self::new(coeffs.iter().map(|&c| C::from_i64(c)).collect())
    }

    // -- Accessors --

    /// The coefficient vector (ascending order, no trailing zeros).
    pub fn coeffs(&self) -> &[C] {
        &self.coeffs
    }

    /// Coefficient of t^k. Returns 0 if k > degree.
    pub fn coeff(&self, k: usize) -> C {
        if k < self.coeffs.len() {
            self.coeffs[k].clone()
        } else {
            C::zero()
        }
    }

    /// The degree, or None for the zero polynomial.
    pub fn degree(&self) -> Option<usize> {
        if self.coeffs.is_empty() {
            None
        } else {
            Some(self.coeffs.len() - 1)
        }
    }

    /// Whether this is the zero polynomial.
    pub fn is_zero(&self) -> bool {
        self.coeffs.is_empty()
    }

    /// The leading coefficient, or None for zero.
    pub fn leading_coefficient(&self) -> Option<C> {
        self.coeffs.last().cloned()
    }

    // -- Operations --

    /// Multiply by a scalar.
    pub fn scale(&self, c: &C) -> Self {
        Self::new(vec_poly::scale(&self.coeffs, c))
    }

    /// Formal derivative dp/dt.
    pub fn derivative(&self) -> Self {
        Self::new(vec_poly::derivative(&self.coeffs))
    }

    /// The power `self^e`, by repeated squaring (`p^0 = 1`).
    pub fn pow(&self, e: usize) -> Self {
        Self::new(vec_poly::pow(&self.coeffs, e))
    }

    /// k-th derivative.
    pub fn nth_derivative(&self, k: usize) -> Self {
        let mut p = self.clone();
        for _ in 0..k {
            p = p.derivative();
        }
        p
    }

    /// Evaluate at a point.
    pub fn evaluate(&self, x: &C) -> C {
        vec_poly::evaluate(&self.coeffs, x)
    }

    /// Shift: return p(t + a).
    pub fn shift(&self, a: &C) -> Self {
        if self.is_zero() {
            return Self::zero();
        }
        let n = self.coeffs.len();
        // Taylor expansion: p(t+a) = Σ p^(k)(a)/k! * t^k
        // We use the matrix method: iterate on coefficients.
        let mut result = self.coeffs.clone();
        for i in (0..n - 1).rev() {
            // result[j] += a * result[j+1] for j = i..n-1
            for j in i..n - 1 {
                let upper = result[j + 1].clone();
                result[j] = result[j].clone() + a.clone() * upper;
            }
        }
        Self::new(result)
    }

    /// Reciprocal polynomial: t^d p(1/t), i.e., reverse the coefficient vector.
    ///
    /// If p(t) = a_0 + a_1 t + ... + a_d t^d, then the reciprocal is
    /// a_d + a_{d-1} t + ... + a_0 t^d.
    ///
    /// Note: if p has a zero constant term (divisible by t), the reciprocal
    /// has lower degree since trailing zeros are stripped.
    pub fn reverse(&self) -> Self {
        Self::new(vec_poly::reciprocal(&self.coeffs))
    }

    /// Reciprocal polynomial with respect to a prescribed degree bound.
    ///
    /// For a polynomial of degree at most `n`, this returns
    ///
    /// ```text
    /// I_n(p)(t) = t^n p(1/t).
    /// ```
    ///
    /// Returns `None` if `self` has degree greater than `n`.
    pub fn reverse_with_degree(&self, n: usize) -> Option<Self> {
        vec_poly::reverse(&self.coeffs, n).map(Self::new)
    }

    /// Dilate: return p(c*t), replacing t with c*t.
    ///
    /// The coefficient of t^k becomes `coeffs[k] * c^k`.
    pub fn dilate(&self, c: &C) -> Self {
        if self.is_zero() {
            return Self::zero();
        }
        let mut power = C::one();
        let coeffs: Vec<C> = self
            .coeffs
            .iter()
            .map(|a| {
                let result = a.clone() * power.clone();
                power = power.clone() * c.clone();
                result
            })
            .collect();
        Self::new(coeffs)
    }

    /// The even Hermite--Biehler part `E`, defined by
    ///
    /// ```text
    /// p(t) = E(t^2) + t O(t^2).
    /// ```
    ///
    /// This returns the polynomial `E`.
    pub fn even_part(&self) -> Self {
        Self::new(self.coeffs.iter().step_by(2).cloned().collect())
    }

    /// The odd Hermite--Biehler part `O`, defined by
    ///
    /// ```text
    /// p(t) = E(t^2) + t O(t^2).
    /// ```
    ///
    /// This returns the polynomial `O`.
    pub fn odd_part(&self) -> Self {
        Self::new(self.coeffs.iter().skip(1).step_by(2).cloned().collect())
    }

    /// Return the Hermite--Biehler decomposition `(E, O)` where
    ///
    /// ```text
    /// p(t) = E(t^2) + t O(t^2).
    /// ```
    pub fn hermite_biehler_decomposition(&self) -> (Self, Self) {
        (self.even_part(), self.odd_part())
    }

    /// Reconstruct a polynomial from its even and odd Hermite--Biehler parts.
    ///
    /// If `p(t) = E(t^2) + t O(t^2)`, then
    /// `Polynomial::from_even_odd_parts(&E, &O) = p`.
    pub fn from_even_odd_parts(even: &Self, odd: &Self) -> Self {
        if even.is_zero() && odd.is_zero() {
            return Self::zero();
        }

        let len = even
            .coeffs
            .len()
            .saturating_mul(2)
            .max(odd.coeffs.len() * 2 + 1);
        let mut coeffs = vec![C::zero(); len];
        for (i, c) in even.coeffs.iter().enumerate() {
            coeffs[2 * i] = c.clone();
        }
        for (i, c) in odd.coeffs.iter().enumerate() {
            coeffs[2 * i + 1] = c.clone();
        }
        Self::new(coeffs)
    }

    /// Check if the polynomial is palindromic (symmetric coefficients).
    ///
    /// A polynomial p(t) of degree d is palindromic if a_i = a_{d-i} for all i,
    /// equivalently t^d p(1/t) = p(t).
    pub fn is_palindromic(&self) -> bool {
        let n = self.coeffs.len();
        for i in 0..n / 2 {
            if self.coeffs[i] != self.coeffs[n - 1 - i] {
                return false;
            }
        }
        true
    }

    /// Stapledon decomposition with respect to a degree bound `n`.
    ///
    /// For a polynomial `p(t)` of degree at most `n`, this returns the unique pair
    /// `(a(t), b(t))` such that
    ///
    /// ```text
    /// p(t) = a(t) + t b(t),
    /// ```
    ///
    /// where `a(t)` is symmetric with center `n/2` and `b(t)` is symmetric with
    /// center `(n-1)/2`.
    ///
    /// Returns `None` if `self` has degree greater than `n`.
    pub fn stapledon_decomposition(&self, n: usize) -> Option<(Self, Self)> {
        let reciprocal = self.reverse_with_degree(n)?;

        let mut a_numerator = vec![C::zero(); n + 2];
        a_numerator[0] = self.coeff(0);
        for (i, coeff) in a_numerator.iter_mut().enumerate().take(n + 1).skip(1) {
            *coeff = self.coeff(i) - reciprocal.coeff(i - 1);
        }
        a_numerator[n + 1] = -self.coeff(0);

        let a = Self::exact_divide_by_one_minus_x(&a_numerator);
        let b = Self::exact_divide_by_one_minus_x(
            &(0..=n)
                .map(|i| reciprocal.coeff(i) - self.coeff(i))
                .collect::<Vec<_>>(),
        );

        Some((a, b))
    }

    fn exact_divide_by_one_minus_x(numerator: &[C]) -> Self {
        if numerator.len() <= 1 {
            return Self::zero();
        }

        let mut quotient = Vec::with_capacity(numerator.len() - 1);
        let mut running = numerator[0].clone();
        quotient.push(running.clone());

        for coeff in numerator
            .iter()
            .skip(1)
            .take(numerator.len().saturating_sub(2))
        {
            running = running + coeff.clone();
            quotient.push(running.clone());
        }

        debug_assert_eq!(
            numerator.last().unwrap().clone() + running,
            C::zero(),
            "numerator should be divisible by 1 - t",
        );

        Self::new(quotient)
    }

    fn strip_trailing_zeros(&mut self) {
        while self.coeffs.last().is_some_and(|c| c.is_zero()) {
            self.coeffs.pop();
        }
    }
}

// ---------------------------------------------------------------------------
// Arithmetic
// ---------------------------------------------------------------------------

impl<C: CoeffRing> PartialEq for Polynomial<C> {
    fn eq(&self, other: &Self) -> bool {
        self.coeffs == other.coeffs
    }
}

impl<C: CoeffRing> Add for Polynomial<C> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self::new(vec_poly::add(&self.coeffs, &rhs.coeffs))
    }
}

impl<C: CoeffRing> Sub for Polynomial<C> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self::new(vec_poly::sub(&self.coeffs, &rhs.coeffs))
    }
}

impl<C: CoeffRing> Neg for Polynomial<C> {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(vec_poly::neg(&self.coeffs))
    }
}

impl<C: CoeffRing> Mul for Polynomial<C> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Self::new(vec_poly::mul(&self.coeffs, &rhs.coeffs))
    }
}

// Implement CoeffRing for Polynomial<C> so polynomials can be used as
// symmetric function coefficients (e.g. for Hall-Littlewood).
impl<C: CoeffRing> CoeffRing for Polynomial<C> {
    fn zero() -> Self {
        Polynomial::zero()
    }
    fn one() -> Self {
        Polynomial::one()
    }
    fn is_zero(&self) -> bool {
        self.coeffs.is_empty()
    }
    fn from_i64(n: i64) -> Self {
        Polynomial::constant(C::from_i64(n))
    }
}

// ---------------------------------------------------------------------------
// Field-dependent operations (GCD, division, monic)
// ---------------------------------------------------------------------------

impl<C: FieldRing> Polynomial<C> {
    /// Polynomial long division. Returns `(quotient, remainder)` such that
    /// `self = quotient * divisor + remainder` with `deg(remainder) < deg(divisor)`.
    ///
    /// Panics if `divisor` is zero.
    pub fn div_rem(&self, divisor: &Self) -> (Self, Self) {
        assert!(!divisor.is_zero(), "division by zero polynomial");
        if self.is_zero() {
            return (Self::zero(), Self::zero());
        }
        let dd = divisor.degree().unwrap();
        let mut rem = self.coeffs.clone();
        let lc_d = divisor.leading_coefficient().unwrap();
        let max_quot_deg = if self.coeffs.len() > dd {
            self.coeffs.len() - 1 - dd
        } else {
            return (Self::zero(), self.clone());
        };
        let mut quot = vec![C::zero(); max_quot_deg + 1];

        while rem.len() > dd {
            let lc_rem = rem.last().cloned().unwrap_or_else(C::zero);
            if lc_rem.is_zero() {
                rem.pop();
                continue;
            }
            let shift = rem.len() - 1 - dd;
            let q = lc_rem.field_div(lc_d.clone());
            quot[shift] = q.clone();
            for (j, c) in divisor.coeffs.iter().enumerate() {
                rem[shift + j] = rem[shift + j].clone() - q.clone() * c.clone();
            }
            rem.pop(); // leading term is now zero
        }
        // Handle rem.len() == dd + 1 (same degree as divisor)
        // Actually the loop condition rem.len() > dd handles this: if rem has
        // dd+1 elements, its degree is dd, and we can still divide once.
        // But rem.len() > dd means rem.len() >= dd+1, so degree >= dd. Correct.
        (Self::new(quot), Self::new(rem))
    }

    /// Polynomial GCD via the Euclidean algorithm. Returns a monic polynomial.
    ///
    /// Returns zero only if both inputs are zero.
    pub fn gcd(&self, other: &Self) -> Self {
        let mut a = self.clone();
        let mut b = other.clone();
        while !b.is_zero() {
            let (_, rem) = a.div_rem(&b);
            a = b;
            b = rem;
        }
        a.make_monic()
    }

    /// Return the squarefree part `p / gcd(p, p')`.
    ///
    /// This removes repeated roots while preserving the overall scaling of `p`.
    /// For nonzero constants, this returns the constant itself.
    pub fn squarefree_part(&self) -> Self {
        if self.is_zero() {
            return Self::zero();
        }
        let dp = self.derivative();
        if dp.is_zero() {
            return self.clone();
        }
        let g = self.gcd(&dp);
        if g.degree() == Some(0) || g.is_zero() {
            self.clone()
        } else {
            self.exact_div(&g)
        }
    }

    /// Alias for [`Self::squarefree_part`].
    pub fn make_squarefree(&self) -> Self {
        self.squarefree_part()
    }

    /// Check whether all roots are simple.
    ///
    /// Returns `false` for the zero polynomial, and `true` for nonzero constants.
    pub fn has_simple_roots(&self) -> bool {
        if self.is_zero() {
            return false;
        }
        let dp = self.derivative();
        if dp.is_zero() {
            return true;
        }
        self.gcd(&dp).degree() == Some(0)
    }

    /// Check whether the polynomial has a repeated root.
    pub fn has_repeated_roots(&self) -> bool {
        !self.is_zero() && !self.has_simple_roots()
    }

    /// Exact division: `self / divisor`, assuming the division is exact.
    ///
    /// Panics if `divisor` does not divide `self` exactly.
    pub fn exact_div(&self, divisor: &Self) -> Self {
        let (q, r) = self.div_rem(divisor);
        assert!(r.is_zero(), "exact_div: nonzero remainder");
        q
    }

    /// Return a monic version of this polynomial (leading coefficient = 1).
    ///
    /// Returns zero for the zero polynomial.
    pub fn make_monic(&self) -> Self {
        match self.leading_coefficient() {
            None => Self::zero(),
            Some(lc) => {
                if lc == C::one() {
                    return self.clone();
                }
                Self::new(
                    self.coeffs
                        .iter()
                        .map(|c| c.clone().field_div(lc.clone()))
                        .collect(),
                )
            }
        }
    }

    /// Lagrange interpolation: construct the unique polynomial of degree < n
    /// passing through the points (x_0, y_0), ..., (x_{n-1}, y_{n-1}).
    ///
    /// Uses Newton's divided differences for numerical stability, then converts
    /// to standard (monomial) form.
    ///
    /// Panics if `points` and `values` have different lengths, or if any two
    /// points coincide.
    pub fn lagrange_interpolation(points: &[C], values: &[C]) -> Self {
        let n = points.len();
        assert_eq!(n, values.len(), "points and values must have same length");
        if n == 0 {
            return Self::zero();
        }
        if n == 1 {
            return Self::constant(values[0].clone());
        }

        // Newton's divided differences
        let mut dd = values.to_vec();
        for j in 1..n {
            for i in (j..n).rev() {
                let num = dd[i].clone() - dd[i - 1].clone();
                let den = points[i].clone() - points[i - j].clone();
                dd[i] = num.field_div(den);
            }
        }

        // Convert Newton form to standard form.
        // Newton form: dd[0] + dd[1]*(t - p[0]) + dd[2]*(t - p[0])*(t - p[1]) + ...
        // Build from highest coefficient down using Horner-like expansion.
        let mut poly = vec![C::zero(); n];
        poly[0] = dd[n - 1].clone();
        for i in (0..n - 1).rev() {
            // poly ← poly * (t - points[i]) + dd[i]
            let mut new_poly = vec![C::zero(); n];
            for j in (0..n).rev() {
                if j > 0 {
                    new_poly[j] = new_poly[j].clone() + poly[j - 1].clone();
                }
                new_poly[j] = new_poly[j].clone() - points[i].clone() * poly[j].clone();
            }
            new_poly[0] = new_poly[0].clone() + dd[i].clone();
            poly = new_poly;
        }

        Self::new(poly)
    }
}

// ---------------------------------------------------------------------------
// Display
// ---------------------------------------------------------------------------

impl<C: CoeffRing> fmt::Display for Polynomial<C> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        }
        let mut terms = Vec::new();
        for (i, c) in self.coeffs.iter().enumerate() {
            if c.is_zero() {
                continue;
            }
            let one = C::one();
            let minus_one = -C::one();
            let term = match (c == &one, c == &minus_one, i) {
                (_, _, 0) => format!("{}", c),
                (true, _, 1) => "t".to_string(),
                (_, true, 1) => "-t".to_string(),
                (_, _, 1) => format!("{}t", c),
                (true, _, e) => format!("t^{}", e),
                (_, true, e) => format!("-t^{}", e),
                (_, _, e) => format!("{}t^{}", c, e),
            };
            terms.push(term);
        }
        if terms.is_empty() {
            return write!(f, "0");
        }
        let mut result = terms[0].clone();
        for term in &terms[1..] {
            if let Some(rest) = term.strip_prefix('-') {
                result.push_str(" - ");
                result.push_str(rest);
            } else {
                result.push_str(" + ");
                result.push_str(term);
            }
        }
        write!(f, "{}", result)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_ops() {
        let p = Polynomial::<i64>::from_i64_coeffs(&[1, 2, 1]); // 1 + 2t + t^2
        let q = Polynomial::<i64>::from_i64_coeffs(&[1, 1]); // 1 + t
        assert_eq!(p.degree(), Some(2));
        assert_eq!(q.degree(), Some(1));
        assert_eq!(p.coeff(0), 1);
        assert_eq!(p.coeff(1), 2);
        assert_eq!(p.coeff(5), 0);
    }

    #[test]
    fn test_add_sub() {
        let p = Polynomial::<i64>::from_i64_coeffs(&[1, 2, 3]);
        let q = Polynomial::<i64>::from_i64_coeffs(&[0, 1, -3, 4]);
        let sum = p.clone() + q.clone();
        assert_eq!(sum, Polynomial::from_i64_coeffs(&[1, 3, 0, 4]));
        let diff = p - q;
        assert_eq!(diff, Polynomial::from_i64_coeffs(&[1, 1, 6, -4]));
    }

    #[test]
    fn test_mul() {
        let p = Polynomial::<i64>::from_i64_coeffs(&[1, 1]); // 1 + t
        let q = p.clone() * p.clone(); // (1+t)^2 = 1 + 2t + t^2
        assert_eq!(q, Polynomial::from_i64_coeffs(&[1, 2, 1]));
    }

    #[test]
    fn test_derivative() {
        let p = Polynomial::<i64>::from_i64_coeffs(&[3, 0, 5, 2]); // 3 + 5t^2 + 2t^3
        let dp = p.derivative(); // 10t + 6t^2
        assert_eq!(dp, Polynomial::from_i64_coeffs(&[0, 10, 6]));
    }

    #[test]
    fn test_nth_derivative() {
        let p = Polynomial::<i64>::from_i64_coeffs(&[1, 3, 3, 1]); // (1+t)^3
        assert_eq!(p.nth_derivative(3), Polynomial::from_i64_coeffs(&[6]));
        assert_eq!(p.nth_derivative(4), Polynomial::zero());
    }

    #[test]
    fn test_evaluate() {
        let p = Polynomial::<i64>::from_i64_coeffs(&[1, 2, 1]); // (1+t)^2
        assert_eq!(p.evaluate(&0), 1);
        assert_eq!(p.evaluate(&1), 4);
        assert_eq!(p.evaluate(&(-1)), 0);
    }

    #[test]
    fn test_display() {
        assert_eq!(
            format!("{}", Polynomial::<i64>::from_i64_coeffs(&[1, 2, 1])),
            "1 + 2t + t^2"
        );
        assert_eq!(
            format!("{}", Polynomial::<i64>::from_i64_coeffs(&[0, -1, 0, 3])),
            "-t + 3t^3"
        );
        assert_eq!(format!("{}", Polynomial::<i64>::zero()), "0");
    }

    #[test]
    fn test_zero_arithmetic() {
        let z = Polynomial::<i64>::zero();
        let p = Polynomial::<i64>::from_i64_coeffs(&[1, 2]);
        assert_eq!(z.clone() + p.clone(), p);
        assert_eq!(z.clone() * p.clone(), z);
        assert!(z.is_zero());
    }

    #[test]
    fn test_polynomial_is_coeff_ring() {
        // Polynomial<i64> implements CoeffRing, so Polynomial<Polynomial<i64>> works
        let inner = Polynomial::<i64>::from_i64_coeffs(&[1, 1]); // 1 + t
        let outer = Polynomial::constant(inner); // constant polynomial whose value is (1+t)
        assert_eq!(outer.degree(), Some(0));
    }

    #[test]
    fn test_shift() {
        // p(t) = t^2, p(t+1) = (t+1)^2 = 1 + 2t + t^2
        let p = Polynomial::<i64>::from_i64_coeffs(&[0, 0, 1]);
        let shifted = p.shift(&1);
        assert_eq!(shifted, Polynomial::from_i64_coeffs(&[1, 2, 1]));
    }

    #[test]
    fn test_reverse() {
        // 1 + 2t + 3t^2 -> 3 + 2t + t^2
        let p = Polynomial::<i64>::from_i64_coeffs(&[1, 2, 3]);
        assert_eq!(p.reverse(), Polynomial::from_i64_coeffs(&[3, 2, 1]));

        // t + 2t^2 -> 2 + t (trailing zero stripped)
        let p = Polynomial::<i64>::from_i64_coeffs(&[0, 1, 2]);
        assert_eq!(p.reverse(), Polynomial::from_i64_coeffs(&[2, 1]));

        assert_eq!(Polynomial::<i64>::zero().reverse(), Polynomial::zero());
    }

    #[test]
    fn test_reverse_with_degree() {
        let p = Polynomial::<i64>::from_i64_coeffs(&[1, 2, 3]);
        assert_eq!(
            p.reverse_with_degree(4),
            Some(Polynomial::from_i64_coeffs(&[0, 0, 3, 2, 1])),
        );
        assert_eq!(p.reverse_with_degree(1), None);
    }

    #[test]
    fn test_dilate() {
        // (1 + t) dilated by 2 -> 1 + 2t
        let p = Polynomial::<i64>::from_i64_coeffs(&[1, 1]);
        assert_eq!(p.dilate(&2), Polynomial::from_i64_coeffs(&[1, 2]));

        // (1 + t^2) dilated by 3 -> 1 + 9t^2
        let p = Polynomial::<i64>::from_i64_coeffs(&[1, 0, 1]);
        assert_eq!(p.dilate(&3), Polynomial::from_i64_coeffs(&[1, 0, 9]));
    }

    #[test]
    fn test_hermite_biehler_decomposition() {
        let p = Polynomial::<i64>::from_i64_coeffs(&[1, 2, 3, 4, 5]);
        let (even, odd) = p.hermite_biehler_decomposition();
        assert_eq!(even, Polynomial::from_i64_coeffs(&[1, 3, 5]));
        assert_eq!(odd, Polynomial::from_i64_coeffs(&[2, 4]));
        assert_eq!(Polynomial::from_even_odd_parts(&even, &odd), p);
    }

    #[test]
    fn test_hermite_biehler_decomposition_zero_and_odd_only() {
        let z = Polynomial::<i64>::zero();
        let (even, odd) = z.hermite_biehler_decomposition();
        assert_eq!(even, Polynomial::zero());
        assert_eq!(odd, Polynomial::zero());

        let p = Polynomial::<i64>::from_i64_coeffs(&[0, 2, 0, 4]);
        let (even, odd) = p.hermite_biehler_decomposition();
        assert_eq!(even, Polynomial::zero());
        assert_eq!(odd, Polynomial::from_i64_coeffs(&[2, 4]));
        assert_eq!(Polynomial::from_even_odd_parts(&even, &odd), p);
    }

    #[test]
    fn test_is_palindromic() {
        assert!(Polynomial::<i64>::from_i64_coeffs(&[1, 2, 1]).is_palindromic());
        assert!(Polynomial::<i64>::from_i64_coeffs(&[1, 11, 11, 1]).is_palindromic());
        assert!(Polynomial::<i64>::from_i64_coeffs(&[1]).is_palindromic());
        assert!(Polynomial::<i64>::zero().is_palindromic());
        assert!(!Polynomial::<i64>::from_i64_coeffs(&[1, 2, 3]).is_palindromic());
    }

    #[test]
    fn test_stapledon_decomposition() {
        let p = Polynomial::<i64>::from_i64_coeffs(&[1, 2, 3]);
        let (a, b) = p.stapledon_decomposition(2).unwrap();
        assert_eq!(a, Polynomial::from_i64_coeffs(&[1, 0, 1]));
        assert_eq!(b, Polynomial::from_i64_coeffs(&[2, 2]));
    }

    #[test]
    fn test_stapledon_decomposition_zero_and_palindromic() {
        let zero = Polynomial::<i64>::zero();
        assert_eq!(
            zero.stapledon_decomposition(4),
            Some((Polynomial::zero(), Polynomial::zero())),
        );

        let p = Polynomial::<i64>::from_i64_coeffs(&[1, 3, 3, 1]);
        let (a, b) = p.stapledon_decomposition(3).unwrap();
        assert_eq!(a, p);
        assert_eq!(b, Polynomial::zero());
    }

    // -- FieldRing tests (using Ratio<i64>) --

    type Q = num_rational::Ratio<i64>;

    fn q_poly(coeffs: &[i64]) -> Polynomial<Q> {
        Polynomial::new(coeffs.iter().map(|&c| Q::from_integer(c)).collect())
    }

    #[test]
    fn test_div_rem() {
        // (t^2 - 1) / (t - 1) = (t + 1) remainder 0
        let f = q_poly(&[-1, 0, 1]);
        let g = q_poly(&[-1, 1]);
        let (q, r) = f.div_rem(&g);
        assert_eq!(q, q_poly(&[1, 1]));
        assert!(r.is_zero());

        // (t^2 + t + 1) / (t + 1) = t, remainder 1
        let f = q_poly(&[1, 1, 1]);
        let g = q_poly(&[1, 1]);
        let (q, r) = f.div_rem(&g);
        assert_eq!(q, q_poly(&[0, 1]));
        assert_eq!(r, q_poly(&[1]));
    }

    #[test]
    fn test_gcd() {
        // gcd((t-1)(t-2), (t-1)(t-3)) = (t-1) (monic)
        let f = q_poly(&[2, -3, 1]); // (t-1)(t-2)
        let g = q_poly(&[3, -4, 1]); // (t-1)(t-3)
        let d = f.gcd(&g);
        assert_eq!(d, q_poly(&[-1, 1])); // t - 1
    }

    #[test]
    fn test_gcd_coprime() {
        // gcd(t-1, t-2) = 1
        let f = q_poly(&[-1, 1]);
        let g = q_poly(&[-2, 1]);
        let d = f.gcd(&g);
        assert_eq!(d, q_poly(&[1])); // constant 1
    }

    #[test]
    fn test_exact_div() {
        // (t-1)(t-2)(t-3) / (t-2) = (t-1)(t-3)
        let f = q_poly(&[-6, 11, -6, 1]);
        let g = q_poly(&[-2, 1]);
        let q = f.exact_div(&g);
        assert_eq!(q, q_poly(&[3, -4, 1]));
    }

    #[test]
    fn test_make_monic() {
        // 2 + 4t -> 1/2 + t  (divide all by lc=4)
        let p = q_poly(&[2, 4]);
        let m = p.make_monic();
        assert_eq!(m.leading_coefficient(), Some(Q::from_integer(1)));
        assert_eq!(m.coeff(0), Q::new(1, 2));

        // already monic
        let p = q_poly(&[-1, 1]);
        assert_eq!(p.make_monic(), p);
    }

    #[test]
    fn test_squarefree_part_and_simple_roots() {
        let repeated = q_poly(&[1, 2, 1]); // (1+t)^2
        assert!(!repeated.has_simple_roots());
        assert!(repeated.has_repeated_roots());
        assert_eq!(repeated.squarefree_part(), q_poly(&[1, 1]));
        assert_eq!(repeated.make_squarefree(), q_poly(&[1, 1]));

        let simple = q_poly(&[1, 0, -1]); // 1 - t^2
        assert!(simple.has_simple_roots());
        assert!(!simple.has_repeated_roots());
        assert_eq!(simple.squarefree_part(), simple);
    }

    #[test]
    fn test_squarefree_part_zero_and_constant() {
        let zero = Polynomial::<Q>::zero();
        assert!(!zero.has_simple_roots());
        assert_eq!(zero.squarefree_part(), zero);

        let constant = q_poly(&[7]);
        assert!(constant.has_simple_roots());
        assert_eq!(constant.squarefree_part(), constant);
    }

    #[test]
    fn test_lagrange_interpolation() {
        // Interpolate through (0,1), (1,3), (2,7) → should give 1 + t + t^2
        let pts: Vec<Q> = vec![Q::from_integer(0), Q::from_integer(1), Q::from_integer(2)];
        let vals: Vec<Q> = vec![Q::from_integer(1), Q::from_integer(3), Q::from_integer(7)];
        let p = Polynomial::lagrange_interpolation(&pts, &vals);
        assert_eq!(p.evaluate(&Q::from_integer(0)), Q::from_integer(1));
        assert_eq!(p.evaluate(&Q::from_integer(1)), Q::from_integer(3));
        assert_eq!(p.evaluate(&Q::from_integer(2)), Q::from_integer(7));
        assert_eq!(p.degree(), Some(2));
        // 1 + t + t^2
        assert_eq!(p.coeff(0), Q::from_integer(1));
        assert_eq!(p.coeff(1), Q::from_integer(1));
        assert_eq!(p.coeff(2), Q::from_integer(1));
    }

    #[test]
    fn test_lagrange_linear() {
        // Two points: (1, 5), (3, 11) → 2 + 3t
        let pts: Vec<Q> = vec![Q::from_integer(1), Q::from_integer(3)];
        let vals: Vec<Q> = vec![Q::from_integer(5), Q::from_integer(11)];
        let p = Polynomial::lagrange_interpolation(&pts, &vals);
        assert_eq!(p, q_poly(&[2, 3]));
    }

    #[test]
    fn test_lagrange_constant() {
        let pts: Vec<Q> = vec![Q::from_integer(42)];
        let vals: Vec<Q> = vec![Q::from_integer(7)];
        let p = Polynomial::lagrange_interpolation(&pts, &vals);
        assert_eq!(p, q_poly(&[7]));
    }
}
