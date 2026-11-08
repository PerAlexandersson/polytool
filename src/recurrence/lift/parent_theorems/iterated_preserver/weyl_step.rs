//! Brändén's Weyl step: the second-order operator
//! `T = A(x) + B(x) D + C(x) D^2` as a real-rootedness preserver.
//!
//! # The theorem
//!
//! P. Brändén, *Linear transformations preserving the Pólya frequency
//! property*, Trans. Amer. Math. Soc. 358 (2006), Corollary 3.5, for
//! `F(x, z) = A(x) + B(x) z + C(x) z^2` and `phi_F = A + B D + C D^2`:
//! if `A != 0` and
//!
//! 1. for every real `xi`, `F(xi, z)` is real-rooted in `z`, that is
//!    `B^2 - 4 A C >= 0` on the real line;
//! 2. `A` strictly interlaces `B` (`deg B = deg A + 1`, the roots of `A`
//!    strictly between those of `B`) or strictly alternates left of `B`
//!    (`deg B = deg A`, `a_1 < b_1 < a_2 < ... < a_k < b_k` with the roots in
//!    increasing order), and `deg A = 0` or `A` and `B` have leading
//!    coefficients of the same sign;
//! 3. `phi_F(x^k) = A x^k + k B x^(k-1) + k (k - 1) C x^(k-2)` has degree
//!    `deg A + k` for `0 <= k <= d`, with leading coefficients of one sign;
//!
//! then `phi_F(f)` is real-rooted for every real-rooted `f` of degree at
//! most `d`.  Here 3 is required for every `k >= 0`, so `T` maps every
//! real-rooted polynomial to a real-rooted polynomial (of degree
//! `deg A + deg f`, so never to zero).  By Corollary 2.5 of the same paper,
//! a linear operator that preserves real-rootedness also preserves the
//! alternation of two polynomials.
//!
//! Conversely, by J. Borcea and P. Brändén, *Multivariate Pólya--Schur
//! classification problems in the Weyl algebra*, Proc. London Math. Soc.
//! 101 (2010), Theorem 1.3, an operator `T` of this form that preserves
//! real-rootedness has a real stable symbol `F(z, -w)`; in particular
//! `B^2 - 4 A C >= 0` on the real line.  A point where the discriminant is
//! negative therefore shows that `T` is not a global real-rootedness
//! preserver (it may still preserve a half line, as Brändén's Euler step
//! does).
//!
//! # What is checked, exactly
//!
//! For a fixed operator ([`check_weyl_step_preserver`], and the iterated
//! route when the coefficients are free of `n`): hypothesis 1 by the
//! square-free decomposition of the discriminant (every factor of odd
//! multiplicity has no real root, Sturm count, and the leading coefficient
//! is positive; a failure gives a rational point with a negative
//! discriminant), hypothesis 2 by the Bézout test of strict interlacing (for
//! alternation, `A` strictly interlaces `(x + M) B` with `-M` below every
//! root of `A`), hypothesis 3 by the sign checker in `k` on
//! `A_a + k B_(a+1) + k (k - 1) C_(a+2)` (`a = deg A`).
//!
//! For index-dependent coefficients (the iterated route): hypothesis 1 when
//! the discriminant is free of `x` (its sign in `n`), quadratic in `x`
//! (`q_2, q_0 >= 0` and `4 q_2 q_0 - q_1^2 >= 0`) or `gamma(n) w(x)` (the sign
//! of `gamma` and the fixed `w`), each by the sign checker; hypothesis 2 for
//! `deg A = 0` (`A`, the leading coefficient of `B` nonzero), `deg A = 1`
//! (the leading coefficients' product positive and the root comparison
//! `A_0 B_1 - A_1 B_0 > 0`, respectively `B_2 (B_2 A_0^2 - B_1 A_0 A_1 + B_0
//! A_1^2) < 0` for `deg B = 2`, by the sign checker) or `A = alpha(n) a(x)`,
//! `B = beta(n) b(x)` (the fixed pair exactly, `alpha beta` by the sign
//! checker); hypothesis 3 by `A_a B_(a+1) >= 0` and `A_a C_(a+2) >= 0` (then
//! `sgn(A_a) phi(x^k)` has leading coefficient `|A_a| + k (...) + k (k - 1)
//! (...) > 0`).  Other shapes stay unresolved, unless the exact check at the
//! first index refutes the hypothesis there.
//!
//! # Root regions
//!
//! When every coefficient of `A`, `B` and `C` has one sign (all `>= 0` or
//! all `<= 0` for every `n`), `T` maps polynomials with coefficients of one
//! sign to such polynomials; a real-rooted polynomial with coefficients of
//! one sign has every root in `(-inf, 0]`.  So the operator then maps every
//! polynomial that is zero or real-rooted with every root in `(-inf, 0]` to
//! such a polynomial.  For `[0, inf)` the same test runs on the reflected
//! operator `A(-x) - B(-x) D + C(-x) D^2`.

use super::*;
use crate::integer_range_sign::IntegerRange;
use crate::recurrence::implicit::algebra::biv_text;
use crate::sturm::SturmChain;

pub(super) const BRANDEN_WEYL_STEP: &str = "Branden's Weyl step: T = A + B D + C D^2 maps every \
     real-rooted polynomial to a real-rooted polynomial when (i) B^2 - 4 A C >= 0 on the real \
     line, (ii) A strictly interlaces B (deg B = deg A + 1) or strictly alternates left of B \
     (deg B = deg A), with leading coefficients of the same sign when deg A >= 1, and (iii) \
     T(x^k) has degree deg A + k with leading coefficients of one sign for every k >= 0 (P. \
     Branden, Linear transformations preserving the Polya frequency property, Trans. Amer. \
     Math. Soc. 358 (2006), Corollary 3.5); (i) is necessary for a real-rootedness preserver \
     (J. Borcea and P. Branden, Proc. London Math. Soc. 101 (2010), Theorem 1.3); with \
     coefficients of one sign T also keeps every root in (-inf, 0]";

/// The verdict of [`check_weyl_step_preserver`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WeylStepPreserverVerdict {
    /// Every hypothesis of Brändén's Corollary 3.5 holds (for every
    /// `k >= 0`): the operator maps every real-rooted polynomial to a
    /// real-rooted polynomial.
    GlobalPreserver,
    /// `B^2 - 4 A C < 0` at `point`: by Borcea and Brändén (2010),
    /// Theorem 1.3, the operator does not preserve real-rootedness on the
    /// whole line.
    NotGlobalPreserver {
        /// A rational point with a negative discriminant.
        point: BigRational,
        /// `B^2 - 4 A C` at `point`.
        discriminant: BigRational,
    },
    /// The discriminant is nonnegative, but Corollary 3.5 does not apply
    /// (hypothesis (ii) or (iii) fails), so nothing is claimed.
    Undecided {
        /// The hypothesis that fails.
        reason: String,
    },
}

/// Decide whether `T = A + B D + C D^2` (coefficients ascending in `x`:
/// `constant = A`, `first = B`, `second = C`) is a real-rootedness preserver
/// on the whole real line by Brändén's Corollary 3.5, or refute it by the
/// necessary condition `B^2 - 4 A C >= 0` of Borcea and Brändén's
/// Theorem 1.3 (module documentation of the iterated-preserver theorem,
/// "Brändén's Weyl step").  Every test is exact.
///
/// # Example
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::recurrence::{check_weyl_step_preserver, BigRational, WeylStepPreserverVerdict};
///
/// let q = |v: i64| BigRational::from_integer(BigInt::from(v));
/// // The unsigned 2-Lah step (4 + x) + (4 + 2x) D + x D^2: discriminant 16.
/// let verdict = check_weyl_step_preserver(&[q(4), q(1)], &[q(4), q(2)], &[q(0), q(1)]);
/// assert_eq!(verdict, WeylStepPreserverVerdict::GlobalPreserver);
/// // The central factorial step (1 + x) + 3x D + x^2 D^2: discriminant x^2 (5 - 4x).
/// let verdict = check_weyl_step_preserver(&[q(1), q(1)], &[q(0), q(3)], &[q(0), q(0), q(1)]);
/// assert!(matches!(verdict, WeylStepPreserverVerdict::NotGlobalPreserver { .. }));
/// ```
pub fn check_weyl_step_preserver(
    constant: &[BigRational],
    first: &[BigRational],
    second: &[BigRational],
) -> WeylStepPreserverVerdict {
    let (a, b, c) = (
        vec_poly::trim(constant),
        vec_poly::trim(first),
        vec_poly::trim(second),
    );
    let discriminant = fixed_discriminant(&a, &b, &c);
    if let Err(point) = nonnegative_on_line(&discriminant) {
        return match point {
            Some(point) => WeylStepPreserverVerdict::NotGlobalPreserver {
                discriminant: vec_poly::evaluate(&discriminant, &point),
                point,
            },
            None => WeylStepPreserverVerdict::Undecided {
                reason: "the discriminant B^2 - 4AC is negative somewhere, but no rational \
                         witness point was isolated"
                    .into(),
            },
        };
    }
    if let Err(reason) = fixed_interlacing(&a, &b) {
        return WeylStepPreserverVerdict::Undecided { reason };
    }
    match fixed_leading_coefficients(&a, &b, &c) {
        Ok(decision) if decision.holds() => WeylStepPreserverVerdict::GlobalPreserver,
        Ok(decision) => WeylStepPreserverVerdict::Undecided {
            reason: decision.describe("the leading coefficient of T(x^k)"),
        },
        Err(reason) => WeylStepPreserverVerdict::Undecided { reason },
    }
}

// ---------------------------------------------------------------------------
// Fixed polynomials
// ---------------------------------------------------------------------------

/// `B^2 - 4 A C`.
fn fixed_discriminant(a: &[Q], b: &[Q], c: &[Q]) -> Vec<Q> {
    vec_poly::trim(&vec_poly::sub(
        &vec_poly::mul(b, b),
        &vec_poly::scale(&vec_poly::mul(a, c), &Q::from_integer(BigInt::from(4))),
    ))
}

fn polynomial(p: &[Q]) -> Polynomial<Q> {
    Polynomial::new(vec_poly::trim(p))
}

/// Whether `p >= 0` on the real line (exact); `Err(Some(x))` gives a
/// rational point with `p(x) < 0`.
fn nonnegative_on_line(p: &[Q]) -> Result<(), Option<Q>> {
    let p = vec_poly::trim(p);
    if vec_poly::is_zero(&p) {
        return Ok(());
    }
    let leading_positive = p.last().is_some_and(Signed::is_positive);
    // Square-free decomposition: the product of the factors of odd
    // multiplicity.
    let whole = polynomial(&p);
    let mut rest = whole.gcd(&polynomial(&vec_poly::derivative(&p)));
    let square_free = whole.div_rem(&rest).0;
    let mut current = square_free.clone();
    let mut odd = Polynomial::new(vec![Q::one()]);
    let mut multiplicity = 1usize;
    while current.coeffs().len() > 1 {
        let next = current.gcd(&rest);
        let factor = current.div_rem(&next).0;
        if multiplicity % 2 == 1 {
            odd = odd * factor;
        }
        rest = rest.div_rem(&next).0;
        current = next;
        multiplicity += 1;
    }
    let odd_roots = odd.coeffs().len() > 1
        && SturmChain::from_bigint_coeffs(&integer_coefficients(odd.coeffs())).count_real_roots()
            > 0;
    if leading_positive && !odd_roots {
        return Ok(());
    }
    Err(negative_point(&p, square_free.coeffs()))
}

/// A rational point where `p < 0`: sample points outside the isolating
/// intervals of the real roots of the square-free part.
fn negative_point(p: &[Q], square_free: &[Q]) -> Option<Q> {
    let integer = integer_coefficients(square_free);
    let intervals = if integer.len() > 1 {
        SturmChain::from_bigint_coeffs(&integer)
            .isolate_roots(&Q::new(BigInt::one(), BigInt::from(1024)))
    } else {
        Vec::new()
    };
    let mut samples = Vec::new();
    match (intervals.first(), intervals.last()) {
        (Some((low, _)), Some((_, high))) => {
            samples.push(low - Q::one());
            samples.push(high + Q::one());
            for pair in intervals.windows(2) {
                samples.push((&pair[0].1 + &pair[1].0) / Q::from_integer(BigInt::from(2)));
            }
        }
        _ => samples.push(Q::zero()),
    }
    samples
        .into_iter()
        .find(|x| vec_poly::evaluate(p, x).is_negative())
}

/// Brändén's hypothesis (ii) for fixed `A`, `B`.
fn fixed_interlacing(a: &[Q], b: &[Q]) -> Result<(), String> {
    let (Some(da), Some(db)) = (vec_poly::degree(a), vec_poly::degree(b)) else {
        return Err("A or B is zero".into());
    };
    if da == 0 {
        return if db <= 1 {
            Ok(())
        } else {
            Err(format!(
                "deg A = 0 and deg B = {db}: A neither interlaces nor alternates left of B"
            ))
        };
    }
    if db != da && db != da + 1 {
        return Err(format!(
            "deg A = {da}, deg B = {db}: neither deg B = deg A nor deg B = deg A + 1"
        ));
    }
    if (&a[da] * &b[db]).is_negative() {
        return Err("A and B have leading coefficients of opposite signs".into());
    }
    let target = if db == da + 1 {
        b.to_vec()
    } else {
        // A strictly alternates left of B iff A strictly interlaces
        // (x + M) B, with -M below every root of A (Cauchy bound + 1).
        let bound = a[..da]
            .iter()
            .map(|c| (c / &a[da]).abs())
            .fold(Q::zero(), |m, c| if c > m { c } else { m });
        vec_poly::mul(b, &[bound + Q::from_integer(BigInt::from(2)), Q::one()])
    };
    let (ia, it) = (integer_coefficients(a), integer_coefficients(&target));
    let holds = crate::is_real_rooted_bigint_coeffs(&ia)
        && crate::is_real_rooted_bigint_coeffs(&it)
        && crate::check_interlacing_bigint_coeffs(&ia, &it) == Some(true);
    if holds {
        Ok(())
    } else if db == da + 1 {
        Err("A does not strictly interlace B (exact Bezout test)".into())
    } else {
        Err("A does not strictly alternate left of B (exact Bezout test)".into())
    }
}

/// Brändén's hypothesis (iii) for fixed `A`, `B`, `C`: the sign decision of
/// `sgn(A_a) (A_a + k B_(a+1) + k (k - 1) C_(a+2)) > 0` for every `k >= 0`,
/// or why the degrees exceed `deg A + k`.
fn fixed_leading_coefficients(a: &[Q], b: &[Q], c: &[Q]) -> Result<IndexSignDecision, String> {
    let da = vec_poly::degree(a).ok_or("A is zero")?;
    if vec_poly::degree(b).is_some_and(|db| db > da + 1)
        || vec_poly::degree(c).is_some_and(|dc| dc > da + 2)
    {
        return Err(
            "deg B > deg A + 1 or deg C > deg A + 2: T(x^k) has degree above deg A + k".into(),
        );
    }
    let lead_a = a[da].clone();
    let b_top = b.get(da + 1).cloned().unwrap_or_else(Q::zero);
    let c_top = c.get(da + 2).cloned().unwrap_or_else(Q::zero);
    let sign = if lead_a.is_negative() {
        -Q::one()
    } else {
        Q::one()
    };
    let lead_k = vec_poly::scale(&[lead_a, &b_top - &c_top, c_top], &sign);
    Ok(decide_index_rational_function_sign(
        &lead_k,
        &[Q::one()],
        &IntegerRange::starting_at(0),
        IndexSignCondition::Positive,
    ))
}

// ---------------------------------------------------------------------------
// The hypotheses of the iterated route
// ---------------------------------------------------------------------------

/// The coefficients `(A, B, C)` of `numerator = A + B D + C D^2` as
/// bivariate matrices, or `None` when the numerator has another shape.
pub(super) fn weyl_step_parts(numerator: &DOp) -> Option<[Bivariate; 3]> {
    let mut parts: [Bivariate; 3] = Default::default();
    for ((shift, order), c) in &numerator.0 {
        if *shift != 0 || *order > 2 {
            return None;
        }
        parts[*order] = c.rows().to_vec();
    }
    (numerator.max_derivative() >= 1).then_some(parts)
}

fn is_index_free(b: &Bivariate) -> bool {
    b.iter().skip(1).all(|row| row.iter().all(Zero::is_zero))
}

fn fixed_row(b: &Bivariate) -> Vec<Q> {
    vec_poly::trim(b.first().map_or(&[][..], Vec::as_slice))
}

fn biv(b: &Bivariate) -> Biv {
    Biv::new(b.clone())
}

fn function(p: Vec<Q>) -> IndexFunction {
    IndexFunction::polynomial(if p.is_empty() { vec![Q::zero()] } else { p })
}

fn product(a: &[Q], b: &[Q]) -> IndexFunction {
    function(vec_poly::mul(a, b))
}

/// The text of `A + B D + C D^2` with its coefficients.
pub(super) fn weyl_step_text(parts: &[Bivariate; 3]) -> String {
    let text = |b: &Bivariate| {
        let text = biv_text(&biv(b), "x");
        if text.is_empty() {
            "0".to_string()
        } else {
            text
        }
    };
    format!(
        "A + B D + C D^2 with A = {}, B = {}, C = {}",
        text(&parts[0]),
        text(&parts[1]),
        text(&parts[2])
    )
}

/// Brändén's hypotheses (i)--(iii) for `S_n = numerator / delta_S(n)` and
/// every `n >= n_first` (the common factor `1/delta_S(n)` changes no sign
/// condition).
pub(super) fn weyl_step_hypotheses(parts: &[Bivariate; 3], n_first: usize) -> Vec<LiftHypothesis> {
    if parts.iter().all(is_index_free) {
        let [a, b, c] = parts.clone().map(|p| fixed_row(&p));
        return fixed_hypotheses(&a, &b, &c);
    }
    let mut out = discriminant_hypotheses(parts, n_first);
    out.extend(interlacing_hypotheses(parts, n_first));
    out.extend(leading_coefficient_hypotheses(parts, n_first));
    out
}

const DISCRIMINANT_STATEMENT: &str = "(i) B^2 - 4 A C >= 0 on the real line (the symbol A + B z \
     + C z^2 is real-rooted in z at every real x)";

fn not_a_global_preserver(point: &Q, value: &Q) -> String {
    format!(
        "B^2 - 4 A C = {} < 0 at x = {}: T is not a real-rootedness preserver on the whole line \
         (Borcea--Branden 2010, Theorem 1.3, necessary condition)",
        format_rational_coeff(value),
        format_rational_coeff(point)
    )
}

fn fixed_hypotheses(a: &[Q], b: &[Q], c: &[Q]) -> Vec<LiftHypothesis> {
    let discriminant = fixed_discriminant(a, b, c);
    let mut out = vec![exact_hypothesis(
        "weyl_step_discriminant_nonnegative",
        DISCRIMINANT_STATEMENT,
        match nonnegative_on_line(&discriminant) {
            Ok(()) => Ok(format!(
                "B^2 - 4 A C = {}: every factor of odd multiplicity has no real root and the \
                 leading coefficient is positive (exact)",
                format_index_polynomial(&discriminant, "x")
            )),
            Err(Some(point)) => Err((
                None,
                not_a_global_preserver(&point, &vec_poly::evaluate(&discriminant, &point)),
            )),
            Err(None) => Err((
                None,
                format!(
                    "B^2 - 4 A C = {} has a real root of odd multiplicity or a negative leading \
                     coefficient: T is not a real-rootedness preserver on the whole line",
                    format_index_polynomial(&discriminant, "x")
                ),
            )),
        },
    )];
    out.push(exact_hypothesis(
        "weyl_step_interlacing",
        INTERLACING_STATEMENT,
        fixed_interlacing(a, b)
            .map(|()| "exact (Bezout test of strict interlacing)".to_string())
            .map_err(|reason| (None, reason)),
    ));
    let name = "weyl_step_leading_coefficients_one_signed";
    let function = "sgn(A_a) (A_a + k B_(a+1) + k (k - 1) C_(a+2)) in k";
    out.push(match fixed_leading_coefficients(a, b, c) {
        // The failing value is an exponent k, not an index n: no witness.
        Ok(decision @ IndexSignDecision::FailsAt { .. }) => exact_hypothesis(
            name,
            LEADING_STATEMENT,
            Err((None, decision.describe(function).replace(" n = ", " k = "))),
        ),
        Ok(decision) => {
            LiftHypothesis::from_sign_decision(name, LEADING_STATEMENT, decision, function)
        }
        Err(reason) => exact_hypothesis(name, LEADING_STATEMENT, Err((None, reason))),
    });
    out
}

const INTERLACING_STATEMENT: &str = "(ii) A strictly interlaces B (deg B = deg A + 1) or \
     strictly alternates left of B (deg B = deg A), with leading coefficients of the same sign \
     when deg A >= 1";

const LEADING_STATEMENT: &str = "(iii) T(x^k) = A x^k + k B x^(k-1) + k (k - 1) C x^(k-2) has \
     degree deg A + k and leading coefficients of one sign for every k >= 0";

/// Hypothesis (i) for index-dependent coefficients.
fn discriminant_hypotheses(parts: &[Bivariate; 3], n_first: usize) -> Vec<LiftHypothesis> {
    let [a, b, c] = parts.clone().map(|p| biv(&p));
    let discriminant = b
        .multiply(&b)
        .subtract(&a.multiply(&c).scale(&Q::from_integer(BigInt::from(4))));
    let rows: Bivariate = discriminant.rows().to_vec();
    let name = "weyl_step_discriminant_nonnegative";
    if c.is_zero() {
        return vec![exact_hypothesis(
            name,
            DISCRIMINANT_STATEMENT,
            Ok("C = 0, so B^2 - 4 A C = B^2 >= 0".into()),
        )];
    }
    match x_degree(&rows) {
        None => vec![exact_hypothesis(
            name,
            DISCRIMINANT_STATEMENT,
            Ok("B^2 - 4 A C = 0".into()),
        )],
        Some(0) => vec![sign_hypothesis(
            name,
            DISCRIMINANT_STATEMENT,
            "B^2 - 4 A C",
            &function(x_column(&rows, 0)),
            n_first,
            IndexSignCondition::Nonnegative,
        )],
        Some(2) => {
            let (q0, q1, q2) = (x_column(&rows, 0), x_column(&rows, 1), x_column(&rows, 2));
            let inner = vec_poly::sub(
                &vec_poly::scale(&vec_poly::mul(&q2, &q0), &Q::from_integer(BigInt::from(4))),
                &vec_poly::mul(&q1, &q1),
            );
            vec![
                sign_hypothesis(
                    "weyl_step_discriminant_leading_nonnegative",
                    &format!(
                        "{DISCRIMINANT_STATEMENT}: with B^2 - 4 A C = q_2 x^2 + q_1 x + q_0, \
                         q_2(n) >= 0"
                    ),
                    "q_2(n)",
                    &function(q2),
                    n_first,
                    IndexSignCondition::Nonnegative,
                ),
                sign_hypothesis(
                    "weyl_step_discriminant_constant_nonnegative",
                    &format!("{DISCRIMINANT_STATEMENT}: q_0(n) >= 0"),
                    "q_0(n)",
                    &function(q0),
                    n_first,
                    IndexSignCondition::Nonnegative,
                ),
                sign_hypothesis(
                    name,
                    &format!(
                        "{DISCRIMINANT_STATEMENT}: 4 q_2(n) q_0(n) - q_1(n)^2 >= 0 (with q_2, \
                         q_0 >= 0, the quadratic is nonnegative on the line)"
                    ),
                    "4 q_2 q_0 - q_1^2",
                    &function(inner),
                    n_first,
                    IndexSignCondition::Nonnegative,
                ),
            ]
        }
        Some(_) => {
            if let Some((gamma, w)) = separate(&rows) {
                if nonnegative_on_line(&w).is_ok() {
                    return vec![signed_fixed_factor(name, &gamma, &w, n_first, true)];
                }
                if nonnegative_on_line(&vec_poly::neg(&w)).is_ok() {
                    return vec![signed_fixed_factor(name, &gamma, &w, n_first, false)];
                }
            }
            vec![refuted_at_first_index(
                name,
                DISCRIMINANT_STATEMENT,
                n_first,
                |n| {
                    let value = discriminant.evaluate_index(&Q::from_integer(BigInt::from(n)));
                    match nonnegative_on_line(&value) {
                        Ok(()) => None,
                        Err(Some(point)) => Some(not_a_global_preserver(
                            &point,
                            &vec_poly::evaluate(&value, &point),
                        )),
                        Err(None) => Some(
                            "B^2 - 4 A C has a real root of odd multiplicity or a negative \
                             leading coefficient there"
                                .into(),
                        ),
                    }
                },
                "the discriminant is not free of x, quadratic in x or gamma(n) w(x) with a \
                 one-signed w",
            )]
        }
    }
}

/// `gamma(n) w(x) >= 0` with `w` of one sign on the line.
fn signed_fixed_factor(
    name: &str,
    gamma: &[Q],
    w: &[Q],
    n_first: usize,
    nonnegative: bool,
) -> LiftHypothesis {
    let mut hypothesis = sign_hypothesis(
        name,
        &format!(
            "{DISCRIMINANT_STATEMENT}: B^2 - 4 A C = gamma(n) w(x) with w(x) = {} {} 0 on the \
             line (exact)",
            format_index_polynomial(w, "x"),
            if nonnegative { ">=" } else { "<=" }
        ),
        "gamma(n)",
        &function(gamma.to_vec()),
        n_first,
        if nonnegative {
            IndexSignCondition::Nonnegative
        } else {
            IndexSignCondition::Nonpositive
        },
    );
    hypothesis.evidence = format!(
        "w(x) one-signed on the line (exact); {}",
        hypothesis.evidence
    );
    hypothesis
}

/// An unresolved hypothesis, or failed when the exact check `fails` refutes
/// it at the first index.
fn refuted_at_first_index(
    name: &str,
    statement: &str,
    n_first: usize,
    fails: impl Fn(usize) -> Option<String>,
    unresolved: &str,
) -> LiftHypothesis {
    match fails(n_first) {
        Some(evidence) => exact_hypothesis(
            name,
            statement,
            Err((Some(n_first), format!("at n = {n_first}: {evidence}"))),
        ),
        None => unresolved_hypothesis(name, statement, unresolved.into()),
    }
}

/// Hypothesis (ii) for index-dependent coefficients.
fn interlacing_hypotheses(parts: &[Bivariate; 3], n_first: usize) -> Vec<LiftHypothesis> {
    let name = "weyl_step_interlacing";
    let (a, b) = (&parts[0], &parts[1]);
    let (Some(da), Some(db)) = (x_degree(a), x_degree(b)) else {
        return vec![exact_hypothesis(
            name,
            INTERLACING_STATEMENT,
            Err((None, "A or B is zero".into())),
        )];
    };
    if db != da && db != da + 1 && !(da == 0 && db <= 1) {
        return vec![exact_hypothesis(
            name,
            INTERLACING_STATEMENT,
            Err((
                None,
                format!("deg_x A = {da}, deg_x B = {db}: neither deg B = deg A nor deg A + 1"),
            )),
        )];
    }
    let column = |p: &Bivariate, j: usize| x_column(p, j);
    match (da, db) {
        (0, _) => vec![
            sign_hypothesis(
                "weyl_step_constant_term_nonzero",
                &format!("{INTERLACING_STATEMENT}: here deg A = 0, A(n) != 0"),
                "A(n)",
                &function(column(a, 0)),
                n_first,
                IndexSignCondition::Nonzero,
            ),
            sign_hypothesis(
                name,
                &format!(
                    "{INTERLACING_STATEMENT}: here deg A = 0 and deg B <= 1 (vacuous), B != 0"
                ),
                "the leading coefficient of B",
                &function(column(b, db)),
                n_first,
                IndexSignCondition::Nonzero,
            ),
        ],
        (1, 1 | 2) => {
            let (a0, a1) = (column(a, 0), column(a, 1));
            let lead = sign_hypothesis(
                "weyl_step_leading_coefficients_same_sign",
                &format!("{INTERLACING_STATEMENT}: A_1(n) B_{db}(n) > 0"),
                "A_1 B_lead",
                &product(&a1, &column(b, db)),
                n_first,
                IndexSignCondition::Positive,
            );
            let comparison = if db == 1 {
                let (b0, b1) = (column(b, 0), column(b, 1));
                sign_hypothesis(
                    name,
                    &format!(
                        "{INTERLACING_STATEMENT}: deg A = deg B = 1, the root -A_0/A_1 of A is \
                         left of the root -B_0/B_1 of B: A_0 B_1 - A_1 B_0 > 0 (with A_1 B_1 > 0)"
                    ),
                    "A_0 B_1 - A_1 B_0",
                    &product(&a0, &b1).sub(&product(&a1, &b0)),
                    n_first,
                    IndexSignCondition::Positive,
                )
            } else {
                let (b0, b1, b2) = (column(b, 0), column(b, 1), column(b, 2));
                // B at the root -A_0/A_1 of A, times A_1^2.
                let value = vec_poly::add(
                    &vec_poly::sub(
                        &vec_poly::mul(&b2, &vec_poly::mul(&a0, &a0)),
                        &vec_poly::mul(&b1, &vec_poly::mul(&a0, &a1)),
                    ),
                    &vec_poly::mul(&b0, &vec_poly::mul(&a1, &a1)),
                );
                sign_hypothesis(
                    name,
                    &format!(
                        "{INTERLACING_STATEMENT}: deg A = 1, deg B = 2, B has opposite signs at \
                         the root of A and at infinity: B_2 (B_2 A_0^2 - B_1 A_0 A_1 + B_0 A_1^2) \
                         < 0"
                    ),
                    "B_2 (B_2 A_0^2 - B_1 A_0 A_1 + B_0 A_1^2)",
                    &product(&b2, &value),
                    n_first,
                    IndexSignCondition::Negative,
                )
            };
            vec![lead, comparison]
        }
        _ => {
            if let (Some((alpha, fixed_a)), Some((beta, fixed_b))) = (separate(a), separate(b)) {
                let fixed = fixed_interlacing_up_to_sign(&fixed_a, &fixed_b);
                let mut out = vec![exact_hypothesis(
                    name,
                    &format!(
                        "{INTERLACING_STATEMENT}: A = alpha(n) a(x), B = beta(n) b(x) with the \
                         roots of a = {} and b = {} interlacing (exact)",
                        format_index_polynomial(&fixed_a, "x"),
                        format_index_polynomial(&fixed_b, "x")
                    ),
                    fixed
                        .clone()
                        .map(|_| "exact (Bezout test of strict interlacing)".to_string())
                        .map_err(|reason| (None, reason)),
                )];
                if let Ok(sign) = fixed {
                    out.push(sign_hypothesis(
                        "weyl_step_leading_coefficients_same_sign",
                        &format!("{INTERLACING_STATEMENT}: the leading coefficients of A and B"),
                        "sgn(lead a lead b) alpha(n) beta(n)",
                        &product(&alpha, &beta).scale(&sign),
                        n_first,
                        IndexSignCondition::Positive,
                    ));
                }
                return out;
            }
            vec![refuted_at_first_index(
                name,
                INTERLACING_STATEMENT,
                n_first,
                |n| {
                    let at = Q::from_integer(BigInt::from(n));
                    fixed_interlacing(&biv(a).evaluate_index(&at), &biv(b).evaluate_index(&at))
                        .err()
                },
                "deg_x A >= 2 with index-dependent coefficients that do not separate",
            )]
        }
    }
}

/// The root pattern of (ii) for fixed `a`, `b`, ignoring the signs of the
/// leading coefficients; `Ok(sign of lead(a) lead(b))`.
fn fixed_interlacing_up_to_sign(a: &[Q], b: &[Q]) -> Result<Q, String> {
    let (Some(da), Some(db)) = (vec_poly::degree(a), vec_poly::degree(b)) else {
        return Err("A or B is zero".into());
    };
    let sign = if (&a[da] * &b[db]).is_negative() {
        -Q::one()
    } else {
        Q::one()
    };
    fixed_interlacing(a, &vec_poly::scale(b, &sign)).map(|()| sign)
}

/// Hypothesis (iii) for index-dependent coefficients.
fn leading_coefficient_hypotheses(parts: &[Bivariate; 3], n_first: usize) -> Vec<LiftHypothesis> {
    let name = "weyl_step_leading_coefficients_one_signed";
    let [a, b, c] = parts;
    let Some(da) = x_degree(a) else {
        return vec![exact_hypothesis(
            name,
            LEADING_STATEMENT,
            Err((None, "A is zero".into())),
        )];
    };
    if x_degree(b).is_some_and(|db| db > da + 1) || x_degree(c).is_some_and(|dc| dc > da + 2) {
        return vec![exact_hypothesis(
            name,
            LEADING_STATEMENT,
            Err((
                None,
                "deg_x B > deg_x A + 1 or deg_x C > deg_x A + 2: T(x^k) has degree above deg A + \
                 k"
                .into(),
            )),
        )];
    }
    let lead_a = x_column(a, da);
    let b_top = x_column(b, da + 1);
    let c_top = x_column(c, da + 2);
    let mut out = Vec::new();
    if b_top.is_empty() && c_top.is_empty() {
        out.push(exact_hypothesis(
            name,
            LEADING_STATEMENT,
            Ok(format!(
                "T(x^k) has leading coefficient A_{da}(n), nonzero by (ii)"
            )),
        ));
        return out;
    }
    for (top, which) in [(&b_top, "B_(a+1)"), (&c_top, "C_(a+2)")] {
        if top.is_empty() {
            continue;
        }
        out.push(sign_hypothesis(
            name,
            &format!(
                "{LEADING_STATEMENT}: with a = deg A, A_a(n) {which}(n) >= 0, so sgn(A_a) times \
                 the leading coefficient A_a + k B_(a+1) + k (k - 1) C_(a+2) is >= |A_a| > 0"
            ),
            &format!("A_a {which}"),
            &product(&lead_a, top),
            n_first,
            IndexSignCondition::Nonnegative,
        ));
    }
    out
}

/// The region hypotheses: every coefficient of `A`, `B`, `C` (of the
/// reflected operator for `[0, inf)`) of one sign for every `n >= n_first`.
pub(super) fn weyl_step_region_hypotheses(
    parts: &[Bivariate; 3],
    reflect: bool,
    n_first: usize,
) -> Vec<LiftHypothesis> {
    // The coefficient columns, as polynomials in n.
    let mut columns: Vec<(String, Vec<Q>)> = Vec::new();
    for (order, (letter, part)) in ["A", "B", "C"].iter().zip(parts.iter()).enumerate() {
        for j in 0..x_width(part) {
            let column = x_column(part, j);
            if column.is_empty() {
                continue;
            }
            let flip = reflect && (j + order) % 2 == 1;
            let column = if flip { vec_poly::neg(&column) } else { column };
            columns.push((format!("{letter}_{j}"), column));
        }
    }
    let name = if reflect {
        "weyl_step_reflected_coefficients_one_signed"
    } else {
        "weyl_step_coefficients_one_signed"
    };
    let operator = if reflect {
        "the reflected operator A(-x) - B(-x) D + C(-x) D^2"
    } else {
        "A + B D + C D^2"
    };
    // The orientation: the sign of the first coefficient for large n.
    let positive = columns
        .first()
        .and_then(|(_, column)| column.last())
        .is_none_or(Signed::is_positive);
    let condition = if positive {
        IndexSignCondition::Nonnegative
    } else {
        IndexSignCondition::Nonpositive
    };
    let relation = if positive { ">= 0" } else { "<= 0" };
    let statement = format!(
        "every coefficient of {operator} is {relation} for every n >= {n_first}, so the operator \
         maps polynomials with coefficients of one sign to such polynomials"
    );
    let (fixed, varying): (Vec<_>, Vec<_>) = columns
        .into_iter()
        .partition(|(_, column)| column.len() <= 1);
    let mut out = Vec::new();
    let failing = fixed.iter().find(|(_, column)| {
        let value = column.first().cloned().unwrap_or_else(Q::zero);
        if positive {
            value.is_negative()
        } else {
            value.is_positive()
        }
    });
    out.push(exact_hypothesis(
        name,
        &statement,
        match failing {
            None => Ok(format!(
                "the coefficients free of n are {relation} (exact); {} coefficient(s) depend on n",
                varying.len()
            )),
            Some((label, _)) => Err((
                None,
                format!("the coefficient {label} free of n has the other sign (exact)"),
            )),
        },
    ));
    for (label, column) in varying {
        out.push(sign_hypothesis(
            &format!("{name}_{label}"),
            &statement,
            &format!("the coefficient {label}(n)"),
            &function(column),
            n_first,
            condition,
        ));
    }
    out
}
