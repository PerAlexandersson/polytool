//! Exact bivariate arithmetic and degree budgets, reusing shared kernels.

use super::*;
use crate::recurrence::generating_function::{exact_divide, gcd_bivariate};

impl super::super::certified::algebra::ExactDeterminantCoefficient for B {
    fn zero() -> Self {
        B::zero()
    }
    fn one() -> Self {
        B::one()
    }
    fn is_zero(&self) -> bool {
        self.is_zero()
    }
    fn multiply(&self, other: &Self) -> Self {
        self.multiply(other)
    }
    fn subtract(&self, other: &Self) -> Self {
        self.subtract(other)
    }
    fn divide_exact(&self, other: &Self) -> Self {
        exact_divide(self, other).expect("Bareiss division is exact over Q[t,b]")
    }
    fn negative(&self) -> Self {
        self.scale(&q(-1))
    }
}

pub(super) fn t(p: &P) -> B {
    B::from_variable_polynomial(p.coeffs())
}
pub(super) fn c(n: i64) -> B {
    B::constant(q(n))
}
pub(super) fn variable() -> B {
    B::monomial(1, 0, Q::one())
}
pub(super) fn rem(a: &B, monic: &B) -> B {
    let degree = monic.index_degree().expect("branch polynomial nonzero");
    let mut a = a.clone();
    while a.index_degree().is_some_and(|d| d >= degree) {
        let shift = a.index_degree().unwrap() - degree;
        let lead = B::from_variable_polynomial(a.rows().last().unwrap());
        let monomial = B::monomial(shift, 0, Q::one());
        a = a.subtract(&monic.multiply(&lead).multiply(&monomial));
    }
    a
}
pub(super) fn resultant(p: &B, f: &B) -> P {
    let rows = |a: &B| {
        a.rows()
            .iter()
            .map(|r| P::new(r.clone()))
            .collect::<Vec<_>>()
    };
    super::super::certified::algebra::resultant(&rows(p), &rows(f))
}
pub(super) fn derivative(p: &B) -> B {
    B::new(
        p.rows()
            .iter()
            .enumerate()
            .skip(1)
            .map(|(i, r)| vec_poly::scale(r, &Q::from_integer(i.into())))
            .collect(),
    )
}
pub(super) fn discriminant(p: &B) -> P {
    let d = p.index_degree().unwrap_or(0);
    resultant(p, &derivative(p)).scale(&q(if (d * (d - 1) / 2) % 2 == 0 { 1 } else { -1 }))
}

/// Rational expressions in t and a branch coordinate; gcd reduction avoids
/// exponential artificial factors when expanding amplitude norms.
#[derive(Clone)]
pub(super) struct R {
    pub n: B,
    pub d: B,
    relation: Option<std::sync::Arc<B>>,
}
impl R {
    pub fn wire(&self) -> EffectiveBkwBranchRational {
        EffectiveBkwBranchRational {
            numerator: encode_b(&self.n),
            denominator: encode_b(&self.d),
        }
    }
    pub fn new(n: B, d: B) -> Self {
        assert!(!d.is_zero(), "internal nonzero symbolic denominator");
        let g = gcd_bivariate(&n, &d);
        Self {
            n: exact_divide(&n, &g).expect("gcd divides"),
            d: exact_divide(&d, &g).expect("gcd divides"),
            relation: None,
        }
    }
    pub fn of(n: B) -> Self {
        Self {
            n,
            d: B::one(),
            relation: None,
        }
    }
    pub fn modulo(mut self, p: &B) -> Self {
        self.relation = Some(std::sync::Arc::new(p.clone()));
        self.combine(self.n.clone(), self.d.clone(), None)
    }
    fn combine(&self, mut n: B, mut d: B, other: Option<&Self>) -> Self {
        let relation = self
            .relation
            .clone()
            .or_else(|| other.and_then(|r| r.relation.clone()));
        if let Some(p) = &relation {
            let reduced_d = rem(&d, p);
            if !reduced_d.is_zero() {
                n = rem(&n, p);
                d = reduced_d;
            }
            // Quotient reduction keeps the branch degree bounded. A full
            // bivariate gcd here can dominate the entire higher-order
            // calculation. The compact supplier explicitly checks the stored
            // pair/cofactor denominators with check_pair_denominators before
            // using compound expressions that could cancel their factors.
            if n.is_zero() {
                d = B::one();
            }
            if n == d {
                n = B::one();
                d = B::one();
            }
            return Self { n, d, relation };
        }
        let mut value = Self::new(n, d);
        value.relation = relation;
        value
    }
    pub fn add(&self, r: &Self) -> Self {
        self.combine(
            self.n.multiply(&r.d).add(&r.n.multiply(&self.d)),
            self.d.multiply(&r.d),
            Some(r),
        )
    }
    pub fn neg(&self) -> Self {
        self.combine(self.n.scale(&q(-1)), self.d.clone(), None)
    }
    pub fn sub(&self, r: &Self) -> Self {
        self.add(&r.neg())
    }
    pub fn mul(&self, r: &Self) -> Self {
        self.combine(self.n.multiply(&r.n), self.d.multiply(&r.d), Some(r))
    }
    pub fn div(&self, r: &Self) -> Self {
        self.combine(self.n.multiply(&r.d), self.d.multiply(&r.n), Some(r))
    }
    pub fn scale(&self, q: &Q) -> Self {
        self.combine(self.n.scale(q), self.d.clone(), None)
    }
    pub fn square(&self) -> Self {
        self.mul(self)
    }
    pub fn pow(&self, n: usize) -> Self {
        self.combine(self.n.power(n), self.d.power(n), None)
    }
    pub fn signed_polynomial(&self, branch: &B) -> B {
        rem(&self.n.multiply(&self.d), branch)
    }
}

pub(super) struct DegreeBudget {
    pub count: EffectiveBkwZeroCount,
    pub explanation: Vec<String>,
    pub data: EffectiveBkwDegreeCertificate,
}
pub(super) fn degree_budget(n: &B, d: &B) -> Result<DegreeBudget, RecurrenceError> {
    // s^n P_n(1/s) is polynomial iff deg_t D_j <= j and deg_t N_j <= j.
    for a in [n, d] {
        if a.rows().iter().enumerate().any(|(j, r)| r.len() > j + 1) {
            return super::degree::newton_budget(n, d);
        }
    }
    let top = |p: &B| {
        P::new(
            p.rows()
                .iter()
                .enumerate()
                .map(|(j, r)| r.get(j).cloned().unwrap_or_else(Q::zero))
                .collect(),
        )
    };
    let nt = top(n);
    let dt = top(d);
    let g = nt.gcd(&dt);
    let nt = nt.exact_div(&g);
    let dt = dt.exact_div(&g);
    let k = dt.coeff(0);
    let nt = nt.scale(&k.recip());
    let dt = dt.scale(&k.recip());
    let (offsets, leading) =
        if nt == P::one() && dt.degree() == Some(1) && dt.coeff(1).is_negative() {
            (vec![0], format!("[t^n]P_n = ({})^n > 0", -dt.coeff(1)))
        } else if nt == P::one() && dt == P::new(vec![q(1), q(0), q(-1)]) {
            (
                vec![0, -1],
                "[t^n]P_n = [n even]; hence degree <= n-1 for odd n and monic degree n for even n"
                    .into(),
            )
        } else {
            return super::degree::newton_budget(n, d);
        };
    let count = EffectiveBkwZeroCount {
        slope: 1,
        offsets,
        divisor: 1,
    };
    Ok(DegreeBudget {
        data: EffectiveBkwDegreeCertificate {
            upper_bound: count.clone(),
            leading_numerator: encode(&nt),
            leading_denominator: encode(&dt),
            positive_top_residues: vec![0],
            newton: None,
        },
        count,
        explanation: vec![
            "deg_t D_j <= j and deg_t N_j <= j imply deg P_n <= n by induction".into(),
            format!(
                "s=1/t: top coefficient GF = {:?}/{:?}",
                encode(&nt),
                encode(&dt)
            ),
            leading,
        ],
    })
}

pub(super) fn empty_certificate(
    d: &[P],
    branch: &B,
    domain: &Domain,
    degree: &DegreeBudget,
) -> EffectiveBkwCertificate {
    EffectiveBkwCertificate {
        coordinate_transport: EffectiveBkwCoordinateTransport::Identity,
        root_region: EffectiveBkwRootRegion::WholeRealLine,
        interval: domain.wire.clone(),
        characteristic: d.iter().rev().map(encode).collect(),
        branch_equation: encode_b(branch),
        branch_selection: String::new(),
        expansion: vec![],
        expansion_data: EffectiveBkwExpansionData {
            pair_real_part: R::of(c(0)).wire(),
            pair_norm_squared: R::of(c(0)).wire(),
            cofactor: vec![],
            remainder_root: None,
            remainder_amplitude: None,
        },
        signs: vec![],
        domination: String::new(),
        cutoff: 1,
        tail_at_cutoff: "0".into(),
        automatic_uniform_constants: None,
        phase: vec![],
        phase_limits: EffectiveBkwPhaseLimits {
            left_eigenvalue_argument: "1".into(),
            right_eigenvalue_argument: "0".into(),
            left_amplitude_argument: String::new(),
            right_amplitude_argument: String::new(),
        },
        exterior: degree.explanation.clone(),
        degree_certificate: degree.data.clone(),
        zero_count: degree.count.clone(),
        finite_rows: vec![],
        lean_declarations: vec![
            "RealRooted.EffectiveBKW.PairExpansion".into(),
            "RealRooted.EffectiveBKW.PairExpansion.exists_alternating".into(),
            "RealRooted.EffectiveBKW.splits_nodup_of_strictAnti_alternating".into(),
        ],
    }
}
