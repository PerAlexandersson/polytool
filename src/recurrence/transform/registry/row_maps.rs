//! Row maps for recurrence transformations.

use super::*;

impl RecurrenceTransformation {
    /// The `kind` tag.
    pub fn kind_name(&self) -> &'static str {
        match self {
            Self::MultiplyPower { .. } => "multiply_power",
            Self::DividePower { .. } => "divide_power",
            Self::Reversal { .. } => "reversal",
            Self::Affine { .. } => "affine",
            Self::Mobius { .. } => "mobius",
            Self::FPolynomial { .. } => "f_polynomial",
            Self::Gamma { .. } => "gamma",
            Self::SymmetricDecomposition { .. } => "symmetric_decomposition",
            Self::Diagonal { .. } => "diagonal",
            Self::DerivativeAntecedent { .. } => "derivative_antecedent",
            Self::Derivative { .. } => "derivative",
        }
    }

    /// The degree schedule, if the map uses one.
    pub fn schedule(&self) -> Option<&DegreeSchedule> {
        match self {
            Self::Reversal { schedule }
            | Self::Mobius { schedule, .. }
            | Self::FPolynomial { schedule }
            | Self::Gamma { schedule }
            | Self::SymmetricDecomposition { schedule } => Some(schedule),
            _ => None,
        }
    }

    /// Number of components of a transformed state.
    pub fn component_count(&self) -> usize {
        match self {
            Self::SymmetricDecomposition { .. } => 2,
            _ => 1,
        }
    }

    /// Whether `index` lies in the map's index domain.
    pub fn contains_index(&self, index: usize) -> bool {
        self.schedule().is_none_or(|s| s.contains(index))
    }

    pub(in crate::recurrence::transform) fn mobius_parameters(
        &self,
    ) -> Result<Option<MobiusParameters<'_>>, TransformError> {
        Ok(match self {
            Self::Reversal { schedule } => Some(MobiusParameters {
                a: q(0),
                b: q(1),
                c: q(1),
                e: q(0),
                schedule,
            }),
            Self::FPolynomial { schedule } => Some(MobiusParameters {
                a: q(1),
                b: q(0),
                c: q(1),
                e: q(1),
                schedule,
            }),
            Self::Mobius {
                a,
                b,
                c,
                e,
                schedule,
            } => Some(MobiusParameters {
                a: parse_q(a, "a")?,
                b: parse_q(b, "b")?,
                c: parse_q(c, "c")?,
                e: parse_q(e, "e")?,
                schedule,
            }),
            _ => None,
        })
    }

    pub(in crate::recurrence::transform) fn affine_parameters(
        &self,
    ) -> Result<Option<(Q, Q)>, TransformError> {
        match self {
            Self::Affine { a, b } => Ok(Some((parse_q(a, "a")?, parse_q(b, "b")?))),
            _ => Ok(None),
        }
    }

    /// The coordinate map as an `ore-algebra` [`Substitution`] (validated),
    /// with the exponent `d_n = (d_0 - h n_0) + h n` of the schedule in the
    /// original index: `x^r` multiplication and removal as power gauges
    /// `x^(+-r)`, the affine map, the reversal, the `f`-polynomial and the
    /// Möbius map (scale `1`).  `None` for the other maps (gamma, symmetric
    /// decomposition, diagonal, derivative antecedent, derivative), which
    /// are not substitutions.
    ///
    /// # Errors
    ///
    /// [`TransformError::SingularMobius`] for a zero affine slope or
    /// `a e - b c = 0`, and the schedule and parameter errors.
    pub(in crate::recurrence) fn substitution(
        &self,
    ) -> Result<Option<Substitution>, TransformError> {
        let power = |power: usize, sign: i64| -> Result<Vec<Q>, TransformError> {
            let power = i64::try_from(power).map_err(|_| TransformError::Limit("power".into()))?;
            Ok(vec![q(sign * power)])
        };
        let degree = |schedule: &DegreeSchedule| -> Result<Vec<Q>, TransformError> {
            Ok(vec![q(schedule.weight_constant()?), q(schedule.step)])
        };
        let substitution = match self {
            Self::MultiplyPower { power: r } => Substitution::PowerGauge {
                power: power(*r, 1)?,
            },
            Self::DividePower { power: r } => Substitution::PowerGauge {
                power: power(*r, -1)?,
            },
            Self::Affine { .. } => {
                let (slope, intercept) = self.affine_parameters()?.expect("affine");
                Substitution::Affine { slope, intercept }
            }
            Self::Reversal { schedule } => Substitution::Reversal {
                degree: degree(schedule)?,
            },
            Self::FPolynomial { schedule } => Substitution::FPolynomial {
                degree: degree(schedule)?,
            },
            Self::Mobius { .. } => {
                let m = self.mobius_parameters()?.expect("mobius");
                Substitution::Mobius(Box::new(MobiusMap {
                    degree: degree(m.schedule)?,
                    a: m.a,
                    b: m.b,
                    c: m.c,
                    e: m.e,
                    scale: Q::one(),
                }))
            }
            _ => return Ok(None),
        };
        substitution.validate().map_err(substitution_error)?;
        Ok(Some(substitution))
    }

    /// Check parameters and schedule; singular maps and invalid schedules are rejected.
    pub fn validate(&self) -> Result<(), TransformError> {
        if let Some(schedule) = self.schedule() {
            schedule.validate()?;
        }
        self.substitution()?;
        match self {
            Self::Diagonal { weight } => weight.validate()?,
            Self::DerivativeAntecedent { g } => {
                derivative_antecedent_parameters(g)?;
            }
            Self::Derivative { order } => derivative_order_check(*order)?,
            _ => {}
        }
        Ok(())
    }

    /// Apply the map to the row `P_index`, returning one coefficient vector per
    /// component (ascending powers, zero as `[0]`).
    ///
    /// Preconditions (degree bound, divisibility, symmetry) are checked exactly
    /// on this row and reported as errors; they are never assumed.
    pub fn map_row(
        &self,
        index: usize,
        row: &[BigRational],
    ) -> Result<Vec<Vec<BigRational>>, TransformError> {
        self.validate()?;
        let row = vec_poly::trim(row);
        match self {
            Self::MultiplyPower { power } => Ok(vec![vec_poly::shift_by(&row, *power)]),
            Self::DividePower { power } => Ok(vec![divide_power(index, &row, *power)?]),
            Self::Affine { .. } => {
                let (a, b) = self.affine_parameters()?.expect("affine");
                Ok(vec![vec_poly::compose(&row, &[b, a])])
            }
            Self::Reversal { schedule } => {
                let d = schedule.degree_at(index)?;
                check_degree(index, &row, d)?;
                Ok(vec![vec_poly::reverse(&row, d).expect("degree checked")])
            }
            Self::FPolynomial { schedule } => {
                let d = schedule.degree_at(index)?;
                check_degree(index, &row, d)?;
                // Integer rows use BigInt arithmetic; both routes are exact.
                let mapped = match integer_coefficients(&row) {
                    Some(integers) => f_polynomial(&Polynomial::new(integers), d)
                        .expect("degree checked")
                        .coeffs()
                        .iter()
                        .cloned()
                        .map(Q::from_integer)
                        .collect::<Vec<_>>(),
                    None => f_polynomial(&Polynomial::new(row), d)
                        .expect("degree checked")
                        .coeffs()
                        .to_vec(),
                };
                Ok(vec![vec_poly::trim(&mapped)])
            }
            Self::Mobius { .. } => {
                let m = self.mobius_parameters()?.expect("mobius");
                let d = m.schedule.degree_at(index)?;
                check_degree(index, &row, d)?;
                Ok(vec![mobius_substitute(&row, &m.a, &m.b, &m.c, &m.e, d)])
            }
            Self::Gamma { schedule } => {
                let d = schedule.degree_at(index)?;
                Ok(vec![gamma_at_index(index, &row, d)?])
            }
            Self::SymmetricDecomposition { schedule } => {
                let d = schedule.degree_at(index)?;
                check_degree(index, &row, d)?;
                let (a, b) = Polynomial::new(row)
                    .stapledon_decomposition(d)
                    .expect("degree checked");
                Ok(vec![vec_poly::trim(a.coeffs()), vec_poly::trim(b.coeffs())])
            }
            Self::Diagonal { weight } => Ok(vec![weight.map_row(index, &row)?]),
            Self::DerivativeAntecedent { g } => {
                let (g, _) = derivative_antecedent_parameters(g)?;
                Ok(vec![derivative_antecedent_row(&g, &row)])
            }
            Self::Derivative { order } => Ok(vec![derivative_row(*order, &row)]),
        }
    }

    /// Reconstruct `P_index` from its transformed components, checking the
    /// image conditions of the inverse exactly.
    pub fn inverse_row(
        &self,
        index: usize,
        components: &[Vec<BigRational>],
    ) -> Result<Vec<BigRational>, TransformError> {
        self.validate()?;
        if components.len() != self.component_count() {
            return Err(TransformError::ComponentCount {
                expected: self.component_count(),
                found: components.len(),
            });
        }
        let image = vec_poly::trim(&components[0]);
        match self {
            Self::MultiplyPower { power } => {
                divide_power(index, &image, *power).map_err(|_| TransformError::ImageCondition {
                    index,
                    reason: format!("x^{power} must divide Q_n"),
                })
            }
            Self::DividePower { power } => Ok(vec_poly::shift_by(&image, *power)),
            Self::Affine { .. } => {
                let (a, b) = self.affine_parameters()?.expect("affine");
                let inverse_a = a.recip();
                let inverse_b = -(&b * &inverse_a);
                Ok(vec_poly::compose(&image, &[inverse_b, inverse_a]))
            }
            Self::Reversal { .. } | Self::FPolynomial { .. } | Self::Mobius { .. } => {
                let m = self.mobius_parameters()?.expect("mobius");
                let d = m.schedule.degree_at(index)?;
                if vec_poly::degree(&image).is_some_and(|degree| degree > d) {
                    return Err(TransformError::ImageCondition {
                        index,
                        reason: format!("deg Q_n must be at most d_n = {d}"),
                    });
                }
                if let Self::Reversal { .. } = self {
                    return Ok(vec_poly::reverse(&image, d).expect("degree checked"));
                }
                let delta = m.delta();
                let unscaled =
                    mobius_substitute(&image, &m.e, &-m.b.clone(), &-m.c.clone(), &m.a, d);
                Ok(vec_poly::scale(
                    &unscaled,
                    &rational_power(&delta, -(d as i64)),
                ))
            }
            Self::Gamma { schedule } => {
                let d = schedule.degree_at(index)?;
                gamma_inverse(index, &image, d)
            }
            Self::SymmetricDecomposition { schedule } => {
                let d = schedule.degree_at(index)?;
                let b = vec_poly::trim(&components[1]);
                if !symmetric_at(&image, Some(d)) {
                    return Err(TransformError::ImageCondition {
                        index,
                        reason: format!("a_n must be symmetric at index d_n = {d}"),
                    });
                }
                if !symmetric_at(&b, d.checked_sub(1)) {
                    return Err(TransformError::ImageCondition {
                        index,
                        reason: format!(
                            "b_n must be symmetric at index d_n - 1 = {}",
                            d as i64 - 1
                        ),
                    });
                }
                Ok(vec_poly::add(&image, &vec_poly::shift(&b)))
            }
            Self::Diagonal { weight } => weight.inverse_row(index, &image),
            Self::DerivativeAntecedent { g } => {
                let (g, _) = derivative_antecedent_parameters(g)?;
                Ok(vec_poly::trim(&derivative_lift_row(&g, &image)))
            }
            Self::Derivative { order } => Err(lossy_inverse_error(index, *order)),
        }
    }

    /// Transport a rule `f P_m = sum s c D^k P_(m-r) + g` (target index `m`)
    /// into the new coordinates.
    ///
    /// The result is an exact consequence of the source rule at every index
    /// where the source rule holds and the transformed rows are defined.
    ///
    /// ```
    /// use polytool::recurrence::*;
    ///
    /// // Eulerian A008292: E_1 = 1, degree n - 1, and its f-polynomials.
    /// let (rule, _seeds, _first) = polytool::oeis::by_id("A008292")
    ///     .unwrap()
    ///     .recurrence_export_parts()
    ///     .unwrap();
    /// let f = RecurrenceTransformation::FPolynomial {
    ///     schedule: DegreeSchedule::new(1, 0, 1),
    /// };
    /// let TransportedRecurrence::Scalar(q) = f.transport(&rule).unwrap() else {
    ///     unreachable!()
    /// };
    /// assert_eq!(q.to_string(), "P(n) = (1 + 2t) P(n-1) + (t + t^2) P'(n-1)");
    /// ```
    pub fn transport(&self, rule: &Recurrence) -> Result<TransportedRecurrence, TransformError> {
        self.validate()?;
        match self {
            Self::MultiplyPower { .. }
            | Self::DividePower { .. }
            | Self::Affine { .. }
            | Self::Reversal { .. }
            | Self::FPolynomial { .. }
            | Self::Mobius { .. } => {
                let substitution = self.substitution()?.expect("a coordinate map");
                transport_scalar(rule, &substitution, &self.gauge()?)
                    .map(TransportedRecurrence::Scalar)
            }
            Self::Gamma { schedule } => {
                let gamma = transport_gamma(rule, schedule)?;
                match gamma.scalar {
                    Some(scalar) => Ok(TransportedRecurrence::Scalar(scalar)),
                    None => Err(TransformError::TransportUnsupported(
                        gamma
                            .obstruction
                            .unwrap_or_else(|| "no scalar gamma rule".into()),
                    )),
                }
            }
            Self::SymmetricDecomposition { schedule } => {
                transport_symmetric_decomposition(rule, schedule).map(TransportedRecurrence::Vector)
            }
            Self::Diagonal { .. } | Self::DerivativeAntecedent { .. } | Self::Derivative { .. } => {
                Err(TransformError::TransportUnsupported(
                    "this map transports implicit relations: use transport_implicit or \
                     implicit_transform_recurrence"
                        .into(),
                ))
            }
        }
    }

    /// The exact row map, inverse, image conditions and conventions.
    pub fn record(&self) -> TransformationRecord {
        let scalar = || {
            vec![TransformComponent {
                name: "Q".into(),
                symmetry_index: None,
            }]
        };
        let schedule_text = self.schedule().map(DegreeSchedule::formula);
        let domain = self
            .schedule()
            .map_or_else(|| "all n >= 0".to_string(), DegreeSchedule::domain_text);
        let mobius_transport = "scalar: exact Möbius-gauge transport; a common power of ell is \
                                cleared, never evaluated at the pole";
        let (row_map, inverse, conditions, components, transport) = match self {
            Self::MultiplyPower { power } => (
                format!("Q_n(x) = x^{power} P_n(x)"),
                format!("P_n(x) = Q_n(x) / x^{power}"),
                vec![format!("inverse: x^{power} divides Q_n")],
                scalar(),
                "scalar: exact monomial-gauge transport".to_string(),
            ),
            Self::DividePower { power } => (
                format!("Q_n(x) = P_n(x) / x^{power}"),
                format!("P_n(x) = x^{power} Q_n(x)"),
                vec![format!("map: x^{power} divides P_n")],
                scalar(),
                "scalar: exact monomial-gauge transport (an identity of rational functions)"
                    .to_string(),
            ),
            Self::Reversal { .. } => (
                "Q_n(x) = x^(d_n) P_n(1/x)".to_string(),
                "P_n(x) = x^(d_n) Q_n(1/x)".to_string(),
                vec![
                    "map: deg P_n <= d_n".to_string(),
                    "inverse: deg Q_n <= d_n".to_string(),
                ],
                scalar(),
                mobius_transport.to_string(),
            ),
            Self::Affine { a, b } => (
                format!("Q_n(x) = P_n(({a}) x + ({b}))"),
                format!("P_n(y) = Q_n((y - ({b}))/({a}))"),
                Vec::new(),
                scalar(),
                "scalar: exact affine transport".to_string(),
            ),
            Self::Mobius { a, b, c, e, .. } => (
                format!(
                    "Q_n(x) = (({c}) x + ({e}))^(d_n) P_n((({a}) x + ({b}))/(({c}) x + ({e})))"
                ),
                format!(
                    "P_n(y) = (({a})({e}) - ({b})({c}))^(-d_n) (({a}) - ({c}) y)^(d_n) \
                     Q_n((({e}) y - ({b}))/(({a}) - ({c}) y))"
                ),
                vec![
                    "map: deg P_n <= d_n".to_string(),
                    "inverse: deg Q_n <= d_n".to_string(),
                ],
                scalar(),
                mobius_transport.to_string(),
            ),
            Self::FPolynomial { .. } => (
                "Q_n(x) = (1 + x)^(d_n) P_n(x/(1 + x))  (decomposition::f_polynomial at bound d_n)"
                    .to_string(),
                "P_n(y) = (1 - y)^(d_n) Q_n(y/(1 - y))".to_string(),
                vec![
                    "map: deg P_n <= d_n".to_string(),
                    "inverse: deg Q_n <= d_n".to_string(),
                ],
                scalar(),
                mobius_transport.to_string(),
            ),
            Self::Gamma { .. } => (
                "P_n(x) = sum_i Q_(n,i) x^i (1 + x)^(d_n - 2i)".to_string(),
                "P_n(x) = sum_i Q_(n,i) x^i (1 + x)^(d_n - 2i)".to_string(),
                vec![
                    "map: x^(d_n) P_n(1/x) = P_n (symmetry at index d_n)".to_string(),
                    "inverse: deg Q_n <= floor(d_n/2)".to_string(),
                ],
                scalar(),
                "scalar: the identity is cleared of powers of 1 + x and split into its \
                 symmetric and antisymmetric parts (w = (1 - x)/(1 + x), w^2 = 1 - 4z); \
                 proportional parts give the scalar rule in G, independent parts are reported \
                 as unsupported; valid where the rows involved are symmetric"
                    .to_string(),
            ),
            Self::SymmetricDecomposition { .. } => (
                "P_n(x) = a_n(x) + x b_n(x), (a_n, b_n) = ((P_n - x Pbar_n)/(1 - x), \
                 (Pbar_n - P_n)/(1 - x)) with Pbar_n = x^(d_n) P_n(1/x)"
                    .to_string(),
                "P_n(x) = a_n(x) + x b_n(x)".to_string(),
                vec![
                    "map: deg P_n <= d_n".to_string(),
                    "inverse: x^(d_n) a_n(1/x) = a_n and x^(d_n - 1) b_n(1/x) = b_n".to_string(),
                ],
                vec![
                    TransformComponent {
                        name: "a".into(),
                        symmetry_index: Some("d_n".into()),
                    },
                    TransformComponent {
                        name: "b".into(),
                        symmetry_index: Some("d_n - 1".into()),
                    },
                ],
                "vector: first-order rules without alternating signs become the complete \
                 two-component system in (a_n, b_n); no component is selected"
                    .to_string(),
            ),
            Self::Diagonal { .. } | Self::DerivativeAntecedent { .. } => {
                let (row_map, inverse, conditions, transport) = implicit_map_record(self);
                (row_map, inverse, conditions, scalar(), transport)
            }
            Self::Derivative { order } => {
                let (row_map, inverse, conditions, transport) = lossy_map_record(*order);
                (row_map, inverse, conditions, scalar(), transport)
            }
        };
        TransformationRecord {
            kind: self.kind_name().to_string(),
            row_map,
            inverse,
            image_conditions: conditions,
            degree_schedule: schedule_text,
            index_domain: domain,
            components,
            transport,
        }
    }
}
