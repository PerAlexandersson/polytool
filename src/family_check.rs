//! Select the first failed requirement of a finite polynomial family.
//!
//! Adapters supply already computed property and interlacing verdicts to
//! [`first_family_failure`]. This module does not recompute them or claim an
//! all-index theorem. Failures are typed so adapters can retain their wording.

/// Required properties of each member and each consecutive pair.
///
/// All requirements default to false; adapters choose their own defaults.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct FamilyRequirements {
    /// Require every member to be real-rooted.
    pub real_rooted: bool,
    /// Require every member to have simple roots.
    pub simple_roots: bool,
    /// Require every member to be palindromic.
    pub palindromic: bool,
    /// Require every member to be gamma-positive.
    pub gamma_positive: bool,
    /// Require every member's coefficients to be unimodal.
    pub unimodal: bool,
    /// Require every member's coefficients to be log-concave.
    pub log_concave: bool,
    /// Require every member's coefficients to be ultra-log-concave.
    pub ultra_log_concave: bool,
    /// Require each supplied consecutive pair to weakly interlace.
    pub consecutive_weak_interlacing: bool,
}

/// Already computed verdicts for one polynomial.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct FamilyPolynomialProperties {
    /// The caller's polynomial index, carried unchanged into a failure.
    pub index: usize,
    /// Whether the polynomial is real-rooted.
    pub real_rooted: bool,
    /// Whether the polynomial has simple roots.
    pub simple_roots: bool,
    /// Whether the polynomial is palindromic.
    pub palindromic: bool,
    /// Whether the polynomial is gamma-positive.
    pub gamma_positive: bool,
    /// Whether its coefficients are unimodal.
    pub unimodal: bool,
    /// Whether its coefficients are log-concave.
    pub log_concave: bool,
    /// Whether its coefficients are ultra-log-concave.
    pub ultra_log_concave: bool,
}

/// A required property of a polynomial, in failure-priority order.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FamilyProperty {
    /// Real-rootedness.
    RealRooted,
    /// Simplicity of roots.
    SimpleRoots,
    /// Palindromicity.
    Palindromic,
    /// Gamma positivity.
    GammaPositive,
    /// Unimodality of coefficients.
    Unimodal,
    /// Log-concavity of coefficients.
    LogConcave,
    /// Ultra-log-concavity of coefficients.
    UltraLogConcave,
}

/// The first failed requirement of a finite family.
///
/// The display wording is suitable for a report; adapters may format variants
/// themselves to preserve an established output schema.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FamilyFailure {
    /// A polynomial failed a required property.
    Polynomial {
        /// The caller's polynomial index.
        index: usize,
        /// The failed property.
        property: FamilyProperty,
    },
    /// A supplied consecutive pair did not have verdict `Some(true)`.
    WeakInterlacing {
        /// The caller's index of the left polynomial.
        left_index: usize,
        /// The caller's index of the right polynomial.
        right_index: usize,
    },
    /// A supplied finite Lace check failed.
    Lace {
        /// Largest minor size checked in the finite truncation.
        max_minor_size: usize,
    },
    /// A supplied recurrence search did not find a recurrence.
    Recurrence,
}

impl std::fmt::Display for FamilyFailure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Polynomial { index, property } => {
                let description = match property {
                    FamilyProperty::RealRooted => "is not real-rooted",
                    FamilyProperty::SimpleRoots => "does not have simple roots",
                    FamilyProperty::Palindromic => "is not palindromic",
                    FamilyProperty::GammaPositive => "is not gamma-positive",
                    FamilyProperty::Unimodal => "is not unimodal",
                    FamilyProperty::LogConcave => "is not log-concave",
                    FamilyProperty::UltraLogConcave => "is not ultra-log-concave",
                };
                write!(f, "polynomial {index} {description}")
            }
            Self::WeakInterlacing {
                left_index,
                right_index,
            } => write!(
                f,
                "consecutive pair {left_index}-{right_index} is not weakly interlacing"
            ),
            Self::Lace { max_minor_size } => write!(
                f,
                "finite Lace truncation is not TNN up to minors of size {max_minor_size}"
            ),
            Self::Recurrence => f.write_str("recurrence search did not find a recurrence"),
        }
    }
}

/// Select the first failed requirement without recomputing any verdict.
///
/// Members are visited in input order, with properties tested in the order
/// [`FamilyProperty`] lists them. Next come supplied consecutive pairs, the
/// optional finite Lace check, and the optional recurrence search.
/// Pair tuples are `(left_index, right_index, weak_interlacing)`; `None` fails
/// a required pair check, just like `Some(false)`. Missing pairs are not
/// synthesized. `lace` is `(totally_nonnegative, max_minor_size)` and
/// `recurrence` is whether a recurrence was found; absent checks are skipped.
///
/// This only reports on the supplied finite data.
///
/// ```
/// use polytool::family_check::{
///     first_family_failure, FamilyFailure, FamilyPolynomialProperties,
///     FamilyProperty, FamilyRequirements,
/// };
///
/// let requirements = FamilyRequirements { real_rooted: true, ..Default::default() };
/// let members = [FamilyPolynomialProperties { index: 3, ..Default::default() }];
/// assert_eq!(
///     first_family_failure(&requirements, members, [], None, None),
///     Some(FamilyFailure::Polynomial { index: 3, property: FamilyProperty::RealRooted }),
/// );
/// ```
pub fn first_family_failure(
    requirements: &FamilyRequirements,
    members: impl IntoIterator<Item = FamilyPolynomialProperties>,
    pairs: impl IntoIterator<Item = (usize, usize, Option<bool>)>,
    lace: Option<(bool, usize)>,
    recurrence: Option<bool>,
) -> Option<FamilyFailure> {
    use FamilyProperty::*;
    for member in members {
        let checks = [
            (requirements.real_rooted, member.real_rooted, RealRooted),
            (requirements.simple_roots, member.simple_roots, SimpleRoots),
            (requirements.palindromic, member.palindromic, Palindromic),
            (
                requirements.gamma_positive,
                member.gamma_positive,
                GammaPositive,
            ),
            (requirements.unimodal, member.unimodal, Unimodal),
            (requirements.log_concave, member.log_concave, LogConcave),
            (
                requirements.ultra_log_concave,
                member.ultra_log_concave,
                UltraLogConcave,
            ),
        ];
        for (required, passed, property) in checks {
            if required && !passed {
                return Some(FamilyFailure::Polynomial {
                    index: member.index,
                    property,
                });
            }
        }
    }
    if requirements.consecutive_weak_interlacing {
        for (left_index, right_index, weak) in pairs {
            if weak != Some(true) {
                return Some(FamilyFailure::WeakInterlacing {
                    left_index,
                    right_index,
                });
            }
        }
    }
    if let Some((false, max_minor_size)) = lace {
        return Some(FamilyFailure::Lace { max_minor_size });
    }
    if recurrence == Some(false) {
        return Some(FamilyFailure::Recurrence);
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn failures_follow_member_then_property_order() {
        use FamilyProperty::*;
        let requirements = FamilyRequirements {
            real_rooted: true,
            simple_roots: true,
            palindromic: true,
            gamma_positive: true,
            unimodal: true,
            log_concave: true,
            ultra_log_concave: true,
            consecutive_weak_interlacing: true,
        };
        let mut member = FamilyPolynomialProperties {
            index: 7,
            ..Default::default()
        };
        for property in [
            RealRooted,
            SimpleRoots,
            Palindromic,
            GammaPositive,
            Unimodal,
            LogConcave,
            UltraLogConcave,
        ] {
            assert_eq!(
                first_family_failure(
                    &requirements,
                    [member],
                    [(7, 9, None)],
                    Some((false, 3)),
                    Some(false)
                ),
                Some(FamilyFailure::Polynomial { index: 7, property }),
            );
            match property {
                RealRooted => member.real_rooted = true,
                SimpleRoots => member.simple_roots = true,
                Palindromic => member.palindromic = true,
                GammaPositive => member.gamma_positive = true,
                Unimodal => member.unimodal = true,
                LogConcave => member.log_concave = true,
                UltraLogConcave => member.ultra_log_concave = true,
            }
        }
        let later = FamilyPolynomialProperties {
            index: 9,
            ..Default::default()
        };
        assert_eq!(
            first_family_failure(&requirements, [member, later], [], None, None),
            Some(FamilyFailure::Polynomial {
                index: 9,
                property: RealRooted
            }),
        );
        // An earlier member's lower-priority failure precedes a later member.
        member.ultra_log_concave = false;
        assert_eq!(
            first_family_failure(&requirements, [member, later], [], None, None),
            Some(FamilyFailure::Polynomial {
                index: 7,
                property: UltraLogConcave
            }),
        );
    }

    #[test]
    fn optional_checks_follow_pairs_lace_recurrence_order() {
        let requirements = FamilyRequirements {
            consecutive_weak_interlacing: true,
            ..Default::default()
        };
        for weak in [None, Some(false)] {
            assert_eq!(
                first_family_failure(
                    &requirements,
                    [],
                    [(2, 4, Some(true)), (4, 8, weak)],
                    Some((false, 3)),
                    Some(false)
                ),
                Some(FamilyFailure::WeakInterlacing {
                    left_index: 4,
                    right_index: 8
                }),
            );
        }
        assert_eq!(
            first_family_failure(&requirements, [], [], Some((false, 3)), Some(false)),
            Some(FamilyFailure::Lace { max_minor_size: 3 }),
        );
        assert_eq!(
            first_family_failure(&requirements, [], [], Some((true, 3)), Some(false)),
            Some(FamilyFailure::Recurrence),
        );
        assert_eq!(
            first_family_failure(&requirements, [], [], None, None),
            None
        );
        assert_eq!(
            first_family_failure(
                &FamilyRequirements::default(),
                [FamilyPolynomialProperties::default()],
                [(0, 1, None)],
                None,
                Some(true)
            ),
            None,
        );
    }
}
