//! Exact identities, precision, expression grammar, and resource limits.
use super::*;
use num_rational::BigRational as Q;
use num_traits::{One, Zero};
fn eval(s: &str, n: usize) -> TruncatedSeries {
    SeriesExpression::parse(s).unwrap().evaluate(n).unwrap()
}
#[test]
fn formal_identities_and_precision() {
    let f = eval("z-z^2/2+z^4/3", 10);
    assert_eq!(f.exponential().unwrap().logarithm().unwrap(), f);
    let unit = eval("2+z-z^2", 10);
    assert_eq!(unit.multiply(&unit.inverse().unwrap()), eval("1", 10));
    assert_eq!(eval("exp(log(1+z))", 10), eval("1+z", 10));
    assert_eq!(eval("log(exp(z-z^3))", 10), eval("z-z^3", 10));
    assert_eq!(eval("derivative(exp(z^2))", 10), eval("2*z*exp(z^2)", 10));
    assert_eq!(
        eval("compose(exp(x)-1-x,z^2)", 10),
        eval("exp(z^2)-1-z^2", 10)
    );
    assert_eq!(eval("trunc(exp(z),2)", 8), eval("1+z+z^2/2", 8));
    assert_eq!(eval("derivative(trunc(exp(z),2))", 8), eval("1+z", 8));
    assert_eq!(eval("trunc(derivative(exp(z)),2)", 8), eval("1+z+z^2/2", 8));
    assert_eq!(eval("-z^2+2*z^2", 3), eval("z^2", 3));
    assert_eq!(f.scale(&Q::zero()).order(), 11);
    assert_eq!(f.scale(&Q::zero()), eval("0", 10));
    assert_eq!(f.truncate(4).multiply(&f).order(), 4);
    assert_eq!(f.truncate(4).derivative().order(), 3);
    assert_eq!(f.truncate(4).truncate(20).order(), 4);
    assert_eq!(f.power(0), eval("1", 10));
    assert_eq!(&f + &f, f.scale(&Q::from_integer(2.into())));
    assert_eq!(&f - &f, eval("0", 10));
    assert_eq!(&f * &Q::zero(), eval("0", 10));
    assert_eq!(&f * &f, f.power(2));
    assert_eq!(-&f, f.scale(&-Q::one()));
    assert!(TruncatedSeries::new(vec![]).inverse().is_err());
    assert!(TruncatedSeries::new(vec![Q::one()]).exponential().is_err());
}
#[test]
fn lists_normalizations_and_rationals() {
    assert_eq!(
        eval("coeffs:[0,1,1/2,1/6]", 6),
        eval("coeffs:exponential:[0,1,1,1]", 6)
    );
    assert_eq!(
        eval("coeffs:ordinary:[0,-1/3,2/7]", 3),
        eval("-z/3+2*z^2/7", 3)
    );
    assert_eq!(eval("coeffs:[]", 4), eval("0", 4));
    assert_eq!(eval("z^1/2", 4), eval("z/2", 4));
    assert_eq!(
        eval("log(1/(1-z))", 6),
        eval("coeffs:[0,1,1/2,1/3,1/4,1/5,1/6]", 6)
    );
    assert_eq!(eval("inv(1-z)", 6), eval("1/(1-z)", 6));
    assert_eq!(
        eval("1/2/3", 0).coefficients(),
        &[Q::new(1.into(), 6.into())]
    );
}
#[test]
fn errors_and_bounds() {
    for s in [
        "",
        "sin(z)",
        "exp(z",
        "z**2",
        "2z",
        "z^-1",
        "z^2^3",
        "coeffs:[1,]",
        "coeffs:[1/0]",
        "coeffs:[\"1]",
        "compose(z)",
        "trunc(z,-1)",
        "z💥",
    ] {
        let e = SeriesExpression::parse(s).unwrap_err();
        assert!(e.to_string().contains("byte"), "{s}: {e}");
    }
    for s in [
        "exp(1)",
        "log(z)",
        "log(2+z)",
        "1/z",
        "1/0",
        "compose(z,1+z)",
    ] {
        assert!(
            SeriesExpression::parse(s).unwrap().evaluate(0).is_err(),
            "{s}"
        );
    }
    assert!(SeriesExpression::parse(&"-".repeat(100)).is_err());
    assert!(SeriesExpression::parse("z^9999999999999999999999999999").is_err());
    let limits = SeriesExpressionLimits {
        coefficient_bits: 8,
        ..SeriesExpressionLimits::default()
    };
    assert!(SeriesExpression::parse_with_limits("256*z", &limits).is_err());
    assert!(SeriesExpression::parse_with_limits("128^8*z", &limits)
        .unwrap()
        .evaluate_with_limits(2, &limits)
        .is_err());
    let limits = SeriesExpressionLimits {
        degree: 3,
        ..limits
    };
    assert!(SeriesExpression::parse("derivative(exp(z))")
        .unwrap()
        .evaluate_with_limits(3, &limits)
        .is_err());
    assert!(SeriesExpression::parse("z")
        .unwrap()
        .evaluate(usize::MAX)
        .is_err());
}
