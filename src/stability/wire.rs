//! Canonical JSON and exact replay shared by the CLI and MCP adapters.

use super::{
    parse::{parse_gaussian, parse_rational, ParseError, ParseLimits},
    *,
};
use exact_algebra::inertia::DiagonalBlock;
use serde_json::{json, Value};

/// Versioned full-report schema.
pub const SCHEMA: &str = "polytool.stability.v1";

/// Parse named/circular region selectors, rejecting ambiguous combinations.
pub fn parse_region(
    region: Option<&str>,
    mobius: Option<&[GaussianRational]>,
    center: Option<&str>,
    radius: Option<&str>,
    limits: ParseLimits,
) -> Result<StabilityRegion, ParseError> {
    let fail = |s: &str| ParseError(s.into());
    if mobius.is_some() && (region.is_some() || center.is_some() || radius.is_some())
        || region.is_some() && (center.is_some() || radius.is_some())
    {
        return Err(fail(
            "region, mobius and centered-disk selectors are exclusive",
        ));
    }
    let result = if let Some(m) = mobius {
        if m.len() != 4 {
            return Err(fail("mobius requires exactly four Gaussian coefficients"));
        }
        for coefficient in m {
            limits.validate(std::slice::from_ref(coefficient))?;
        }
        StabilityRegion::from(
            MobiusTransform::new(m[0].clone(), m[1].clone(), m[2].clone(), m[3].clone())
                .map_err(|e| fail(&e.to_string()))?,
        )
    } else if center.is_some() || radius.is_some() {
        let (Some(center), Some(radius)) = (center, radius) else {
            return Err(fail("disk-center and radius are required together"));
        };
        StabilityRegion::CenteredDisk {
            center: parse_gaussian(center, limits)?,
            radius: parse_rational(radius, limits)?,
        }
    } else {
        match region.unwrap_or("upper") {
            "upper" => StabilityRegion::Upper,
            "lower" => StabilityRegion::Lower,
            "left" => StabilityRegion::Left,
            "right" => StabilityRegion::Right,
            "disk" => StabilityRegion::Disk,
            "exterior" => StabilityRegion::Exterior,
            _ => return Err(fail("unknown stability region")),
        }
    };
    result
        .mobius_transform()
        .map_err(|e| fail(&e.to_string()))?;
    Ok(result)
}

/// Resolve an optional convention using the region's default.
pub fn parse_convention(
    value: Option<&str>,
    region: &StabilityRegion,
) -> Result<StabilityConvention, ParseError> {
    match value {
        None => Ok(region.default_convention()),
        Some("avoid") => Ok(StabilityConvention::Avoid),
        Some("inside") => Ok(StabilityConvention::Inside),
        _ => Err(ParseError("unknown stability convention".into())),
    }
}

/// Exact Gaussian number with rational parts encoded as strings.
pub fn gaussian_json(value: &GaussianRational) -> Value {
    json!({"real": value.real.to_string(), "imaginary":value.imaginary.to_string()})
}

/// Parse the common Gaussian JSON input (integer, exact literal, or parts).
pub fn parse_gaussian_json(
    value: &Value,
    limits: ParseLimits,
) -> Result<GaussianRational, ParseError> {
    match value {
        Value::Number(n) if n.is_i64() || n.is_u64() => parse_gaussian(&n.to_string(), limits),
        Value::String(s) => parse_gaussian(s, limits),
        Value::Object(parts) if parts.len() == 2 => {
            let part = |key| {
                parts
                    .get(key)
                    .and_then(Value::as_str)
                    .ok_or_else(|| ParseError("Gaussian parts must be rational strings".into()))
            };
            Ok(GaussianRational::new(
                parse_rational(part("real")?, limits)?,
                parse_rational(part("imaginary")?, limits)?,
            ))
        }
        _ => Err(ParseError(
            "expected exact Gaussian integer, string or {real,imaginary}; JSON floats are rejected"
                .into(),
        )),
    }
}

/// Compact per-polynomial result; no input echo or certificate.
pub fn summary(index: usize, report: &StabilityReport) -> Value {
    let counts = report
        .counts
        .map(|c| json!({"inside":c.inside,"boundary":c.boundary,"outside":c.outside}));
    let meaning = match report.convention {
        StabilityConvention::Avoid => "no root in the open region",
        StabilityConvention::Inside => "all roots in the open region",
    };
    let method = match report.method {
        StabilityMethod::ZeroPolynomial => "zero_polynomial",
        StabilityMethod::RealSturm => "sturm",
        StabilityMethod::BezoutianSturm => "bezoutian_sturm",
    };
    json!({"index":index,"degree":report.degree,"counts":counts,"verdict":report.verdict,"meaning":meaning,"method":method})
}

fn options(report: &StabilityReport) -> Value {
    let mut value = match &report.region {
        StabilityRegion::CenteredDisk { center, radius } => {
            json!({"disk_center":gaussian_json(center),"radius":radius.to_string()})
        }
        StabilityRegion::Mobius(m) => {
            json!({"mobius":m.coefficients().into_iter().map(gaussian_json).collect::<Vec<_>>()})
        }
        StabilityRegion::Upper => json!({"region":"upper"}),
        StabilityRegion::Lower => json!({"region":"lower"}),
        StabilityRegion::Left => json!({"region":"left"}),
        StabilityRegion::Right => json!({"region":"right"}),
        StabilityRegion::Disk => json!({"region":"disk"}),
        StabilityRegion::Exterior => json!({"region":"exterior"}),
    };
    value["convention"] = json!(match report.convention {
        StabilityConvention::Avoid => "avoid",
        StabilityConvention::Inside => "inside",
    });
    value
}

/// Full canonical item including every exact certificate entry.
pub fn full_item(
    index: usize,
    coefficients: &[GaussianRational],
    report: &StabilityReport,
) -> Value {
    let strings = |v: &[num_bigint::BigInt]| v.iter().map(ToString::to_string).collect::<Vec<_>>();
    let certificate = report.certificate.as_ref().map(|c| json!({
        "transformed": {"coefficients":c.transformed.coefficients.iter().map(gaussian_json).collect::<Vec<_>>(),"original_degree":c.transformed.original_degree,"boundary_at_infinity":c.transformed.boundary_at_infinity},
        "real_part":strings(&c.real_part),"imaginary_part":strings(&c.imaginary_part),"common_factor":strings(&c.common_factor),
        "real_root_layers":c.real_root_layers.iter().map(|l|json!({"polynomial":strings(&l.polynomial),"distinct_real_roots":l.distinct_real_roots})).collect::<Vec<_>>(),
        "bezoutian":c.bezoutian.iter().map(|row|strings(row)).collect::<Vec<_>>(),
        "decomposition":{"permutation":c.decomposition.permutation,"lower":c.decomposition.lower.iter().map(|row|row.iter().map(ToString::to_string).collect::<Vec<_>>()).collect::<Vec<_>>(),
        "blocks":c.decomposition.blocks.iter().map(|b|match b {DiagonalBlock::Single(q)=>json!({"single":q.to_string()}),DiagonalBlock::Pair(q)=>json!({"pair":q.to_string()})}).collect::<Vec<_>>()}
    }));
    let mut value = summary(index, report);
    let end = coefficients
        .iter()
        .rposition(|c| !c.is_zero())
        .map_or(0, |i| i + 1);
    value["input"] = json!(coefficients[..end]
        .iter()
        .map(gaussian_json)
        .collect::<Vec<_>>());
    value["options"] = options(report);
    value["certificate"] = json!(certificate);
    value
}

/// Wrap canonical full items in their versioned replay envelope.
pub fn full_report(items: Vec<Value>) -> Value {
    json!({"schema":SCHEMA,"version":env!("CARGO_PKG_VERSION"),"items":items})
}

/// Recompute a bounded untrusted full report and compare all mathematical fields.
/// The schema selects the replay contract. `version` is a nonempty producer
/// version string, not a compatibility gate; it need not match this executable.
/// Malformed input, extra fields and altered evidence are rejected.
pub fn verify_report(value: &Value, limits: ParseLimits) -> Result<bool, ParseError> {
    let bad = || ParseError("malformed stability report".into());
    if value["schema"] != SCHEMA {
        return Err(bad());
    }
    let version = value["version"]
        .as_str()
        .filter(|s| !s.is_empty())
        .ok_or_else(bad)?;
    let items = value["items"].as_array().ok_or_else(bad)?;
    if items.is_empty() {
        return Err(bad());
    }
    let mut expected = Vec::new();
    for (index, item) in items.iter().enumerate() {
        let input = item["input"].as_array().ok_or_else(bad)?;
        if input.len() > limits.max_degree.saturating_add(1) {
            return Err(bad());
        }
        let coefficients = input
            .iter()
            .map(|c| parse_gaussian_json(c, limits))
            .collect::<Result<Vec<_>, _>>()?;
        limits.validate(&coefficients)?;
        let opts = &item["options"];
        let mobius = opts
            .get("mobius")
            .map(|m| {
                m.as_array().ok_or_else(bad).and_then(|a| {
                    a.iter()
                        .map(|c| parse_gaussian_json(c, limits))
                        .collect::<Result<Vec<_>, _>>()
                })
            })
            .transpose()?;
        let center = opts
            .get("disk_center")
            .map(|c| parse_gaussian_json(c, limits))
            .transpose()?;
        let center_text = center
            .as_ref()
            .map(|c| format!("({})+({})*I", c.real, c.imaginary));
        let region = parse_region(
            opts.get("region").and_then(Value::as_str),
            mobius.as_deref(),
            center_text.as_deref(),
            opts.get("radius").and_then(Value::as_str),
            limits,
        )?;
        let convention = parse_convention(opts.get("convention").and_then(Value::as_str), &region)?;
        let report = check_univariate_stability_with_convention(&coefficients, &region, convention)
            .map_err(|e| ParseError(e.to_string()))?;
        expected.push(full_item(index, &coefficients, &report));
    }
    let mut expected = full_report(expected);
    expected["version"] = json!(version);
    Ok(expected == *value)
}
