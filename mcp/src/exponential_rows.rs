//! Bounded adapter for exact exponential-family generation and finite checks.
use crate::parsing::invalid_params;
use polytool::exponential_rows::*;
use rmcp::ErrorData;
use schemars::JsonSchema;
use serde::Deserialize;

#[derive(Clone, Debug, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ExponentialRowsToolRequest {
    /// Request: {family:{unmarked,marked},rows:[n,...],
    /// check?:none|real_rooted|interlacing|first_failure,summary?:bool,finite_rank?:[]}.
    /// Block specs: `{kind:sizes,sizes:[1,2]}`, `{kind:finite,weights:["1","1/2"]}`,
    /// `{kind:tail,minimum:2}`, `{kind:exponential,scale:"1",rate:"1",minimum:2}`,
    /// `{kind:polynomial_substitution,seed:SPEC,coefficients:["2","3","1"]}`,
    /// or `{kind:ode_substitution,seed:SPEC,a:"1",b:"0",lambda:"1"}`.
    /// Expressions: `{kind:series_expression,expression:"exp(z)-1"}`;
    /// +,-,*,/,^,exp,log,inv,compose,trunc,derivative,coeffs:[...].
    /// Expressions: rank <=128, bytes <=2048, nodes <=128,
    /// depth <=16, exponents/list degrees <=64, literal and intermediate
    /// numerator/denominator <=4096 bits; derivative lookahead <=144.
    /// Bounds: JSON <=65536 bytes; 1..512 rows, rank <=1000; generated degree
    /// <=256, checked degree <=64; output coefficients <=65536. Structured-weight
    /// rationals: text <=80 bytes, numerator/denominator <=256 bits. Composition depth
    /// <=4 and composed ranks <=128; polynomial composition degree <=64.
    /// <=16 finite-rank references.
    pub request: serde_json::Value,
}
fn coefficient(s: &str) -> Result<(), ErrorData> {
    if s.len() > 80 {
        return Err(invalid_params("rational text exceeds 80 bytes"));
    }
    let q =
        polytool::recurrence::parse_rational_coeff(s).map_err(|e| invalid_params(e.to_string()))?;
    if q.numer().bits() > 256 || q.denom().bits() > 256 {
        return Err(invalid_params("rational coefficient exceeds 256 bits"));
    }
    Ok(())
}
fn spec(s: &BlockWeightSpec, depth: usize, max_n: usize) -> Result<(), ErrorData> {
    if depth > 4 {
        return Err(invalid_params("composition depth exceeds 4"));
    }
    match s {
        BlockWeightSpec::SeriesExpression { expression } => {
            if max_n > 128 {
                return Err(invalid_params("expression rank exceeds 128"));
            }
            let limits = polytool::power_series::SeriesExpressionLimits {
                bytes: 2048,
                nodes: 128,
                depth: 16,
                exponent: 64,
                degree: 144,
                coefficient_bits: 4096,
            };
            polytool::power_series::SeriesExpression::parse_with_limits(expression, &limits)
                .and_then(|e| e.evaluate_with_limits(max_n, &limits))
                .map_err(|e| invalid_params(e.to_string()))?;
        }
        BlockWeightSpec::Finite { weights } => {
            if weights.len() > 1000 {
                return Err(invalid_params("weight count exceeds 1000"));
            }
            for w in weights {
                coefficient(w)?;
            }
        }
        BlockWeightSpec::Sizes { sizes } => {
            if sizes.len() > 1000 || sizes.iter().any(|&r| r == 0 || r > 1000) {
                return Err(invalid_params(
                    "sizes must be in 1..=1000, at most 1000 entries",
                ));
            }
        }
        BlockWeightSpec::Tail { minimum } | BlockWeightSpec::Exponential { minimum, .. } => {
            if *minimum == 0 || *minimum > 1000 {
                return Err(invalid_params("minimum must be in 1..=1000"));
            }
        }
        BlockWeightSpec::PolynomialSubstitution { seed, coefficients } => {
            if max_n > 128 || coefficients.len() > 65 {
                return Err(invalid_params(
                    "polynomial composition requires rank <=128 and degree <=64",
                ));
            }
            spec(seed, depth + 1, max_n)?;
            for c in coefficients {
                coefficient(c)?;
            }
        }
        BlockWeightSpec::OdeSubstitution { seed, a, b, lambda }
        | BlockWeightSpec::RationalOdeSubstitution {
            seed, a, b, lambda, ..
        } => {
            if max_n > 128 {
                return Err(invalid_params("composed rank exceeds 128"));
            }
            spec(seed, depth + 1, max_n)?;
            for s in [a, b, lambda] {
                coefficient(s)?;
            }
            if let BlockWeightSpec::RationalOdeSubstitution { d, .. } = s {
                coefficient(d)?;
            }
        }
    }
    if let BlockWeightSpec::Exponential { scale, rate, .. } = s {
        coefficient(scale)?;
        coefficient(rate)?;
    }
    Ok(())
}
pub(crate) fn run(request: ExponentialRowsToolRequest) -> Result<serde_json::Value, ErrorData> {
    if request.request.to_string().len() > 65536 {
        return Err(invalid_params("request JSON exceeds 65536 bytes"));
    }
    let request: ExponentialRowsRequest =
        serde_json::from_value(request.request).map_err(|e| invalid_params(e.to_string()))?;
    if request.rows.is_empty() || request.rows.len() > 512 || request.rows.iter().any(|&n| n > 1000)
    {
        return Err(invalid_params("require 1..512 rows with ranks <=1000"));
    }
    let max_n = *request.rows.iter().max().expect("nonempty rows");
    spec(&request.family.unmarked, 0, max_n)?;
    spec(&request.family.marked, 0, max_n)?;
    let weights = request
        .family
        .marked
        .block_weights(max_n)
        .map_err(|e| invalid_params(e.to_string()))?;
    use num_bigint::BigInt;
    let minimum = weights
        .iter()
        .position(|w| w.numer() != &BigInt::from(0))
        .unwrap_or(max_n + 1);
    let limit = if request.check == ExponentialRowCheck::None {
        256
    } else {
        64
    };
    let maximum_degree = max_n / minimum.max(1);
    if maximum_degree > limit {
        return Err(invalid_params(format!("degree bound exceeds {limit}")));
    }
    if request
        .rows
        .iter()
        .map(|n| n / minimum.max(1) + 1)
        .sum::<usize>()
        > 65536
    {
        return Err(invalid_params("output coefficient count exceeds 65536"));
    }
    if request.finite_rank.len() > 16 {
        return Err(invalid_params("finite-rank references exceed 16"));
    }
    for reference in &request.finite_rank {
        if reference.reference.len() > 65 {
            return Err(invalid_params("reference degree exceeds 64"));
        }
        for w in &reference.reference {
            coefficient(w)?;
        }
        if let Some(p) = &reference.pencil {
            if p.factor.len() > 64 {
                return Err(invalid_params("pencil degree exceeds 63"));
            }
            for w in p.factor.iter().chain([&p.a, &p.b]) {
                coefficient(w)?;
            }
        }
    }
    let report = check_exponential_rows(&request).map_err(|e| invalid_params(e.to_string()))?;
    serde_json::to_value(report).map_err(|e| invalid_params(e.to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::PolynomialToolsServer;
    use rmcp::handler::server::wrapper::Parameters;
    use serde_json::json;
    #[test]
    fn adapter_and_schema_reject_excessive_inputs() {
        let request = json!({"family":{"unmarked":{"kind":"sizes","sizes":[1,2]},"marked":{"kind":"sizes","sizes":[7]}},"rows":[304,305],"check":"first_failure"});
        let server = PolynomialToolsServer::new();
        let rmcp::Json(result) = server
            .generate_exponential_rows(Parameters(ExponentialRowsToolRequest {
                request: request.clone(),
            }))
            .unwrap();
        assert_eq!(result["first_failure"]["index"], 305);
        for patch in [
            json!({"rows":[1001]}),
            json!({"rows":[500],"family":{"unmarked":{"kind":"sizes","sizes":[]},"marked":{"kind":"tail","minimum":1}}}),
            json!({"rows":[1],"family":{"unmarked":{"kind":"finite","weights":["1/0"]},"marked":{"kind":"tail","minimum":1}}}),
        ] {
            let mut bad = request.clone();
            for (k, v) in patch.as_object().unwrap() {
                bad[k] = v.clone();
            }
            assert!(server
                .generate_exponential_rows(Parameters(ExponentialRowsToolRequest { request: bad }))
                .is_err());
        }
    }
}

#[cfg(test)]
mod rational_ode_tests {
    use super::*;
    use crate::PolynomialToolsServer;
    use rmcp::{handler::server::wrapper::Parameters, Json};
    use serde_json::json;
    #[test]
    fn rational_ode_is_bounded_and_never_labelled_as_a_lean_theorem() {
        let server = PolynomialToolsServer::new();
        let mut request = json!({"family":{"unmarked":{"kind":"sizes","sizes":[1]},"marked":{"kind":"rational_ode_substitution","seed":{"kind":"sizes","sizes":[1]},"a":"1","b":"0","lambda":"0","d":"1"}},"rows":[6,8],"check":"real_rooted"});
        let Json(report) = server
            .generate_exponential_rows(Parameters(ExponentialRowsToolRequest {
                request: request.clone(),
            }))
            .unwrap();
        assert_eq!(report["routes"][0]["evidence"], "paper_proof");
        assert!(report["first_failure"].is_null());
        request["rows"] = json!([129]);
        request["check"] = json!("none");
        assert!(server
            .generate_exponential_rows(Parameters(ExponentialRowsToolRequest { request }))
            .is_err());
    }
}

#[cfg(test)]
mod expression_tests {
    use super::*;
    use serde_json::json;
    #[test]
    fn exact_expressions_and_resource_bounds() {
        let request = json!({"family":{"unmarked":{"kind":"series_expression","expression":"z"},"marked":{"kind":"series_expression","expression":"log(1/(1-z))"}},"rows":[6,8],"check":"real_rooted"});
        let report = run(ExponentialRowsToolRequest {
            request: request.clone(),
        })
        .unwrap();
        assert_eq!(report["routes"][0]["evidence"], "paper_proof");
        assert!(report["first_failure"].is_null());
        for b in [
            "z^65",
            "exp(1)",
            "1+z",
            "compose(z,1+z)",
            "exp(z",
            &"-".repeat(129),
            &"1+".repeat(1100),
        ] {
            let mut bad = request.clone();
            bad["family"]["marked"]["expression"] = json!(b);
            assert!(
                run(ExponentialRowsToolRequest { request: bad }).is_err(),
                "{b}"
            );
        }
        let mut bad = request.clone();
        bad["rows"] = json!([129]);
        assert!(run(ExponentialRowsToolRequest { request: bad }).is_err());
        let mut bad = request.clone();
        bad["family"]["marked"]["expression"] = json!(format!("({})^64*z", "9".repeat(100)));
        assert!(run(ExponentialRowsToolRequest { request: bad }).is_err());
        let mut counterexample = request;
        counterexample["family"]["marked"]["expression"] = json!("z^2*(1+z^4)");
        let report = run(ExponentialRowsToolRequest {
            request: counterexample,
        })
        .unwrap();
        assert_eq!(report["first_failure"]["index"], 6);
    }
}
