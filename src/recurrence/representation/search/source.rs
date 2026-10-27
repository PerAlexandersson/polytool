//! Source resolution: validation of a request and construction of the source rows.

use super::*;

pub(super) struct Source {
    pub(super) kind: &'static str,
    pub(super) rule: Option<(Recurrence, Vec<Vec<Q>>, RecurrenceJson)>,
    pub(super) origin: usize,
    pub(super) rows: Vec<Vec<Q>>,
}

pub(super) fn validate(request: &RepresentationSearchRequest) -> Result<(), RecurrenceError> {
    let limits = &request.limits;
    if request.transformations.len() > MAX_REPRESENTATION_TRANSFORMATIONS {
        return Err(RecurrenceError::LimitExceeded(format!(
            "at most {MAX_REPRESENTATION_TRANSFORMATIONS} transformations are supported"
        )));
    }
    for (i, t) in request.transformations.iter().enumerate() {
        t.validate()
            .map_err(|e| RecurrenceError::InvalidInput(format!("transformations[{i}]: {e}")))?;
    }
    let bounds = &request.search;
    let degree = [
        bounds.max_rec_len,
        bounds.max_var_deg,
        bounds.max_idx_deg,
        bounds.max_diff_deg,
        bounds.vector_denominator_var_deg,
        bounds.vector_denominator_idx_deg,
    ]
    .into_iter()
    .max()
    .unwrap_or(0);
    if degree > limits.max_search_degree {
        return Err(RecurrenceError::LimitExceeded(format!(
            "search bound {degree} exceeds max_search_degree = {}",
            limits.max_search_degree
        )));
    }
    if bounds.max_candidates > limits.max_search_candidates {
        return Err(RecurrenceError::LimitExceeded(format!(
            "max_candidates {} exceeds max_search_candidates = {}",
            bounds.max_candidates, limits.max_search_candidates
        )));
    }
    if limits.max_depth > 16 || limits.max_representations > 4096 {
        return Err(RecurrenceError::LimitExceeded(
            "max_depth is at most 16 and max_representations at most 4096".into(),
        ));
    }
    if let Some(root) = &request.root_domain {
        root.validate()?;
    }
    if request.space_shapes.len() > limits.max_space_shapes {
        return Err(RecurrenceError::InvalidInput(
            "more space_shapes than max_space_shapes".into(),
        ));
    }
    Ok(())
}

pub(super) fn resolve(request: &RepresentationSearchRequest) -> Result<Source, RecurrenceError> {
    let limits = &request.limits;
    let defining = match (&request.source_recurrence, &request.oeis_id) {
        (Some(_), Some(_)) => {
            return Err(RecurrenceError::InvalidInput(
                "give at most one of source_recurrence and oeis_id".into(),
            ))
        }
        (Some(json), None) => {
            let (rule, first, seeds) = json.to_recurrence_parts()?;
            Some(("defining_rule", rule, first, seeds))
        }
        (None, Some(id)) => {
            let definition = crate::oeis::by_id(id).ok_or_else(|| {
                RecurrenceError::InvalidInput(format!("unknown OEIS family `{id}`"))
            })?;
            let (rule, seeds, first_row) = definition.recurrence_export_parts()?;
            let first = usize::try_from(first_row).map_err(|_| {
                RecurrenceError::InvalidInput(format!(
                    "OEIS family `{id}` starts at a negative row index"
                ))
            })?;
            Some(("oeis", rule, first, seeds))
        }
        (None, None) => None,
    };
    let explicit = request
        .rows
        .as_ref()
        .map(|rows| {
            rows.iter()
                .map(|r| decode_poly(r))
                .collect::<Result<Vec<_>, _>>()
        })
        .transpose()?;
    let (kind, rule, origin, rows) = match defining {
        Some((kind, rule, first, seeds)) => {
            if rule.terms.iter().any(|t| t.offset == 0) {
                return Err(RecurrenceError::InvalidInput(
                    "source rule terms need positive offsets".into(),
                ));
            }
            if seeds.is_empty() || rule.max_offset() > seeds.len() {
                return Err(RecurrenceError::InvalidInput(
                    "source rule needs at least max-offset seeds".into(),
                ));
            }
            if rule
                .terms
                .iter()
                .any(|t| t.deriv_order > TransformLimits::default().max_rule_derivative_order)
            {
                return Err(RecurrenceError::LimitExceeded(
                    "source rule derivative order is too large".into(),
                ));
            }
            if request.row_origin.is_some_and(|o| o != first) {
                return Err(RecurrenceError::InvalidInput(format!(
                    "row_origin differs from the source first index {first}"
                )));
            }
            let count = explicit.as_ref().map_or(
                request
                    .row_count
                    .unwrap_or(DEFAULT_REPRESENTATION_ROW_COUNT),
                Vec::len,
            );
            if count > limits.max_rows {
                return Err(RecurrenceError::LimitExceeded(format!(
                    "{count} rows exceed max_rows = {}",
                    limits.max_rows
                )));
            }
            let generated = rule
                .generate_rows_rational(&seeds, first, count.max(seeds.len()))
                .map_err(|e| {
                    RecurrenceError::InvalidInput(format!("source rule generation failed: {e}"))
                })?
                .into_iter()
                .map(|r| vec_poly::trim(&r))
                .collect::<Vec<_>>();
            if let Some(rows) = &explicit {
                if rows[..] != generated[..rows.len()] {
                    return Err(RecurrenceError::InvalidInput(
                        "rows disagree with the rows generated by the source rule".into(),
                    ));
                }
            }
            let json = RecurrenceJson::from_recurrence_rational(&rule, first, &seeds, None);
            (kind, Some((rule, seeds, json)), first, generated)
        }
        None => {
            let rows = explicit.ok_or(RecurrenceError::InvalidInput(
                "give source_recurrence, oeis_id or rows".into(),
            ))?;
            ("rows_only", None, request.row_origin.unwrap_or(0), rows)
        }
    };
    if rows.is_empty() {
        return Err(RecurrenceError::InvalidInput(
            "at least one row is required".into(),
        ));
    }
    if rows.len() > limits.max_rows {
        return Err(RecurrenceError::LimitExceeded(format!(
            "{} rows exceed max_rows = {}",
            rows.len(),
            limits.max_rows
        )));
    }
    if let Some(degree) = rows.iter().filter_map(|r| vec_poly::degree(r)).max() {
        if degree > limits.max_degree {
            return Err(RecurrenceError::LimitExceeded(format!(
                "row degree {degree} exceeds max_degree"
            )));
        }
    }
    Ok(Source {
        kind,
        rule,
        origin,
        rows,
    })
}
