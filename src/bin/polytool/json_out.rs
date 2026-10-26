//! Hand-written JSON fragments shared by the text/JSON report printers.

use num_bigint::BigInt;

fn json_escape(input: &str) -> String {
    let mut out = String::with_capacity(input.len());
    for ch in input.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c.is_control() => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out
}

pub(crate) fn json_string(input: &str) -> String {
    format!("\"{}\"", json_escape(input))
}

pub(crate) fn json_bigint_vec(values: &[BigInt]) -> String {
    format!(
        "[{}]",
        values
            .iter()
            .map(|value| json_string(&value.to_string()))
            .collect::<Vec<_>>()
            .join(",")
    )
}

pub(crate) fn json_usize_vec(values: &[usize]) -> String {
    format!(
        "[{}]",
        values
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>()
            .join(",")
    )
}

pub(crate) fn json_string_vec(values: &[String]) -> String {
    format!(
        "[{}]",
        values
            .iter()
            .map(|value| json_string(value))
            .collect::<Vec<_>>()
            .join(",")
    )
}

pub(crate) fn json_bool_option(value: Option<bool>) -> String {
    value
        .map(|v| if v { "true" } else { "false" }.to_string())
        .unwrap_or_else(|| "null".to_string())
}
