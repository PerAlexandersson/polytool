//! JSON fragments and result printing shared by the CLI report adapters.

use num_bigint::BigInt;

/// Print a complete JSON result, or fail before writing data on an error.
pub(crate) fn print_json_result<T: serde::Serialize>(result: Result<T, String>, compact: bool) {
    let text = result.and_then(|value| {
        if compact {
            serde_json::to_string(&value)
        } else {
            serde_json::to_string_pretty(&value)
        }
        .map_err(|error| error.to_string())
    });
    match text {
        Ok(text) => println!("{text}"),
        Err(error) => {
            eprintln!("{error}");
            std::process::exit(2);
        }
    }
}

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
