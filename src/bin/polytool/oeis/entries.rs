//! List local OEIS records with the shared keyword reader.

use polytool::oeis::list_local_oeis_entries;
use std::path::Path;

pub(super) fn run(args: &[String]) {
    let mut data = None;
    let mut keywords = Vec::new();
    let mut json = false;
    let mut iter = args.iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--data" => data = iter.next(),
            "--keyword" => {
                let Some(value) = iter.next() else {
                    crate::input::usage_error(
                        "oeis entries",
                        "--keyword requires a comma-separated list",
                    );
                };
                for keyword in value.split(',').map(str::trim) {
                    if keyword.is_empty() || keyword.starts_with('-') {
                        crate::input::usage_error("oeis entries", "keywords must be nonempty");
                    }
                    keywords.push(keyword.to_string());
                }
            }
            "--json" => json = true,
            _ => crate::input::usage_error("oeis entries", "invalid argument"),
        }
    }
    let Some(data) = data else {
        crate::input::usage_error("oeis entries", "--data DIR is required");
    };
    let entries = list_local_oeis_entries(Path::new(data), &keywords).unwrap_or_else(|error| {
        eprintln!("cannot list local OEIS entries: {error}");
        std::process::exit(2);
    });
    if json {
        println!("{}", serde_json::to_string_pretty(&entries).unwrap());
    } else {
        for id in entries {
            println!("{id}");
        }
    }
}
