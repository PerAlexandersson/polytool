use polytool_mcp::{PolynomialToolsServer, ReportStoreConfig};
use rmcp::{transport::stdio, ServiceExt};

const BIN_NAME: &str = "polytool-mcp";
const CONTACT: &str = "Per Alexandersson (@PerAlexandersson, <per.w.alexandersson@gmail.com>)";

#[derive(Debug, Eq, PartialEq)]
enum StartupMode {
    Serve(ReportStoreConfig),
    Help,
    Version,
}

fn parse_args(args: &[String]) -> Result<StartupMode, String> {
    if args.len() == 2 {
        match args[1].as_str() {
            "--help" | "-h" => return Ok(StartupMode::Help),
            "--version" | "-V" => return Ok(StartupMode::Version),
            _ => {}
        }
    }
    let mut config = ReportStoreConfig::from_env()?;
    let mut no_store = false;
    let mut arguments = args.iter().skip(1);
    while let Some(flag) = arguments.next() {
        if flag == "--no-report-store" {
            no_store = true;
            continue;
        }
        if ![
            "--report-dir",
            "--report-max-bytes",
            "--report-max-age-seconds",
        ]
        .contains(&flag.as_str())
        {
            return Err(format!("unknown option: {flag}"));
        }
        let value = arguments
            .next()
            .ok_or_else(|| format!("missing value for {flag}"))?;
        match flag.as_str() {
            "--report-dir" if value.is_empty() => {
                return Err("report directory must not be empty".into())
            }
            "--report-dir" => config.directory = Some(value.into()),
            "--report-max-bytes" => {
                config.max_bytes = value
                    .parse()
                    .map_err(|_| format!("invalid value for {flag}"))?
            }
            "--report-max-age-seconds" => {
                config.max_age = std::time::Duration::from_secs(
                    value
                        .parse()
                        .map_err(|_| format!("invalid value for {flag}"))?,
                )
            }
            _ => unreachable!(),
        }
    }
    if no_store {
        config.directory = None;
    }
    Ok(StartupMode::Serve(config))
}

fn version_text() -> String {
    format!("{} {}", BIN_NAME, env!("CARGO_PKG_VERSION"))
}

fn help_text() -> String {
    format!(
        "\
{name} {version}
{description}

Usage:
  {name} [--report-dir DIR] [--no-report-store]
  {name} --help
  {name} --version

Options:
  -h, --help       Print this help text and exit
  -V, --version    Print version information and exit
  --report-dir DIR              Store full reports in DIR (overrides environment)
  --no-report-store             Disable storage, even with a configured directory
  --report-max-bytes N          Store capacity in bytes (default 536870912)
  --report-max-age-seconds N    Report lifetime in seconds (default 1209600)

Environment:
  POLYTOOL_MCP_REPORT_DIR, POLYTOOL_MCP_REPORT_MAX_BYTES,
  POLYTOOL_MCP_REPORT_MAX_AGE_SECONDS. No directory is chosen implicitly.

MCP usage:
  Run without arguments from an MCP client. The server uses stdio transport and
  exposes tools only; it does not provide resources, prompts, sampling, HTTP, or
  arbitrary filesystem access. With a report directory configured it stores
  full reports and accepts only content-addressed report ids for reads.

find_recurrence arguments:
  Provide exactly one of polynomials, coefficients, expressions, or text.
  Search controls are top-level: skip_prefix; min_recurrence_length;
  max_recurrence_length; min_variable_degree; max_variable_degree;
  min_index_degree; max_index_degree; min_derivative_order;
  max_derivative_order; try_inhomogeneous; min_inhomogeneous_variable_degree;
  max_inhomogeneous_variable_degree; min_inhomogeneous_index_degree;
  max_inhomogeneous_index_degree; try_denominator; try_alternating_sign;
  try_laurent; max_laurent_power; max_denominator_variable_degree;
  max_denominator_index_degree; min_equation_margin; no_held_out_rows;
  fit_extra_rows; modular_prefilter; and max_candidates. The names used up to 0.6.0
  (min_rec_len, max_var_deg, max_idx_deg, ...) are accepted as aliases.
  Legacy nested options remain accepted; top-level values take precedence.
  Set include_code=true to include Mathematica, Python, Sage, and recurrence JSON.
  Omitting include_code returns the compact response.

Contact:
  Report issues in the Git repository, or contact {contact}.
",
        name = BIN_NAME,
        version = env!("CARGO_PKG_VERSION"),
        description = env!("CARGO_PKG_DESCRIPTION"),
        contact = CONTACT,
    )
}

#[tokio::main]
async fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let config = match parse_args(&args) {
        Ok(StartupMode::Help) => {
            print!("{}", help_text());
            return Ok(());
        }
        Ok(StartupMode::Version) => {
            println!("{}", version_text());
            return Ok(());
        }
        Ok(StartupMode::Serve(config)) => config,
        Err(message) => {
            eprintln!("error: {message}");
            eprintln!("Try '{BIN_NAME} --help' for usage.");
            std::process::exit(2);
        }
    };

    let service = PolynomialToolsServer::with_report_store(config)
        .serve(stdio())
        .await?;
    service.waiting().await?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(values: &[&str]) -> Vec<String> {
        values.iter().map(|value| value.to_string()).collect()
    }

    #[test]
    fn parses_cli_flags() {
        assert!(matches!(
            parse_args(&args(&["polytool-mcp"])),
            Ok(StartupMode::Serve(_))
        ));
        assert_eq!(
            parse_args(&args(&["polytool-mcp", "--help"])),
            Ok(StartupMode::Help)
        );
        assert_eq!(
            parse_args(&args(&["polytool-mcp", "-h"])),
            Ok(StartupMode::Help)
        );
        assert_eq!(
            parse_args(&args(&["polytool-mcp", "--version"])),
            Ok(StartupMode::Version)
        );
        assert!(parse_args(&args(&["polytool-mcp", "--bad"])).is_err());
        assert!(parse_args(&args(&["polytool-mcp", "--help", "extra"])).is_err());
    }

    #[test]
    fn help_mentions_transport_and_contact() {
        let help = help_text();
        assert!(help.contains("stdio transport"));
        assert!(help.contains("Search controls are top-level"));
        assert!(help.contains("max_candidates"));
        assert!(help.contains("include_code=true"));
        assert!(help.contains(CONTACT));
    }

    #[test]
    fn report_flags_override_defaults_and_no_store_wins() {
        let Ok(StartupMode::Serve(config)) = parse_args(&args(&[
            "polytool-mcp",
            "--report-dir",
            "reports",
            "--report-max-bytes",
            "1024",
            "--report-max-age-seconds",
            "60",
        ])) else {
            panic!("valid report flags");
        };
        assert_eq!(config.directory, Some("reports".into()));
        assert_eq!(config.max_bytes, 1024);
        assert_eq!(config.max_age, std::time::Duration::from_secs(60));
        let Ok(StartupMode::Serve(config)) = parse_args(&args(&[
            "polytool-mcp",
            "--no-report-store",
            "--report-dir",
            "reports",
        ])) else {
            panic!("valid disable flag");
        };
        assert!(config.directory.is_none());
        for flags in [
            vec!["--report-dir"],
            vec!["--report-dir", ""],
            vec!["--report-max-bytes", "bad"],
            vec!["--report-max-age-seconds", "-1"],
        ] {
            let mut values = vec!["polytool-mcp"];
            values.extend(flags);
            assert!(parse_args(&args(&values)).is_err());
        }
    }
}
