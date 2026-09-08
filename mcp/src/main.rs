use polytool_mcp::PolynomialToolsServer;
use rmcp::{transport::stdio, ServiceExt};

const BIN_NAME: &str = "polytool-mcp";
const CONTACT: &str = "Per Alexandersson (@PerAlexandersson, <per.w.alexandersson@gmail.com>)";

#[derive(Debug, Eq, PartialEq)]
enum StartupMode {
    Serve,
    Help,
    Version,
}

fn parse_args(args: &[String]) -> Result<StartupMode, String> {
    match args {
        [_program] => Ok(StartupMode::Serve),
        [_program, flag] if flag == "--help" || flag == "-h" => Ok(StartupMode::Help),
        [_program, flag] if flag == "--version" || flag == "-V" => Ok(StartupMode::Version),
        [_program, flag] => Err(format!("unknown option: {flag}")),
        [_program, flag, ..] => Err(format!("unexpected extra arguments after: {flag}")),
        [] => Ok(StartupMode::Serve),
    }
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
  {name}
  {name} --help
  {name} --version

Options:
  -h, --help       Print this help text and exit
  -V, --version    Print version information and exit

MCP usage:
  Run without arguments from an MCP client. The server uses stdio transport and
  exposes tools only; it does not provide resources, prompts, sampling, HTTP, or
  filesystem access.

find_recurrence arguments:
  Provide exactly one of polynomials, coefficients, expressions, or text.
  Search controls are top-level: skip_prefix; min_rec_len; max_rec_len;
  min_var_deg; max_var_deg; min_idx_deg; max_idx_deg; min_diff_deg;
  max_diff_deg; try_inhomogeneous; min_inhomo_var_deg; max_inhomo_var_deg;
  min_inhomo_idx_deg; max_inhomo_idx_deg; try_denominator;
  try_alternating_sign; max_denom_var_deg; max_denom_idx_deg; min_margin;
  no_verify; fit_extra_rows; modular_prefilter; and max_candidates.
  Legacy nested options remain accepted; top-level values take precedence.
  Set include_code=false to omit Mathematica, Python, Sage, and recurrence JSON.
  Omitting include_code preserves the full response.

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
    match parse_args(&args) {
        Ok(StartupMode::Help) => {
            print!("{}", help_text());
            return Ok(());
        }
        Ok(StartupMode::Version) => {
            println!("{}", version_text());
            return Ok(());
        }
        Ok(StartupMode::Serve) => {}
        Err(message) => {
            eprintln!("error: {message}");
            eprintln!("Try '{BIN_NAME} --help' for usage.");
            std::process::exit(2);
        }
    }

    let service = PolynomialToolsServer::new().serve(stdio()).await?;
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
        assert_eq!(parse_args(&args(&["polytool-mcp"])), Ok(StartupMode::Serve));
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
        assert!(help.contains("include_code=false"));
        assert!(help.contains(CONTACT));
    }
}
