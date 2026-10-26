//! `polytool tnn-network`: total non-negativity of the coefficient triangle of
//! a monic polynomial sequence (or of a lower-unitriangular matrix), its
//! canonical planar network, and optional chain polynomials.

use crate::help::print_coefficient_input_help;
use crate::input::read_polys_bigint;
use crate::json_out::{json_bigint_vec, json_string, json_usize_vec};
use num_bigint::BigInt;
use polytool::*;

/// Fixed guardrail printed with every report.
pub(crate) const TNN_NETWORK_NOTE: &str = "TNN of the coefficient matrix does not imply \
     that the row polynomials are real-rooted: [[1,0,0],[1,1,0],[1,1,1]] is TNN and has \
     the row 1 + t + t^2. Brändén--Saud Maia Leite (arXiv:2412.06595, Theorem 3.7) \
     concerns the chain polynomials, not the rows. A finite triangle certifies only \
     this finite prefix.";

const CHAIN_THEOREM: &str = "Brändén--Saud Maia Leite, arXiv:2412.06595, Theorem 3.7";

pub(crate) fn print_tnn_network_help() {
    println!("Usage:");
    println!("  polytool tnn-network [--matrix] [--chain] [--json]");
    println!("  polytool tnn-network --help");
    println!();
    println!("Decide exactly whether the coefficient triangle R = (r_{{n,k}}) of the input");
    println!("rows is a lower-unitriangular totally nonnegative (TNN) matrix.  On success,");
    println!("print the weights lambda[n][k] of its canonical planar network Gamma_N(lambda)");
    println!("(lambda[n][k] is the weight of the vertical edge (n+1,k) -> (n,k); horizontal");
    println!("edges have weight 1) and re-evaluate the network's path matrix against R");
    println!("(reconstruct -> evaluate -> compare).  On failure, say why: a shape problem");
    println!("(not lower-unitriangular), a negative entry, or an explicit negative minor");
    println!("(rows, columns, exact value; zero-based indices).");
    println!();
    print_coefficient_input_help();
    println!();
    println!("  Default (polynomial rows): row n is a monic polynomial P_n of degree n,");
    println!("  r_{{n,k}} = [t^k] P_n.  The rows may start at P_0 = 1 or at P_1, in which");
    println!("  case P_0 = 1 is prepended.");
    println!("  With --matrix: each line is a matrix row, either a square lower-triangular");
    println!("  matrix (every row as long as the number of rows) or triangular rows (row n");
    println!("  has n + 1 entries).  Shape errors then name the offending entry.");
    println!();
    println!("Options:");
    println!("  --matrix  Read matrix rows instead of monic polynomial rows");
    println!("  --chain   Also compute the chain polynomials p_0 = 1,");
    println!("            p_n = t * sum_{{k<n}} r_{{n,k}} p_k, and check that each is real-rooted");
    println!("            with zeros in [-1, 0] and that p_n interlaces p_{{n+1}}");
    println!("  --json    Emit machine-readable JSON (schema polytool.tnn-network.v1)");
    println!();
    println!("Guardrails:");
    println!("  TNN of R does NOT make the row polynomials real-rooted:");
    println!("  [[1,0,0],[1,1,0],[1,1,1]] is TNN, yet its last row is 1 + t + t^2.");
    println!("  Brändén--Saud Maia Leite (arXiv:2412.06595, Thm 3.7): if R is TNN then the");
    println!("  chain polynomials are real-rooted in [-1, 0] and interlace.  The --chain");
    println!("  report keeps the TNN verdict separate from these polynomial checks.");
    println!("  A finite triangle certifies only the supplied prefix.");
    println!();
    println!("Exit status:");
    println!("  0  R is TNN and the network certificate re-evaluates to R");
    println!("  1  R is not lower-unitriangular TNN (the report says why), or a check that");
    println!("     the theory guarantees failed");
    println!("  2  invalid input or options");
    println!();
    println!("Examples:");
    println!("  printf '1,1\\n1,2,1\\n1,3,3,1\\n' | polytool tnn-network --chain");
    println!("  printf '1\\n1,1\\n1,1,1\\n' | polytool tnn-network");
    println!("  printf '1,0,0\\n1,1,0\\n2,0,1\\n' | polytool tnn-network --matrix --json");
}

#[derive(Default)]
struct TnnNetworkOptions {
    json: bool,
    matrix: bool,
    chain: bool,
}

fn parse_tnn_network_args(args: &[String]) -> Result<TnnNetworkOptions, String> {
    let mut options = TnnNetworkOptions::default();
    for arg in args {
        match arg.as_str() {
            "--json" => options.json = true,
            "--matrix" => options.matrix = true,
            "--chain" => options.chain = true,
            other => return Err(format!("unknown tnn-network option: {other}")),
        }
    }
    Ok(options)
}

/// Exact analysis shared by the text and JSON printers.
struct TnnNetworkAnalysis {
    input_mode: &'static str,
    input_rows: usize,
    prepended_p0: bool,
    coefficient_matrix: Option<Vec<Vec<BigInt>>>,
    verdict: Result<(), NetworkError>,
    lambda: Option<Vec<Vec<BigRational>>>,
    certificate_verified: Option<bool>,
    chain: Option<Result<ChainPolynomialReport, NetworkError>>,
}

impl TnnNetworkAnalysis {
    fn tnn(&self) -> bool {
        self.verdict.is_ok()
    }

    fn succeeded(&self) -> bool {
        let chain_ok = match &self.chain {
            Some(Ok(report)) => !report.contradicts_theorem(),
            _ => true,
        };
        self.tnn() && self.certificate_verified == Some(true) && chain_ok
    }
}

/// Triangular rows from the input lines in matrix mode: a square matrix goes
/// through `lower_triangular_rows_from_square`; other shapes are taken as
/// triangular rows and validated later.
fn matrix_rows(lines: Vec<Vec<BigInt>>) -> Result<Vec<Vec<BigInt>>, NetworkError> {
    let size = lines.len();
    if lines.iter().all(|row| row.len() == size) {
        lower_triangular_rows_from_square(&lines)
    } else {
        Ok(lines)
    }
}

fn analyze(lines: Vec<Vec<BigInt>>, options: &TnnNetworkOptions) -> TnnNetworkAnalysis {
    let input_rows = lines.len();
    let (input_mode, rows) = if options.matrix {
        ("matrix", matrix_rows(lines))
    } else {
        let polys: Vec<Polynomial<BigInt>> = lines.into_iter().map(Polynomial::new).collect();
        (
            "polynomials",
            coefficient_matrix_from_monic_polynomials(&polys),
        )
    };
    let prepended_p0 = !options.matrix && rows.as_ref().is_ok_and(|r| r.len() == input_rows + 1);
    let verdict = rows
        .as_ref()
        .map_err(Clone::clone)
        .and_then(|r| check_lower_unitriangular_tnn(r));
    let (lambda, certificate_verified) = match (&rows, &verdict) {
        (Ok(r), Ok(())) => match reconstruct_canonical_tnn_network(r) {
            Ok(network) => {
                let verified = verify_path_matrix_certificate(r, &network);
                (Some(network.lambda), Some(verified))
            }
            Err(_) => (None, Some(false)),
        },
        _ => (None, None),
    };
    let chain = options.chain.then(|| {
        rows.as_ref()
            .map_err(Clone::clone)
            .and_then(|r| check_chain_polynomials_bigint(r))
    });
    TnnNetworkAnalysis {
        input_mode,
        input_rows,
        prepended_p0,
        coefficient_matrix: rows.ok(),
        verdict,
        lambda,
        certificate_verified,
        chain,
    }
}

/// Stable snake-case name of a [`NetworkError`] variant.
pub(crate) fn network_error_kind(error: &NetworkError) -> &'static str {
    match error {
        NetworkError::EmptyMatrix => "empty_matrix",
        NetworkError::EmptyPolynomialSequence => "empty_polynomial_sequence",
        NetworkError::NotSquare { .. } => "not_square",
        NetworkError::NonzeroAboveDiagonal { .. } => "nonzero_above_diagonal",
        NetworkError::UnsupportedStartingDegree { .. } => "unsupported_starting_degree",
        NetworkError::UnexpectedPolynomialDegree { .. } => "unexpected_polynomial_degree",
        NetworkError::NonMonicPolynomial { .. } => "non_monic_polynomial",
        NetworkError::InvalidRowLength { .. } => "invalid_row_length",
        NetworkError::NegativeEntry { .. } => "negative_entry",
        NetworkError::NonUnitDiagonal { .. } => "non_unit_diagonal",
        NetworkError::NotTnn { .. } => "not_tnn",
        NetworkError::ZeroTailViolation { .. } => "zero_tail_violation",
        NetworkError::NegativeReducedEntry { .. } => "negative_reduced_entry",
        NetworkError::ReconstructionMismatch => "reconstruction_mismatch",
    }
}

/// The offending position of a [`NetworkError`]: `(row, col, value)`.
fn network_error_location(error: &NetworkError) -> (Option<usize>, Option<usize>, Option<String>) {
    match error {
        NetworkError::NotSquare { row, .. } | NetworkError::InvalidRowLength { row, .. } => {
            (Some(*row), None, None)
        }
        NetworkError::UnexpectedPolynomialDegree { index, .. } => (Some(*index), None, None),
        NetworkError::NonMonicPolynomial { index, leading } => (
            Some(*index),
            None,
            Some(
                leading
                    .as_ref()
                    .map_or("0".to_string(), ToString::to_string),
            ),
        ),
        NetworkError::NonzeroAboveDiagonal { row, col, value }
        | NetworkError::NegativeEntry { row, col, value } => {
            (Some(*row), Some(*col), Some(value.to_string()))
        }
        NetworkError::NonUnitDiagonal { row, value } => {
            (Some(*row), Some(*row), Some(value.to_string()))
        }
        NetworkError::ZeroTailViolation { row, value } => {
            (Some(*row), Some(0), Some(value.to_string()))
        }
        NetworkError::NegativeReducedEntry { row, col, value } => {
            (Some(*row), Some(*col), Some(value.to_string()))
        }
        _ => (None, None, None),
    }
}

fn json_option_usize(value: Option<usize>) -> String {
    value.map_or("null".to_string(), |v| v.to_string())
}

fn json_option_string(value: Option<&str>) -> String {
    value.map_or("null".to_string(), json_string)
}

fn json_rational_rows(rows: &[Vec<BigRational>]) -> String {
    format!(
        "[{}]",
        rows.iter()
            .map(|row| {
                format!(
                    "[{}]",
                    row.iter()
                        .map(|q| json_string(&q.to_string()))
                        .collect::<Vec<_>>()
                        .join(",")
                )
            })
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn json_bigint_rows(rows: &[Vec<BigInt>]) -> String {
    format!(
        "[{}]",
        rows.iter()
            .map(|row| json_bigint_vec(row))
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn network_error_json(error: &NetworkError) -> String {
    let (row, col, value) = network_error_location(error);
    let witness = match error {
        NetworkError::NotTnn {
            witness: Some(minor),
            ..
        } => format!(
            "{{\"rows\":{},\"cols\":{},\"value\":{}}}",
            json_usize_vec(&minor.rows),
            json_usize_vec(&minor.cols),
            json_string(&minor.value.to_string())
        ),
        _ => "null".to_string(),
    };
    format!(
        "{{\"kind\":{},\"message\":{},\"row\":{},\"col\":{},\"value\":{},\"witness\":{}}}",
        json_string(network_error_kind(error)),
        json_string(&error.to_string()),
        json_option_usize(row),
        json_option_usize(col),
        json_option_string(value.as_deref()),
        witness
    )
}

fn chain_json(chain: &Result<ChainPolynomialReport, NetworkError>) -> String {
    let report = match chain {
        Ok(report) => report,
        Err(error) => {
            return format!(
                "{{\"theorem\":{},\"error\":{},\"polynomials\":[],\"interlacing\":[],\
                 \"conclusion_holds\":null,\"contradicts_theorem\":null}}",
                json_string(CHAIN_THEOREM),
                json_string(&error.to_string())
            )
        }
    };
    let polynomials = report
        .chain_polynomials
        .iter()
        .enumerate()
        .map(|(n, p)| {
            format!(
                "{{\"n\":{},\"polynomial\":{},\"coefficients\":{},\"real_rooted\":{},\
                 \"roots_in_minus_one_zero\":{}}}",
                n,
                json_string(&format_poly_bigint_coeffs(p)),
                json_bigint_vec(p),
                report.real_rooted[n],
                report.roots_in_minus_one_zero[n]
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    let interlacing = report
        .interlacing
        .iter()
        .enumerate()
        .map(|(n, ok)| {
            format!(
                "{{\"left\":{},\"right\":{},\"interlaces\":{}}}",
                n,
                n + 1,
                ok
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    format!(
        "{{\"theorem\":{},\"error\":null,\"polynomials\":[{}],\"interlacing\":[{}],\
         \"conclusion_holds\":{},\"contradicts_theorem\":{}}}",
        json_string(CHAIN_THEOREM),
        polynomials,
        interlacing,
        report.chain_conclusion_holds(),
        report.contradicts_theorem()
    )
}

fn analysis_json(analysis: &TnnNetworkAnalysis) -> String {
    format!(
        "{{\"schema\":\"polytool.tnn-network.v1\",\"input_mode\":{},\"input_rows\":{},\
         \"prepended_p0\":{},\"coefficient_matrix\":{},\"tnn\":{},\"failure\":{},\
         \"lambda\":{},\"certificate_verified\":{},\"chain\":{},\"note\":{}}}",
        json_string(analysis.input_mode),
        analysis.input_rows,
        analysis.prepended_p0,
        analysis
            .coefficient_matrix
            .as_deref()
            .map_or("null".to_string(), json_bigint_rows),
        analysis.tnn(),
        analysis
            .verdict
            .as_ref()
            .err()
            .map_or("null".to_string(), network_error_json),
        analysis
            .lambda
            .as_deref()
            .map_or("null".to_string(), json_rational_rows),
        analysis
            .certificate_verified
            .map_or("null".to_string(), |v| v.to_string()),
        analysis
            .chain
            .as_ref()
            .map_or("null".to_string(), chain_json),
        json_string(TNN_NETWORK_NOTE)
    )
}

fn join_display<T: ToString>(values: &[T]) -> String {
    values
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>()
        .join(", ")
}

fn print_chain_text(chain: &Result<ChainPolynomialReport, NetworkError>) {
    println!();
    println!("Chain polynomials ({CHAIN_THEOREM}):");
    let report = match chain {
        Ok(report) => report,
        Err(error) => {
            println!("- Not available: {error}");
            return;
        }
    };
    for (n, p) in report.chain_polynomials.iter().enumerate() {
        let status = if report.roots_in_minus_one_zero[n] {
            "real-rooted, zeros in [-1, 0]"
        } else if report.real_rooted[n] {
            "real-rooted, some zero outside [-1, 0]"
        } else {
            "not real-rooted"
        };
        println!("  p_{n} = {}: {status}", format_poly_bigint_coeffs(p));
    }
    let failed: Vec<String> = report
        .interlacing
        .iter()
        .enumerate()
        .filter(|(_, ok)| !**ok)
        .map(|(n, _)| format!("p_{n} / p_{}", n + 1))
        .collect();
    if failed.is_empty() {
        println!("- Consecutive interlacing p_n << p_(n+1): all hold");
    } else {
        println!(
            "- Consecutive interlacing p_n << p_(n+1): fails for {}",
            failed.join(", ")
        );
    }
    println!(
        "- Theorem conclusion (real-rooted in [-1, 0] and interlacing): {}",
        if report.chain_conclusion_holds() {
            "holds"
        } else {
            "fails"
        }
    );
    if report.contradicts_theorem() {
        println!("- WARNING: the matrix is TNN but the conclusion fails (contradicts the theorem)");
    } else if !report.is_tnn() {
        println!("- The matrix is not TNN, so the theorem does not apply to these polynomials.");
    }
}

fn print_analysis_text(analysis: &TnnNetworkAnalysis) {
    println!("TNN network check");
    println!(
        "- Input: {} {} row(s){}",
        analysis.input_rows,
        if analysis.input_mode == "matrix" {
            "matrix"
        } else {
            "polynomial"
        },
        if analysis.prepended_p0 {
            "; P_0 = 1 prepended"
        } else {
            ""
        }
    );
    if let Some(rows) = &analysis.coefficient_matrix {
        println!("- Coefficient matrix (triangular rows):");
        for (n, row) in rows.iter().enumerate() {
            println!("  {n}: {}", join_display(row));
        }
    }
    match &analysis.verdict {
        Ok(()) => println!("- Verdict: lower-unitriangular and totally nonnegative (TNN)"),
        Err(error) => {
            println!("- Verdict: not lower-unitriangular TNN");
            println!("- Reason ({}): {error}", network_error_kind(error));
            if let NetworkError::NotTnn {
                witness: Some(minor),
                ..
            } = error
            {
                println!(
                    "- Witness: rows {:?}, columns {:?} (zero-based) give the minor {}",
                    minor.rows, minor.cols, minor.value
                );
            }
        }
    }
    if let Some(lambda) = &analysis.lambda {
        println!("- Canonical network weights lambda[n][k] (edge (n+1,k) -> (n,k)):");
        if lambda.is_empty() {
            println!("  (none: a single row)");
        }
        for (n, row) in lambda.iter().enumerate() {
            println!("  n={n}: {}", join_display(row));
        }
    }
    if let Some(verified) = analysis.certificate_verified {
        println!(
            "- Certificate (reconstruct -> evaluate -> compare): {}",
            if verified {
                "the network's path matrix equals the input"
            } else {
                "FAILED"
            }
        );
    }
    if let Some(chain) = &analysis.chain {
        print_chain_text(chain);
    }
    println!();
    println!("Note: {TNN_NETWORK_NOTE}");
}

pub(crate) fn cmd_tnn_network(args: &[String]) {
    let options = match parse_tnn_network_args(args) {
        Ok(options) => options,
        Err(error) => {
            eprintln!("{error}");
            eprintln!("Run `polytool help tnn-network` for usage.");
            std::process::exit(2);
        }
    };
    let lines = read_polys_bigint();
    if lines.is_empty() {
        eprintln!("tnn-network: no input rows");
        std::process::exit(2);
    }
    let analysis = analyze(lines, &options);
    if options.json {
        println!("{}", analysis_json(&analysis));
    } else {
        print_analysis_text(&analysis);
    }
    if !analysis.succeeded() {
        std::process::exit(1);
    }
}
