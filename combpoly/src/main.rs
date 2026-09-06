use clap::{Args, Parser, Subcommand, ValueEnum};
use num_bigint::BigInt;
use std::io::{self, Write};

use combpoly::{
    catalan, lattice_path_matroid, order, parking, permutation, polynomial_builder, statistics,
    word,
};
use polytool::real_rootedness::{format_poly, is_log_concave, is_real_rooted};
use polytool::recurrence;
use statistics::Stat;

#[derive(Debug, Clone, Copy, ValueEnum)]
enum OrderType {
    Bruhat,
    Weak,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum OutputFormat {
    Text,
    Csv,
    Json,
}

#[derive(Parser)]
#[command(name = "combpoly")]
#[command(about = "Explore combinatorial polynomials from permutations and words")]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Compute a generating polynomial over a set of objects
    Poly {
        #[command(flatten)]
        source: Source,

        /// Statistic to use
        #[arg(long, value_enum, default_value_t = Stat::Des)]
        stat: Stat,

        /// Check real-rootedness
        #[arg(long)]
        real_rooted: bool,

        /// Check log-concavity
        #[arg(long)]
        log_concave: bool,
    },

    /// For each permutation, compute polynomial over its order ideal
    Scan {
        /// Size of permutations (generates S_n)
        #[arg(long)]
        size: u8,

        #[command(flatten)]
        filters: PermFilters,

        /// Order type for ideal computation
        #[arg(long, value_enum, default_value_t = OrderType::Bruhat)]
        ideal: OrderType,

        /// Statistic to use
        #[arg(long, value_enum, default_value_t = Stat::Des)]
        stat: Stat,

        /// Check real-rootedness of each polynomial
        #[arg(long)]
        real_rooted: bool,

        /// Stop on first non-real-rooted polynomial
        #[arg(long)]
        halt: bool,
    },

    /// List combinatorial objects
    List {
        #[command(flatten)]
        source: Source,
    },

    /// Find a polynomial recurrence P_n(t) for n = min..=max
    Recurrence {
        /// Minimum n (1-based)
        #[arg(long, default_value_t = 1)]
        min_n: u8,

        /// Maximum n
        #[arg(long)]
        max_n: u8,

        /// Statistic to use
        #[arg(long, value_enum, default_value_t = Stat::Des)]
        stat: Stat,

        #[command(flatten)]
        filters: PermFilters,

        /// Automatically search for the simplest recurrence
        #[arg(long)]
        auto: bool,

        /// Show search progress (with --auto)
        #[arg(long)]
        verbose: bool,

        // --- Manual parameters (ignored when --auto is set) ---
        /// Max degree of coefficients in t
        #[arg(long, default_value_t = 1)]
        var_deg: usize,

        /// Max degree of coefficients in n
        #[arg(long, default_value_t = 1)]
        idx_deg: usize,

        /// Max derivative order
        #[arg(long, default_value_t = 0)]
        diff_deg: usize,

        /// Recurrence length (how many previous terms)
        #[arg(long, default_value_t = 2)]
        rec_len: usize,

        /// Allow inhomogeneous recurrence
        #[arg(long)]
        inhomogeneous: bool,

        /// Degree in t of LHS denominator
        #[arg(long, default_value_t = 0)]
        denom_var_deg: usize,

        /// Degree in n of LHS denominator
        #[arg(long, default_value_t = 0)]
        denom_idx_deg: usize,

        // --- Auto-mode upper bounds ---
        /// (--auto) Max recurrence length to search
        #[arg(long, default_value_t = 5)]
        max_rec_len: usize,

        /// (--auto) Max var_deg to search
        #[arg(long, default_value_t = 3)]
        max_var_deg: usize,

        /// (--auto) Max idx_deg to search
        #[arg(long, default_value_t = 3)]
        max_idx_deg: usize,

        /// (--auto) Max diff_deg to search
        #[arg(long, default_value_t = 2)]
        max_diff_deg: usize,

        /// (--auto) Minimum equation surplus
        #[arg(long, default_value_t = 1)]
        min_margin: usize,
    },

    /// Enumerate Dyck paths and compute exact lattice-path-matroid h*-vectors
    LpmHstarTable {
        /// Dyck semilength / area-sequence length
        #[arg(long)]
        semilength: Option<usize>,

        /// Enumerate all semilengths up to this value.  If --semilength is
        /// also present, enumerate semilength..=max-semilength.
        #[arg(long)]
        max_semilength: Option<usize>,

        /// Output format
        #[arg(long, value_enum, default_value_t = OutputFormat::Text)]
        format: OutputFormat,

        /// Include the full basis list in text/json output
        #[arg(long)]
        bases: bool,

        /// Check real-rootedness of each h*-polynomial when coefficients fit in i64
        #[arg(long)]
        real_rooted: bool,
    },

    /// Print the cyclic-interval hyperplanes for one LPM area sequence
    LpmHyperplanes {
        /// Dyck area sequence, e.g. 0,1,1,2
        #[arg(long)]
        area: String,

        /// Output format
        #[arg(long, value_enum, default_value_t = OutputFormat::Text)]
        format: OutputFormat,
    },

    /// Verify the Dyck snake-contact volume formula for LPMs
    LpmSnakeVolume {
        /// Single Dyck area sequence to check, e.g. 0,1,1,2
        #[arg(long)]
        area: Option<String>,

        /// Dyck semilength / area-sequence length
        #[arg(long)]
        semilength: Option<usize>,

        /// Enumerate all semilengths up to this value.  If --semilength is
        /// also present, enumerate semilength..=max-semilength.
        #[arg(long)]
        max_semilength: Option<usize>,

        /// Output format
        #[arg(long, value_enum, default_value_t = OutputFormat::Text)]
        format: OutputFormat,

        /// Print every checked row in text output
        #[arg(long)]
        all: bool,
    },
}

#[derive(Args)]
struct Source {
    /// Generate all permutations of S_n
    #[arg(long)]
    perms: Option<u8>,

    /// Generate all parking functions of size n
    #[arg(long)]
    parking: Option<u8>,

    /// Keep only words with no equal adjacent letters
    #[arg(long)]
    tieless: bool,

    /// Keep only words with no strict peaks w_{i-1} < w_i > w_{i+1}
    #[arg(long)]
    strict_peakless: bool,

    /// Keep only words with no weak-left peaks w_{i-1} <= w_i > w_{i+1}
    #[arg(long)]
    weak_peakless: bool,

    /// Bruhat lower ideal of a permutation (e.g., 321 or 3,2,1)
    #[arg(long, value_name = "PERM")]
    bruhat_ideal: Option<String>,

    /// Weak lower ideal of a permutation (e.g., 321 or 3,2,1)
    #[arg(long, value_name = "PERM")]
    weak_ideal: Option<String>,

    /// Generate words on a board (use with --content and --board)
    #[arg(long)]
    words: bool,

    /// Content vector for words (comma-separated, e.g., 2,2,1)
    #[arg(long, value_name = "ALPHA")]
    content: Option<String>,

    /// Board/partition for words (e.g., 22233 or 2,2,2,3,3)
    #[arg(long, value_name = "LAMBDA")]
    board: Option<String>,

    /// Skew lower board (e.g., 11 or 1,1)
    #[arg(long, value_name = "MU")]
    skew: Option<String>,

    #[command(flatten)]
    filters: PermFilters,
}

#[derive(Args)]
struct PermFilters {
    /// Pattern(s) to avoid (e.g., 312). Repeatable.
    #[arg(long, value_name = "PATTERN")]
    avoiding: Vec<String>,

    /// Arrow pattern(s) to avoid (e.g., "12;1->3"). Repeatable.
    #[arg(long, value_name = "ARROW_PATTERN")]
    arrow_avoiding: Vec<String>,

    /// Only alternating permutations (up-down: p1 < p2 > p3 < ...)
    #[arg(long)]
    alternating: bool,

    /// Only derangements (no fixed points)
    #[arg(long)]
    derangement: bool,

    /// Only permutations starting with this value
    #[arg(long, value_name = "VAL")]
    starts_with: Option<u8>,

    /// Only permutations ending with this value
    #[arg(long, value_name = "VAL")]
    ends_with: Option<u8>,
}

fn build_constraints(filters: &PermFilters) -> permutation::PermConstraints {
    permutation::PermConstraints {
        avoiding: filters
            .avoiding
            .iter()
            .map(|s| permutation::parse_sequence(s))
            .collect(),
        avoiding_arrow: filters
            .arrow_avoiding
            .iter()
            .map(|s| permutation::parse_arrow_pattern(s))
            .collect(),
        derangement: filters.derangement,
        alternating: filters.alternating,
        involution: false, // TODO: add --involution CLI flag
        starts_with: filters.starts_with,
        ends_with: filters.ends_with,
    }
}

fn get_objects(source: &Source) -> Vec<Vec<u8>> {
    if let Some(n) = source.perms {
        let constraints = build_constraints(&source.filters);
        let has_constraints = !constraints.avoiding.is_empty()
            || !constraints.avoiding_arrow.is_empty()
            || constraints.derangement
            || constraints.alternating
            || constraints.starts_with.is_some()
            || constraints.ends_with.is_some();
        let perms = if has_constraints {
            permutation::filtered_permutations(n, &constraints)
        } else {
            permutation::all_permutations(n)
        };
        perms
    } else if let Some(n) = source.parking {
        let mut parking_functions = parking::all_parking_functions(n);
        parking_functions.retain(|word| {
            (!source.tieless || word.windows(2).all(|pair| pair[0] != pair[1]))
                && (!source.strict_peakless || statistics::peaks(word) == 0)
                && (!source.weak_peakless || statistics::weak_left_peaks(word) == 0)
        });
        parking_functions
    } else if let Some(ref s) = source.bruhat_ideal {
        let perm = permutation::parse_sequence(s);
        order::bruhat_lower_ideal(&perm)
    } else if let Some(ref s) = source.weak_ideal {
        let perm = permutation::parse_sequence(s);
        order::weak_lower_ideal(&perm)
    } else if source.words || source.content.is_some() {
        let content_str = source
            .content
            .as_ref()
            .expect("--content required for words");
        let content = word::parse_content(content_str);
        let n: usize = content.iter().sum();
        let k = content.len() as u8;
        let board = if let Some(ref b) = source.board {
            word::parse_board(b)
        } else {
            vec![k; n]
        };
        let skew = source.skew.as_ref().map(|s| word::parse_board(s));
        word::words_on_board(&content, &board, skew.as_deref())
    } else {
        eprintln!(
            "No source specified. Use --perms, --parking, --bruhat-ideal, --weak-ideal, or --content."
        );
        std::process::exit(1);
    }
}

fn format_obj(p: &[u8]) -> String {
    if p.iter().all(|&x| x < 10) {
        p.iter().map(|x| x.to_string()).collect::<String>()
    } else {
        let strs: Vec<String> = p.iter().map(|x| x.to_string()).collect();
        format!("[{}]", strs.join(","))
    }
}

fn main() {
    let cli = Cli::parse();

    match cli.command {
        Commands::Poly {
            source,
            stat,
            real_rooted,
            log_concave,
        } => {
            let objects = get_objects(&source);
            let coeffs = polynomial_builder::build_generating_polynomial(&objects, stat);

            println!("Objects: {}", objects.len());
            println!("Stat: {}", stat);
            println!("Polynomial: {}", format_poly(&coeffs));
            println!("Coefficients: {:?}", coeffs);

            if real_rooted {
                let rr = is_real_rooted(&coeffs);
                println!("Real-rooted: {}", if rr { "yes" } else { "NO" });
            }
            if log_concave {
                let lc = is_log_concave(&coeffs);
                println!("Log-concave: {}", if lc { "yes" } else { "NO" });
            }
        }

        Commands::Scan {
            size,
            filters,
            ideal,
            stat,
            real_rooted,
            halt,
        } => {
            let constraints = build_constraints(&filters);
            let perms = permutation::filtered_permutations(size, &constraints);

            let total = perms.len();
            eprintln!(
                "Scanning {} permutations (ideal: {:?}, stat: {})",
                total, ideal, stat
            );

            let mut count = 0;
            let mut counterexamples = 0;

            for pi in &perms {
                let ideal_set = match ideal {
                    OrderType::Bruhat => order::bruhat_lower_ideal(pi),
                    OrderType::Weak => order::weak_lower_ideal(pi),
                };
                let coeffs = polynomial_builder::build_generating_polynomial(&ideal_set, stat);

                let mut line = format!("{}\t{:?}", format_obj(pi), coeffs);

                if real_rooted {
                    let rr = is_real_rooted(&coeffs);
                    line.push_str(if rr { "\tRR" } else { "\tNOT-RR" });
                    if !rr {
                        counterexamples += 1;
                        if halt {
                            println!("{}", line);
                            eprintln!("COUNTEREXAMPLE found! Stopping.");
                            std::process::exit(0);
                        }
                    }
                }

                println!("{}", line);
                count += 1;

                if total >= 500 && count % 500 == 0 {
                    eprint!("\r[{}/{}]", count, total);
                    io::stderr().flush().ok();
                }
            }

            if total >= 500 {
                eprintln!();
            }

            if real_rooted {
                eprintln!("---");
                if counterexamples == 0 {
                    eprintln!("All {} polynomials are real-rooted.", total);
                } else {
                    eprintln!(
                        "{} of {} polynomials are NOT real-rooted.",
                        counterexamples, total
                    );
                }
            }
        }

        Commands::List { source } => {
            let objects = get_objects(&source);
            eprintln!("Count: {}", objects.len());
            for obj in &objects {
                println!("{}", format_obj(obj));
            }
        }

        Commands::Recurrence {
            min_n,
            max_n,
            stat,
            filters,
            auto,
            verbose,
            var_deg,
            idx_deg,
            diff_deg,
            rec_len,
            inhomogeneous,
            denom_var_deg,
            denom_idx_deg,
            max_rec_len,
            max_var_deg,
            max_idx_deg,
            max_diff_deg,
            min_margin,
        } => {
            // Compute P_n(t) for each n in the range.
            let mut polys: Vec<Vec<i64>> = Vec::new();
            for n in min_n..=max_n {
                let constraints = build_constraints(&filters);
                let perms = permutation::filtered_permutations(n, &constraints);
                let coeffs = polynomial_builder::build_generating_polynomial(&perms, stat);
                eprintln!(
                    "P_{}(t) = {}  (coeffs: {:?})",
                    n,
                    format_poly(&coeffs),
                    coeffs
                );
                polys.push(coeffs);
            }

            eprintln!("---");

            if auto {
                let search = recurrence::AdaptiveSearchOptions {
                    max_rec_len,
                    max_var_deg,
                    max_idx_deg,
                    max_diff_deg,
                    try_inhomogeneous: inhomogeneous,
                    try_denominator: denom_var_deg > 0 || denom_idx_deg > 0,
                    min_margin,
                    verbose,
                    ..Default::default()
                };

                match recurrence::find_recurrence_adaptive(&polys, &search) {
                    Some(result) => {
                        eprintln!(
                            "Found after {} candidates (unknowns={}, equations={}, \
                             rec_len={}, var_deg={}, idx_deg={}, diff_deg={})",
                            result.candidates_tried,
                            result.num_unknowns,
                            result.num_equations,
                            result.opts.rec_len,
                            result.opts.var_deg,
                            result.opts.idx_deg,
                            result.opts.diff_deg,
                        );
                        println!("{}", result.recurrence);
                    }
                    None => println!("No recurrence found within search bounds."),
                }
            } else {
                let opts = recurrence::RecurrenceOptions {
                    var_deg,
                    idx_deg,
                    diff_deg,
                    rec_len,
                    homogeneous: !inhomogeneous,
                    denom_var_deg,
                    denom_idx_deg,
                    ..Default::default()
                };

                match recurrence::find_polynomial_recurrence(&polys, &opts) {
                    Some(rec) => println!("{}", rec),
                    None => println!("No recurrence found with the given parameters."),
                }
            }
        }

        Commands::LpmHstarTable {
            semilength,
            max_semilength,
            format,
            bases,
            real_rooted,
        } => {
            let Some((start, end)) = lpm_semilength_range(semilength, max_semilength) else {
                eprintln!("error: use --semilength N, --max-semilength N, or both");
                std::process::exit(1);
            };
            if start > end {
                eprintln!("error: --semilength must be <= --max-semilength");
                std::process::exit(1);
            }

            let mut rows = Vec::new();
            for n in start..=end {
                match lattice_path_matroid::lpm_hstar_table(n) {
                    Ok(mut table) => rows.append(&mut table),
                    Err(err) => {
                        eprintln!("error at semilength {n}: {err}");
                        std::process::exit(1);
                    }
                }
            }
            print_lpm_hstar_table(&rows, format, bases, real_rooted);
        }

        Commands::LpmHyperplanes { area, format } => {
            let area = parse_area_sequence_arg(&area).unwrap_or_else(|err| {
                eprintln!("error: {err}");
                std::process::exit(1);
            });
            if !catalan::is_area_sequence(&area) {
                eprintln!("error: not a valid Dyck area sequence: {area:?}");
                std::process::exit(1);
            }
            let intervals = catalan::peak_cliques(&area);
            let inequalities =
                lattice_path_matroid::lpm_cyclic_interval_inequalities_from_area_sequence(&area)
                    .expect("area sequence was already validated");
            print_lpm_hyperplanes(&area, &intervals, &inequalities, format);
        }

        Commands::LpmSnakeVolume {
            area,
            semilength,
            max_semilength,
            format,
            all,
        } => {
            let rows = if let Some(area) = area {
                let area = parse_area_sequence_arg(&area).unwrap_or_else(|err| {
                    eprintln!("error: {err}");
                    std::process::exit(1);
                });
                if !catalan::is_area_sequence(&area) {
                    eprintln!("error: not a valid Dyck area sequence: {area:?}");
                    std::process::exit(1);
                }
                let hstar = lattice_path_matroid::lpm_hstar_from_area_sequence(&area)
                    .unwrap_or_else(|err| {
                        eprintln!("error: {err}");
                        std::process::exit(1);
                    });
                let hstar_volume = hstar.iter().fold(BigInt::from(0), |acc, coeff| acc + coeff);
                let snake_volume =
                    lattice_path_matroid::lpm_snake_contact_volume_from_area_sequence(&area)
                        .unwrap_or_else(|err| {
                            eprintln!("error: {err}");
                            std::process::exit(1);
                        });
                vec![lattice_path_matroid::LpmSnakeContactVolumeRow {
                    area_sequence: area.clone(),
                    dyck_word: catalan::area_sequence_to_dyck_word(&area),
                    intervals: catalan::peak_cliques(&area),
                    rank: catalan::peak_cliques(&area).len(),
                    hstar,
                    hstar_volume,
                    snake_volume,
                }]
            } else {
                let Some((start, end)) = lpm_semilength_range(semilength, max_semilength) else {
                    eprintln!("error: use --area, --semilength N, or --max-semilength N");
                    std::process::exit(1);
                };
                if start > end {
                    eprintln!("error: --semilength must be <= --max-semilength");
                    std::process::exit(1);
                }

                let mut rows = Vec::new();
                for n in start..=end {
                    match lattice_path_matroid::lpm_snake_contact_volume_table(n) {
                        Ok(mut table) => rows.append(&mut table),
                        Err(err) => {
                            eprintln!("error at semilength {n}: {err}");
                            std::process::exit(1);
                        }
                    }
                }
                rows
            };
            print_lpm_snake_volume_check(&rows, format, all);
        }
    }
}

fn parse_area_sequence_arg(input: &str) -> Result<Vec<u8>, String> {
    let trimmed = input.trim();
    let trimmed = trimmed
        .strip_prefix('[')
        .unwrap_or(trimmed)
        .strip_suffix(']')
        .unwrap_or(trimmed);
    if trimmed.trim().is_empty() {
        return Ok(Vec::new());
    }

    let has_separator = trimmed
        .chars()
        .any(|ch| ch == ',' || ch == ';' || ch.is_ascii_whitespace());
    if !has_separator && trimmed.chars().all(|ch| ch.is_ascii_digit()) {
        return trimmed
            .chars()
            .map(|ch| {
                ch.to_digit(10)
                    .and_then(|digit| u8::try_from(digit).ok())
                    .ok_or_else(|| format!("invalid area entry `{ch}`"))
            })
            .collect();
    }

    trimmed
        .split(|ch: char| ch == ',' || ch == ';' || ch.is_ascii_whitespace())
        .filter(|part| !part.is_empty())
        .map(|part| {
            part.parse::<u8>()
                .map_err(|_| format!("invalid area entry `{part}`"))
        })
        .collect()
}

fn lpm_semilength_range(
    semilength: Option<usize>,
    max_semilength: Option<usize>,
) -> Option<(usize, usize)> {
    match (semilength, max_semilength) {
        (Some(n), None) => Some((n, n)),
        (None, Some(max_n)) => Some((0, max_n)),
        (Some(n), Some(max_n)) => Some((n, max_n)),
        (None, None) => None,
    }
}

fn print_lpm_hyperplanes(
    area: &[u8],
    intervals: &[(usize, usize)],
    inequalities: &[lattice_path_matroid::CyclicIntervalInequality],
    format: OutputFormat,
) {
    match format {
        OutputFormat::Text => {
            println!(
                "n={} area={} intervals={} rank={} cyclic_inequalities={}",
                area.len(),
                format_u8_list(area),
                format_intervals(intervals),
                intervals.len(),
                inequalities.len(),
            );
            for inequality in inequalities {
                println!(
                    "  {} <= {}",
                    format_cyclic_interval_sum(inequality, area.len()),
                    inequality.rank
                );
            }
        }
        OutputFormat::Csv => {
            println!("area,start,end,mask,rank");
            for inequality in inequalities {
                println!(
                    "\"{}\",{},{},{},{}",
                    format_u8_list(area),
                    inequality.start,
                    inequality.end,
                    inequality.mask,
                    inequality.rank,
                );
            }
        }
        OutputFormat::Json => {
            println!("[");
            for (i, inequality) in inequalities.iter().enumerate() {
                let comma = if i + 1 == inequalities.len() { "" } else { "," };
                println!(
                    "  {{\"area\":{},\"start\":{},\"end\":{},\"mask\":{},\"rank\":{}}}{}",
                    format_u8_json(area),
                    inequality.start,
                    inequality.end,
                    inequality.mask,
                    inequality.rank,
                    comma,
                );
            }
            println!("]");
        }
    }
}

fn print_lpm_snake_volume_check(
    rows: &[lattice_path_matroid::LpmSnakeContactVolumeRow],
    format: OutputFormat,
    show_all: bool,
) {
    match format {
        OutputFormat::Text => {
            let mismatches: Vec<_> = rows
                .iter()
                .filter(|row| row.snake_volume != row.hstar_volume)
                .collect();
            println!(
                "checked_rows={} mismatches={}",
                rows.len(),
                mismatches.len()
            );
            for row in rows
                .iter()
                .filter(|row| show_all || row.snake_volume != row.hstar_volume)
            {
                println!(
                    "n={} area={} dyck={} intervals={} rank={} h*={} h*(1)={} \
                     snake_volume={} match={}",
                    row.area_sequence.len(),
                    format_u8_list(&row.area_sequence),
                    row.dyck_word,
                    format_intervals(&row.intervals),
                    row.rank,
                    format_bigint_list(&row.hstar),
                    row.hstar_volume,
                    row.snake_volume,
                    row.snake_volume == row.hstar_volume,
                );
            }
        }
        OutputFormat::Csv => {
            println!("semilength,area,dyck,intervals,rank,hstar,hstar_volume,snake_volume,match");
            for row in rows {
                println!(
                    "{},\"{}\",\"{}\",\"{}\",{},\"{}\",{},{},{}",
                    row.area_sequence.len(),
                    format_u8_list(&row.area_sequence),
                    row.dyck_word,
                    format_intervals(&row.intervals),
                    row.rank,
                    format_bigint_list(&row.hstar),
                    row.hstar_volume,
                    row.snake_volume,
                    row.snake_volume == row.hstar_volume,
                );
            }
        }
        OutputFormat::Json => {
            println!("[");
            for (i, row) in rows.iter().enumerate() {
                let comma = if i + 1 == rows.len() { "" } else { "," };
                println!(
                    "  {{\"semilength\":{},\"area\":{},\"dyck\":\"{}\",\
                     \"intervals\":{},\"rank\":{},\"hstar\":{},\
                     \"hstar_volume\":{},\"snake_volume\":{},\"match\":{}}}{}",
                    row.area_sequence.len(),
                    format_u8_json(&row.area_sequence),
                    row.dyck_word,
                    format_intervals_json(&row.intervals),
                    row.rank,
                    format_bigint_json(&row.hstar),
                    row.hstar_volume,
                    row.snake_volume,
                    row.snake_volume == row.hstar_volume,
                    comma,
                );
            }
            println!("]");
        }
    }
}

fn format_cyclic_interval_sum(
    inequality: &lattice_path_matroid::CyclicIntervalInequality,
    ground_size: usize,
) -> String {
    let terms = (0..ground_size)
        .filter(|&i| (inequality.mask & (1usize << i)) != 0)
        .map(|i| format!("x_{}", i + 1))
        .collect::<Vec<_>>();
    terms.join(" + ")
}

fn print_lpm_hstar_table(
    rows: &[lattice_path_matroid::LpmHstarRow],
    format: OutputFormat,
    include_bases: bool,
    real_rooted: bool,
) {
    match format {
        OutputFormat::Text => {
            eprintln!("Rows: {}", rows.len());
            for row in rows {
                let rr = real_rooted.then(|| hstar_real_rooted_label(&row.hstar));
                println!(
                    "n={} area={} dyck={} intervals={} rank={} bases={} dim={} h*={}{}",
                    row.area_sequence.len(),
                    format_u8_list(&row.area_sequence),
                    row.dyck_word,
                    format_intervals(&row.intervals),
                    row.rank,
                    row.num_bases,
                    row.dimension,
                    format_bigint_list(&row.hstar),
                    rr.map(|label| format!(" rr={label}")).unwrap_or_default(),
                );
                if include_bases {
                    let matroid = lattice_path_matroid::LatticePathMatroid::from_peak_intervals(
                        row.area_sequence.len(),
                        row.intervals.clone(),
                    )
                    .expect("table row should have valid intervals");
                    println!("  bases={}", format_basis_list(&matroid.bases()));
                }
            }
        }
        OutputFormat::Csv => {
            println!("semilength,area,dyck,intervals,rank,num_bases,dimension,hstar,real_rooted");
            for row in rows {
                let rr = if real_rooted {
                    hstar_real_rooted_label(&row.hstar).to_string()
                } else {
                    String::new()
                };
                println!(
                    "{},\"{}\",\"{}\",\"{}\",{},{},{},\"{}\",\"{}\"",
                    row.area_sequence.len(),
                    format_u8_list(&row.area_sequence),
                    row.dyck_word,
                    format_intervals(&row.intervals),
                    row.rank,
                    row.num_bases,
                    row.dimension,
                    format_bigint_list(&row.hstar),
                    rr,
                );
            }
        }
        OutputFormat::Json => {
            println!("[");
            for (i, row) in rows.iter().enumerate() {
                let comma = if i + 1 == rows.len() { "" } else { "," };
                let rr_field = if real_rooted {
                    format!(
                        ",\"real_rooted\":\"{}\"",
                        hstar_real_rooted_label(&row.hstar)
                    )
                } else {
                    String::new()
                };
                let bases_field = if include_bases {
                    let matroid = lattice_path_matroid::LatticePathMatroid::from_peak_intervals(
                        row.area_sequence.len(),
                        row.intervals.clone(),
                    )
                    .expect("table row should have valid intervals");
                    format!(",\"bases\":{}", format_basis_json(&matroid.bases()))
                } else {
                    String::new()
                };
                println!(
                    "  {{\"semilength\":{},\"area\":{},\"dyck\":\"{}\",\"intervals\":{},\"rank\":{},\
                     \"num_bases\":{},\"dimension\":{},\"hstar\":{}{}{}}}{}",
                    row.area_sequence.len(),
                    format_u8_json(&row.area_sequence),
                    row.dyck_word,
                    format_intervals_json(&row.intervals),
                    row.rank,
                    row.num_bases,
                    row.dimension,
                    format_bigint_json(&row.hstar),
                    rr_field,
                    bases_field,
                    comma,
                );
            }
            println!("]");
        }
    }
}

fn hstar_real_rooted_label(hstar: &[BigInt]) -> &'static str {
    match lattice_path_matroid::hstar_to_i64(hstar) {
        Some(coeffs) => {
            if is_real_rooted(&coeffs) {
                "yes"
            } else {
                "NO"
            }
        }
        None => "too-large",
    }
}

fn format_u8_list(values: &[u8]) -> String {
    format!(
        "[{}]",
        values
            .iter()
            .map(u8::to_string)
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn format_u8_json(values: &[u8]) -> String {
    format_u8_list(values)
}

fn format_bigint_list(values: &[BigInt]) -> String {
    format!(
        "[{}]",
        values
            .iter()
            .map(BigInt::to_string)
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn format_bigint_json(values: &[BigInt]) -> String {
    format_bigint_list(values)
}

fn format_intervals(intervals: &[(usize, usize)]) -> String {
    format!(
        "[{}]",
        intervals
            .iter()
            .map(|(start, end)| format!("{start}-{end}"))
            .collect::<Vec<_>>()
            .join(";")
    )
}

fn format_intervals_json(intervals: &[(usize, usize)]) -> String {
    format!(
        "[{}]",
        intervals
            .iter()
            .map(|(start, end)| format!("[{start},{end}]"))
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn format_basis_list(bases: &[Vec<usize>]) -> String {
    format!(
        "[{}]",
        bases
            .iter()
            .map(|basis| {
                format!(
                    "{{{}}}",
                    basis
                        .iter()
                        .map(usize::to_string)
                        .collect::<Vec<_>>()
                        .join(",")
                )
            })
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn format_basis_json(bases: &[Vec<usize>]) -> String {
    format!(
        "[{}]",
        bases
            .iter()
            .map(|basis| {
                format!(
                    "[{}]",
                    basis
                        .iter()
                        .map(usize::to_string)
                        .collect::<Vec<_>>()
                        .join(",")
                )
            })
            .collect::<Vec<_>>()
            .join(",")
    )
}
