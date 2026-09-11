//! Export exact literature families as polytool-compatible coefficient rows.
use polytool::sequences::literature::{
    delannoy_polynomials_bigint, delannoy_square_polynomials_bigint,
    eulerian_square_polynomials_bigint, hoggatt_polynomials_bigint,
};

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() < 2 || args.iter().any(|s| s == "--help") {
        eprintln!("Usage: literature_sequence_rows delannoy|delannoy-square|eulerian-square MAX_N");
        eprintln!("       literature_sequence_rows hoggatt MAX_N M Q");
        std::process::exit(if args.iter().any(|s| s == "--help") {
            0
        } else {
            2
        });
    }
    let n: usize = args[1]
        .parse()
        .expect("MAX_N must be a nonnegative integer");
    let rows = match args[0].as_str() {
        "delannoy" if args.len() == 2 => delannoy_polynomials_bigint(n),
        "delannoy-square" if args.len() == 2 => delannoy_square_polynomials_bigint(n),
        "eulerian-square" if args.len() == 2 => eulerian_square_polynomials_bigint(n),
        "hoggatt" if args.len() == 4 => hoggatt_polynomials_bigint(
            n,
            args[2].parse().expect("M must be positive"),
            args[3].parse().expect("Q must be positive"),
        ),
        _ => {
            eprintln!("Unknown family or wrong argument count; use --help");
            std::process::exit(2);
        }
    };
    for row in rows {
        println!(
            "{}",
            row.iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(",")
        );
    }
}
