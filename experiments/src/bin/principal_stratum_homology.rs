use clap::{Parser, ValueEnum};
use combinatoric_core::{cancel_units, integral_homology, UnitReductionOptions};
use experiments::principal_stratum::{modular_field_betti, rational_field_betti, PrincipalStratumModel};
use polytool::SparseModRowOrder;

#[derive(Clone, Copy, Debug, ValueEnum)]
enum PivotOrder { Input, Nonzeros, Markowitz, DynamicMarkowitz }
impl From<PivotOrder> for SparseModRowOrder {
    fn from(value: PivotOrder) -> Self { match value {
        PivotOrder::Input => Self::Input,
        PivotOrder::Nonzeros => Self::IncreasingNonzeros,
        PivotOrder::Markowitz => Self::IncreasingMarkowitzCost,
        PivotOrder::DynamicMarkowitz => Self::DynamicMarkowitz,
    }}
}

#[derive(Parser, Debug)]
#[command(about = "Bounded exact homology computation for principal-stratum composition complexes")]
struct Args {
    #[arg(long, default_value = "1,1", help = "Comma-separated positive parts of omega")]
    omega: String,
    #[arg(long, default_value_t = 4)] d: u32,
    #[arg(long, default_value_t = 251)] prime: u64,
    #[arg(long, value_enum, default_value_t = PivotOrder::Input)] pivot_order: PivotOrder,
    #[arg(long, default_value_t = 250_000)] max_cells: usize,
    #[arg(long, default_value_t = 5_000_000)] max_nnz: usize,
    #[arg(long, default_value_t = 50_000_000)] max_reduction_nnz: usize,
    #[arg(long, default_value_t = 1_000_000)] max_pivots: usize,
    #[arg(long, default_value_t = 0, help = "Run a bounded rational oracle only when nonzero")] rational_dense_budget: usize,
    #[arg(long, help = "Skip unit cancellation and integral assembly")] field_only: bool,
}

fn main() -> Result<(), String> {
    let args = Args::parse();
    let omega = args.omega.split(',').map(|part| part.trim().parse::<u32>().map_err(|error| format!("invalid omega part {part:?}: {error}"))).collect::<Result<Vec<_>, _>>()?;
    let (model, build) = PrincipalStratumModel::build_bounded(omega, args.d, args.max_cells).map_err(|error| error.to_string())?;
    if model.initial_nnz() > args.max_nnz { return Err(format!("boundary NNZ {} exceeds --max-nnz {}", model.initial_nnz(), args.max_nnz)); }
    println!("input omega={:?} d={} cells={} euler={} initial_nnz={}", model.omega(), model.d(), model.cell_count(), model.euler_characteristic(), model.initial_nnz());
    println!("timing enumeration_ms={} assembly_validation_ms={}", build.enumeration_time.as_millis(), build.assembly_and_validation_time.as_millis());
    println!("cells_by_degree={:?}", model.cells_by_degree().iter().map(|(degree, cells)| (*degree, cells.len())).collect::<Vec<_>>());
    let modular_started = std::time::Instant::now();
    let modular = modular_field_betti(&model, args.prime, args.pivot_order.into()).map_err(|error| error.to_string())?;
    println!("field {} betti={:?} ranks={:?} time_ms={}", modular.evidence_label, modular.betti_numbers, modular.ranks, modular_started.elapsed().as_millis());
    if args.rational_dense_budget > 0 { println!("rational_oracle={:?}", rational_field_betti(&model, args.rational_dense_budget).map_err(|error| error.to_string())?); }
    if args.field_only { return Ok(()); }
    let reduction_started = std::time::Instant::now();
    let reduction = cancel_units(model.complex(), UnitReductionOptions { max_pivots: args.max_pivots, max_nnz: args.max_reduction_nnz, record_certificate: false, ..UnitReductionOptions::default() }).map_err(|error| error.to_string())?;
    println!("unit_reduction_ms={} pivots={} initial_peak_final_nnz=({},{},{}) residual_counts={:?}", reduction_started.elapsed().as_millis(), reduction.stats.pivot_updates, reduction.stats.initial_nnz, reduction.stats.peak_nnz, reduction.stats.final_nnz, reduction.reduced.generator_counts());
    let smith_started = std::time::Instant::now();
    let groups = integral_homology(&reduction.reduced, Default::default()).map_err(|error| error.to_string())?;
    println!("integral_groups={groups:?} smith_ms={}", smith_started.elapsed().as_millis());
    Ok(())
}
