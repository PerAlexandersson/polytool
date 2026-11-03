//! `polytool atlas`: the real-rooted atlas engine (see `polytool::atlas`).
//!
//! Subcommands: `seed`, `grow`, `oeis-index`, `match-oeis`, `audit`,
//! `report`, and the internal `job` (one child-process job).

use crate::atlas_jobs::{cmd_atlas_job, run_child, run_pool, JobItem, JobRequest, JobResponse};
use polytool::atlas::grow::{plan_forward_round, ForwardBudget};
use polytool::atlas::oeis_match::{
    build_oeis_terms_cache, match_atlas, OeisMatchOptions, OeisTermsIndex,
};
use polytool::atlas::seed::{catalog_sources, named_family_sources, AnalyzedDerivation};
use polytool::atlas::{
    atlas_report, audit_atlas_with, integrate_derivations, AtlasAuditOptions, AtlasCertificate,
    AtlasDerivation, AtlasStore, CertificateAnalysis, IntegrationOptions, IntegrationSummary,
};
use polytool::recurrence::RecurrenceJson;
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::time::{Duration, Instant};

pub(crate) fn print_atlas_help() {
    println!("Usage:");
    println!("  polytool atlas seed DIR [options]");
    println!("  polytool atlas grow DIR [options]");
    println!("  polytool atlas oeis-index --oeis-data OEIS_DIR --cache FILE");
    println!("  polytool atlas match-oeis DIR --cache FILE [--oeis-data OEIS_DIR]");
    println!("  polytool atlas audit DIR [--json] [--sample K [--seed S]]");
    println!("  polytool atlas report DIR [--json]");
    println!();
    println!("The real-rooted atlas: a proof DAG of integer polynomial sequences, each");
    println!("proved real-rooted (with a root region) by polytool's automatic routes.");
    println!("DIR holds nodes.jsonl, edges.jsonl and links.jsonl (append-only JSONL; the");
    println!("last revision of a node id wins).  Every edge carries a replayable polytool");
    println!("certificate; `audit` replays all of them from the roots.");
    println!();
    println!("seed: try every polytool named family and every bundled OEIS catalog entry.");
    println!("  A node is added only when a theorem check (three-term, Liu--Wang,");
    println!("  derivative lag, Ma--Wang), the parent provers (periodic product, paper");
    println!("  theorem, root domain) or a lift from an existing node (preserver image,");
    println!("  gamma, basis map) proves it.  Unproved entries stay out of the atlas.");
    println!("  --jobs N            Parallel child processes (default 4)");
    println!("  --item-timeout S    Seconds per source (default 60); a source over the");
    println!("                      limit is reported as timeout");
    println!("  --only LIST         named, catalog, or comma-separated A-numbers / names");
    println!("  --limit N           First N sources only");
    println!("  --no-lifts          Root provers only");
    println!("  --lift-batch N      Sources per lift child (default 40)");
    println!("  --lift-passes N     Lift passes (default 1)");
    println!();
    println!("grow: one forward round over the existing nodes (preserver images, basis");
    println!("  maps, gamma maps, Veronese sections r <= 3, Obreschkoff combinations");
    println!("  P_n + P_(n+1)), with exact image rules, integer rows and de-duplication");
    println!("  by fingerprint (a hit adds an edge, a merged definition or a link).");
    println!("  --jobs N, --item-timeout S (default 60), --batch N (tasks per child,");
    println!("  default 40), --max-tasks N (default 3000), --max-basis-maps N (default 6),");
    println!("  --operations LIST   Comma-separated subset of preserver_image,basis_map,");
    println!("                      gamma_polynomial,veronese_section,obreschkoff_combination");
    println!("  --per-operation N   At most N tasks of each operation (after --operations)");
    println!("  --parents FILE      Only tasks on the node ids listed in FILE (one per line),");
    println!("                      e.g. the nodes before the round when it is run in parts");
    println!();
    println!("oeis-index: build the OEIS terms cache from a local OEIS copy (DIR/seq/A123/");
    println!("  A123456.seq); keep the cache outside the data directory.");
    println!("match-oeis: match every node's triangle, reversed triangle and row sums");
    println!("  against the cache (built first when --oeis-data is given and the cache is");
    println!("  missing); hits are links, never identifications.");
    println!("audit: replay every edge, check every property and every stored row; exit");
    println!("  status 1 on any problem.  --sample K replays only the incoming edges of K");
    println!("  pseudo-random nodes, of nodes with several incoming edges and of nodes in");
    println!("  candidate-equal links (other edges are taken as recorded); the row checks");
    println!("  of every node's claims always run.");
    println!("report: counts of nodes, edges per operation, links, OEIS hits.");
}

struct Options {
    positional: Vec<String>,
    values: BTreeMap<String, String>,
    flags: Vec<String>,
}

fn parse(args: &[String], with_value: &[&str]) -> Result<Options, String> {
    let mut options = Options {
        positional: Vec::new(),
        values: BTreeMap::new(),
        flags: Vec::new(),
    };
    let mut i = 0;
    while i < args.len() {
        let arg = &args[i];
        if let Some(name) = arg.strip_prefix("--") {
            if with_value.contains(&name) {
                let value = args
                    .get(i + 1)
                    .ok_or_else(|| format!("--{name} needs a value"))?;
                options.values.insert(name.to_string(), value.clone());
                i += 2;
                continue;
            }
            options.flags.push(name.to_string());
        } else {
            options.positional.push(arg.clone());
        }
        i += 1;
    }
    Ok(options)
}

impl Options {
    fn number(&self, name: &str, default: usize) -> Result<usize, String> {
        match self.values.get(name) {
            None => Ok(default),
            Some(value) => value
                .parse()
                .map_err(|_| format!("--{name} expects a nonnegative integer")),
        }
    }

    fn directory(&self) -> Result<PathBuf, String> {
        self.positional
            .first()
            .map(PathBuf::from)
            .ok_or_else(|| "missing data directory".to_string())
    }

    fn flag(&self, name: &str) -> bool {
        self.flags.iter().any(|f| f == name)
    }
}

fn fail(message: impl std::fmt::Display) -> ! {
    eprintln!("polytool atlas: {message}");
    std::process::exit(2);
}

pub(crate) fn cmd_atlas(args: &[String]) {
    let Some(sub) = args.first() else {
        print_atlas_help();
        std::process::exit(1);
    };
    let rest = &args[1..];
    let result = match sub.as_str() {
        "job" => {
            cmd_atlas_job();
            Ok(())
        }
        "seed" => cmd_seed(rest),
        "grow" => cmd_grow(rest),
        "oeis-index" => cmd_oeis_index(rest),
        "match-oeis" => cmd_match_oeis(rest),
        "audit" => cmd_audit(rest),
        "report" => cmd_report(rest),
        other => Err(format!("unknown atlas subcommand `{other}`")),
    };
    if let Err(message) = result {
        fail(message);
    }
}

/// Integrate job items in order and print a summary line.
fn integrate_items(
    store: &mut AtlasStore,
    derivations: Vec<AnalyzedDerivation>,
    equality_timeout: Duration,
) -> Result<IntegrationSummary, String> {
    let list: Vec<(AtlasDerivation, Option<CertificateAnalysis>)> =
        derivations.into_iter().map(|(d, a)| (d, Some(a))).collect();
    let mut equality =
        |existing: &RecurrenceJson, new: &RecurrenceJson| -> Option<AtlasCertificate> {
            let request = JobRequest::Equality {
                existing: existing.clone(),
                new: new.clone(),
            };
            run_child(&request, equality_timeout)
                .ok()
                .and_then(|r| r.items.into_iter().next())
                .and_then(|item| item.equality)
        };
    let summary =
        integrate_derivations(store, &list, &IntegrationOptions::default(), &mut equality)
            .map_err(|e| e.to_string())?;
    store.commit().map_err(|e| e.to_string())?;
    Ok(summary)
}

fn print_summary(label: &str, summary: &IntegrationSummary) {
    eprintln!(
        "{label}: {} derivations, {} edges added ({} present), {} nodes added, {} definitions \
         merged, {} links added",
        summary.derivations,
        summary.edges_added,
        summary.edges_present,
        summary.nodes_added,
        summary.definitions_merged,
        summary.links_added
    );
    for (reason, count) in &summary.rejected {
        eprintln!("  rejected {reason}: {count}");
        for example in summary.examples.get(reason).into_iter().flatten() {
            eprintln!("    e.g. {example}");
        }
    }
    for contradiction in &summary.contradictions {
        eprintln!("  CONTRADICTION: {contradiction}");
    }
}

fn status_counts(items: &[(String, JobItem)]) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for (_, item) in items {
        let key = if item.status.starts_with("error") {
            "error".to_string()
        } else {
            item.status.clone()
        };
        *counts.entry(key).or_default() += 1;
    }
    counts
}

fn flatten_responses(
    labels: Vec<Vec<String>>,
    responses: Vec<Result<JobResponse, String>>,
) -> Vec<(String, JobItem)> {
    let mut out = Vec::new();
    for (labels, response) in labels.into_iter().zip(responses) {
        match response {
            Ok(response) => {
                for (label, item) in labels.into_iter().zip(response.items) {
                    out.push((label, item));
                }
            }
            Err(error) => {
                for label in labels {
                    out.push((label, JobItem::failed(error.clone(), 0.0)));
                }
            }
        }
    }
    out
}

fn cmd_seed(args: &[String]) -> Result<(), String> {
    let options = parse(
        args,
        &[
            "jobs",
            "item-timeout",
            "only",
            "limit",
            "lift-batch",
            "lift-passes",
        ],
    )?;
    let directory = options.directory()?;
    let jobs = options.number("jobs", 4)?;
    let item_timeout = options.number("item-timeout", 60)? as u64;
    let limit = options.number("limit", usize::MAX)?;
    let lift_batch = options.number("lift-batch", 40)?.max(1);
    let lift_passes = options.number("lift-passes", 1)?;
    let started = Instant::now();
    let mut sources = named_family_sources();
    let (catalog, skipped) = catalog_sources();
    sources.extend(catalog);
    if let Some(only) = options.values.get("only") {
        let wanted: Vec<&str> = only.split(',').collect();
        sources.retain(|s| {
            wanted.iter().any(|w| match *w {
                "named" => s.source.starts_with("named_family:"),
                "catalog" => s.source.starts_with("oeis_catalog:"),
                other => s.oeis_id.as_deref() == Some(other) || s.name == other,
            })
        });
    }
    sources.truncate(limit);
    eprintln!(
        "seed: {} sources ({} catalog entries skipped: rule not runnable)",
        sources.len(),
        skipped.len()
    );
    let mut store = AtlasStore::open(&directory).map_err(|e| e.to_string())?;
    // Root provers, one child per source.
    let phase = Instant::now();
    let requests: Vec<JobRequest> = sources
        .iter()
        .map(|spec| JobRequest::Roots { spec: spec.clone() })
        .collect();
    let labels: Vec<Vec<String>> = sources.iter().map(|s| vec![s.source.clone()]).collect();
    let responses = run_pool(
        requests,
        jobs,
        |_| Duration::from_secs(item_timeout),
        "roots",
    );
    let items = flatten_responses(labels, responses);
    eprintln!(
        "roots: {:?} in {:.0} s",
        status_counts(&items),
        phase.elapsed().as_secs_f64()
    );
    let mut proved_sources = vec![false; sources.len()];
    let mut derivations = Vec::new();
    let mut owner = Vec::new();
    for (i, (_, item)) in items.iter().enumerate() {
        for derivation in &item.derivations {
            derivations.push(derivation.clone());
            owner.push(i);
        }
    }
    let summary = integrate_items(&mut store, derivations, Duration::from_secs(item_timeout))?;
    print_summary("roots integrated", &summary);
    for (outcome, &i) in summary.outcomes.iter().zip(&owner) {
        if outcome.is_some() {
            proved_sources[i] = true;
        }
    }
    let timeouts: Vec<&String> = items
        .iter()
        .filter(|(_, item)| item.status == "timeout")
        .map(|(label, _)| label)
        .collect();
    if !timeouts.is_empty() {
        eprintln!("roots timeouts: {timeouts:?}");
    }
    if !options.flag("no-lifts") {
        for pass in 0..lift_passes {
            let phase = Instant::now();
            let open: Vec<usize> = (0..sources.len()).filter(|&i| !proved_sources[i]).collect();
            if open.is_empty() {
                break;
            }
            let mut requests = Vec::new();
            let mut labels = Vec::new();
            let mut owners = Vec::new();
            for chunk in open.chunks(lift_batch) {
                requests.push(JobRequest::Lifts {
                    directory: directory.clone(),
                    specs: chunk.iter().map(|&i| sources[i].clone()).collect(),
                    item_timeout_seconds: item_timeout,
                });
                labels.push(chunk.iter().map(|&i| sources[i].source.clone()).collect());
                owners.extend(chunk.iter().copied());
            }
            let responses = run_pool(
                requests,
                jobs,
                |request| match request {
                    JobRequest::Lifts { specs, .. } => {
                        Duration::from_secs(item_timeout * specs.len() as u64 + 1800)
                    }
                    _ => Duration::from_secs(item_timeout),
                },
                "lifts",
            );
            let items = flatten_responses(labels, responses);
            eprintln!(
                "lifts pass {}: {} open sources, {:?} in {:.0} s",
                pass + 1,
                open.len(),
                status_counts(&items),
                phase.elapsed().as_secs_f64()
            );
            let mut derivations = Vec::new();
            let mut owner = Vec::new();
            for ((_, item), &i) in items.iter().zip(&owners) {
                for derivation in &item.derivations {
                    derivations.push(derivation.clone());
                    owner.push(i);
                }
            }
            let summary =
                integrate_items(&mut store, derivations, Duration::from_secs(item_timeout))?;
            print_summary("lifts integrated", &summary);
            let mut progress = false;
            for (outcome, &i) in summary.outcomes.iter().zip(&owner) {
                if let Some(node) = outcome {
                    let own = store.node(node).is_some_and(|n| {
                        n.definitions.iter().any(|d| d.source == sources[i].source)
                    });
                    if own && !proved_sources[i] {
                        proved_sources[i] = true;
                        progress = true;
                    }
                }
            }
            if !progress {
                break;
            }
        }
    }
    let proved = proved_sources.iter().filter(|p| **p).count();
    eprintln!(
        "seed: {proved} of {} sources proved; {} nodes, {} edges, {} links; {:.0} s",
        sources.len(),
        store.nodes().len(),
        store.edges().len(),
        store.links().len(),
        started.elapsed().as_secs_f64()
    );
    Ok(())
}

fn cmd_grow(args: &[String]) -> Result<(), String> {
    let options = parse(
        args,
        &[
            "jobs",
            "item-timeout",
            "batch",
            "max-tasks",
            "max-basis-maps",
            "operations",
            "per-operation",
            "parents",
        ],
    )?;
    let directory = options.directory()?;
    let jobs = options.number("jobs", 4)?;
    let item_timeout = options.number("item-timeout", 60)? as u64;
    let batch = options.number("batch", 40)?.max(1);
    let mut budget = ForwardBudget {
        max_tasks: usize::MAX,
        max_basis_maps: options
            .number("max-basis-maps", ForwardBudget::default().max_basis_maps)?,
        ..ForwardBudget::default()
    };
    if options.values.contains_key("parents") || options.values.contains_key("per-operation") {
        // The caps are then applied after the filters, per operation.
        budget.max_preserver_tasks = usize::MAX;
        budget.max_basis_map_tasks = usize::MAX;
        budget.max_gamma_tasks = usize::MAX;
        budget.max_veronese_tasks = usize::MAX;
        budget.max_obreschkoff_tasks = usize::MAX;
    }
    let started = Instant::now();
    let mut store = AtlasStore::open(&directory).map_err(|e| e.to_string())?;
    let mut tasks = plan_forward_round(&store, &budget).map_err(|e| e.to_string())?;
    if let Some(operations) = options.values.get("operations") {
        let wanted: Vec<&str> = operations.split(',').collect();
        tasks.retain(|t| wanted.contains(&t.operation()));
    }
    if let Some(file) = options.values.get("parents") {
        let text = std::fs::read_to_string(file).map_err(|e| format!("--parents: {e}"))?;
        let parents: std::collections::BTreeSet<&str> = text
            .lines()
            .map(str::trim)
            .filter(|l| !l.is_empty())
            .collect();
        tasks.retain(|t| parents.contains(t.node()));
    }
    if let Some(cap) = options.values.get("per-operation") {
        let cap: usize = cap
            .parse()
            .map_err(|_| "--per-operation expects a nonnegative integer")?;
        let mut seen: BTreeMap<&str, usize> = BTreeMap::new();
        tasks.retain(|t| {
            let count = seen.entry(t.operation()).or_default();
            *count += 1;
            *count <= cap
        });
    }
    tasks.truncate(options.number("max-tasks", ForwardBudget::default().max_tasks)?);
    let mut planned: BTreeMap<&str, usize> = BTreeMap::new();
    for task in &tasks {
        *planned.entry(task.operation()).or_default() += 1;
    }
    eprintln!("grow: {} tasks planned {planned:?}", tasks.len());
    let mut requests = Vec::new();
    let mut labels = Vec::new();
    let mut operations = Vec::new();
    for chunk in tasks.chunks(batch) {
        requests.push(JobRequest::Forward {
            tasks: chunk.to_vec(),
            item_timeout_seconds: item_timeout,
        });
        labels.push(chunk.iter().map(|t| t.label()).collect());
        operations.extend(chunk.iter().map(|t| t.operation()));
    }
    let responses = run_pool(
        requests,
        jobs,
        |request| match request {
            JobRequest::Forward { tasks, .. } => {
                Duration::from_secs(item_timeout * tasks.len() as u64 + 300)
            }
            _ => Duration::from_secs(item_timeout),
        },
        "grow",
    );
    let items = flatten_responses(labels, responses);
    let timed_out: Vec<&String> = items
        .iter()
        .filter(|(_, item)| item.status == "timeout")
        .map(|(label, _)| label)
        .collect();
    if !timed_out.is_empty() {
        eprintln!("grow: {} timed-out tasks (not rerun):", timed_out.len());
        for label in &timed_out {
            eprintln!("  timeout: {label}");
        }
    }
    let mut by_operation: BTreeMap<String, BTreeMap<String, usize>> = BTreeMap::new();
    let mut derivations = Vec::new();
    for ((_, item), operation) in items.iter().zip(&operations) {
        let status = if item.status.starts_with("error")
            || (item.derivations.is_empty() && item.status == "ok")
        {
            "no_result".to_string()
        } else {
            item.status.clone()
        };
        *by_operation
            .entry(operation.to_string())
            .or_default()
            .entry(status)
            .or_default() += 1;
        derivations.extend(item.derivations.iter().cloned());
    }
    eprintln!(
        "grow: task outcomes {by_operation:?} in {:.0} s",
        started.elapsed().as_secs_f64()
    );
    let summary = integrate_items(&mut store, derivations, Duration::from_secs(item_timeout))?;
    print_summary("grow integrated", &summary);
    eprintln!(
        "grow: {} nodes, {} edges, {} links; {:.0} s",
        store.nodes().len(),
        store.edges().len(),
        store.links().len(),
        started.elapsed().as_secs_f64()
    );
    Ok(())
}

fn cmd_oeis_index(args: &[String]) -> Result<(), String> {
    let options = parse(args, &["oeis-data", "cache"])?;
    let data = options
        .values
        .get("oeis-data")
        .ok_or("--oeis-data is required")?;
    let cache = options.values.get("cache").ok_or("--cache is required")?;
    let started = Instant::now();
    let count = build_oeis_terms_cache(&PathBuf::from(data), &PathBuf::from(cache))
        .map_err(|e| e.to_string())?;
    eprintln!(
        "oeis-index: {count} entries written to {cache} in {:.0} s",
        started.elapsed().as_secs_f64()
    );
    Ok(())
}

fn cmd_match_oeis(args: &[String]) -> Result<(), String> {
    let options = parse(args, &["oeis-data", "cache"])?;
    let directory = options.directory()?;
    let cache = PathBuf::from(options.values.get("cache").ok_or("--cache is required")?);
    let started = Instant::now();
    if !cache.exists() {
        let data = options
            .values
            .get("oeis-data")
            .ok_or("the cache is missing and --oeis-data is not given")?;
        let count =
            build_oeis_terms_cache(&PathBuf::from(data), &cache).map_err(|e| e.to_string())?;
        eprintln!("match-oeis: cache built ({count} entries)");
    }
    let match_options = OeisMatchOptions::default();
    let index = OeisTermsIndex::load(&cache, &match_options).map_err(|e| e.to_string())?;
    eprintln!(
        "match-oeis: {} entries indexed in {:.0} s",
        index.len(),
        started.elapsed().as_secs_f64()
    );
    let mut store = AtlasStore::open(&directory).map_err(|e| e.to_string())?;
    let summary = match_atlas(&mut store, &index, &match_options);
    store.commit().map_err(|e| e.to_string())?;
    eprintln!(
        "match-oeis: {} nodes, hits {:?}, {} links added, {} nodes matched outside the \
         catalog; {:.0} s",
        summary.nodes,
        summary.hits_by_target,
        summary.links_added,
        summary.new_triangle_matches.len(),
        started.elapsed().as_secs_f64()
    );
    Ok(())
}

fn cmd_audit(args: &[String]) -> Result<(), String> {
    let options = parse(args, &["sample", "seed"])?;
    let directory = options.directory()?;
    let started = Instant::now();
    let store = AtlasStore::open(&directory).map_err(|e| e.to_string())?;
    let audit_options = AtlasAuditOptions {
        sample_nodes: match options.values.get("sample") {
            None => None,
            Some(_) => Some(options.number("sample", 0)?),
        },
        sample_seed: options.number("seed", 1)? as u64,
    };
    let report = audit_atlas_with(&store, &audit_options);
    if options.flag("json") {
        let value = serde_json::json!({
            "schema": "polytool.atlas-audit.v1",
            "nodes": report.nodes,
            "edges": report.edges,
            "links": report.links,
            "edges_accepted": report.edges_accepted,
            "edges_trusted": report.edges_trusted,
            "nodes_replayed": report.nodes_replayed,
            "properties_checked": report.properties_checked,
            "row_claims_checked": report.row_claims_checked,
            "edges_by_operation": report.edges_by_operation,
            "problems": report.problems,
            "passed": report.passed(),
            "seconds": started.elapsed().as_secs_f64(),
        });
        println!(
            "{}",
            serde_json::to_string_pretty(&value).unwrap_or_default()
        );
    } else {
        println!(
            "audit: {} nodes ({} replayed), {} edges ({} replayed and accepted, {} taken as \
             recorded), {} links; {} properties and {} row claims checked; {} problems; {:.0} s",
            report.nodes,
            report.nodes_replayed,
            report.edges,
            report.edges_accepted,
            report.edges_trusted,
            report.links,
            report.properties_checked,
            report.row_claims_checked,
            report.problems.len(),
            started.elapsed().as_secs_f64()
        );
        for problem in &report.problems {
            println!("  {problem}");
        }
    }
    if !report.passed() {
        std::process::exit(1);
    }
    Ok(())
}

fn cmd_report(args: &[String]) -> Result<(), String> {
    let options = parse(args, &[])?;
    let directory = options.directory()?;
    let store = AtlasStore::open(&directory).map_err(|e| e.to_string())?;
    let report = atlas_report(&store);
    if options.flag("json") {
        println!(
            "{}",
            serde_json::to_string_pretty(&report).map_err(|e| e.to_string())?
        );
        return Ok(());
    }
    println!("nodes: {} {:?}", report.nodes, report.nodes_by_source);
    println!("catalog ids identified: {}", report.catalog_ids_identified);
    println!(
        "catalog tails identified (rows after prefix rows): {}",
        report.catalog_tail_ids_identified
    );
    println!(
        "lean-flagged catalog entries reached: {} of {} (and {} as tails only)",
        report.lean_flagged_reached, report.lean_flagged_total, report.lean_flagged_reached_as_tail
    );
    println!(
        "tails of entries flagged known_not_real_rooted: {:?}",
        report.tail_identified_known_not_real_rooted
    );
    println!(
        "identified, not lean-flagged: {} {:?}",
        report.identified_not_lean_flagged.len(),
        report.identified_not_lean_flagged
    );
    println!(
        "identified but flagged known_not_real_rooted: {:?}",
        report.identified_known_not_real_rooted
    );
    println!("edges: {} {:?}", report.edges, report.edges_by_operation);
    println!("links: {:?}", report.links_by_kind);
    println!("oeis links: {:?}", report.oeis_links);
    println!("nodes by tightest region: {:?}", report.nodes_by_region);
    println!(
        "nodes with interlacing: {}; defining basis maps: {}",
        report.nodes_with_interlacing, report.nodes_defining_basis_maps
    );
    println!(
        "nodes with several incoming proof edges: {} (several distinct routes: {})",
        report.nodes_with_several_proofs.len(),
        report.nodes_with_several_routes
    );
    for (node, name, operations) in &report.nodes_with_several_proofs {
        println!("  {node} {name}: {operations:?}");
    }
    println!(
        "nodes matched to OEIS entries outside the catalog: {}",
        report.nodes_matched_outside_catalog.len()
    );
    for (node, name, ids) in &report.nodes_matched_outside_catalog {
        println!("  {node} {name}: {ids:?}");
    }
    Ok(())
}
