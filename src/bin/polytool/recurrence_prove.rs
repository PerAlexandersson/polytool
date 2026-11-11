//! Discover and replay general recurrence proofs through the library registry.

use crate::input::{error_text, parse_usize_option, read_stdin_or_exit};
use crate::recurrence_input::read_recurrence_json;
use polytool::recurrence::{
    discover_recurrence_proofs_with_stage_options, AffineCoordinateOptions,
    AlternativeRepresentationOptions, CatalogLiftOptions, CatalogLiftSearch,
    RealRootednessTheoremOptions, RealRootednessTheoremReport, RecurrenceProofDiscoveryReport,
    RecurrenceProofStageOptions,
};

pub(crate) fn cmd_recurrence_prove(args: &[String]) {
    match run(args) {
        Ok(proved) => {
            if !proved {
                std::process::exit(1);
            }
        }
        Err(error) => {
            eprintln!("{error}");
            std::process::exit(2);
        }
    }
}

fn run(args: &[String]) -> Result<bool, String> {
    let (mut path, mut candidate, mut delay, mut representations) = (None, None, None, None);
    let (mut lifts, mut affine) = (None, None);
    let (mut json, mut verify) = (false, false);
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--recurrence" => {
                i += 1;
                path = Some(
                    args.get(i)
                        .ok_or("--recurrence needs a file or -")?
                        .as_str(),
                );
            }
            "--candidate" => candidate = Some(parse_usize_option(args, &mut i, "--candidate")?),
            "--max-start-delay" => {
                delay = Some(parse_usize_option(args, &mut i, "--max-start-delay")?)
            }
            "--max-representation-candidates" => {
                representations = Some(parse_usize_option(
                    args,
                    &mut i,
                    "--max-representation-candidates",
                )?)
            }
            "--max-lift-candidates" => {
                lifts = Some(parse_usize_option(args, &mut i, "--max-lift-candidates")?)
            }
            "--max-affine-coordinates" => {
                affine = Some(parse_usize_option(
                    args,
                    &mut i,
                    "--max-affine-coordinates",
                )?)
            }
            "--json" => json = true,
            "--verify" => verify = true,
            other => return Err(format!("unknown recurrence-prove option: {other}")),
        }
        i += 1;
    }
    if verify {
        if path.is_some()
            || candidate.is_some()
            || delay.is_some()
            || representations.is_some()
            || lifts.is_some()
            || affine.is_some()
        {
            return Err("--verify cannot be combined with search flags".into());
        }
        let report: RecurrenceProofDiscoveryReport =
            serde_json::from_str(&read_stdin_or_exit()).map_err(error_text)?;
        if report.maximum_start_delay > 16 {
            return Err("maximum_start_delay exceeds 16".into());
        }
        report.verify().map_err(error_text)?;
        if json {
            println!(
                "{}",
                serde_json::json!({"verified": true, "status": report.status, "lean_replay": report.lean_replay})
            );
        } else {
            println!(
                "Report verified; discovery status: {}. Lean replay: {}.",
                report.status, report.lean_replay
            );
        }
        return Ok(true);
    }
    let maximum_start_delay = delay.unwrap_or(3);
    if maximum_start_delay > 16 {
        return Err("--max-start-delay must be at most 16".into());
    }
    // Without --recurrence the record is the command input (stdin or
    // --input FILE).
    let record = read_recurrence_json(path.unwrap_or("-"), candidate)?;
    let alternative_options = AlternativeRepresentationOptions {
        max_candidates: representations
            .unwrap_or(AlternativeRepresentationOptions::default().max_candidates),
        ..Default::default()
    };
    if alternative_options.max_candidates > 100_000 {
        return Err("--max-representation-candidates must be at most 100000".into());
    }
    let lift_options = CatalogLiftOptions {
        max_candidates: lifts.unwrap_or(CatalogLiftOptions::default().max_candidates),
        ..Default::default()
    };
    if lift_options.max_candidates > 1000 {
        return Err("--max-lift-candidates must be at most 1000".into());
    }
    let affine_options = AffineCoordinateOptions {
        max_coordinates: affine.unwrap_or(AffineCoordinateOptions::default().max_coordinates),
    };
    if affine_options.max_coordinates > 64 {
        return Err("--max-affine-coordinates must be at most 64".into());
    }
    let report = discover_recurrence_proofs_with_stage_options(
        &record,
        &RealRootednessTheoremOptions {
            maximum_start_delay,
            ..Default::default()
        },
        &RecurrenceProofStageOptions {
            alternative_representations: alternative_options,
            catalog_lifts: lift_options,
            affine_coordinates: affine_options,
        },
    )
    .map_err(error_text)?;
    if json {
        println!(
            "{}",
            serde_json::to_string_pretty(&report).map_err(error_text)?
        );
    } else {
        print_report_text(&report);
    }
    Ok(report.status == "proved")
}

fn print_attempts_text(attempts: &[RealRootednessTheoremReport], indent: &str) {
    for attempt in attempts {
        println!("\n{indent}{:?}: {:?}", attempt.theorem, attempt.status);
        if let Some(conclusion) = &attempt.conclusion {
            println!("{indent}  {}", conclusion.statement);
        } else {
            println!("{indent}  {}", attempt.normal_form);
        }
        for hypothesis in &attempt.hypotheses {
            println!(
                "{indent}  {} [{:?}]: {}",
                hypothesis.name, hypothesis.status, hypothesis.evidence
            );
        }
        println!("{indent}  Reference: {}", attempt.reference);
        if !attempt.lean_declarations.is_empty() {
            println!(
                "{indent}  Related Lean: {} (replay not attempted)",
                attempt.lean_declarations.join(", ")
            );
        }
        if let Some(certificate) = attempt.degree_law_reversal_certificate() {
            println!(
                "{indent}  Closure: polynomial reversal, {:?}",
                certificate.transform
            );
            print_attempts_text(
                std::slice::from_ref(certificate.transformed_proof.as_ref()),
                &format!("{indent}    "),
            );
        }
        if let Some(certificate) = attempt.index_section_certificate() {
            if let Some(proof) = certificate
                .selected
                .and_then(|i| certificate.parent_attempts.get(i))
            {
                println!(
                    "{indent}  Closure: Cayley--Hamilton index section, period {}; {:?}",
                    certificate.period, certificate.transform
                );
                print_attempts_text(std::slice::from_ref(proof), &format!("{indent}    "));
            }
        }
    }
}

fn print_report_text(report: &RecurrenceProofDiscoveryReport) {
    println!("Proof discovery: {}", report.status);
    if let Some(index) = report.real_rooted_from_index {
        println!("Every row P_n with n >= {index} is nonzero and real-rooted.");
    }
    if let Some(rows) = &report.real_rooted_outside_zero_residue_classes {
        let zero = &rows.zero_residue_classes;
        println!(
            "Zero residue classes: every row P_n with n >= {} and n mod {} in {:?} is zero; \
             every other row P_n with n >= {} is nonzero and real-rooted.",
            rows.from_index, zero.period, zero.residues, rows.from_index
        );
    }
    if let Some(reason) = &report.not_applicable_reason {
        println!("Not applicable: {reason}");
    }
    println!("{}", report.statement);
    if report.trivially_real_rooted {
        println!("Trivially real-rooted: product seed class, proved by product induction.");
    }
    if let Some(conditional) = report.trivially_real_rooted_if_rule_holds {
        if conditional != report.trivially_real_rooted {
            println!(
                "Trivially real-rooted if the rule holds: {}.",
                if conditional { "yes" } else { "no" }
            );
        }
    }
    print_attempts_text(&report.attempts, "");
    if let Some(tail) = &report.tail {
        println!("\nTail from index {}: {}", tail.first_index, tail.reason);
        print_attempts_text(&tail.attempts, "  ");
    }
    if let Some(search) = &report.alternative_representations {
        println!(
            "\nAlternative representations: {} fitted rule(s), {} finder configuration(s); {} \
             examined",
            search.fitted_rules,
            search.configurations_considered,
            search.alternatives.len()
        );
        for alternative in &search.alternatives {
            let equality = &alternative.equality;
            let rule = equality
                .alternative
                .to_recurrence_parts()
                .map(|(rule, _, _)| rule.to_string())
                .unwrap_or_default();
            println!("\n  Alternative rule: {rule}");
            println!(
                "  Equality with the family: {} {}",
                equality.status,
                equality.method.as_deref().unwrap_or("")
            );
            println!("  {}", equality.statement);
            for reason in &equality.reasons {
                println!("  - {reason}");
            }
            print_attempts_text(&alternative.attempts, "  ");
        }
    }
    if let Some(search) = &report.catalog_lifts {
        print_catalog_lifts_text(search);
    }
    if let Some(search) = &report.affine_coordinates {
        println!(
            "\nAffine coordinates: {} of {} candidate coordinate(s) tried",
            search.coordinates.len(),
            search.candidates_available
        );
        for candidate in &search.coordinates {
            let coordinate = &candidate.coordinate;
            println!(
                "\n  x = ({}) u + ({}) [{}]: {}",
                coordinate.slope,
                coordinate.intercept,
                coordinate.reason,
                if candidate.proved {
                    "proved"
                } else {
                    "not proved"
                }
            );
            if candidate.proved {
                print_attempts_text(&candidate.attempts, "  ");
                if let Some(tail) = &candidate.tail {
                    print_attempts_text(&tail.attempts, "    ");
                }
            }
        }
        match &search.conclusion {
            Some(conclusion) => println!("  Conclusion: {}", conclusion.statement),
            None => println!("  No coordinate proves the record."),
        }
    }
}

fn print_catalog_lifts_text(search: &CatalogLiftSearch) {
    println!(
        "\nCatalog lifts: {} catalog famil{} certified",
        search.candidates_tried(),
        if search.candidates_tried() == 1 {
            "y"
        } else {
            "ies"
        }
    );
    for identification in &search.identifications {
        println!(
            "  {}: equality {}; {}",
            identification.catalog_id,
            identification.equality.status,
            identification
                .reason
                .as_deref()
                .unwrap_or(match identification.parent_proof {
                    Some(polytool::recurrence::CatalogParentProofKind::Polytool) => {
                        "parent proved by polytool"
                    }
                    _ => "parent attested by the catalog flag",
                })
        );
    }
    if let Some(stage) = &search.basis_maps {
        match (&stage.search, &stage.error) {
            (Some(lifts), _) => println!(
                "  Basis maps: {} tried, {} certified candidate(s), proved: {}",
                lifts.maps_tried,
                lifts.lifts.len(),
                lifts
                    .proved()
                    .map_or("none", |lift| lift.candidate.basis.name())
            ),
            (None, Some(error)) => println!("  Basis maps: not run ({error})"),
            (None, None) => {}
        }
        for coefficient in &stage.coefficient_identifications {
            println!(
                "  Coefficient family of lift {}: {} (equality {})",
                coefficient.lift,
                coefficient.identification.catalog_id,
                coefficient.identification.equality.status
            );
        }
    }
    match &search.conclusion {
        Some(conclusion) => println!("  Conclusion: {}", conclusion.statement),
        None => println!("  No lift concludes."),
    }
}
