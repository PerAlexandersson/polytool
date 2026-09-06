use std::collections::BTreeMap;

use combinatoric_core::next_permutation;
use sym_poly_core::UnivariatePolynomial;
use sym_poly_multipoly::beta_grothendieck_to_lascoux;

type BetaPolynomial = UnivariatePolynomial<i64>;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.get(1).map(String::as_str) == Some("--fit") {
        let permutation = parse_permutation(args.get(2).expect("missing permutation"));
        let max_dilation = args
            .get(3)
            .and_then(|value| value.parse().ok())
            .unwrap_or(6);
        print_affine_onset(&permutation, max_dilation);
        return;
    }
    if args.get(1).map(String::as_str) == Some("--scan") {
        let rank = args
            .get(2)
            .and_then(|value| value.parse().ok())
            .unwrap_or(4);
        let max_dilation = args
            .get(3)
            .and_then(|value| value.parse().ok())
            .unwrap_or(5);
        scan_rank(rank, max_dilation);
        return;
    }
    if args.len() >= 2 {
        let permutation = parse_permutation(&args[1]);
        let max_dilation = args
            .get(2)
            .and_then(|value| value.parse().ok())
            .unwrap_or(4);
        print_stretched_expansions(&permutation, max_dilation);
        return;
    }

    for (permutation, max_dilation) in [
        (vec![2, 1, 4, 3], 4),
        (vec![2, 1, 3, 5, 4], 3),
        (vec![3, 2, 1, 5, 4], 2),
    ] {
        print_stretched_expansions(&permutation, max_dilation);
    }
}

fn print_affine_onset(permutation: &[usize], max_dilation: usize) {
    assert!(max_dilation >= 3, "an affine fit needs three dilations");
    let code = lehmer_code(permutation);
    let active_rank = code
        .iter()
        .rposition(|&part| part > 0)
        .map(|index| index + 1)
        .unwrap_or(0);
    let beta = BetaPolynomial::variable();
    let expansions = (1..=max_dilation)
        .map(|dilation| normalized_expansion(permutation, dilation, active_rank, &beta))
        .collect::<Vec<_>>();
    let onset = (0..=max_dilation - 3)
        .find(|&start| affine_tail(&expansions[start..]))
        .map(|start| start + 1);
    let term_counts = expansions.iter().map(Vec::len).collect::<Vec<_>>();
    println!(
        "u={permutation:?}, code={:?}, term_counts={term_counts:?}, affine_onset={onset:?}",
        trim_trailing_zeroes(&code)
    );
}

fn scan_rank(rank: usize, max_dilation: usize) {
    assert!(max_dilation >= 3, "a linearity scan needs three dilations");
    let beta = BetaPolynomial::variable();
    let mut permutation = (1..=rank).collect::<Vec<_>>();
    let mut checked = 0usize;
    let mut failures = Vec::new();
    let mut onset_counts = BTreeMap::<usize, usize>::new();
    let mut onset_examples = BTreeMap::<usize, Vec<Vec<usize>>>::new();
    let mut largest_expansion = (0usize, Vec::new());

    loop {
        let code = lehmer_code(&permutation);
        let active_rank = code
            .iter()
            .rposition(|&part| part > 0)
            .map(|index| index + 1)
            .unwrap_or(0);
        let expansions = (1..=max_dilation)
            .map(|dilation| normalized_expansion(&permutation, dilation, active_rank, &beta))
            .collect::<Vec<_>>();

        let max_terms = expansions.iter().map(Vec::len).max().unwrap_or(0);
        if max_terms > largest_expansion.0 {
            largest_expansion = (max_terms, permutation.clone());
        }

        let onset = (0..=max_dilation - 3).find(|&start| affine_tail(&expansions[start..]));
        match onset {
            Some(start) => {
                *onset_counts.entry(start + 1).or_insert(0) += 1;
                onset_examples
                    .entry(start + 1)
                    .or_default()
                    .push(permutation.clone());
            }
            None => failures.push(permutation.clone()),
        }
        checked += 1;

        if !next_permutation(&mut permutation) {
            break;
        }
    }

    println!("S_{rank} scan through N={max_dilation}: checked {checked}");
    println!("  affine onsets: {onset_counts:?}");
    for (onset, examples) in onset_examples {
        if onset > 1 {
            let preview = examples.iter().take(10).collect::<Vec<_>>();
            println!(
                "  onset N={onset}: {} permutations, first examples {preview:?}",
                examples.len()
            );
        }
    }
    println!("  failures: {} {failures:?}", failures.len());
    println!(
        "  largest expansion: {} terms for {:?}",
        largest_expansion.0, largest_expansion.1
    );
}

fn normalized_expansion(
    permutation: &[usize],
    dilation: usize,
    active_rank: usize,
    beta: &BetaPolynomial,
) -> Vec<(Vec<u32>, BetaPolynomial)> {
    let stretched = stretch_permutation(permutation, dilation);
    beta_grothendieck_to_lascoux(&stretched, beta)
        .into_iter()
        .map(|(index, coefficient)| {
            assert!(
                index[active_rank..].iter().all(|&part| part == 0),
                "Lascoux index has support beyond the fixed active rank"
            );
            (index[..active_rank].to_vec(), coefficient)
        })
        .collect()
}

fn affine_tail(expansions: &[Vec<(Vec<u32>, BetaPolynomial)>]) -> bool {
    if expansions.len() < 3 {
        return false;
    }
    let term_count = expansions[0].len();
    if expansions
        .iter()
        .any(|expansion| expansion.len() != term_count)
    {
        return false;
    }

    for term in 0..term_count {
        let coefficient = &expansions[0][term].1;
        if expansions
            .iter()
            .any(|expansion| &expansion[term].1 != coefficient)
        {
            return false;
        }

        let first_difference = vector_difference(&expansions[1][term].0, &expansions[0][term].0);
        for window in expansions.windows(2).skip(1) {
            if vector_difference(&window[1][term].0, &window[0][term].0) != first_difference {
                return false;
            }
        }
    }
    true
}

fn vector_difference(left: &[u32], right: &[u32]) -> Vec<i64> {
    left.iter()
        .zip(right)
        .map(|(&left_part, &right_part)| i64::from(left_part) - i64::from(right_part))
        .collect()
}

fn print_stretched_expansions(permutation: &[usize], max_dilation: usize) {
    println!(
        "u={permutation:?}, code={:?}",
        trim_trailing_zeroes(&lehmer_code(permutation))
    );
    let beta = BetaPolynomial::variable();

    for dilation in 1..=max_dilation {
        let stretched = stretch_permutation(permutation, dilation);
        let expansion = beta_grothendieck_to_lascoux(&stretched, &beta);
        let trimmed = expansion
            .into_iter()
            .map(|(index, coefficient)| (trim_trailing_zeroes_u32(&index), coefficient))
            .collect::<BTreeMap<_, _>>();

        println!(
            "  N={dilation}: ambient_rank={}, terms={}",
            stretched.len(),
            trimmed.len()
        );
        for (index, coefficient) in trimmed {
            println!("    {} * L_{index:?}", format_beta_polynomial(&coefficient));
        }
    }
}

fn format_beta_polynomial(polynomial: &BetaPolynomial) -> String {
    let terms = polynomial
        .coeffs()
        .iter()
        .enumerate()
        .filter_map(|(degree, &coefficient)| {
            if coefficient == 0 {
                return None;
            }
            Some(match (coefficient, degree) {
                (1, 0) => "1".to_string(),
                (value, 0) => value.to_string(),
                (1, 1) => "beta".to_string(),
                (value, 1) => format!("{value}*beta"),
                (1, power) => format!("beta^{power}"),
                (value, power) => format!("{value}*beta^{power}"),
            })
        })
        .collect::<Vec<_>>();
    if terms.is_empty() {
        "0".to_string()
    } else {
        terms.join(" + ")
    }
}

fn parse_permutation(text: &str) -> Vec<usize> {
    if text.contains(',') {
        return text
            .split(',')
            .map(|part| part.trim().parse().expect("invalid permutation entry"))
            .collect();
    }
    text.chars()
        .map(|character| {
            character
                .to_digit(10)
                .expect("compact permutations must contain digits") as usize
        })
        .collect()
}

fn lehmer_code(permutation: &[usize]) -> Vec<usize> {
    let mut code = vec![0usize; permutation.len()];
    for i in 0..permutation.len() {
        for j in i + 1..permutation.len() {
            code[i] += usize::from(permutation[i] > permutation[j]);
        }
    }
    code
}

fn stretch_permutation(permutation: &[usize], dilation: usize) -> Vec<usize> {
    let code = lehmer_code(permutation);
    let needed_len = code
        .iter()
        .enumerate()
        .filter(|(_, part)| **part > 0)
        .map(|(i, part)| i + 1 + dilation * part)
        .max()
        .unwrap_or(permutation.len())
        .max(permutation.len());
    let mut stretched_code = vec![0usize; needed_len];
    for (target, part) in stretched_code.iter_mut().zip(code) {
        *target = dilation * part;
    }
    from_lehmer_code(&stretched_code)
}

fn from_lehmer_code(code: &[usize]) -> Vec<usize> {
    let mut available = (1..=code.len()).collect::<Vec<_>>();
    code.iter().map(|&part| available.remove(part)).collect()
}

fn trim_trailing_zeroes(values: &[usize]) -> Vec<usize> {
    let mut trimmed = values.to_vec();
    while trimmed.last() == Some(&0) {
        trimmed.pop();
    }
    trimmed
}

fn trim_trailing_zeroes_u32(values: &[u32]) -> Vec<u32> {
    let mut trimmed = values.to_vec();
    while trimmed.last() == Some(&0) {
        trimmed.pop();
    }
    trimmed
}
