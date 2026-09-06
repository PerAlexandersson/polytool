use combinatoric_core::next_permutation;
use sym_poly_multipoly::{schubert_polynomial, MultiPoly};

fn main() {
    let mut checks = 0usize;
    checks += check_rank(2, 5);
    checks += check_rank(3, 3);
    println!("verified {checks} Schubert dual-basis pairings");
}

fn check_rank(rank: usize, max_degree: usize) -> usize {
    let rho: Vec<usize> = (0..rank).map(|i| rank - 1 - i).collect();
    let mut checks = 0usize;

    for degree in 0..=max_degree {
        let compositions = weak_compositions(degree, rank);
        for lambda in &compositions {
            for mu in &compositions {
                let shift = mu
                    .iter()
                    .zip(&rho)
                    .map(|(&part, &rho_part)| part.saturating_sub(rho_part))
                    .max()
                    .unwrap_or(0);
                let dual_code: Vec<usize> = rho
                    .iter()
                    .zip(mu)
                    .map(|(&rho_part, &part)| rho_part + shift - part)
                    .collect();
                let ambient_rank = required_rank(lambda)
                    .max(required_rank(&dual_code))
                    .max(rank);
                let schubert_lambda = schubert_from_code(lambda, ambient_rank);
                let schubert_dual = schubert_from_code(&dual_code, ambient_rank);
                let product = schubert_lambda * schubert_dual;

                let pairing = alternating_coefficient_sum(&product, &rho, shift, ambient_rank);
                let expected = i64::from(lambda == mu);
                assert_eq!(
                    pairing, expected,
                    "duality failed in rank {rank} for lambda={lambda:?}, mu={mu:?}"
                );
                checks += 1;
            }
        }
    }

    checks
}

fn alternating_coefficient_sum(
    product: &MultiPoly<i64>,
    rho: &[usize],
    shift: usize,
    ambient_rank: usize,
) -> i64 {
    let mut sigma: Vec<usize> = (0..rho.len()).collect();
    let mut result = 0i64;

    loop {
        let mut exponent = vec![0u32; ambient_rank];
        for i in 0..rho.len() {
            exponent[i] = (shift + rho[sigma[i]]) as u32;
        }
        let sign = if inversion_count(&sigma).is_multiple_of(2) {
            1
        } else {
            -1
        };
        result += sign * product.coefficient(&exponent);

        if !next_permutation(&mut sigma) {
            break;
        }
    }

    result
}

fn schubert_from_code(code: &[usize], ambient_rank: usize) -> MultiPoly<i64> {
    let mut padded_code = vec![0usize; ambient_rank];
    padded_code[..code.len()].copy_from_slice(code);
    schubert_polynomial(&from_lehmer_code(&padded_code))
}

fn required_rank(code: &[usize]) -> usize {
    code.iter()
        .enumerate()
        .map(|(i, &part)| i + 1 + part)
        .max()
        .unwrap_or(code.len())
        .max(code.len())
}

fn from_lehmer_code(code: &[usize]) -> Vec<usize> {
    let mut available: Vec<usize> = (1..=code.len()).collect();
    let mut permutation = Vec::with_capacity(code.len());
    for &part in code {
        permutation.push(available.remove(part));
    }
    permutation
}

fn weak_compositions(total: usize, length: usize) -> Vec<Vec<usize>> {
    fn generate(
        remaining: usize,
        slots: usize,
        current: &mut Vec<usize>,
        output: &mut Vec<Vec<usize>>,
    ) {
        if slots == 1 {
            current.push(remaining);
            output.push(current.clone());
            current.pop();
            return;
        }

        for part in 0..=remaining {
            current.push(part);
            generate(remaining - part, slots - 1, current, output);
            current.pop();
        }
    }

    let mut output = Vec::new();
    generate(total, length, &mut Vec::new(), &mut output);
    output
}

fn inversion_count(permutation: &[usize]) -> usize {
    let mut inversions = 0usize;
    for i in 0..permutation.len() {
        for j in i + 1..permutation.len() {
            inversions += usize::from(permutation[i] > permutation[j]);
        }
    }
    inversions
}
