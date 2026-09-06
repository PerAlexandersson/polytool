use sym_poly_multipoly::{pi_i, rothe_diagram, schubert_polynomial, MultiPoly};

#[derive(Debug)]
struct Orthodontia {
    initial_intervals: Vec<usize>,
    word: Vec<usize>,
    multiplicities: Vec<usize>,
}

fn main() {
    let cases = [
        (vec![2, 1, 4, 3], 1usize, 7usize),
        (vec![3, 2, 1, 6, 5, 4], 2usize, 5usize),
    ];

    for (perm, first_dilation, last_dilation) in cases {
        println!("u={perm:?}, code={:?}", lehmer_code(&perm));
        for dilation in first_dilation..=last_dilation {
            let stretched = stretch_perm(&perm, dilation);
            let columns = rothe_columns(&stretched);
            let data = orthodontia(columns, perm.len());
            let reconstructed = orthodontic_polynomial(stretched.len(), &data);
            let schubert = schubert_polynomial::<i64>(&stretched);
            println!(
                "  N={dilation}: k={:?}, word={:?}, m={:?}, match={}",
                data.initial_intervals,
                data.word,
                data.multiplicities,
                reconstructed == schubert
            );
            assert_eq!(reconstructed, schubert);
        }
    }
}

fn orthodontia(mut columns: Vec<Vec<usize>>, height: usize) -> Orthodontia {
    let mut initial_intervals = Vec::with_capacity(height);
    for size in 1..=height {
        let interval = (1..=size).collect::<Vec<_>>();
        initial_intervals.push(columns.iter().filter(|column| **column == interval).count());
        columns.retain(|column| *column != interval);
    }

    let mut word = Vec::new();
    let mut multiplicities = Vec::new();
    while let Some(first) = columns.first() {
        let tooth = (1..height)
            .find(|&row| !first.contains(&row) && first.contains(&(row + 1)))
            .expect("a non-interval column has a missing tooth");

        for column in &mut columns {
            for row in column.iter_mut() {
                if *row == tooth {
                    *row += 1;
                } else if *row == tooth + 1 {
                    *row -= 1;
                }
            }
            column.sort_unstable();
        }

        let interval = (1..=tooth).collect::<Vec<_>>();
        let multiplicity = columns.iter().filter(|column| **column == interval).count();
        columns.retain(|column| *column != interval);
        word.push(tooth);
        multiplicities.push(multiplicity);
    }

    Orthodontia {
        initial_intervals,
        word,
        multiplicities,
    }
}

fn orthodontic_polynomial(num_vars: usize, data: &Orthodontia) -> MultiPoly<i64> {
    let mut result = MultiPoly::constant(num_vars, 1);
    for (&tooth, &multiplicity) in data.word.iter().zip(&data.multiplicities).rev() {
        result = interval_monomial(num_vars, tooth, multiplicity) * result;
        result = pi_i(&result, tooth - 1);
    }

    for (idx, &multiplicity) in data.initial_intervals.iter().enumerate() {
        result = interval_monomial(num_vars, idx + 1, multiplicity) * result;
    }
    result
}

fn interval_monomial(num_vars: usize, size: usize, power: usize) -> MultiPoly<i64> {
    let mut exponent = vec![0; num_vars];
    exponent[..size].fill(power as u32);
    MultiPoly::x_power(num_vars, exponent)
}

fn rothe_columns(perm: &[usize]) -> Vec<Vec<usize>> {
    let diagram = rothe_diagram(perm);
    let max_col = diagram.iter().map(|cell| cell.col).max().unwrap_or(0);
    (1..=max_col)
        .filter_map(|col| {
            let rows = diagram
                .iter()
                .filter(|cell| cell.col == col)
                .map(|cell| cell.row)
                .collect::<Vec<_>>();
            (!rows.is_empty()).then_some(rows)
        })
        .collect()
}

fn lehmer_code(perm: &[usize]) -> Vec<usize> {
    (0..perm.len())
        .map(|i| (i + 1..perm.len()).filter(|&j| perm[j] < perm[i]).count())
        .collect()
}

fn stretch_perm(perm: &[usize], dilation: usize) -> Vec<usize> {
    let code = lehmer_code(perm);
    let needed_len = code
        .iter()
        .enumerate()
        .filter(|(_, entry)| **entry > 0)
        .map(|(i, entry)| i + 1 + dilation * entry)
        .max()
        .unwrap_or(perm.len())
        .max(perm.len());
    let mut stretched_code = vec![0; needed_len];
    for (i, entry) in code.into_iter().enumerate() {
        stretched_code[i] = dilation * entry;
    }
    from_lehmer_code(&stretched_code)
}

fn from_lehmer_code(code: &[usize]) -> Vec<usize> {
    let mut available = (1..=code.len()).collect::<Vec<_>>();
    code.iter().map(|&entry| available.remove(entry)).collect()
}
