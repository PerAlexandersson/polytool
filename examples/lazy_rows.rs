//! Lazy row generation: pause, continue and collect, for the OEIS catalog and a
//! plain recurrence.
use num_bigint::BigInt;
use polytool::oeis::by_id;
use polytool::recurrence::{BigRational, BivarPoly, Recurrence, RecurrenceTerm};
use polytool::sequences::chebyshev_polynomials_t_iter;

fn main() {
    // OEIS A008292 (Eulerian triangle): rows are computed on demand.
    let eulerian = by_id("A008292").expect("A008292 is bundled");
    let mut rows = eulerian
        .generate_rows_iter()
        .expect("the bundled prefix decodes");
    let first_three: Vec<Vec<BigInt>> = rows
        .by_ref()
        .take(3)
        .collect::<Result<_, _>>()
        .expect("the recurrence is exact");
    println!("first three rows: {first_three:?}");

    // The same iterator continues with row 4 without recomputing rows 1-3.
    let fourth = rows.next().expect("endless stream").expect("exact row");
    println!("row {}: {fourth:?}", rows.next_displayed_row().unwrap() - 1);

    // Collecting a prefix of a fresh iterator is exactly the batch form.
    let batch = eulerian.generate_rows(4).expect("batch rows");
    let lazy: Vec<Vec<BigInt>> = eulerian
        .generate_rows_iter()
        .unwrap()
        .take(4)
        .collect::<Result<_, _>>()
        .unwrap();
    assert_eq!(batch, lazy);
    assert_eq!(batch[3], fourth);

    // Start at a displayed row; earlier rows are traversed but not retained.
    let from_row_five = eulerian
        .generate_rows_from_iter(5)
        .unwrap()
        .next()
        .unwrap()
        .unwrap();
    println!("row 5: {from_row_five:?}");

    // The shared kernel for an arbitrary recurrence: P_n = n P_{n-1}, P_0 = 1.
    let n_coefficient = BivarPoly::from_coefficients_normalized(vec![
        vec![BigRational::from_integer(0.into())],
        vec![BigRational::from_integer(1.into())],
    ]);
    let factorials = Recurrence {
        terms: vec![RecurrenceTerm::new(1, 0, n_coefficient)],
        denominator: None,
        inhomogeneous: None,
    };
    let seed = vec![vec![BigRational::from_integer(1.into())]];
    let mut iterator = factorials.generate_rows_rational_iter(&seed, 0);
    let _ = iterator.by_ref().take(5).count();
    println!(
        "after 5 rows: {} steps evaluated, {} rows retained, next index {:?}",
        iterator.steps_evaluated(),
        iterator.retained_row_count(),
        iterator.next_index()
    );
    assert_eq!(
        iterator.next().unwrap().unwrap(),
        vec![BigRational::from_integer(120.into())]
    );

    // A named family: Chebyshev polynomials T_n, endless and fused.
    let t5 = chebyshev_polynomials_t_iter().nth(5).unwrap();
    println!("T_5 coefficients: {t5:?}");
}
