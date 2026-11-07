//! A rectangular Schur deformation, evaluated by exact dual Jacobi–Trudi minors.
use polytool::sequences::schur::rectangular_elementary_substitution_polynomials_bigint;

fn main() {
    let alphabet = vec![1.into(); 4];
    let rows = rectangular_elementary_substitution_polynomials_bigint(&alphabet, 2, 1, 6)
        .expect("positive alphabet and valid indices");
    for (n, row) in rows.iter().enumerate() {
        println!("F_{n}: {row:?}");
    }
    assert_eq!(rows[3], vec![50.into(), 40.into(), 1.into()]);
}
