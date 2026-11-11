//! A replayable all-index PF Riordan proof, with exact generating-function output.

use polytool::oeis;
use polytool::recurrence::{
    check_rational_pf_riordan_theorem, rational_generating_function, riordan_decomposition,
};

fn main() {
    let entry = oeis::by_id("A216182").unwrap();
    let (rule, seeds) = entry.recurrence_parts().unwrap();
    let first = entry.recurrence_first_index;
    let gf = rational_generating_function(&rule, first, &seeds).unwrap();
    let riordan = riordan_decomposition(&gf).unwrap();
    println!("F(x,z) = {gf}");
    println!("d(z) = {}; K(z) = {}", riordan.d, riordan.k);
    let proof =
        check_rational_pf_riordan_theorem(&rule, first, &seeds, &Default::default()).unwrap();
    assert!(proof.is_proved());
    proof.verify().unwrap();
    println!("{}", proof.conclusion.unwrap().statement);
}
