//! Small checked walk-through of `combinatoric_core::colored_exchange`.
//!
//! Run with `cargo run -p combinatoric-core --example colored_exchange`.
//! Every printed claim is asserted, so the example doubles as a smoke test.

use combinatoric_core::colored_exchange::{word_to_letters, CertificateError};
use combinatoric_core::FixedContentFamily;
use num_bigint::BigUint;

fn letters(words: &[Vec<usize>]) -> String {
    words
        .iter()
        .map(|w| word_to_letters(w))
        .collect::<Vec<_>>()
        .join(" ")
}

fn main() {
    // The cyclic three-word family: every pair exchanges by a full 3-cycle.
    let cyclic = FixedContentFamily::from_letter_words(3, &["ABC", "BCA", "CAB"]).unwrap();
    cyclic.check_symmetric_cycle_exchange().unwrap();
    let witness = cyclic.exchange_witness(0, 1, 0).unwrap().unwrap();
    cyclic.verify_exchange(0, 1, 0, &witness).unwrap();
    println!(
        "{{{}}} satisfies SCE; witness for words 0,1 at position 0: {:?}",
        letters(cyclic.words()),
        witness.positions
    );

    // The exact two-color failure and its repair.
    let bad = FixedContentFamily::from_letter_words(2, &["ABAB", "BABA"]).unwrap();
    let obstruction = bad.check_symmetric_cycle_exchange().unwrap_err();
    println!("{{ABAB, BABA}} fails: {obstruction}");
    let good = FixedContentFamily::from_letter_words(2, &["ABAB", "ABBA", "BAAB", "BABA"]).unwrap();
    good.check_symmetric_cycle_exchange().unwrap();
    println!(
        "adding both offspring repairs it: {{{}}} satisfies SCE",
        letters(good.words())
    );

    // Exact weighted A-marginal of the cyclic family with weights 10^(3p + c).
    let weights: Vec<Vec<BigUint>> = (0..3)
        .map(|p| {
            (0..3)
                .map(|c| BigUint::from(10u32).pow(3 * p + c))
                .collect()
        })
        .collect();
    let marginal = cyclic.weighted_color_marginal(0, &weights).unwrap();
    for (positions, coefficient) in &marginal {
        println!("  A-marginal coefficient of x_{positions:?}: {coefficient}");
    }
    assert_eq!(marginal[&vec![0]], BigUint::from(10u32).pow(4 + 8));

    // Binary-port decomposition of {AACB, ABCA, BACA} across the cut {0, 1}.
    let glued = FixedContentFamily::from_letter_words(3, &["AACB", "ABCA", "BACA"]).unwrap();
    let binary = glued.binary_port_certificate(&[0, 1]).unwrap();
    binary.verify(&glued).unwrap();
    println!(
        "{{{}}} = binary gluing of {{{}}} and {{{}}} (ports last, counted color {}, low count {})",
        letters(glued.words()),
        letters(binary.left.words()),
        letters(binary.right.words()),
        word_to_letters(&[binary.counted_color]),
        binary.low_count
    );

    // Private-block decomposition of the cloned cyclic family.
    let cloned =
        FixedContentFamily::from_letter_words(4, &["ADBC", "DABC", "BDCA", "DBCA", "CDAB", "DCAB"])
            .unwrap();
    let block = cloned.private_block_certificate(&[0, 1], &[3]).unwrap();
    block.verify(&cloned).unwrap();
    println!(
        "{{{}}} = substitution of {{{}}} into {{{}}} with allowed positions {:?}",
        letters(cloned.words()),
        letters(block.internal.words()),
        letters(block.compressed.words()),
        block.allowed_positions
    );

    // A corrupted certificate is rejected with a precise reason.
    let mut corrupted = block.clone();
    corrupted.allowed_positions.insert(0, vec![0]);
    match corrupted.verify(&cloned) {
        Err(CertificateError::ReconstructionMismatch { missing, .. }) => {
            println!("narrowing the A-list loses {}", letters(&missing));
        }
        other => panic!("unexpected verification result {other:?}"),
    }
}
