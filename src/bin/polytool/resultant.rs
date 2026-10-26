//! `polytool resultant` and `polytool discriminant`.

use crate::input::{read_polys_bigint, strip_trailing_zeros_bigint};
use polytool::*;

pub(crate) fn cmd_resultant() {
    let polys = read_polys_bigint();
    if polys.len() < 2 {
        eprintln!("Need exactly two polynomials for resultant.");
        return;
    }
    let p = strip_trailing_zeros_bigint(&polys[0]);
    let q = strip_trailing_zeros_bigint(&polys[1]);
    let r = resultant_bigint_coeffs(p, q);
    println!(
        "Res({}, {}) = {}",
        format_poly_bigint_coeffs(p),
        format_poly_bigint_coeffs(q),
        r
    );
}

pub(crate) fn cmd_discriminant() {
    for coeffs in read_polys_bigint() {
        let c = strip_trailing_zeros_bigint(&coeffs);
        let d = discriminant_bigint_coeffs(c);
        println!("disc({}) = {}", format_poly_bigint_coeffs(c), d);
    }
}
