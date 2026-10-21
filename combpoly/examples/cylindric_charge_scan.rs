//! Scan the cocharge cyclic sieving conjectures for cylindric staircases `(1 0^z)^y`.
//!
//! For each shape: values of `CC` at primitive `d`-th roots (`d | n`) against the conjectured
//! counts, and the Alexandersson--Amini test for `C_n`. For rectangles (`z = 1`) also the
//! group `C_{2a}` and the divisibility of `X` by `prod_{j <= n/2} (1 + q^j)`.
//!
//! Usage: `cylindric_charge_scan Z_MAX Y_MAX N_MAX` (shapes with `n = yw <= N_MAX`).
use combpoly::cylindric::*;
use num_bigint::BigInt;
use num_traits::{Signed, Zero};
use polytool::cyclic_sieving::{
    cyclic_sieving_orbit_decomposition_bigint, root_of_unity_evaluation_bigint,
};

fn divide_by_one_plus_q_pow(p: &[BigInt], j: usize) -> Option<Vec<BigInt>> {
    // exact division by 1 + q^j
    let mut r = p.to_vec();
    if r.len() <= j {
        return None;
    }
    let mut quo = vec![BigInt::zero(); r.len() - j];
    for i in 0..quo.len() {
        quo[i] = r[i].clone();
        let c = r[i].clone();
        r[i] -= &c;
        r[i + j] -= &c;
    }
    if r.iter().all(|c| c.is_zero()) {
        Some(quo)
    } else {
        None
    }
}

fn main() {
    let args: Vec<usize> = std::env::args()
        .skip(1)
        .map(|a| a.parse().unwrap())
        .collect();
    let (zmax, ymax, nmax) = (args[0], args[1], args[2]);
    let (mut shapes, mut fails) = (0, 0);
    for z in 1..=zmax {
        for y in 2..=ymax {
            for w in 1.. {
                let n = y * w;
                if n > nmax {
                    break;
                }
                shapes += 1;
                let cc = staircase_cocharge_polynomial_bigint(z, y, w);
                let mut bad = vec![];
                for d in (2..=n).filter(|d| n % d == 0) {
                    let ev = root_of_unity_evaluation_bigint(&cc, d, 1).integer_value;
                    if ev != staircase_cocharge_root_of_unity_prediction(z, y, w, d) {
                        bad.push(format!("value d={d}"));
                    }
                }
                let aa = cyclic_sieving_orbit_decomposition_bigint(&cc, n).unwrap();
                if !aa.realizable {
                    bad.push("A-A at N=n".into());
                }
                if z == 1 {
                    let aa2 = cyclic_sieving_orbit_decomposition_bigint(&cc, 2 * w).unwrap();
                    if !aa2.realizable {
                        bad.push("A-A at N=2a".into());
                    }
                    let mut quo = staircase_charge_polynomial_bigint(1, y, w);
                    let mut ok = true;
                    for j in 1..=n / 2 {
                        match divide_by_one_plus_q_pow(&quo, j) {
                            Some(qq) => quo = qq,
                            None => {
                                ok = false;
                                break;
                            }
                        }
                    }
                    let nonneg = ok && quo.iter().all(|c| !c.is_negative());
                    let trimmed: Vec<&BigInt> = {
                        let lo = quo.iter().position(|c| !c.is_zero()).unwrap_or(0);
                        let hi = quo.iter().rposition(|c| !c.is_zero()).unwrap_or(0);
                        quo[lo..=hi].iter().collect()
                    };
                    let pal = trimmed.iter().eq(trimmed.iter().rev());
                    if !(ok && nonneg && pal) {
                        bad.push(format!(
                            "(-q;q) factor: divides={ok} nonneg={nonneg} palindromic={pal}"
                        ));
                    }
                }
                if !bad.is_empty() {
                    fails += 1;
                    println!("z={z} y={y} w={w} n={n}: FAIL {bad:?}");
                }
            }
        }
    }
    println!("shapes checked: {shapes}, failures: {fails}");
}
