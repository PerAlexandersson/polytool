//! Generate benchmark fixtures for adaptive recurrence search.
//!
//! Run from the Rust workspace root:
//!
//! ```text
//! cargo run -p polytool --example generate_recurrence_benchmarks
//! ```

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Zero};
use std::fs;
use std::path::Path;

type Rat = BigRational;
type Poly = Vec<Rat>;
type NextFn = dyn Fn(usize, &[Poly]) -> Poly;

struct Fixture {
    slug: &'static str,
    title: &'static str,
    features: &'static str,
    suggested_args: &'static str,
    recurrence: &'static str,
    initial: Vec<Poly>,
    next: Box<NextFn>,
}

fn rat(n: i64) -> Rat {
    Rat::from_integer(BigInt::from(n))
}

fn rat_usize(n: usize) -> Rat {
    Rat::from_integer(BigInt::from(n))
}

fn trim(mut p: Poly) -> Poly {
    while p.len() > 1 && p.last().is_some_and(|c| c.is_zero()) {
        p.pop();
    }
    if p.is_empty() {
        vec![Rat::zero()]
    } else {
        p
    }
}

fn poly(coeffs: &[i64]) -> Poly {
    trim(coeffs.iter().map(|&c| rat(c)).collect())
}

fn prev(rows: &[Poly], offset: usize) -> &Poly {
    &rows[rows.len() - offset]
}

fn add(a: &[Rat], b: &[Rat]) -> Poly {
    let len = a.len().max(b.len());
    let mut out = vec![Rat::zero(); len];
    for i in 0..len {
        if i < a.len() {
            out[i] += a[i].clone();
        }
        if i < b.len() {
            out[i] += b[i].clone();
        }
    }
    trim(out)
}

fn add3(a: &[Rat], b: &[Rat], c: &[Rat]) -> Poly {
    add(&add(a, b), c)
}

fn sum_polys(polys: &[Poly]) -> Poly {
    let mut out = vec![Rat::zero()];
    for p in polys {
        out = add(&out, p);
    }
    out
}

fn scale(a: &[Rat], c: Rat) -> Poly {
    if c.is_zero() {
        return vec![Rat::zero()];
    }
    trim(a.iter().map(|x| x.clone() * c.clone()).collect())
}

fn scale_i(a: &[Rat], c: i64) -> Poly {
    scale(a, rat(c))
}

fn shift_t(a: &[Rat], amount: usize) -> Poly {
    if amount == 0 {
        return a.to_vec();
    }
    if a.iter().all(|c| c.is_zero()) {
        return vec![Rat::zero()];
    }
    let mut out = vec![Rat::zero(); amount];
    out.extend_from_slice(a);
    trim(out)
}

fn mul_linear(constant: Rat, t_coeff: Rat, p: &[Rat]) -> Poly {
    add(&scale(p, constant), &scale(&shift_t(p, 1), t_coeff))
}

fn mul_t_poly(coeffs: &[Rat], p: &[Rat]) -> Poly {
    if coeffs.iter().all(|c| c.is_zero()) || p.iter().all(|c| c.is_zero()) {
        return vec![Rat::zero()];
    }
    let mut out = vec![Rat::zero(); coeffs.len() + p.len() - 1];
    for (i, coeff) in coeffs.iter().enumerate() {
        if coeff.is_zero() {
            continue;
        }
        for (j, value) in p.iter().enumerate() {
            if value.is_zero() {
                continue;
            }
            out[i + j] += coeff.clone() * value.clone();
        }
    }
    trim(out)
}

fn derivative(p: &[Rat]) -> Poly {
    if p.len() <= 1 {
        return vec![Rat::zero()];
    }
    trim(
        p.iter()
            .enumerate()
            .skip(1)
            .map(|(i, c)| c.clone() * rat_usize(i))
            .collect(),
    )
}

fn second_derivative(p: &[Rat]) -> Poly {
    nth_derivative(p, 2)
}

fn nth_derivative(p: &[Rat], order: usize) -> Poly {
    let mut out = p.to_vec();
    for _ in 0..order {
        out = derivative(&out);
    }
    out
}

fn alternating_sign(n: usize) -> Rat {
    if n & 1 == 0 {
        Rat::one()
    } else {
        -Rat::one()
    }
}

fn generate_rows(fixture: &Fixture, rows: usize) -> Vec<Poly> {
    let mut out = fixture.initial.clone();
    while out.len() < rows {
        let n = out.len() + 1;
        let next = (fixture.next)(n, &out);
        out.push(trim(next));
    }
    out
}

fn format_rat(x: &Rat) -> String {
    if x.denom() == &BigInt::one() {
        x.numer().to_string()
    } else {
        format!("{}/{}", x.numer(), x.denom())
    }
}

fn format_row(p: &[Rat]) -> String {
    trim(p.to_vec())
        .iter()
        .map(format_rat)
        .collect::<Vec<_>>()
        .join(", ")
}

fn fixtures() -> Vec<Fixture> {
    vec![
        Fixture {
            slug: "01_scalar_geometric",
            title: "Scalar geometric",
            features: "constant coefficient, first order",
            suggested_args: "--max-rec-len 1 --max-var-deg 0 --max-idx-deg 0 --max-diff-deg 0",
            recurrence: "P_n = 2 P_{n-1}",
            initial: vec![poly(&[1])],
            next: Box::new(|_, rows| scale_i(prev(rows, 1), 2)),
        },
        Fixture {
            slug: "02_scalar_fibonacci",
            title: "Scalar Fibonacci",
            features: "constant coefficient, second order",
            suggested_args: "--max-rec-len 2 --max-var-deg 0 --max-idx-deg 0 --max-diff-deg 0",
            recurrence: "P_n = P_{n-1} + P_{n-2}",
            initial: vec![poly(&[1]), poly(&[1])],
            next: Box::new(|_, rows| add(prev(rows, 1), prev(rows, 2))),
        },
        Fixture {
            slug: "03_binomial_powers",
            title: "Binomial powers",
            features: "t-dependent coefficient",
            suggested_args: "--max-rec-len 1 --max-var-deg 1 --max-idx-deg 0 --max-diff-deg 0",
            recurrence: "P_n = (1+t) P_{n-1}",
            initial: vec![poly(&[1])],
            next: Box::new(|_, rows| mul_linear(Rat::one(), Rat::one(), prev(rows, 1))),
        },
        Fixture {
            slug: "04_chebyshev_t",
            title: "Chebyshev T",
            features: "t-dependent coefficient, second order",
            suggested_args: "--max-rec-len 2 --max-var-deg 1 --max-idx-deg 0 --max-diff-deg 0",
            recurrence: "P_n = 2t P_{n-1} - P_{n-2}",
            initial: vec![poly(&[1]), poly(&[0, 1])],
            next: Box::new(|_, rows| add(&scale_i(&shift_t(prev(rows, 1), 1), 2), &scale_i(prev(rows, 2), -1))),
        },
        Fixture {
            slug: "05_factorial_index",
            title: "Factorial index coefficient",
            features: "n-dependent coefficient",
            suggested_args: "--max-rec-len 1 --max-var-deg 0 --max-idx-deg 1 --max-diff-deg 0",
            recurrence: "P_n = n P_{n-1}",
            initial: vec![poly(&[1])],
            next: Box::new(|n, rows| scale(prev(rows, 1), rat_usize(n))),
        },
        Fixture {
            slug: "06_affine_product",
            title: "Affine product",
            features: "n- and t-dependent coefficient",
            suggested_args: "--max-rec-len 1 --max-var-deg 1 --max-idx-deg 1 --max-diff-deg 0",
            recurrence: "P_n = (n+t) P_{n-1}",
            initial: vec![poly(&[1])],
            next: Box::new(|n, rows| mul_linear(rat_usize(n), Rat::one(), prev(rows, 1))),
        },
        Fixture {
            slug: "07_hermite_indexed_second_order",
            title: "Hermite indexed second order",
            features: "n-dependent coefficient, second order",
            suggested_args: "--max-rec-len 2 --max-var-deg 1 --max-idx-deg 1 --max-diff-deg 0",
            recurrence: "P_n = 2t P_{n-1} - 2(n-2) P_{n-2}",
            initial: vec![poly(&[1]), poly(&[0, 2])],
            next: Box::new(|n, rows| {
                add(
                    &scale_i(&shift_t(prev(rows, 1), 1), 2),
                    &scale(prev(rows, 2), -rat_usize(2 * (n - 2))),
                )
            }),
        },
        Fixture {
            slug: "08_inhomogeneous_linear",
            title: "Inhomogeneous linear",
            features: "inhomogeneous, degree one in n and t",
            suggested_args: "--inhomogeneous --max-rec-len 1 --max-var-deg 0 --max-idx-deg 0 --max-diff-deg 0 --max-inhomo-var-deg 1 --max-inhomo-idx-deg 1",
            recurrence: "P_n = P_{n-1} + n + t",
            initial: vec![poly(&[1])],
            next: Box::new(|n, rows| add(prev(rows, 1), &[rat_usize(n), Rat::one()])),
        },
        Fixture {
            slug: "09_inhomogeneous_quadratic",
            title: "Inhomogeneous quadratic",
            features: "inhomogeneous, degree two in n and t",
            suggested_args: "--inhomogeneous --max-rec-len 1 --max-var-deg 0 --max-idx-deg 0 --max-diff-deg 0 --max-inhomo-var-deg 2 --max-inhomo-idx-deg 2",
            recurrence: "P_n = P_{n-1} + n^2 + nt + t^2",
            initial: vec![poly(&[0])],
            next: Box::new(|n, rows| {
                let n_rat = rat_usize(n);
                let n_sq = n_rat.clone() * n_rat.clone();
                add(prev(rows, 1), &[n_sq, n_rat, Rat::one()])
            }),
        },
        Fixture {
            slug: "10_alternating_scalar",
            title: "Alternating scalar",
            features: "alternating sign",
            suggested_args: "--alternating-sign --max-rec-len 1 --max-var-deg 0 --max-idx-deg 0 --max-diff-deg 0",
            recurrence: "P_n = (-1)^n P_{n-1}",
            initial: vec![poly(&[1])],
            next: Box::new(|n, rows| scale(prev(rows, 1), alternating_sign(n))),
        },
        Fixture {
            slug: "11_alternating_fibonacci",
            title: "Alternating Fibonacci",
            features: "alternating sign, second order",
            suggested_args: "--alternating-sign --max-rec-len 2 --max-var-deg 0 --max-idx-deg 0 --max-diff-deg 0",
            recurrence: "P_n = P_{n-1} + (-1)^n P_{n-2}",
            initial: vec![poly(&[1]), poly(&[1])],
            next: Box::new(|n, rows| add(prev(rows, 1), &scale(prev(rows, 2), alternating_sign(n)))),
        },
        Fixture {
            slug: "12_eulerian_derivative",
            title: "Eulerian derivative",
            features: "first derivative, n- and t-dependent coefficients",
            suggested_args: "--max-rec-len 1 --max-var-deg 2 --max-idx-deg 1 --max-diff-deg 1",
            recurrence: "P_n = (1+(n-2)t)P_{n-1} + (t-t^2)P'_{n-1}",
            initial: vec![poly(&[1])],
            next: Box::new(|n, rows| {
                let d = derivative(prev(rows, 1));
                add3(
                    &mul_linear(Rat::one(), rat((n as i64) - 2), prev(rows, 1)),
                    &shift_t(&d, 1),
                    &scale_i(&shift_t(&d, 2), -1),
                )
            }),
        },
        Fixture {
            slug: "13_derivative_appell",
            title: "Derivative Appell-style",
            features: "first derivative, n- and t-dependent coefficient",
            suggested_args: "--max-rec-len 1 --max-var-deg 1 --max-idx-deg 1 --max-diff-deg 1",
            recurrence: "P_n = (n+t)P_{n-1} + P'_{n-1}",
            initial: vec![poly(&[1])],
            next: Box::new(|n, rows| {
                add(
                    &mul_linear(rat_usize(n), Rat::one(), prev(rows, 1)),
                    &derivative(prev(rows, 1)),
                )
            }),
        },
        Fixture {
            slug: "14_second_derivative",
            title: "Second derivative",
            features: "second derivative",
            suggested_args: "--max-rec-len 1 --max-var-deg 2 --max-idx-deg 0 --max-diff-deg 2",
            recurrence: "P_n = (1+t)P_{n-1} + t^2 P''_{n-1}",
            initial: vec![poly(&[1, 1])],
            next: Box::new(|_, rows| {
                add(
                    &mul_linear(Rat::one(), Rat::one(), prev(rows, 1)),
                    &shift_t(&second_derivative(prev(rows, 1)), 2),
                )
            }),
        },
        Fixture {
            slug: "15_mixed_derivative_second_order",
            title: "Mixed derivative second order",
            features: "first derivative, n- and t-dependent coefficient, second order",
            suggested_args: "--max-rec-len 2 --max-var-deg 1 --max-idx-deg 1 --max-diff-deg 1",
            recurrence: "P_n = (1+nt)P_{n-1} + tP'_{n-1} + P_{n-2}",
            initial: vec![poly(&[1]), poly(&[1, 1])],
            next: Box::new(|n, rows| {
                add3(
                    &mul_linear(Rat::one(), rat_usize(n), prev(rows, 1)),
                    &shift_t(&derivative(prev(rows, 1)), 1),
                    prev(rows, 2),
                )
            }),
        },
        Fixture {
            slug: "16_denominator_linear_index",
            title: "Linear index denominator",
            features: "LHS denominator in n",
            suggested_args: "--denominator --max-rec-len 1 --max-var-deg 0 --max-idx-deg 0 --max-diff-deg 0 --max-denom-idx-deg 1",
            recurrence: "(n+1)P_n = P_{n-1}",
            initial: vec![poly(&[1])],
            next: Box::new(|n, rows| scale(prev(rows, 1), Rat::new(BigInt::one(), BigInt::from(n + 1)))),
        },
        Fixture {
            slug: "17_denominator_quadratic_index",
            title: "Quadratic index denominator",
            features: "LHS denominator, quadratic in n",
            suggested_args: "--denominator --max-rec-len 1 --max-var-deg 0 --max-idx-deg 0 --max-diff-deg 0 --max-denom-idx-deg 2",
            recurrence: "(1+n+n^2)P_n = P_{n-1}",
            initial: vec![poly(&[1])],
            next: Box::new(|n, rows| {
                let denom = n * n + n + 1;
                scale(prev(rows, 1), Rat::new(BigInt::one(), BigInt::from(denom)))
            }),
        },
        Fixture {
            slug: "18_denominator_with_t_rhs",
            title: "Denominator with t-dependent RHS",
            features: "LHS denominator, t-dependent RHS",
            suggested_args: "--denominator --max-rec-len 1 --max-var-deg 1 --max-idx-deg 0 --max-diff-deg 0 --max-denom-idx-deg 1",
            recurrence: "(n+1)P_n = (1+t)P_{n-1}",
            initial: vec![poly(&[1])],
            next: Box::new(|n, rows| {
                scale(
                    &mul_linear(Rat::one(), Rat::one(), prev(rows, 1)),
                    Rat::new(BigInt::one(), BigInt::from(n + 1)),
                )
            }),
        },
        Fixture {
            slug: "19_denominator_derivative",
            title: "Denominator derivative",
            features: "LHS denominator, first derivative",
            suggested_args: "--denominator --max-rec-len 1 --max-var-deg 1 --max-idx-deg 0 --max-diff-deg 1 --max-denom-idx-deg 1",
            recurrence: "(n+1)P_n = P_{n-1} + tP'_{n-1}",
            initial: vec![poly(&[1, 1])],
            next: Box::new(|n, rows| {
                scale(
                    &add(prev(rows, 1), &shift_t(&derivative(prev(rows, 1)), 1)),
                    Rat::new(BigInt::one(), BigInt::from(n + 1)),
                )
            }),
        },
        Fixture {
            slug: "20_complex_mixed_alternating_derivative",
            title: "Complex mixed alternating derivative",
            features: "alternating sign, first derivative, n- and t-dependent coefficients, second order",
            suggested_args: "--alternating-sign --max-rec-len 2 --max-var-deg 1 --max-idx-deg 1 --max-diff-deg 1",
            recurrence: "P_n = (n+t)P_{n-1} + (1-t)P'_{n-1} + (-1)^n tP_{n-2}",
            initial: vec![poly(&[1]), poly(&[1, 1])],
            next: Box::new(|n, rows| {
                add3(
                    &mul_linear(rat_usize(n), Rat::one(), prev(rows, 1)),
                    &add(&derivative(prev(rows, 1)), &scale_i(&shift_t(&derivative(prev(rows, 1)), 1), -1)),
                    &scale(&shift_t(prev(rows, 2), 1), alternating_sign(n)),
                )
            }),
        },
        Fixture {
            slug: "21_fifth_order_mixed_coefficients",
            title: "Fifth-order mixed coefficients",
            features: "fifth-order recurrence, n^2- and t^5-dependent coefficients",
            suggested_args: "--max-rec-len 5 --max-var-deg 5 --max-idx-deg 2 --max-diff-deg 0",
            recurrence: "P_n = tP_{n-1} + (t^2-t-n)P_{n-2} + (t^2+tn+1)P_{n-3} + (t^4+t^3n^2+1)P_{n-4} - t^5P_{n-5}",
            initial: vec![
                poly(&[1]),
                poly(&[1, 1]),
                poly(&[2, -1, 1]),
                poly(&[1, 0, 2, 1]),
                poly(&[3, 1, -1, 2, 1]),
            ],
            next: Box::new(|n, rows| {
                let n_rat = rat_usize(n);
                let n_sq = n_rat.clone() * n_rat.clone();
                let term1 = shift_t(prev(rows, 1), 1);
                let term2 = mul_t_poly(&[-n_rat.clone(), rat(-1), Rat::one()], prev(rows, 2));
                let term3 = mul_t_poly(&[Rat::one(), n_rat.clone(), Rat::one()], prev(rows, 3));
                let term4 = mul_t_poly(
                    &[Rat::one(), Rat::zero(), Rat::zero(), n_sq, Rat::one()],
                    prev(rows, 4),
                );
                let term5 = scale_i(&shift_t(prev(rows, 5), 5), -1);
                sum_polys(&[term1, term2, term3, term4, term5])
            }),
        },
        Fixture {
            slug: "22_fifth_derivative_generic",
            title: "Fifth-derivative generic",
            features: "second-order recurrence, derivatives through order five, n- and t-dependent coefficients",
            suggested_args: "--max-rec-len 2 --max-var-deg 2 --max-idx-deg 1 --max-diff-deg 5",
            recurrence: "P_n = (2+t-nt^2)P_{n-1} + (-1+3t)P'_{n-1} + (1+t^2)P''_{n-2} + (4-2n+t)P'''_{n-2} + (-3+t+nt^2)P^{(4)}_{n-2} + (5+t-2nt)P^{(5)}_{n-2}",
            initial: vec![
                poly(&[1, 1, 2, -1, 0, 1]),
                poly(&[2, -1, 3, 0, 1, 2]),
            ],
            next: Box::new(|n, rows| {
                let n_rat = rat_usize(n);
                let term1 = mul_t_poly(&[rat(2), Rat::one(), -n_rat.clone()], prev(rows, 1));
                let term2 = mul_t_poly(&[rat(-1), rat(3)], &derivative(prev(rows, 1)));
                let term3 = mul_t_poly(
                    &[Rat::one(), Rat::zero(), Rat::one()],
                    &nth_derivative(prev(rows, 2), 2),
                );
                let term4 = mul_t_poly(
                    &[rat(4) - rat_usize(2 * n), Rat::one()],
                    &nth_derivative(prev(rows, 2), 3),
                );
                let term5 = mul_t_poly(
                    &[rat(-3), Rat::one(), n_rat],
                    &nth_derivative(prev(rows, 2), 4),
                );
                let term6 = mul_t_poly(
                    &[rat(5), rat(1) - rat_usize(2 * n)],
                    &nth_derivative(prev(rows, 2), 5),
                );
                sum_polys(&[term1, term2, term3, term4, term5, term6])
            }),
        },
        Fixture {
            slug: "23_sparse_second_derivative_lag",
            title: "Sparse second-derivative lag",
            features: "second derivative, third-order recurrence, sparse t^2 coefficients",
            suggested_args: "--max-rec-len 3 --max-var-deg 2 --max-idx-deg 1 --max-diff-deg 2 --fit-extra-rows 2",
            recurrence: "P_n = tP_{n-1} + t^2P''_{n-1} - t^2P_{n-3} + nt^2P''_{n-3}",
            initial: vec![
                poly(&[1, 2, 3, 1]),
                poly(&[2, -1, 1, 4, 1]),
                poly(&[1, 0, 3, -2, 2]),
            ],
            next: Box::new(|n, rows| {
                let term1 = shift_t(prev(rows, 1), 1);
                let term2 = shift_t(&second_derivative(prev(rows, 1)), 2);
                let term3 = scale_i(&shift_t(prev(rows, 3), 2), -1);
                let term4 = scale(
                    &shift_t(&second_derivative(prev(rows, 3)), 2),
                    rat_usize(n),
                );
                sum_polys(&[term1, term2, term3, term4])
            }),
        },
    ]
}

fn write_fixture(base: &Path, fixture: &Fixture, rows: &[Poly]) -> std::io::Result<()> {
    let rows_dir = base.join("rows");
    fs::create_dir_all(&rows_dir)?;
    let mut text = String::new();
    for row in rows {
        text.push_str(&format_row(row));
        text.push('\n');
    }
    fs::write(rows_dir.join(format!("{}.txt", fixture.slug)), text)
}

fn write_manifest(base: &Path, fixtures: &[Fixture]) -> std::io::Result<()> {
    let mut tsv = String::new();
    tsv.push_str("slug\ttitle\tfeatures\tsuggested_args\trecurrence\trows_file\tjson_file\n");
    for fixture in fixtures {
        tsv.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\trows/{}.txt\tjson/{}.json\n",
            fixture.slug,
            fixture.title,
            fixture.features,
            fixture.suggested_args,
            fixture.recurrence,
            fixture.slug,
            fixture.slug
        ));
    }
    fs::write(base.join("manifest.tsv"), tsv)
}

fn write_readme(base: &Path, fixtures: &[Fixture], row_count: usize) -> std::io::Result<()> {
    let mut md = String::new();
    md.push_str("# Recurrence Benchmark Fixtures\n\n");
    md.push_str("Generated exact coefficient rows for adaptive recurrence-search benchmarks.\n");
    md.push_str(&format!(
        "Each `rows/*.txt` file contains {row_count} polynomials, one dense coefficient list per line,\n"
    ));
    md.push_str("with coefficients in ascending powers of `t`. These files contain no headers,\n");
    md.push_str("so they can be piped directly into `polytool recurrence`.\n");
    md.push_str("Metadata lives separately in `manifest.tsv` and in the table below.\n");
    md.push_str("The matching `json/*.json` files are recurrence JSON records emitted by\n");
    md.push_str("`polytool recurrence --json`; they include minimal initial conditions and\n");
    md.push_str("can regenerate or extend the raw row files with `recurrence-generate`.\n\n");
    md.push_str("Regenerate from the Rust workspace root with:\n\n");
    md.push_str("```sh\n");
    md.push_str("cargo run -p polytool --example generate_recurrence_benchmarks\n");
    md.push_str("bash polytool/fixtures/recurrence-benchmarks/regenerate-json.sh\n");
    md.push_str("```\n\n");
    md.push_str("Example timing command:\n\n");
    md.push_str("```sh\n");
    md.push_str("time cargo run -q -p polytool --bin polytool -- recurrence \\\n");
    md.push_str("  --max-rec-len 1 --max-var-deg 1 --max-idx-deg 0 --max-diff-deg 0 \\\n");
    md.push_str("  < polytool/fixtures/recurrence-benchmarks/rows/03_binomial_powers.txt\n");
    md.push_str("```\n\n");
    md.push_str("Example regeneration command:\n\n");
    md.push_str("```sh\n");
    md.push_str("polytool recurrence-generate \\\n");
    md.push_str(
        "  --recurrence polytool/fixtures/recurrence-benchmarks/json/03_binomial_powers.json \\\n",
    );
    md.push_str("  --rows 50\n");
    md.push_str("```\n\n");
    md.push_str("| slug | features | suggested args | recurrence |\n");
    md.push_str("|---|---|---|---|\n");
    for fixture in fixtures {
        md.push_str(&format!(
            "| `{}` | {} | `{}` | `{}` |\n",
            fixture.slug, fixture.features, fixture.suggested_args, fixture.recurrence
        ));
    }
    fs::write(base.join("README.md"), md)
}

fn main() -> std::io::Result<()> {
    let base = Path::new("polytool/fixtures/recurrence-benchmarks");
    fs::create_dir_all(base.join("rows"))?;

    let fixtures = fixtures();
    let row_count = 50;
    for fixture in &fixtures {
        let rows = generate_rows(fixture, row_count);
        write_fixture(base, fixture, &rows)?;
    }
    write_manifest(base, &fixtures)?;
    write_readme(base, &fixtures, row_count)?;

    println!(
        "wrote {} recurrence benchmark fixtures to {}",
        fixtures.len(),
        base.display()
    );
    Ok(())
}
