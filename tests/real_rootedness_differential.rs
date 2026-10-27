//! Randomized differential and ground-truth checks of every exact
//! real-rootedness, root-counting and interlacing method (about 64,000
//! polynomials and 20,000 interlacing pairs). They take about a minute in a
//! release build (15 in debug), so they are ignored by default:
//!
//! ```text
//! cargo test --release -p polytool --test real_rootedness_differential -- --ignored
//! ```

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{Signed, ToPrimitive, Zero};
use polytool::real_rootedness::{isolate_real_roots_rational_coeffs, RealInterval};
use polytool::root_count::{isolate_real_roots_bigint_coeffs, IsolatedRealRoot};
use polytool::sturm::SturmChain;
use polytool::*;
// Shadows two trait names that `polytool::*` would otherwise import
// ambiguously next to the explicit imports above.
#[allow(dead_code, unused_imports)]
mod hide {
    pub trait CoeffRing {}
    pub trait FieldRing {}
}
#[allow(unused_imports)]
use hide::{CoeffRing, FieldRing};
use std::collections::BTreeSet;

type Q = BigRational;
/// A boxed entry-point call checked for panics.
type Probe = Box<dyn Fn() + std::panic::RefUnwindSafe>;

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn range(&mut self, lo: i64, hi: i64) -> i64 {
        lo + (self.next() % ((hi - lo + 1) as u64)) as i64
    }
    fn chance(&mut self, num: u64, den: u64) -> bool {
        self.next() % den < num
    }
}

fn big(v: &[i64]) -> Vec<BigInt> {
    v.iter().map(|&c| BigInt::from(c)).collect()
}

fn mul(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    let mut r = vec![BigInt::from(0); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            r[i + j] += x * y;
        }
    }
    r
}

fn trim(mut p: Vec<BigInt>) -> Vec<BigInt> {
    while p.last().is_some_and(|c| c.is_zero()) {
        p.pop();
    }
    p
}

fn eval_q(p: &[BigInt], x: &Q) -> Q {
    p.iter()
        .rev()
        .fold(Q::from_integer(BigInt::from(0)), |acc, c| {
            acc * x + Q::from_integer(c.clone())
        })
}

fn sign_q(x: &Q) -> i8 {
    if x.is_positive() {
        1
    } else if x.is_negative() {
        -1
    } else {
        0
    }
}

fn is_square(n: i64) -> bool {
    if n < 0 {
        return false;
    }
    let r = (n as f64).sqrt().round() as i64;
    (r - 1..=r + 1).any(|s| s >= 0 && s * s == n)
}

#[derive(Clone, Debug)]
enum Quad {
    Complex(i64, i64, i64),
    RealIrr(i64, i64, i64),
}

/// Ground truth from an explicit factorization.
#[derive(Debug)]
struct Truth {
    coeffs: Vec<BigInt>,
    real_rooted: bool,
    distinct_real: usize,
    distinct_positive: usize,
    real_with_multiplicity: usize,
    rational_roots: BTreeSet<Q>,
    has_real_irrational: bool,
    hurwitz_stable: bool,
}

fn primitive_positive_quad(a: i64, b: i64, c: i64) -> (i64, i64, i64) {
    fn gcd(a: i64, b: i64) -> i64 {
        if b == 0 {
            a.abs()
        } else {
            gcd(b, a % b)
        }
    }
    let g = gcd(gcd(a, b), c);
    let (mut a, mut b, mut c) = (a / g, b / g, c / g);
    if a < 0 {
        a = -a;
        b = -b;
        c = -c;
    }
    (a, b, c)
}

fn random_factored(rng: &mut Rng, max_linear: i64, max_quads: i64) -> Truth {
    let mut coeffs = vec![BigInt::from(
        rng.range(1, 3) * if rng.chance(1, 2) { 1 } else { -1 },
    )];
    let mut rational_roots = BTreeSet::new();
    let mut real_with_multiplicity = 0usize;
    let mut real_quads: BTreeSet<(i64, i64, i64)> = BTreeSet::new();
    let mut real_rooted = true;
    let mut hurwitz_stable = true;
    let mut linear: Vec<(i64, i64)> = Vec::new();
    for _ in 0..rng.range(0, max_linear) {
        if !linear.is_empty() && rng.chance(1, 3) {
            let f = linear[(rng.next() % linear.len() as u64) as usize];
            linear.push(f);
        } else {
            let a = rng.range(1, 3) * if rng.chance(1, 2) { 1 } else { -1 };
            let b = rng.range(-4, 4);
            linear.push((a, b));
        }
    }
    for &(a, b) in &linear {
        coeffs = mul(&coeffs, &big(&[b, a]));
        let root = Q::new(BigInt::from(-b), BigInt::from(a));
        if !root.is_negative() {
            hurwitz_stable = false;
        }
        rational_roots.insert(root);
        real_with_multiplicity += 1;
    }
    let mut quads = Vec::new();
    for _ in 0..rng.range(0, max_quads) {
        loop {
            let a = rng.range(1, 3) * if rng.chance(1, 2) { 1 } else { -1 };
            let b = rng.range(-5, 5);
            let c = rng.range(-5, 5);
            if c == 0 {
                continue;
            }
            let disc = b * b - 4 * a * c;
            if disc < 0 {
                quads.push(Quad::Complex(a, b, c));
                break;
            }
            if disc > 0 && !is_square(disc) {
                quads.push(Quad::RealIrr(a, b, c));
                break;
            }
        }
    }
    let mut extra_positive = 0usize;
    for q in &quads {
        let (a, b, c) = match q {
            Quad::Complex(a, b, c) | Quad::RealIrr(a, b, c) => (*a, *b, *c),
        };
        coeffs = mul(&coeffs, &big(&[c, b, a]));
        let (pa, pb, pc) = primitive_positive_quad(a, b, c);
        if !(pb > 0 && pc > 0) {
            hurwitz_stable = false;
        }
        match q {
            Quad::Complex(..) => real_rooted = false,
            Quad::RealIrr(..) => {
                real_with_multiplicity += 2;
                if real_quads.insert((pa, pb, pc)) {
                    // roots positive count
                    let prod_pos = pc > 0; // c/a with a>0
                    let sum_pos = -pb > 0;
                    extra_positive += if !prod_pos {
                        1
                    } else if sum_pos {
                        2
                    } else {
                        0
                    };
                }
            }
        }
    }
    let distinct_real = rational_roots.len() + 2 * real_quads.len();
    let distinct_positive =
        rational_roots.iter().filter(|r| r.is_positive()).count() + extra_positive;
    let coeffs = trim(coeffs);
    Truth {
        coeffs,
        real_rooted,
        distinct_real,
        distinct_positive,
        real_with_multiplicity,
        rational_roots,
        has_real_irrational: !real_quads.is_empty(),
        hurwitz_stable,
    }
}

fn with_trailing_zeros(rng: &mut Rng, mut p: Vec<BigInt>) -> Vec<BigInt> {
    for _ in 0..rng.range(0, 2) {
        p.push(BigInt::from(0));
    }
    p
}

fn fits_i64(p: &[BigInt]) -> Option<Vec<i64>> {
    p.iter().map(|c| c.to_i64()).collect()
}

/// All real-rootedness verdicts available for one polynomial.
fn all_verdicts(p: &[BigInt]) -> Vec<(&'static str, bool)> {
    let mut v = vec![
        ("fast", is_real_rooted_fast_bigint_coeffs(p)),
        ("default", is_real_rooted_bigint_coeffs(p)),
        ("prs", is_real_rooted_prs_bigint_coeffs(p)),
        ("uspensky", is_real_rooted_uspensky_bigint_coeffs(p)),
        ("sturm_rational", is_real_rooted_sturm_bigint_coeffs(p)),
        ("bezout_psd", is_real_rooted_bezout_bigint_coeffs(p)),
        (
            "bezout_squarefree",
            is_real_rooted_bezout_squarefree_bigint_coeffs(p),
        ),
        ("hermite", is_real_rooted_hermite_bigint_coeffs(p)),
        (
            "enclosure",
            isolate_real_roots_rational_coeffs(
                &p.iter()
                    .map(|c| Q::from_integer(c.clone()))
                    .collect::<Vec<_>>(),
            )
            .is_real_rooted(),
        ),
    ];
    if let Some(rr) = is_real_rooted_one_signed_bigint_coeffs(p) {
        v.push(("one_signed", rr));
    }
    if let Some(small) = fits_i64(p) {
        v.push(("i64_default", is_real_rooted(&small)));
        v.push(("i64_uspensky", is_real_rooted_uspensky_i64(&small)));
        v.push(("i64_bezout", is_real_rooted_bezout(&small)));
        v.push(("i64_bezout_sf", is_real_rooted_bezout_squarefree(&small)));
        v.push(("i64_hermite", is_real_rooted_hermite(&small)));
        v.push(("i64_sturm", is_real_rooted_sturm(&small)));
    }
    v
}

fn all_real_counts(p: &[BigInt]) -> Vec<(&'static str, usize)> {
    let iso = isolate_real_roots_rational_coeffs(
        &p.iter()
            .map(|c| Q::from_integer(c.clone()))
            .collect::<Vec<_>>(),
    );
    vec![
        ("prs", count_real_roots_prs_bigint_coeffs(p)),
        ("uspensky", count_real_roots_uspensky_bigint_coeffs(p)),
        (
            "sturm_rational",
            SturmChain::from_bigint_coeffs(p).count_real_roots(),
        ),
        (
            "dyadic_isolation",
            isolate_real_roots_bigint_coeffs(p).len(),
        ),
        ("enclosure", iso.roots.len()),
    ]
}

fn check_dyadic_isolation(p: &[BigInt], truth: &Truth, failures: &mut Vec<String>) {
    let iso = isolate_real_roots_bigint_coeffs(p);
    let to_q =
        |d: &polytool::root_count::DyadicRational| Q::new(d.numerator().clone(), d.denominator());
    for w in iso.windows(2) {
        let (a, b) = (to_q(w[0].upper()), to_q(w[1].lower()));
        let strict = matches!(w[0], IsolatedRealRoot::Exact(_))
            || matches!(w[1], IsolatedRealRoot::Exact(_));
        if a > b || (strict && a == b) {
            failures.push(format!("isolation order/overlap p={p:?}"));
        }
    }
    for r in &truth.rational_roots {
        let hits = iso
            .iter()
            .filter(|e| match e {
                IsolatedRealRoot::Exact(x) => &to_q(x) == r,
                IsolatedRealRoot::Interval { lower, upper } => &to_q(lower) < r && r < &to_q(upper),
            })
            .count();
        if hits != 1 {
            failures.push(format!(
                "dyadic isolation: root {r} in {hits} enclosures, p={p:?}"
            ));
        }
    }
    for e in &iso {
        if let IsolatedRealRoot::Interval { lower, upper } = e {
            if eval_q(p, &to_q(lower)).is_zero() || eval_q(p, &to_q(upper)).is_zero() {
                failures.push(format!("dyadic isolation endpoint is a root p={p:?}"));
            }
            // exactly one root of the squarefree part inside: check sign change for odd-
            // multiplicity cases is not reliable; rely on count comparisons instead.
        } else if let IsolatedRealRoot::Exact(x) = e {
            if !eval_q(p, &to_q(x)).is_zero() {
                failures.push(format!("dyadic exact root is not a root p={p:?}"));
            }
        }
    }
    // Enclosures (rational API) with multiplicities.
    let qp: Vec<Q> = p.iter().map(|c| Q::from_integer(c.clone())).collect();
    let enc = isolate_real_roots_rational_coeffs(&qp);
    if enc.real_root_count != truth.real_with_multiplicity {
        failures.push(format!(
            "enclosure multiplicity count {} != {} p={p:?}",
            enc.real_root_count, truth.real_with_multiplicity
        ));
    }
    for r in &truth.rational_roots {
        let hits = enc
            .roots
            .iter()
            .filter(|e| enc.root_in(e, &RealInterval::point(r.clone())))
            .count();
        if hits != 1 {
            failures.push(format!("enclosure: root {r} matched {hits} times p={p:?}"));
        }
    }
}

fn report(name: &str, tried: usize, failures: &[String]) {
    eprintln!("[{name}] tried {tried}, failures {}", failures.len());
    for f in failures.iter().take(25) {
        eprintln!("  {f}");
    }
}

#[test]
#[ignore = "differential fuzz, about a minute in release; run with --ignored"]
fn probe_factored_ground_truth() {
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    let mut failures = Vec::new();
    let mut tried = 0;
    for iter in 0..6000 {
        let (ml, mq) = if iter % 3 == 0 { (9, 2) } else { (6, 2) };
        let truth = random_factored(&mut rng, ml, mq);
        let p = with_trailing_zeros(&mut rng, truth.coeffs.clone());
        tried += 1;
        for (name, v) in all_verdicts(&p) {
            if v != truth.real_rooted {
                failures.push(format!(
                    "rr {name}: got {v}, expected {} p={:?}",
                    truth.real_rooted, truth.coeffs
                ));
            }
        }
        for (name, c) in all_real_counts(&p) {
            if c != truth.distinct_real {
                failures.push(format!(
                    "count {name}: got {c}, expected {} p={:?}",
                    truth.distinct_real, truth.coeffs
                ));
            }
        }
        let pos = [
            ("prs", count_positive_roots_prs_bigint_coeffs(&p)),
            ("uspensky", count_positive_roots_uspensky_bigint_coeffs(&p)),
        ];
        for (name, c) in pos {
            if c != truth.distinct_positive {
                failures.push(format!(
                    "positive {name}: got {c}, expected {} p={:?}",
                    truth.distinct_positive, truth.coeffs
                ));
            }
        }
        let sf1 = squarefree_degree_bigint_coeffs(&p);
        let sf2 = SturmChain::from_bigint_coeffs(&p).square_free_degree();
        let sf3 = vec_poly::degree(&squarefree_part_bigint_coeffs(&p)).unwrap_or(0);
        if sf1 != sf2 || sf2 != sf3 {
            failures.push(format!(
                "squarefree degree mismatch {sf1} {sf2} {sf3} p={p:?}"
            ));
        }
        // Hurwitz.
        let hs = is_hurwitz_stable_bigint(&p);
        if hs != truth.hurwitz_stable {
            failures.push(format!(
                "hurwitz: got {hs}, expected {} p={:?}",
                truth.hurwitz_stable, truth.coeffs
            ));
        }
        // Discriminant zero iff repeated root.
        let has_repeat = sf1 != vec_poly::degree(&p).unwrap_or(0);
        let disc = discriminant_bigint_coeffs(&p);
        if vec_poly::degree(&p).unwrap_or(0) >= 1 && disc.is_zero() != has_repeat {
            failures.push(format!(
                "discriminant zero={} repeat={has_repeat} p={p:?}",
                disc.is_zero()
            ));
        }
        // Isolation and enclosure checks.
        check_dyadic_isolation(&p, &truth, &mut failures);
        // Root sign counts when all real roots are rational.
        if !truth.has_real_irrational {
            let test = big(&[rng.range(-3, 3), rng.range(-3, 3), rng.range(-2, 2)]);
            let mut expected = (0usize, 0usize, 0usize);
            for r in &truth.rational_roots {
                match sign_q(&eval_q(&test, r)) {
                    1 => expected.0 += 1,
                    -1 => expected.1 += 1,
                    _ => expected.2 += 1,
                }
            }
            let got = count_root_signs_prs_bigint_coeffs(&p, &test);
            if (got.positive, got.negative, got.zero) != expected
                || got.total != truth.rational_roots.len()
            {
                failures.push(format!(
                    "root signs got {got:?} expected {expected:?} p={p:?} q={test:?}"
                ));
            }
            let taq = tarski_query_prs_bigint_coeffs(&p, &test);
            if taq != expected.0 as isize - expected.1 as isize {
                failures.push(format!("tarski got {taq} p={p:?} q={test:?}"));
            }
        }
    }
    report("factored", tried, &failures);
    assert!(failures.is_empty());
}

#[test]
#[ignore = "differential fuzz, about a minute in release; run with --ignored"]
fn probe_exhaustive_small_coefficients() {
    let mut failures = Vec::new();
    let mut tried = 0usize;
    for degree in 0usize..=8 {
        let total = 3usize.pow((degree + 1) as u32);
        for mut mask in 0..total {
            let mut coeffs = Vec::with_capacity(degree + 1);
            for _ in 0..=degree {
                coeffs.push((mask % 3) as i64 - 1);
                mask /= 3;
            }
            let p = big(&coeffs);
            tried += 1;
            let verdicts = all_verdicts(&p);
            let reference = verdicts[0].1;
            if verdicts.iter().any(|(_, v)| *v != reference) {
                failures.push(format!("rr disagreement {coeffs:?}: {verdicts:?}"));
            }
            let counts = all_real_counts(&p);
            if counts.iter().any(|(_, c)| *c != counts[0].1) {
                failures.push(format!("count disagreement {coeffs:?}: {counts:?}"));
            }
            let a = count_positive_roots_prs_bigint_coeffs(&p);
            let b = count_positive_roots_uspensky_bigint_coeffs(&p);
            if a != b {
                failures.push(format!("positive disagreement {coeffs:?}: {a} {b}"));
            }
        }
    }
    report("exhaustive {-1,0,1} deg<=8", tried, &failures);
    assert!(failures.is_empty());
}

#[test]
#[ignore = "differential fuzz, about a minute in release; run with --ignored"]
fn probe_random_wide_coefficients() {
    let mut rng = Rng(0xDEAD_BEEF_1234_5678);
    let mut failures = Vec::new();
    let mut tried = 0usize;
    for _ in 0..8000 {
        let degree = rng.range(0, 11) as usize;
        let bound = [3, 10, 1000, 1_000_000_000][(rng.next() % 4) as usize];
        let mut coeffs: Vec<i64> = (0..=degree).map(|_| rng.range(-bound, bound)).collect();
        if rng.chance(1, 4) {
            // one-signed
            for c in &mut coeffs {
                *c = c.abs();
            }
        }
        if rng.chance(1, 5) {
            // palindromic
            let n = coeffs.len();
            for i in 0..n / 2 {
                coeffs[n - 1 - i] = coeffs[i];
            }
        }
        let p = big(&coeffs);
        tried += 1;
        let verdicts = all_verdicts(&p);
        let reference = verdicts[0].1;
        if verdicts.iter().any(|(_, v)| *v != reference) {
            failures.push(format!("rr disagreement {coeffs:?}: {verdicts:?}"));
        }
        let counts = all_real_counts(&p);
        if counts.iter().any(|(_, c)| *c != counts[0].1) {
            failures.push(format!("count disagreement {coeffs:?}: {counts:?}"));
        }
        let a = count_positive_roots_prs_bigint_coeffs(&p);
        let b = count_positive_roots_uspensky_bigint_coeffs(&p);
        if a != b {
            failures.push(format!("positive disagreement {coeffs:?}: {a} {b}"));
        }
    }
    report("random wide", tried, &failures);
    assert!(failures.is_empty());
}

/// Square and products of real-rooted families with huge, close roots.
#[test]
#[ignore = "differential fuzz, about a minute in release; run with --ignored"]
fn probe_close_large_roots() {
    let mut rng = Rng(0x1234_5678_9ABC_DEF1);
    let mut failures = Vec::new();
    let mut tried = 0usize;
    for _ in 0..600 {
        let center = BigInt::from(rng.range(-1_000_000_000, 1_000_000_000))
            * BigInt::from(rng.range(1, 1_000_000));
        let k = rng.range(1, 6);
        let mut p = big(&[1]);
        let mut roots = BTreeSet::new();
        for _ in 0..k {
            let r = &center + BigInt::from(rng.range(-3, 3));
            roots.insert(r.clone());
            p = mul(&p, &[-r, BigInt::from(1)]);
        }
        let complex = rng.chance(1, 3);
        if complex {
            // (t - center)^2 + 1, nearly real
            p = mul(
                &p,
                &[
                    &center * &center + 1,
                    -BigInt::from(2) * &center,
                    BigInt::from(1),
                ],
            );
        }
        tried += 1;
        for (name, v) in all_verdicts(&p) {
            if v == complex {
                failures.push(format!("close rr {name}: got {v} p={p:?}"));
            }
        }
        for (name, c) in all_real_counts(&p) {
            if c != roots.len() {
                failures.push(format!(
                    "close count {name}: got {c} expected {} p={p:?}",
                    roots.len()
                ));
            }
        }
    }
    report("close large roots", tried, &failures);
    assert!(failures.is_empty());
}

// ---------------------------------------------------------------- interlacing

fn poly_from_roots(roots: &[Q], scale: i64) -> Vec<BigInt> {
    // multiply by denominators: (den t - num)
    let mut p = big(&[scale]);
    for r in roots {
        p = mul(&p, &[-r.numer().clone(), r.denom().clone()]);
    }
    p
}

fn directed_oracle(p_roots: &[Q], q_roots: &[Q]) -> Option<bool> {
    let mut a = p_roots.to_vec();
    let mut b = q_roots.to_vec();
    a.sort_by(|x, y| y.cmp(x));
    b.sort_by(|x, y| y.cmp(x));
    if b.len() != a.len() && b.len() != a.len() + 1 {
        return None;
    }
    let mut merged = Vec::new();
    for i in 0..b.len() {
        merged.push(b[i].clone());
        if i < a.len() {
            merged.push(a[i].clone());
        }
    }
    Some(merged.windows(2).all(|w| w[0] >= w[1]))
}

fn strict_oracle(p_roots: &[Q], q_roots: &[Q]) -> Option<bool> {
    let mut a = p_roots.to_vec();
    let mut b = q_roots.to_vec();
    a.sort_by(|x, y| y.cmp(x));
    b.sort_by(|x, y| y.cmp(x));
    if b.len() != a.len() + 1 {
        return None;
    }
    let mut merged = Vec::new();
    for i in 0..b.len() {
        merged.push(b[i].clone());
        if i < a.len() {
            merged.push(a[i].clone());
        }
    }
    Some(merged.windows(2).all(|w| w[0] > w[1]))
}

fn random_roots(rng: &mut Rng, n: usize) -> Vec<Q> {
    (0..n)
        .map(|_| {
            let den = if rng.chance(1, 4) { 2 } else { 1 };
            Q::new(BigInt::from(rng.range(-4, 4)), BigInt::from(den))
        })
        .collect()
}

#[test]
#[ignore = "differential fuzz, about a minute in release; run with --ignored"]
fn probe_interlacing_against_root_oracle() {
    let mut rng = Rng(0x0BAD_C0DE_F00D_1234);
    let mut failures = Vec::new();
    let mut tried = 0usize;
    for _ in 0..20000 {
        let m = rng.range(0, 5) as usize;
        let n = m + rng.range(0, 1) as usize;
        let mut a = random_roots(&mut rng, m);
        let mut b = random_roots(&mut rng, n);
        // Encourage shared roots.
        if !a.is_empty() && !b.is_empty() && rng.chance(1, 2) {
            let i = (rng.next() % a.len() as u64) as usize;
            let j = (rng.next() % b.len() as u64) as usize;
            b[j] = a[i].clone();
        }
        if rng.chance(1, 4) {
            a.sort();
            b.sort();
        }
        let sa = rng.range(1, 3) * if rng.chance(1, 2) { 1 } else { -1 };
        let sb = rng.range(1, 3) * if rng.chance(1, 2) { 1 } else { -1 };
        let mut p = poly_from_roots(&a, sa);
        let mut q = poly_from_roots(&b, sb);
        // Sometimes multiply in a complex factor (to one or both).
        let complex_mode = rng.range(0, 6);
        let cfac = big(&[rng.range(1, 3), rng.range(-1, 1), 1]); // disc < 0 since c>=1, b in [-1,1]
        let cfac = if cfac[1].abs() * cfac[1].abs() < BigInt::from(4) * &cfac[0] {
            cfac
        } else {
            big(&[1, 0, 1])
        };
        let (expected_weak, expected_strict) = match complex_mode {
            0 => {
                // complex factor in both: weak/strict both false (if degrees allow)
                p = mul(&p, &cfac);
                q = mul(&q, &cfac);
                (
                    directed_oracle(&a, &b).map(|_| false),
                    strict_oracle(&a, &b).map(|_| false),
                )
            }
            1 => {
                // complex factor only in p; q gets two real roots so degrees stay aligned
                p = mul(&p, &cfac);
                q = mul(&q, &big(&[-1, 0, 1]));
                (
                    directed_oracle(&a, &b).map(|_| false),
                    strict_oracle(&a, &b).map(|_| false),
                )
            }
            2 => {
                // complex factor only in q
                q = mul(&q, &cfac);
                p = mul(&p, &big(&[-1, 0, 1]));
                (
                    directed_oracle(&a, &b).map(|_| false),
                    strict_oracle(&a, &b).map(|_| false),
                )
            }
            _ => (directed_oracle(&a, &b), strict_oracle(&a, &b)),
        };
        if rng.chance(1, 4) {
            p.push(BigInt::from(0));
        }
        tried += 1;
        let weak = check_weak_interlacing_bigint_coeffs(&p, &q);
        if weak != expected_weak {
            failures.push(format!(
                "weak p={p:?} q={q:?} roots a={a:?} b={b:?}: got {weak:?} expected {expected_weak:?}"
            ));
        }
        if let (Some(pi), Some(qi)) = (fits_i64(&p), fits_i64(&q)) {
            let w2 = check_weak_interlacing(&pi, &qi);
            if w2 != weak {
                failures.push(format!("weak i64 vs bigint p={pi:?} q={qi:?}"));
            }
        }
        // Strict (degree difference exactly one).
        let strict = check_interlacing_bigint_coeffs(&p, &q);
        let exp_strict = if vec_poly::degree(&p).is_none() {
            Some(is_real_rooted_bigint_coeffs(&q))
        } else {
            expected_strict
        };
        if strict != exp_strict {
            failures.push(format!(
                "strict p={p:?} q={q:?}: got {strict:?} expected {exp_strict:?}"
            ));
        }
        // Sturm strict fallback (None if not real-rooted, Some(false) on wrong degree).
        if let (Some(pi), Some(qi)) = (fits_i64(&p), fits_i64(&q)) {
            let both_rr = is_real_rooted(&pi) && is_real_rooted(&qi);
            let sturm = check_interlacing_sturm(&pi, &qi);
            let exp = if !both_rr {
                None
            } else {
                Some(exp_strict.unwrap_or(false) && vec_poly::degree(&pi).is_some())
            };
            if sturm != exp && vec_poly::degree(&pi).is_some() {
                failures.push(format!(
                    "sturm strict p={pi:?} q={qi:?}: got {sturm:?} expected {exp:?}"
                ));
            }
        }
        // Undirected consistency.
        let either = interlaces_either_order_bigint_coeffs(&p, &q);
        let either_rev = interlaces_either_order_bigint_coeffs(&q, &p);
        if either != either_rev {
            failures.push(format!("either not symmetric p={p:?} q={q:?}"));
        }
    }
    report("interlacing oracle", tried, &failures);
    assert!(failures.is_empty());
}

// ---------------------------------------------------------------- i64 panics

#[test]
#[ignore = "differential fuzz, about a minute in release; run with --ignored"]
fn probe_i64_extreme_inputs_do_not_panic() {
    let mut rng = Rng(0x5555_AAAA_3333_CCCC);
    let pool = [
        i64::MIN,
        i64::MIN + 1,
        i64::MAX,
        -1,
        0,
        1,
        2,
        -2,
        1 << 40,
        -(1 << 40),
        3_037_000_499,
        -3_037_000_500,
    ];
    let mut failures = Vec::new();
    let mut tried = 0usize;
    let pick = |rng: &mut Rng| pool[(rng.next() % pool.len() as u64) as usize];
    let prev = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));
    for _ in 0..300 {
        let dp = rng.range(0, 4) as usize;
        let p: Vec<i64> = (0..=dp).map(|_| pick(&mut rng)).collect();
        let dq = rng.range(0, 4) as usize;
        let q: Vec<i64> = (0..=dq).map(|_| pick(&mut rng)).collect();
        tried += 1;
        let calls: Vec<(&str, Probe)> = vec![
            (
                "is_real_rooted",
                Box::new({
                    let p = p.clone();
                    move || {
                        is_real_rooted(&p);
                    }
                }),
            ),
            (
                "is_real_rooted_bezout",
                Box::new({
                    let p = p.clone();
                    move || {
                        is_real_rooted_bezout(&p);
                    }
                }),
            ),
            (
                "is_real_rooted_bezout_squarefree",
                Box::new({
                    let p = p.clone();
                    move || {
                        is_real_rooted_bezout_squarefree(&p);
                    }
                }),
            ),
            (
                "is_real_rooted_hermite",
                Box::new({
                    let p = p.clone();
                    move || {
                        is_real_rooted_hermite(&p);
                    }
                }),
            ),
            (
                "is_strictly_real_rooted_bezout",
                Box::new({
                    let p = p.clone();
                    move || {
                        is_strictly_real_rooted_bezout(&p);
                    }
                }),
            ),
            (
                "is_real_rooted_sturm",
                Box::new({
                    let p = p.clone();
                    move || {
                        is_real_rooted_sturm(&p);
                    }
                }),
            ),
            (
                "is_real_rooted_uspensky_i64",
                Box::new({
                    let p = p.clone();
                    move || {
                        is_real_rooted_uspensky_i64(&p);
                    }
                }),
            ),
            (
                "has_simple_roots",
                Box::new({
                    let p = p.clone();
                    move || {
                        has_simple_roots(&p);
                    }
                }),
            ),
            (
                "discriminant",
                Box::new({
                    let p = p.clone();
                    move || {
                        discriminant(&p);
                    }
                }),
            ),
            (
                "resultant",
                Box::new({
                    let (p, q) = (p.clone(), q.clone());
                    move || {
                        resultant(&p, &q);
                    }
                }),
            ),
            (
                "bezout_matrix",
                Box::new({
                    let (p, q) = (p.clone(), q.clone());
                    move || {
                        bezout_matrix(&p, &q);
                    }
                }),
            ),
            (
                "check_interlacing",
                Box::new({
                    let (p, q) = (p.clone(), q.clone());
                    move || {
                        check_interlacing(&p, &q);
                    }
                }),
            ),
            (
                "check_weak_interlacing",
                Box::new({
                    let (p, q) = (p.clone(), q.clone());
                    move || {
                        check_weak_interlacing(&p, &q);
                    }
                }),
            ),
            (
                "interlaces_either_order",
                Box::new({
                    let (p, q) = (p.clone(), q.clone());
                    move || {
                        interlaces_either_order(&p, &q);
                    }
                }),
            ),
            (
                "check_interlacing_sturm",
                Box::new({
                    let (p, q) = (p.clone(), q.clone());
                    move || {
                        check_interlacing_sturm(&p, &q);
                    }
                }),
            ),
            (
                "approximate_real_roots",
                Box::new({
                    let p = p.clone();
                    move || {
                        approximate_real_roots(&p);
                    }
                }),
            ),
            (
                "is_log_concave",
                Box::new({
                    let p = p.clone();
                    move || {
                        is_log_concave(&p);
                    }
                }),
            ),
        ];
        for (name, f) in calls {
            if std::panic::catch_unwind(&f).is_err() {
                failures.push(format!("{name} panicked on p={p:?} q={q:?}"));
            }
        }
    }
    std::panic::set_hook(prev);
    report("i64 extremes", tried, &failures);
    assert!(failures.is_empty());
}

// ---------------------------------------------------------------- high degree

/// High multiplicities and large degree: abnormal (degree-jumping) PRS chains
/// and the modular positive-definiteness path (Bézout dimension >= 30).
#[test]
#[ignore = "differential fuzz, about a minute in release; run with --ignored"]
fn probe_high_degree_and_multiplicities() {
    let mut rng = Rng(0x7777_1111_2222_3333);
    let mut failures = Vec::new();
    let mut tried = 0usize;
    for _ in 0..120 {
        let mut p = big(&[rng.range(1, 2) * if rng.chance(1, 2) { 1 } else { -1 }]);
        let mut roots = BTreeSet::new();
        let nfactors = rng.range(1, 6);
        for _ in 0..nfactors {
            let r = rng.range(-6, 6);
            let mult = rng.range(1, 6);
            roots.insert(r);
            for _ in 0..mult {
                p = mul(&p, &big(&[-r, 1]));
            }
        }
        let complex = rng.chance(1, 3);
        if complex {
            let m = rng.range(1, 3);
            for _ in 0..m {
                p = mul(&p, &big(&[2, 2, 1]));
            }
        }
        tried += 1;
        for (name, v) in all_verdicts(&p) {
            if v == complex {
                failures.push(format!("hd rr {name}: got {v} p={p:?}"));
            }
        }
        for (name, c) in all_real_counts(&p) {
            if c != roots.len() {
                failures.push(format!(
                    "hd count {name}: got {c} exp {} p={p:?}",
                    roots.len()
                ));
            }
        }
        let pos = roots.iter().filter(|r| **r > 0).count();
        let a = count_positive_roots_prs_bigint_coeffs(&p);
        let b = count_positive_roots_uspensky_bigint_coeffs(&p);
        if a != pos || b != pos {
            failures.push(format!("hd positive {a} {b} exp {pos} p={p:?}"));
        }
    }
    // Interlacing with degree 30..34 (modular PD path).
    for _ in 0..40 {
        let n = rng.range(30, 34) as usize;
        let mut all: Vec<i64> = (0..(2 * n + 1) as i64).map(|i| 3 * i - 50).collect();
        // perturb to create violations sometimes
        let violate = rng.chance(1, 2);
        let mut g_roots: Vec<Q> = Vec::new();
        let mut f_roots: Vec<Q> = Vec::new();
        for (i, r) in all.iter_mut().enumerate() {
            let q = Q::from_integer(BigInt::from(*r));
            if i % 2 == 0 {
                g_roots.push(q);
            } else {
                f_roots.push(q);
            }
        }
        if violate {
            let i = (rng.next() % f_roots.len() as u64) as usize;
            f_roots[i] = f_roots[i].clone() + Q::from_integer(BigInt::from(rng.range(2, 9)));
        }
        let f = poly_from_roots(&f_roots, 1);
        let g = poly_from_roots(&g_roots, if rng.chance(1, 2) { 1 } else { -2 });
        tried += 1;
        let exp_strict = strict_oracle(&f_roots, &g_roots);
        let exp_weak = directed_oracle(&f_roots, &g_roots);
        let s = check_interlacing_bigint_coeffs(&f, &g);
        let w = check_weak_interlacing_bigint_coeffs(&f, &g);
        if s != exp_strict {
            failures.push(format!("hd strict got {s:?} exp {exp_strict:?} n={n}"));
        }
        if w != exp_weak {
            failures.push(format!("hd weak got {w:?} exp {exp_weak:?} n={n}"));
        }
        let srr = is_strictly_real_rooted_bezout_bigint_coeffs(&g);
        let hrr = is_strictly_real_rooted_hermite_bigint_coeffs(&g);
        if !srr || !hrr {
            failures.push(format!("hd strict rr bezout={srr} hermite={hrr} n={n}"));
        }
    }
    report("high degree", tried, &failures);
    assert!(failures.is_empty());
}
