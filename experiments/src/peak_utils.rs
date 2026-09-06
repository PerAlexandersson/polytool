use combpoly::statistics::peaks;

/// Trim trailing zero coefficients, keeping at least one entry.
pub fn pt(p: &[i64]) -> Vec<i64> {
    let mut v = p.to_vec();
    while v.len() > 1 && v.last() == Some(&0) {
        v.pop();
    }
    v
}

/// Test whether a coefficient vector is identically zero.
pub fn pz(p: &[i64]) -> bool {
    p.iter().all(|&c| c == 0)
}

/// Add two coefficient vectors and trim trailing zeros.
pub fn pa(a: &[i64], b: &[i64]) -> Vec<i64> {
    let len = a.len().max(b.len());
    let mut result = vec![0i64; len];
    for (i, &coeff) in a.iter().enumerate() {
        result[i] += coeff;
    }
    for (i, &coeff) in b.iter().enumerate() {
        result[i] += coeff;
    }
    pt(&result)
}

/// Multiply a polynomial by `t`.
pub fn pmt(p: &[i64]) -> Vec<i64> {
    let mut result = vec![0i64; p.len() + 1];
    for (i, &coeff) in p.iter().enumerate() {
        result[i + 1] = coeff;
    }
    pt(&result)
}

/// Convert a Ferrers board encoding into the associated 312-avoiding permutation.
pub fn board_to_perm(board: &[u8]) -> Vec<u8> {
    let n = board.len();
    let mut perm = vec![0u8; n];
    let mut used = vec![false; n + 1];
    for i in 0..n {
        for c in (1..=(board[i] as usize).min(n)).rev() {
            if !used[c] {
                perm[i] = c as u8;
                used[c] = true;
                break;
            }
        }
    }
    perm
}

/// Test whether a permutation avoids the classical pattern 312.
pub fn is_312_avoiding(perm: &[u8]) -> bool {
    let n = perm.len();
    for i in 0..n {
        for j in i + 1..n {
            for k in j + 1..n {
                if perm[k] < perm[i] && perm[i] < perm[j] {
                    return false;
                }
            }
        }
    }
    true
}

/// Count peaks in a permutation.
pub fn peak_count(perm: &[u8]) -> usize {
    peaks(perm)
}

/// Enumerate weakly increasing Ferrers boards of size `n`.
pub fn gen_boards(n: usize) -> Vec<Vec<u8>> {
    let mut result = vec![];
    let mut current = vec![];
    gb(n, n, 0, &mut current, &mut result);
    result
}

fn gb(n: usize, mx: usize, d: usize, current: &mut Vec<u8>, result: &mut Vec<Vec<u8>>) {
    if d == n {
        result.push(current.clone());
        return;
    }
    for v in (d + 1).max(if d > 0 { current[d - 1] as usize } else { 1 })..=mx {
        current.push(v as u8);
        gb(n, mx, d + 1, current, result);
        current.pop();
    }
}
