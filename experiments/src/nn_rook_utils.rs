pub use combpoly::rook_placements::{
    add_coefficient_vectors as add, adjacent_non_nesting_packet_global as adjacent_packet_global,
    directed_weak_interlaces as interlaces, divide_by_t as div_t,
    fixed_remaining_non_nesting_packet_lgv as fixed_remaining_packet_lgv,
    fixed_row_non_nesting_count_lgv as fixed_row_count_lgv,
    has_nonnegative_coefficients as has_nonnegative_coeffs, integer_partitions as partitions,
    marked_delta_non_nesting_count_lgv as marked_delta_count_lgv,
    monomial_coefficient_vector as monomial, multiply_by_t as mul_t,
    non_nesting_rook_polynomial as nn_rook, polynomial_degree_zero_convention as degree,
    row_mask_size as subset_size, strip_ferrers_columns as strip_columns,
    subtract_coefficient_vectors as sub, trim_coefficient_vector as trim,
};

pub struct TestResult {
    pub name: &'static str,
    pub pass: usize,
    pub total: usize,
    pub fails: Vec<String>,
}

impl TestResult {
    pub fn new(name: &'static str) -> Self {
        Self {
            name,
            pass: 0,
            total: 0,
            fails: Vec::new(),
        }
    }

    pub fn check(&mut self, ok: bool, msg: impl FnOnce() -> String) {
        self.total += 1;
        if ok {
            self.pass += 1;
        } else if self.fails.len() < 8 {
            self.fails.push(msg());
        }
    }

    pub fn print(&self) {
        println!("{}: {}/{}", self.name, self.pass, self.total);
        if self.fails.is_empty() {
            println!("  All pass!\n");
        } else {
            for fail in &self.fails {
                println!("  {fail}");
            }
            println!();
        }
    }
}
