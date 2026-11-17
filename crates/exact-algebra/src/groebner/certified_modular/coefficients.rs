//! Private word-prime coefficients with context-free field identities.
//!
//! Every input coefficient is tagged with its validated prime. Only the field's
//! zero and one lack a prime; arithmetic binds them to the other operand's field.
//! A rational value is retained for operations involving only these identities.
//! Tagging is never erased by arithmetic, including for residues zero and one.

use super::*;
use crate::modular_arithmetic::{add_modulo_u64, multiply_modulo_u64, subtract_modulo_u64};
use num_traits::{One, Zero};
use std::ops::{Add, Mul, Sub};

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct WordPrimeField;
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) enum WordPrimeCoefficient {
    Identity(BigRational),
    Residue { value: u64, prime: ValidatedPrime },
}
#[derive(Clone, Copy)]
enum Operation {
    Add,
    Subtract,
    Multiply,
    Divide,
}
impl WordPrimeCoefficient {
    fn apply(self, other: Self, operation: Operation) -> Self {
        use WordPrimeCoefficient::*;
        let prime = match (&self, &other) {
            (Residue { prime: a, .. }, Residue { prime: b, .. }) => {
                assert_eq!(a, b, "modular coefficient primes differ");
                *a
            }
            (Residue { prime, .. }, _) | (_, Residue { prime, .. }) => *prime,
            (Identity(_), Identity(_)) => {
                let (Identity(a), Identity(b)) = (self, other) else {
                    unreachable!()
                };
                return Identity(match operation {
                    Operation::Add => a + b,
                    Operation::Subtract => a - b,
                    Operation::Multiply => a * b,
                    Operation::Divide => a / b,
                });
            }
        };
        let residue = |coefficient| match coefficient {
            Residue { value, .. } => value,
            Identity(value) => {
                reduce_coefficient(&value, prime).expect("field identity denominator is invertible")
            }
        };
        let (a, b) = (residue(self), residue(other));
        let value = match operation {
            Operation::Add => add_modulo_u64(a, b, prime.get()),
            Operation::Subtract => subtract_modulo_u64(a, b, prime.get()),
            Operation::Multiply => multiply_modulo_u64(a, b, prime.get()),
            Operation::Divide => multiply_modulo_u64(
                a,
                prime.inverse(b).expect("nonzero field divisor"),
                prime.get(),
            ),
        };
        Residue { value, prime }
    }
    pub(super) fn residue(&self, prime: ValidatedPrime) -> u64 {
        match self {
            Self::Residue {
                value,
                prime: actual,
            } => {
                assert_eq!(*actual, prime, "modular coefficient primes differ");
                *value
            }
            Self::Identity(value) => {
                reduce_coefficient(value, prime).expect("field identity denominator is invertible")
            }
        }
    }
}
impl Add for WordPrimeCoefficient {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        self.apply(rhs, Operation::Add)
    }
}
impl Sub for WordPrimeCoefficient {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        self.apply(rhs, Operation::Subtract)
    }
}
impl Mul for WordPrimeCoefficient {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        self.apply(rhs, Operation::Multiply)
    }
}
impl CoefficientRing for WordPrimeField {
    type Element = WordPrimeCoefficient;
    fn zero() -> Self::Element {
        WordPrimeCoefficient::Identity(BigRational::zero())
    }
    fn one() -> Self::Element {
        WordPrimeCoefficient::Identity(BigRational::one())
    }
    fn is_zero(coefficient: &Self::Element) -> bool {
        match coefficient {
            WordPrimeCoefficient::Identity(value) => value.is_zero(),
            WordPrimeCoefficient::Residue { value, .. } => *value == 0,
        }
    }
}
impl ExactField for WordPrimeField {
    fn divide(a: Self::Element, b: Self::Element) -> Self::Element {
        a.apply(b, Operation::Divide)
    }
}
