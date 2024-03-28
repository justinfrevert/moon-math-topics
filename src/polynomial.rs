use crate::circuit_field::{CircuitField, FieldElementZero};
use std::fmt::Display;
use std::ops::{Add, Mul};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Polynomial<F>(Vec<F>);

impl<F: Add + Mul> Polynomial<F> {
    pub fn new(coefficients: Vec<F>) -> Self {
        Polynomial(coefficients)
    }
}

impl<F: Add<Output = F> + Mul<Output = F> + FieldElementZero + Clone> Mul for Polynomial<F> {
    type Output = Self;
    fn mul(self, other: Polynomial<F>) -> <Self as Mul<Polynomial<F>>>::Output {
        let mut prod = vec![self.0[0].zero(); self.0.len() + other.0.len() - 1];
        for (i, lhs) in self.0.iter().enumerate() {
            for (j, rhs) in other.0.iter().enumerate() {
                // Figure out how we can do it without cloning
                prod[i + j] = prod[i + j].clone() + lhs.clone() * rhs.clone();
            }
        }
        Polynomial::new(prod)
    }
}

impl<F: Display> Display for Polynomial<F> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let mut result = String::new();
        for (i, coeff) in self.0.iter().enumerate() {
            if i > 0 {
                result.push_str(" + ");
                result.push_str(&format!("{}x^{}", coeff, i));
            } else {
                result.push_str(&format!("{}", coeff));
            }
        }

        write!(f, "{}", result.replace("00", ""))
    }
}

#[test]
fn multiplication() {
    let lhs = Polynomial::new(vec![1, 1, 5, 3]);
    let rhs = Polynomial::new(vec![3, 5, 2]);
    let ans = lhs * rhs;
    assert_eq!(ans, Polynomial::new(vec![3, 8, 22, 36, 25, 6]));
}

// #[test]
// fn multiplication2() {
//     let field = CircuitField(13);

//     let lhs = Polynomial::new(vec![-field.element(5), field.element(1)]);
//     let rhs = Polynomial::new(vec![-field.element(7), field.element(1)]);
//     let ans = lhs * rhs;
//     assert_eq!(ans, Polynomial::new(vec![field.element(1)]));
// }
