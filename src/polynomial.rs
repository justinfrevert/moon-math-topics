use crate::circuit_field::{CircuitField, CircuitFieldElement, FieldElementOne, FieldElementZero};
use std::fmt::Display;
use std::ops::{Add, Div, Mul, Neg};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Polynomial<F>(Vec<F>);

impl<F: Add + Mul + Neg<Output = F> + FieldElementZero + FieldElementOne> Polynomial<F> {
    pub fn new(coefficients: Vec<F>) -> Self {
        Polynomial(coefficients)
    }
    fn is_zero(&self) -> bool {
        self.0.is_empty() || self.0.iter().all(|coeff| coeff.is_zero())
    }
    fn leading_coefficient(&self) -> Option<&F> {
        self.0.last().clone()
    }

    // pub fn lagrange_polynomial(points: Vec<(F, F)>, field: CircuitField) -> Self {
    // pub fn lagrange_polynomial(points: Vec<(CircuitFieldElement, CircuitFieldElement)>, field: CircuitField) -> Self {
    //     // let mut x = field.element(1);
    //     let x = points[0].0.one();

    //     let sum = vec![];

    //     for (j, (x_j, y_j)) in points.into_iter().enumerate() {
    //         // let inner_product = field.element(1);
    //         let inner_product = Polynomial::new(vec![field.element(1)]);
    //         for (i, (x_i, y_i)) in points.into_iter().enumerate() {
    //             if i != j {
    //                 // let x_i = ;
    //                 // let x_j = ;
    //                 // let numerator = Polynomial::new(vec![field.element(1)]);
    //                 let numerator = Polynomial::new(vec![-x_i, x]);
    //                 let denominator = Polynomial::new(vec![-x_i, x_j]);

    //                 inner_product = inner_product * (numerator / denominator)
    //             }
    //         }
    //         sum.push(inner_product);
    //     }
    //     // outer_sum
    //     sum
    // }

    // pub fn lagrange_polynomial(points: Vec<(CircuitFieldElement, CircuitFieldElement)>, field: CircuitField) -> Self {
    //     // let mut x = field.element(1);
    //     let x = points[0].0.one();

    //     for (j, (x_j, y_j)) in points.into_iter().enumerate() {
    //         // let inner_product = field.element(1);
    //         // let inner_product = Polynomial::new(vec![field.element(1)]);
    //         for (k, (x_k, y_k)) in points.into_iter().enumerate() {
    //             if j == k {
    //                 // Do nothing
    //             } else {

    //             }
    //         }
    //     }
    //     sum
    // }
}

impl<F: Add<Output = F> + Mul<Output = F> + FieldElementZero + FieldElementOne + Clone + Neg<Output = F>> Mul for Polynomial<F> {
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

impl Div for Polynomial<CircuitFieldElement> {
    type Output = Self;
    fn div(self, divisor: Self) -> Self::Output {
        let zero = self.0[0].zero();
        let mut quotient = Polynomial::new(vec![zero; self.0.len() - divisor.0.len() + 1]);
        let mut remainder: Polynomial<CircuitFieldElement> = self.clone();
        // Can unwrap here because we know self is not zero.
        let divisor_leading_inv = divisor.leading_coefficient().unwrap().invert();
        while !remainder.is_zero() && remainder.0.len() >= divisor.0.len() {
            let cur_q_coeff = remainder.leading_coefficient().unwrap() * &divisor_leading_inv;
            let cur_q_degree = remainder.0.len() - divisor.0.len();
            quotient.0[cur_q_degree] = cur_q_coeff.clone();

            for (i, div_coeff) in divisor.0.iter().enumerate() {
                remainder.0[cur_q_degree + i] -= &(cur_q_coeff.clone() * div_coeff.clone());
                // remainder.0[cur_q_degree + i] = remainder.clone().0.clone()[cur_q_degree.clone() + i.clone()] - (cur_q_coeff.clone() * div_coeff.clone());
            }
            while let Some(true) = remainder.0.last().map(|c| c.is_zero().into()) {
                remainder.0.pop();
            }
        }
        quotient
    }
}

// impl<F: Add + Mul + Div + Neg<Output = F> + FieldElementZero + Clone> Div for Polynomial<F> {
//     type Output = Self;
//     fn div(self, divisor: Self) -> Self::Output {
//         let zero = self.0[0].zero();

//         let mut quotient = Polynomial::new(vec![zero; self.0.len() - divisor.0.len() + 1]);

//         let mut remainder: Polynomial<F> = self.clone();
//         // Can unwrap here because we know self is not zero.
//         let divisor_leading_inv = divisor.leading_coefficient().unwrap().invert();
//         while !remainder.is_zero() && remainder.0.len() >= divisor.0.len() {
//             let cur_q_coeff = remainder.leading_coefficient().unwrap() * &divisor_leading_inv;
//             let cur_q_degree = remainder.0.len() - divisor.0.len();
//             quotient.0[cur_q_degree] = cur_q_coeff;

//             for (i, div_coeff) in divisor.0.iter().enumerate() {
//                 remainder.0[cur_q_degree + i] -= &(cur_q_coeff * div_coeff);
//             }
//             while let Some(true) = remainder.0.last().map(|c| c.is_zero().into()) {
//                 remainder.0.pop();
//             }
//         }
//         quotient
//     }
// }

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

// Helper function to simplify instantiating polynomials with less boilerplate
pub fn poly_from_field_and_integers(coefficients: Vec<u64>, field: CircuitField) -> Polynomial<CircuitFieldElement> {
    let mut elements = vec![];
    for c in coefficients {
        elements.push(field.element(c))
    }
    Polynomial(elements)
}

#[test]
fn multiplication() {
    let field = CircuitField(4200);
    // let lhs = Polynomial::new(vec![1, 1, 5, 3]);
    let lhs = poly_from_field_and_integers(vec![1, 1, 5, 3], field.clone());
    let rhs = poly_from_field_and_integers(vec![3, 5, 2], field.clone());
    let ans = lhs * rhs;
    assert_eq!(ans, poly_from_field_and_integers(vec![3, 8, 22, 36, 25, 6], field.clone()));
}

// #[test]
// fn multiplication2() {
//     let field = CircuitField(13);

//     let lhs = Polynomial::new(vec![-field.element(5), field.element(1)]);
//     let rhs = Polynomial::new(vec![-field.element(7), field.element(1)]);
//     let ans = lhs * rhs;
//     assert_eq!(ans, Polynomial::new(vec![field.element(1)]));
// }

#[test]
fn divides_polynomials() {
    let field = CircuitField(42);
    //  2x^2+5x+3
    let dividend = poly_from_field_and_integers(vec![2, 5, 3], field.clone());
    // x + 1
    let divisor = poly_from_field_and_integers(vec![1, 1], field.clone());
    // 2x+3
    let ans = poly_from_field_and_integers(vec![2, 3], field.clone());
    assert_eq!(dividend / divisor, ans)
}