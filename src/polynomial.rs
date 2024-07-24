use blstrs::{G1Projective, G2Projective, Scalar};
use crypto_bigint::U512;

use crate::circuit_field::{CircuitField, CircuitFieldElement, FieldElementOne, FieldElementZero};
use std::cmp::max;
use std::fmt::Display;
use std::ops::{Add, Div, Mul, Neg, Sub};
use group::{ff::Field, Group};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Polynomial<F>(pub Vec<F>);

// New implementation using blstrs library
impl Polynomial<Scalar> {
    pub fn evaluate(&self, tau: Scalar) -> Scalar {
        let mut result = Scalar::ZERO;
        let mut tau_power = Scalar::ONE;
        for coeff in self.0.iter() {
            result += coeff * tau_power;
            tau_power *= tau;
        }
        result
    }

    fn trim(&self) -> Self {
        let mut new_coeffs = self.0.clone();
        while let Some(true) = new_coeffs.last().map(|c| Field::is_zero(c).into()) {
            new_coeffs.pop();
        }
        Polynomial::new(new_coeffs)
    }

    pub fn div_rem(self, divisor: Self) -> (Self, Self) {
        let zero = Scalar::ZERO;

        if self.is_zero() {
            let zero_poly = Polynomial::new(vec![zero]);
            return (zero_poly.clone(), zero_poly);
        }

        let mut quotient = Polynomial::new(vec![zero; self.0.len() - divisor.0.len() + 1]);
        let mut remainder: Polynomial<Scalar> = self;
        let divisor_leading_inv = divisor.leading_coefficient().unwrap().invert().unwrap();
        while !remainder.is_zero() && remainder.0.len() >= divisor.0.len() {
            let cur_q_coeff = remainder.leading_coefficient().unwrap() * &divisor_leading_inv;
            let cur_q_degree = remainder.0.len() - divisor.0.len();
            quotient.0[cur_q_degree] = cur_q_coeff.clone();
            for (i, div_coeff) in divisor.0.iter().enumerate() {
                if cur_q_degree + i < remainder.0.len() {
                    let old_value = remainder.0[cur_q_degree + i].clone();
                    remainder.0[cur_q_degree + i] -= &(cur_q_coeff.clone() * div_coeff.clone());
                }
            }
            remainder = remainder.trim();
        }
        (quotient, remainder)
    }

    // TODO: these functions aren't that useful. Can probably be replaced by their inner operations
    pub fn evaluate_polynomial_commitment(
        &self,
        tau_powers: &[G1Projective],
        // poly: &Polynomial<Scalar>
    ) -> G1Projective {
        G1Projective::multi_exp(&tau_powers, &self.0)
    }
    
    pub fn evaluate_polynomial_commitment_g2(
        &self,
        tau_powers: &[G2Projective],
        // poly: &Polynomial<Scalar>
    ) -> G2Projective {
        G2Projective::multi_exp(&tau_powers, &self.0)
    }

}

// Old, manual implementation for circuit fields
impl Polynomial<CircuitFieldElement> {
    pub fn evaluate(&self, point: CircuitFieldElement) -> CircuitFieldElement {
        let mut total = self.0[0].zero();
        let field = &self.0[0].field;
        for (i, coefficient) in self.0.iter().enumerate() {
            // total += pow(point, i) * coefficient
            let current = CircuitFieldElement::new(U512::from_u128(i as u128), field.clone());
            let current_field_element = CircuitFieldElement::new(point.pow(current), field.clone());
            total += current_field_element * coefficient.clone()
        }
        total
    }

    pub fn div_rem(self, divisor: Self) -> (Self, Self) {
        let zero = self.0[0].zero();
        let mut quotient = Polynomial::new(vec![zero; self.0.len() - divisor.0.len() + 1]);
        let mut remainder: Polynomial<CircuitFieldElement> = self;

        let divisor_leading_inv = divisor.leading_coefficient().unwrap().invert();

        while !remainder.is_zero() && remainder.0.len() >= divisor.0.len() {
            let cur_q_coeff = remainder.leading_coefficient().unwrap() * &divisor_leading_inv;
            let cur_q_degree = remainder.0.len() - divisor.0.len();
            quotient.0[cur_q_degree] = cur_q_coeff.clone();

            for (i, div_coeff) in divisor.0.iter().enumerate() {
                if cur_q_degree + i < remainder.0.len() {
                    remainder.0[cur_q_degree + i] -= &(cur_q_coeff.clone() * div_coeff.clone());
                }
            }
            remainder = remainder.trim();
        }

        (quotient, remainder)
    }

    fn trim(&self) -> Self {
        let mut new_coeffs = self.0.clone();
        while let Some(true) = new_coeffs.last().map(|c| c.is_zero()) {
            new_coeffs.pop();
        }
        Polynomial::new(new_coeffs)
    }
}

impl<F: Add + Mul + Neg<Output = F> + FieldElementZero + FieldElementOne> Polynomial<F> {
    pub fn new(coefficients: Vec<F>) -> Self {
        Polynomial(coefficients)
    }
    pub fn is_zero(&self) -> bool {
        self.0.is_empty() || self.0.iter().all(|coeff| coeff.is_zero())
    }

    fn leading_coefficient(&self) -> Option<&F> {
        self.0.last().clone()
    }

    // Original lagrange interpolation with assumptions of the learning implementation of the circuit field
    pub fn lagrange_interpolation(
        points: &Vec<(CircuitFieldElement, CircuitFieldElement)>,
        field: CircuitField,
    ) -> Polynomial<CircuitFieldElement> {
        let zero = field.element(U512::ZERO);
        let one = field.element(U512::ONE);
        let mut interpolation = Polynomial::new(vec![zero.clone()]);

        for (i, (x_i, y_i)) in points.iter().enumerate() {
            let mut numerator = Polynomial::new(vec![one.clone()]);
            let mut denominator = one.clone();

            for (j, (x_j, _)) in points.iter().enumerate() {
                if x_i != x_j {
                    numerator = numerator * Polynomial::new(vec![-x_j.clone(), one.clone()]);
                    denominator = denominator * (x_i.clone() - x_j.clone());
                }
            }

            let denominator_poly = Polynomial::new(vec![denominator]);
            let lagrange_basis_poly = numerator / denominator_poly;

            // Represent y as a polynomial to multiply with other values here
            let y_polynomial = Polynomial::new(vec![y_i.clone()]);
            interpolation = interpolation + (lagrange_basis_poly * y_polynomial);
        }
        interpolation
    }

    pub fn lagrange_interpolation_scalar(
        points: &Vec<(Scalar, Scalar)>,
    ) -> Polynomial<Scalar> {
        let zero = Scalar::ZERO;
        let one = Scalar::ONE;
        let mut interpolation = Polynomial::new(vec![zero.clone()]);

        for (i, (x_i, y_i)) in points.iter().enumerate() {
            let mut numerator = Polynomial::new(vec![one.clone()]);
            let mut denominator = one.clone();

            for (j, (x_j, _)) in points.iter().enumerate() {
                if x_i != x_j {
                    numerator = numerator * Polynomial::new(vec![-x_j.clone(), one.clone()]);
                    denominator = denominator * (x_i.clone() - x_j.clone());
                }
            }

            let denominator_poly = Polynomial::new(vec![denominator]);
            let lagrange_basis_poly = numerator / denominator_poly;

            // Represent y as a polynomial to multiply with other values here
            let y_polynomial = Polynomial::new(vec![y_i.clone()]);
            interpolation = interpolation + (lagrange_basis_poly * y_polynomial);
        }
        interpolation
    }
}

impl<
        F: Add<Output = F>
            + Mul<Output = F>
            + FieldElementZero
            + FieldElementOne
            + Clone
            + Neg<Output = F>,
    > Mul for Polynomial<F>
{
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

impl Mul<Scalar> for Polynomial<Scalar> {
    type Output = Polynomial<Scalar>;

    fn mul(self, rhs: Scalar) -> Polynomial<Scalar> {
        let new_coefficients: Vec<Scalar> = self.0.into_iter().map(|coeff| coeff * rhs).collect();
        Polynomial::new(new_coefficients)
    }
}

impl Mul<Scalar> for &Polynomial<Scalar> {
    type Output = Polynomial<Scalar>;

    fn mul(self, rhs: Scalar) -> Polynomial<Scalar> {
        let new_coefficients: Vec<Scalar> = self.0.iter().map(|&coeff| coeff * rhs).collect();
        Polynomial::new(new_coefficients)
    }
}

impl Mul<Polynomial<Scalar>> for Scalar {
    type Output = Polynomial<Scalar>;

    fn mul(self, rhs: Polynomial<Scalar>) -> Polynomial<Scalar> {
        let new_coefficients: Vec<Scalar> = rhs.0.into_iter().map(|coeff| self * coeff).collect();
        Polynomial::new(new_coefficients)
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

impl Div for Polynomial<Scalar> {
    type Output = Self;
    fn div(self, divisor: Self) -> Self::Output {
        let zero = self.0[0].zero();
        let mut quotient = Polynomial::new(vec![zero; self.0.len() - divisor.0.len() + 1]);
        let mut remainder: Polynomial<Scalar> = self.clone();
        // Can unwrap here because we know self is not zero.
        let divisor_leading_inv = divisor.leading_coefficient().unwrap().invert().unwrap();
        while !remainder.is_zero() && remainder.0.len() >= divisor.0.len() {
            let cur_q_coeff = remainder.leading_coefficient().unwrap() * &divisor_leading_inv;
            let cur_q_degree = remainder.0.len() - divisor.0.len();
            quotient.0[cur_q_degree] = cur_q_coeff.clone();

            for (i, div_coeff) in divisor.0.iter().enumerate() {
                remainder.0[cur_q_degree + i] -= &(cur_q_coeff.clone() * div_coeff.clone());
                // remainder.0[cur_q_degree + i] = remainder.clone().0.clone()[cur_q_degree.clone() + i.clone()] - (cur_q_coeff.clone() * div_coeff.clone());
            }
            while let Some(true) = remainder.0.last().map(|c| Field::is_zero(c).into()) {
                remainder.0.pop();
            }
        }
        quotient
    }
}

impl Add for Polynomial<CircuitFieldElement> {
    type Output = Polynomial<CircuitFieldElement>;

    fn add(self, rhs: Polynomial<CircuitFieldElement>) -> Polynomial<CircuitFieldElement> {
        let mut out = vec![];
        let zero = self.0[0].zero();
        for i in 0..(max(self.0.len(), rhs.0.len())) {
            let left = self.0.get(i).unwrap_or(&zero);
            let right = rhs.0.get(i).unwrap_or(&zero);

            out.push(left + right);
        }
        Polynomial::new(out)
    }
}

impl Add for Polynomial<Scalar> {
    type Output = Polynomial<Scalar>;

    fn add(self, rhs: Polynomial<Scalar>) -> Polynomial<Scalar> {
        let mut out = vec![];
        let zero = self.0[0].zero();
        for i in 0..(max(self.0.len(), rhs.0.len())) {
            let left = self.0.get(i).unwrap_or(&zero);
            let right = rhs.0.get(i).unwrap_or(&zero);

            out.push(left + right);
        }
        Polynomial::new(out)
    }
}

impl Sub for Polynomial<CircuitFieldElement> {
    type Output = Polynomial<CircuitFieldElement>;
    fn sub(self, rhs: Polynomial<CircuitFieldElement>) -> Polynomial<CircuitFieldElement> {
        let mut out = vec![];
        let zero = self.0[0].zero();
        for i in 0..max(self.0.len(), rhs.0.len()) {
            let left = self.0.get(i).unwrap_or(&zero);
            let right = rhs.0.get(i).unwrap_or(&zero);

            out.push(left.clone() - right.clone());
        }
        Polynomial::new(out)
    }
}

impl Sub for Polynomial<Scalar> {
    type Output = Polynomial<Scalar>;
    fn sub(self, rhs: Polynomial<Scalar>) -> Polynomial<Scalar> {
        let mut out = vec![];
        let zero = self.0[0].zero();
        for i in 0..max(self.0.len(), rhs.0.len()) {
            let left = self.0.get(i).unwrap_or(&zero);
            let right = rhs.0.get(i).unwrap_or(&zero);

            out.push(left.clone() - right.clone());
        }
        Polynomial::new(out)
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
pub fn poly_from_field_and_integers(
    coefficients: Vec<u64>,
    field: CircuitField,
) -> Polynomial<CircuitFieldElement> {
    let mut elements = vec![];
    for c in coefficients {
        elements.push(field.element(c.into()))
    }
    Polynomial(elements)
}

#[test]
fn multiplication() {
    let field = CircuitField(U512::from_u32(4201));
    // let lhs = Polynomial::new(vec![1, 1, 5, 3]);
    let lhs = poly_from_field_and_integers(vec![1, 1, 5, 3], field.clone());
    let rhs = poly_from_field_and_integers(vec![3, 5, 2], field.clone());
    let ans = lhs * rhs;
    assert_eq!(
        ans,
        poly_from_field_and_integers(vec![3, 8, 22, 36, 25, 6], field.clone())
    );
}

#[test]
fn adds() {
    let field = CircuitField(U512::from_u32(4200));
    let lft = Polynomial::new(vec![-field.element(U512::from_u32(1)), field.element(U512::from_u32(2)), field.element(U512::from_u32(1))]);
    let rhs = Polynomial::new(vec![field.element(U512::from_u32(6)), -field.element(U512::from_u32(3)), field.element(U512::from_u32(2))]);
    let ans = Polynomial::new(vec![field.element(U512::from_u32(5)), -field.element(U512::from_u32(1)), field.element(U512::from_u32(3))]);
    assert_eq!(lft + rhs, ans);
}

#[test]
fn divides_polynomials() {
    let field = CircuitField(U512::from_u32(43));
    //  2x^2+5x+3
    let dividend = poly_from_field_and_integers(vec![2, 5, 3], field.clone());
    // x + 1
    let divisor = poly_from_field_and_integers(vec![1, 1], field.clone());
    // 2x+3
    let ans = poly_from_field_and_integers(vec![2, 3], field.clone());
    assert_eq!(dividend / divisor, ans)
}

#[test]
fn lagrange_inerpolation() {
    let field = CircuitField(U512::from_u32(13));

    let point1 = (field.element(U512::from_u32(2)), field.element(U512::from_u32(4)));
    let point2 = (field.element(U512::from_u32(1)), field.element(U512::from_u32(3)));

    let poly = Polynomial::<CircuitFieldElement>::lagrange_interpolation(
        &vec![point1, point2],
        field.clone(),
    );

    let expected_poly = poly_from_field_and_integers(vec![2, 1], field.clone());

    assert_eq!(poly, expected_poly);
}
#[test]
fn lagrange_inerpolation_larger() {
    let field = CircuitField(U512::from_u32(13));

    let point1 = (field.element(U512::from_u32(2)), field.element(U512::from_u32(4)));
    let point2 = (field.element(U512::from_u32(1)), field.element(U512::from_u32(3)));
    let point3 = (field.element(U512::from_u32(7)), field.element(U512::from_u32(11)));

    let poly = Polynomial::<CircuitFieldElement>::lagrange_interpolation(
        &vec![point1, point2, point3],
        field.clone(),
    );

    let expected = poly_from_field_and_integers(vec![3, 6, 7], field.clone());
    assert_eq!(poly, expected);
}

#[test]
fn lagrange_inerpolation_2() {
    let field = CircuitField(U512::from_u32(13));

    // Points from example in text
    let points = vec![
        (field.clone().element(U512::from_u32(5)), field.clone().element(U512::from_u32(1))),
        (field.clone().element(U512::from_u32(7)), field.clone().element(U512::from_u32(0))),
    ];
    let points2 = vec![
        (field.clone().element(U512::from_u32(5)), field.clone().element(U512::from_u32(0))),
        (field.clone().element(U512::from_u32(7)), field.clone().element(U512::from_u32(1))),
    ];

    let poly = Polynomial::<CircuitFieldElement>::lagrange_interpolation(&points, field.clone());
    let poly2 = Polynomial::<CircuitFieldElement>::lagrange_interpolation(&points2, field.clone());
    let poly3 = Polynomial::<CircuitFieldElement>::lagrange_interpolation(&vec![points, points2].concat(), field.clone());
    let expected = poly_from_field_and_integers(vec![10, 6], field.clone());
    let expected2 = poly_from_field_and_integers(vec![4, 7], field.clone());
    let expected3 = poly_from_field_and_integers(vec![12, 7, 7], field);
    assert_eq!(poly, expected);
    assert_eq!(poly2, expected2);
    assert_eq!(poly3, expected3);
}

#[test]
fn eval_g_3() {
    // 39 == x^3 -4x^2 +3x -1
    let evaluate_at = 5_u64;
    let evaluate_at_commitment = G1Projective::identity() * Scalar::from(5_u64);

    let x3 = G1Projective::generator() * Scalar::from(evaluate_at.pow(3));
    let x2 = G1Projective::generator() * Scalar::from(evaluate_at.pow(2));
    let x = G1Projective::generator() * Scalar::from(evaluate_at);

    // assert_eq!(x, x + (Scalar::from(2) * ));

    // let poly = Polynomial::new(vec![
    //     Scalar::from(1),
    //     -Scalar::from(4),
    //     Scalar::from(3),
    //     -Scalar::from(1)
    // ]);

    let poly = Polynomial::new(vec![
        -Scalar::from(1),
        Scalar::from(3),
        -Scalar::from(4),
        Scalar::from(1)
    ]);

    let lhs = G1Projective::generator() * Scalar::from(39);

    let rhs = x3 * Scalar::from(1)
        + x2 * -Scalar::from(4)
        + x * Scalar::from(3)
        + G1Projective::generator() * -Scalar::from(1);

    // let rhs = x3 * Scalar::from(1)
    //     + x2 * -Scalar::from(4)
    //     + x * Scalar::from(3)
    //     + G1Projective::generator() * -Scalar::from(1);

    let mut computed_rhs = G1Projective::identity();

    let mut first_answer = G1Projective::identity();
    let mut second_answer = G1Projective::identity();
    let mut third_answer = G1Projective::identity();
    let mut fourth_answer = G1Projective::identity();

    for (rev_idx, coef) in poly.0.iter().enumerate().rev() {
        if rev_idx == 3 {
            println!("Three");
            let x = G1Projective::generator() * Scalar::from(evaluate_at.pow(rev_idx as u32));
            assert_eq!(x, x3);
            
            // let x_com = {
            //     let tt = G1Projective::generator() * evaluate_at;
            // };

            assert_eq!(*coef, Scalar::from(1));

            first_answer = x * coef;
            computed_rhs += first_answer;
            assert_eq!(first_answer, x3 * Scalar::from(1));
        } else if rev_idx == 2 {
            println!("Two");

            let x = G1Projective::generator() * Scalar::from(evaluate_at.pow(rev_idx as u32));
            assert_eq!(x, x2);
            assert_eq!(*coef, -Scalar::from(4));

            second_answer = x * coef;
            assert_eq!(second_answer, x2 * -Scalar::from(4));

            assert_eq!(
                x2 * -Scalar::from(4),
                second_answer
                );

                computed_rhs += second_answer;

        } else if rev_idx == 1 {
            println!("One");
            let x_local = G1Projective::generator() * Scalar::from(evaluate_at.pow(rev_idx as u32));
            assert_eq!(x, x_local);
            assert_eq!(*coef, Scalar::from(3));

            third_answer = x_local * coef;
            assert_eq!(third_answer, x * Scalar::from(3));

            assert_eq!(
                x * Scalar::from(3),
                third_answer
            );
            computed_rhs += third_answer;
        } else if rev_idx == 0 {
            assert_eq!(*coef, -Scalar::from(1));

            fourth_answer = G1Projective::generator() * coef;
            computed_rhs += fourth_answer;
        }

        if rev_idx > 0 {
            let x = G1Projective::generator() * Scalar::from(evaluate_at.pow(rev_idx as u32));
            assert_eq!(first_answer, x3 * Scalar::from(1));

            computed_rhs += x * coef;
        } else if rev_idx == 0 {
            computed_rhs += G1Projective::generator() * coef;
        }

        assert_eq!(
            x3 * Scalar::from(1),
            first_answer
        );
    }

    assert_eq!(
        x3 * Scalar::from(1)
        + x2 * -Scalar::from(4)
        + x * Scalar::from(3)
        + G1Projective::generator() * -Scalar::from(1),
        first_answer + second_answer + third_answer + fourth_answer
    );

    assert_eq!(poly.evaluate(Scalar::from(5)), Scalar::from(39));

    assert_eq!(lhs, rhs);

    let tau_powers: Vec<G1Projective> = (0..4)  
        .map(|i| G1Projective::generator() * Scalar::from(evaluate_at.pow(i as u32)))
        .collect();

    // assert_eq!(rhs, evaluate_polynomial_commitment(&tau_powers, &poly));
    assert_eq!(rhs, poly.evaluate_polynomial_commitment(&tau_powers));

}


// #[test]
// fn subtracts() {
//     // 3x^2+8x+10
//     //  2x^2+2+2

//     let field = CircuitField(U512::from_u32(13));
//     let a = poly_from_field_and_integers(vec![10, 8, 3], field.clone());
//     let b = poly_from_field_and_integers(vec![2, 2, 2], field.clone());

//     assert_eq!(a - b, poly_from_field_and_integers(vec![6, 8, 1], field));
// }
