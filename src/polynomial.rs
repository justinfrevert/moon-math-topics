use crypto_bigint::{ConstZero, U512};

use crate::circuit_field::{CircuitField, CircuitFieldElement, FieldElementOne, FieldElementZero};
use std::cmp::max;
use std::fmt::Display;
use std::ops::{Add, Div, Mul, Neg, Sub};
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Polynomial<F>(pub Vec<F>);

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

// #[test]
// fn subtracts() {
//     // 3x^2+8x+10
//     //  2x^2+2+2

//     let field = CircuitField(U512::from_u32(13));
//     let a = poly_from_field_and_integers(vec![10, 8, 3], field.clone());
//     let b = poly_from_field_and_integers(vec![2, 2, 2], field.clone());

//     assert_eq!(a - b, poly_from_field_and_integers(vec![6, 8, 1], field));
// }
