use crypto_bigint::{ConstZero, Constants, U512};
use group::ff::Field;
use rand::RngCore;

use crate::{
    circuit_field::{CircuitFieldElement, FieldElementOne, FieldElementZero},
    polynomial::Polynomial,
    r1cs::R1CS,
    CircuitField,
};
use blstrs::Scalar;

use std::{
    collections::HashSet,
    fmt::Debug,
    ops::{Add, Mul, Neg},
};

struct InterpolationPoints<F> {
    pub a: Vec<Vec<(F, F)>>,
    pub b: Vec<Vec<(F, F)>>,
    pub c: Vec<Vec<(F, F)>>,
}

impl InterpolationPoints<Scalar> {
    fn interpolate_scalars(&self) -> InterpolatedPolynomials<Scalar> {
        let a = self
            .a
            .iter()
            .map(|column_points| Polynomial::<Scalar>::lagrange_interpolation_scalar(column_points))
            .collect();

        let b = self
            .b
            .iter()
            .map(|column_points| Polynomial::<Scalar>::lagrange_interpolation_scalar(column_points))
            .collect();

        let c = self
            .c
            .iter()
            .map(|column_points| Polynomial::<Scalar>::lagrange_interpolation_scalar(column_points))
            .collect();

        InterpolatedPolynomials { a, b, c }
    }
}

// impl<F: Add + Mul + Neg<Output = F> + FieldElementZero + FieldElementOne>InterpolationPoints<CircuitFieldElement> {
impl InterpolationPoints<CircuitFieldElement> {
    fn interpolate(&self, field: CircuitField) -> InterpolatedPolynomials<CircuitFieldElement> {
        let a = self
            .a
            .iter()
            .map(|column_points| {
                Polynomial::<CircuitFieldElement>::lagrange_interpolation(
                    column_points,
                    field.clone(),
                )
            })
            .collect();

        let b = self
            .b
            .iter()
            .map(|column_points| {
                Polynomial::<CircuitFieldElement>::lagrange_interpolation(
                    column_points,
                    field.clone(),
                )
            })
            .collect();

        let c = self
            .c
            .iter()
            .map(|column_points| {
                Polynomial::<CircuitFieldElement>::lagrange_interpolation(
                    column_points,
                    field.clone(),
                )
            })
            .collect();

        InterpolatedPolynomials { a, b, c }
    }
}

pub struct InterpolatedPolynomials<F> {
    pub a: Vec<Polynomial<F>>,
    pub b: Vec<Polynomial<F>>,
    pub c: Vec<Polynomial<F>>,
}

// Definition for quadratic arithmetic program
#[derive(Clone, Debug)]
pub struct QAP<E> {
    pub target_polynomial: Polynomial<E>,
    pub a: Vec<Polynomial<E>>,
    pub b: Vec<Polynomial<E>>,
    pub c: Vec<Polynomial<E>>,
}

impl QAP<Scalar> {
    fn get_target_polynomial_scalar(
        constraint_polynomials: Vec<Polynomial<Scalar>>,
    ) -> Polynomial<Scalar> {
        let one_polynomial = Polynomial::new(vec![Scalar::ONE]);
        constraint_polynomials
            .into_iter()
            .fold(one_polynomial, |acc, x| acc * x.clone())
    }

    // Takes columns of an a, b, or c R1CS matrix, and returns points to interpolate for each column of the matrix, assumes scalars
    fn get_interpolation_points_per_column_scalar(
        columns: Vec<Vec<Scalar>>,
        x_evaluation_points: &Vec<Scalar>,
    ) -> Vec<Vec<(Scalar, Scalar)>> {
        let mut points = vec![];
        for column in columns.iter() {
            let mut column_points = vec![];
            for (column_value, x_eval_point) in column.iter().zip(x_evaluation_points) {
                column_points.push((x_eval_point.clone(), column_value.clone()));
            }
            points.push(column_points);
        }
        points
    }

    // Second implementation of `new` with more realistic Scalar type
    pub fn new_with_scalars(r1cs: R1CS<Scalar>, rng: impl RngCore + Clone) -> Result<Self, QAPError> {
        let constraint_count = r1cs.a.len();
        let mut x_evaluation_points = vec![];
        let mut constraint_polynomials = vec![];

        // Set of the arbitrary element hashes
        let mut random_elements_set: HashSet<[u8; 32]> = HashSet::new();

        for _ in 0..constraint_count {
            let mut arbitrary_element = Scalar::random(rng.clone());
            // Get hashable data of arbitrary element
            let mut arbitrary_element_id = arbitrary_element.to_bytes_be();

            // Ensure uniqueness of the elements, and that they are not zero
            while random_elements_set.get(&arbitrary_element_id).is_some() {
                arbitrary_element = Scalar::random(rng.clone());
                arbitrary_element_id = arbitrary_element.to_bytes_be();
            }

            x_evaluation_points.push(arbitrary_element.clone());

            random_elements_set.insert(arbitrary_element_id.clone());
            // Populate contraint polynomials for target polynomial
            constraint_polynomials.push(Polynomial::new(vec![
                -arbitrary_element.clone(),
                Scalar::ONE,
            ]));
        }

        let target_polynomial = Self::get_target_polynomial_scalar(constraint_polynomials);

        let (columns_a, columns_b, columns_c) = r1cs.extract_all_columns();

        let interpolation_points = InterpolationPoints {
            a: Self::get_interpolation_points_per_column_scalar(columns_a, &x_evaluation_points),
            b: Self::get_interpolation_points_per_column_scalar(columns_b, &x_evaluation_points),
            c: Self::get_interpolation_points_per_column_scalar(columns_c, &x_evaluation_points),
        };

        let interpolated_polynomials: InterpolatedPolynomials<Scalar> =
            interpolation_points.interpolate_scalars();

        let qap = QAP {
            target_polynomial,
            a: interpolated_polynomials.a,
            b: interpolated_polynomials.b,
            c: interpolated_polynomials.c,
        };
        Ok(qap)
    }

    pub fn verify(&self) -> bool {
        // let z = Polynomial::new(vec![field.element(U512::ZERO)]);
        let z = Polynomial::new(vec![Scalar::ZERO]);

        // Turn list of polynomials for a, b, and c into a single polynomial each
        let (a, b, c) = self.a.iter().zip(self.b.iter().zip(self.c.iter())).fold(
            (z.clone(), z.clone(), z),
            |acc, (a, (b, c))| {
                let a = acc.0 + a.clone();
                let b = acc.1 + b.clone();
                let c = acc.2 + c.clone();
                (a, b, c)
            },
        );

        let qap_evaluation = a * b - c;

        let (_, polynomial_remainder) = qap_evaluation.div_rem(self.target_polynomial.clone());
        polynomial_remainder.is_zero()
    }
}

// impl<F: Add + Mul + Clone + Debug + Display + FieldElementZero + Into<CircuitFieldElement>> QAP {
impl QAP<CircuitFieldElement> {
    // Original implementation, for learning only
    pub fn new(r1cs: R1CS<CircuitFieldElement>, field: CircuitField) -> Result<Self, QAPError> {
        let constraint_count = r1cs.a.len();
        let mut x_evaluation_points = vec![];
        let mut constraint_polynomials = vec![];

        let mut random_elements_set: HashSet<CircuitFieldElement> = HashSet::new();

        for _ in 0..constraint_count {
            let mut arbitrary_element = field.random_element();
            // let mut arbitrary_element = i;

            // Ensure uniqueness of the elements, and that they are not zero
            while random_elements_set.get(&arbitrary_element).is_some() {
                arbitrary_element = field.random_element();
            }

            x_evaluation_points.push(arbitrary_element.clone());
            random_elements_set.insert(arbitrary_element.clone());
            // Populate contraint polynomials for target polynomial
            constraint_polynomials.push(Polynomial::new(vec![
                -arbitrary_element.clone(),
                field.element(U512::ONE),
            ]));
        }

        let x_evaluation_points = vec![
            field.element(U512::from_u32(5)),
            field.element(U512::from_u32(7)),
        ];

        let target_polynomial = Self::get_target_polynomial(constraint_polynomials, field.clone());

        let (columns_a, columns_b, columns_c) = r1cs.extract_all_columns();

        let interpolation_points = InterpolationPoints {
            a: Self::get_interpolation_points_per_column(columns_a, &x_evaluation_points),
            b: Self::get_interpolation_points_per_column(columns_b, &x_evaluation_points),
            c: Self::get_interpolation_points_per_column(columns_c, &x_evaluation_points),
        };

        let interpolated_polynomials: InterpolatedPolynomials<CircuitFieldElement> =
            interpolation_points.interpolate(field);

        let qap = QAP {
            target_polynomial,
            a: interpolated_polynomials.a,
            b: interpolated_polynomials.b,
            c: interpolated_polynomials.c,
        };

        Ok(qap)
    }

    // Takes columns of an a, b, or c R1CS matrix, and returns points to interpolate for each column of the matrix
    fn get_interpolation_points_per_column(
        columns: Vec<Vec<CircuitFieldElement>>,
        x_evaluation_points: &Vec<CircuitFieldElement>,
    ) -> Vec<Vec<(CircuitFieldElement, CircuitFieldElement)>> {
        let mut points = vec![];
        for column in columns.iter() {
            let mut column_points = vec![];
            for (column_value, x_eval_point) in column.iter().zip(x_evaluation_points) {
                column_points.push((x_eval_point.clone(), column_value.clone()));
            }
            points.push(column_points);
        }
        points
    }

    fn get_target_polynomial(
        constraint_polynomials: Vec<Polynomial<CircuitFieldElement>>,
        field: CircuitField,
    ) -> Polynomial<CircuitFieldElement> {
        let one_polynomial = Polynomial::new(vec![field.element(U512::ONE)]);
        constraint_polynomials
            .into_iter()
            .fold(one_polynomial, |acc, x| acc * x.clone())
    }

    pub fn verify(&self, field: CircuitField) -> bool {
        let z = Polynomial::new(vec![field.element(U512::ZERO)]);

        // Turn list of polynomials for a, b, and c into a single polynomial each
        let (a, b, c) = self.a.iter().zip(self.b.iter().zip(self.c.iter())).fold(
            (z.clone(), z.clone(), z),
            |acc, (a, (b, c))| {
                let a = acc.0 + a.clone();
                let b = acc.1 + b.clone();
                let c = acc.2 + c.clone();
                (a, b, c)
            },
        );

        let qap_evaluation = a * b - c;
        let (_, polynomial_remainder) = qap_evaluation.div_rem(self.target_polynomial.clone());
        polynomial_remainder.is_zero()
    }
}

#[derive(Debug)]
pub enum QAPError {
    FieldTooSmall,
}

#[test]
fn gets_target_polynomials() {
    let field = CircuitField(U512::from_u32(13));
    let polys = vec![
        Polynomial::new(vec![
            -field.element(U512::from_u32(7)),
            field.element(U512::from_u32(1)),
        ]),
        Polynomial::new(vec![
            -field.element(U512::from_u32(5)),
            field.element(U512::from_u32(1)),
        ]),
    ];

    let target_polynomial = QAP::get_target_polynomial(polys, field.clone());
    assert_eq!(
        target_polynomial,
        Polynomial::new(vec![
            field.element(U512::from_u32(9)),
            field.element(U512::from_u32(1)),
            field.element(U512::from_u32(1))
        ])
    );
}

#[test]
fn determines_correct_interpolation_points_for_columns() {
    let field = CircuitField(U512::from_u32(13));
    let column = vec![
        vec![
            field.element(U512::from_u32(0)),
            field.element(U512::from_u32(0)),
        ],
        vec![
            field.element(U512::from_u32(0)),
            field.element(U512::from_u32(0)),
        ],
        vec![
            field.element(U512::from_u32(1)),
            field.element(U512::from_u32(0)),
        ],
        vec![
            field.element(U512::from_u32(0)),
            field.element(U512::from_u32(0)),
        ],
        vec![
            field.element(U512::from_u32(0)),
            field.element(U512::from_u32(1)),
        ],
        vec![
            field.element(U512::from_u32(0)),
            field.element(U512::from_u32(0)),
        ],
    ];
    // Example numbers and expected result from text
    let x_evaluation_points = vec![
        field.element(U512::from_u32(5)),
        field.element(U512::from_u32(7)),
    ];
    let column_interpolation =
        QAP::get_interpolation_points_per_column(column, &x_evaluation_points);

    let expected: Vec<Vec<(CircuitFieldElement, CircuitFieldElement)>> = vec![
        vec![
            (
                field.element(U512::from_u32(5)),
                field.element(U512::from_u32(0)),
            ),
            (
                field.element(U512::from_u32(7)),
                field.element(U512::from_u32(0)),
            ),
        ],
        vec![
            (
                field.element(U512::from_u32(5)),
                field.element(U512::from_u32(0)),
            ),
            (
                field.element(U512::from_u32(7)),
                field.element(U512::from_u32(0)),
            ),
        ],
        vec![
            (
                field.element(U512::from_u32(5)),
                field.element(U512::from_u32(1)),
            ),
            (
                field.element(U512::from_u32(7)),
                field.element(U512::from_u32(0)),
            ),
        ],
        vec![
            (
                field.element(U512::from_u32(5)),
                field.element(U512::from_u32(0)),
            ),
            (
                field.element(U512::from_u32(7)),
                field.element(U512::from_u32(0)),
            ),
        ],
        vec![
            (
                field.element(U512::from_u32(5)),
                field.element(U512::from_u32(0)),
            ),
            (
                field.element(U512::from_u32(7)),
                field.element(U512::from_u32(1)),
            ),
        ],
        vec![
            (
                field.element(U512::from_u32(5)),
                field.element(U512::from_u32(0)),
            ),
            (
                field.element(U512::from_u32(7)),
                field.element(U512::from_u32(0)),
            ),
        ],
    ];

    assert_eq!(column_interpolation, expected)
}
