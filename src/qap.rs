use crate::{circuit_field::CircuitFieldElement, polynomial::Polynomial, r1cs::R1CS, CircuitField};

use std::{collections::HashSet, fmt::Debug};

struct InterpolationPoints {
    pub a: Vec<Vec<(CircuitFieldElement, CircuitFieldElement)>>,
    pub b: Vec<Vec<(CircuitFieldElement, CircuitFieldElement)>>,
    pub c: Vec<Vec<(CircuitFieldElement, CircuitFieldElement)>>,
}

impl InterpolationPoints {
    fn interpolate(&self, field: CircuitField) -> InterpolatedPolynomials {
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

pub struct InterpolatedPolynomials {
    pub a: Vec<Polynomial<CircuitFieldElement>>,
    pub b: Vec<Polynomial<CircuitFieldElement>>,
    pub c: Vec<Polynomial<CircuitFieldElement>>,
}

// Definition for quadratic arithmetic program
#[derive(Debug)]
pub struct QAP {
    pub target_polynomial: Polynomial<CircuitFieldElement>,
    pub a: Vec<Polynomial<CircuitFieldElement>>,
    pub b: Vec<Polynomial<CircuitFieldElement>>,
    pub c: Vec<Polynomial<CircuitFieldElement>>,
}

// impl<F: Add + Mul + Clone + Debug + Display + FieldElementZero + Into<CircuitFieldElement>> QAP {
impl QAP {
    pub fn new(r1cs: R1CS<CircuitFieldElement>, field: CircuitField) -> Result<Self, QAPError> {
        // TODO: remove unwrap
        if field.0 < r1cs.a.len().try_into().unwrap() {
            return Err(QAPError::FieldTooSmall);
        }

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
                field.element(1),
            ]));
        }

        let x_evaluation_points = vec![field.element(5), field.element(7)];

        let target_polynomial = Self::get_target_polynomial(constraint_polynomials, field.clone());

        let (columns_a, columns_b, columns_c) = r1cs.extract_all_columns();

        let interpolation_points = InterpolationPoints {
            a: Self::get_interpolation_points_per_column(columns_a, &x_evaluation_points),
            b: Self::get_interpolation_points_per_column(columns_b, &x_evaluation_points),
            c: Self::get_interpolation_points_per_column(columns_c, &x_evaluation_points),
        };

        let interpolated_polynomials: InterpolatedPolynomials =
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
        let one_polynomial = Polynomial::new(vec![field.element(1)]);
        constraint_polynomials
            .into_iter()
            .fold(one_polynomial, |acc, x| acc * x.clone())
    }

    pub fn verify(&self, field: CircuitField) -> bool {
        let z = Polynomial::new(vec![field.element(0)]);

        // Turn list of polynomials for a, b, and c into a single polynomial each
        let (a, b, c) = self.a
            .iter()
            .zip(self.b.iter().zip(self.c.iter()))
            .fold((z.clone(), z.clone(), z), |acc, (a, (b, c))| {
                let a = acc.0 + a.clone();
                let b = acc.1 + b.clone();
                let c = acc.2 + c.clone();
                (a, b, c)
            });

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
    let field = CircuitField(13);
    let polys = vec![
        Polynomial::new(vec![-field.element(7), field.element(1)]),
        Polynomial::new(vec![-field.element(5), field.element(1)]),
    ];

    let target_polynomial = QAP::get_target_polynomial(polys, field.clone());
    assert_eq!(
        target_polynomial,
        Polynomial::new(vec![field.element(9), field.element(1), field.element(1)])
    );
}

#[test]
fn determines_correct_interpolation_points_for_columns() {
    let field = CircuitField(13);
    let column = vec![
        vec![field.element(0), field.element(0)],
        vec![field.element(0), field.element(0)],
        vec![field.element(1), field.element(0)],
        vec![field.element(0), field.element(0)],
        vec![field.element(0), field.element(1)],
        vec![field.element(0), field.element(0)],
    ];
    // Example numbers and expected result from text
    let x_evaluation_points = vec![field.element(5), field.element(7)];
    let column_interpolation =
        QAP::get_interpolation_points_per_column(column, &x_evaluation_points);

    let expected: Vec<Vec<(CircuitFieldElement, CircuitFieldElement)>> = vec![
        vec![
            (field.element(5), field.element(0)),
            (field.element(7), field.element(0)),
        ],
        vec![
            (field.element(5), field.element(0)),
            (field.element(7), field.element(0)),
        ],
        vec![
            (field.element(5), field.element(1)),
            (field.element(7), field.element(0)),
        ],
        vec![
            (field.element(5), field.element(0)),
            (field.element(7), field.element(0)),
        ],
        vec![
            (field.element(5), field.element(0)),
            (field.element(7), field.element(1)),
        ],
        vec![
            (field.element(5), field.element(0)),
            (field.element(7), field.element(0)),
        ],
    ];

    assert_eq!(column_interpolation, expected)
}
