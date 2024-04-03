use crate::{
    circuit::{Circuit, Operation},
    circuit_field::{CircuitFieldElement, FieldElementZero},
    polynomial::Polynomial,
    r1cs::R1CS,
    CircuitField, Node,
};
use std::{collections::{HashMap, HashSet}, fmt::Debug};
use std::fmt::Display;
use std::ops::{Add, Mul};

struct InterpolationPoints {
    pub a: Vec<Vec<(CircuitFieldElement, CircuitFieldElement)>>,
    pub b: Vec<Vec<(CircuitFieldElement, CircuitFieldElement)>>,
    pub c: Vec<Vec<(CircuitFieldElement, CircuitFieldElement)>>,
}

impl InterpolationPoints {
    fn interpolate(&self, field: CircuitField) ->  InterpolatedPolynomials {
        let a = self.a.iter().map(|column_points| {
            Polynomial::<CircuitFieldElement>::lagrange_interpolation(column_points, field.clone())
        }).collect();

        let b = self.b.iter().map(|column_points| {
            Polynomial::<CircuitFieldElement>::lagrange_interpolation(column_points, field.clone())
        }).collect();

        let c = self.c.iter().map(|column_points| {
            Polynomial::<CircuitFieldElement>::lagrange_interpolation(column_points, field.clone())
        }).collect();

        InterpolatedPolynomials {
            a,
            b,
            c,
        }
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

        let mut random_elements_set: HashSet<CircuitFieldElement>  = HashSet::new();

        let temporary_elements = vec![
            field.element(5),
            field.element(7)
        ];

        // Populate the polynomials used for finding the t polynomial, and keep their values for later
        // for i in 0..constraint_count {
            // for i in 0..constraint_count {
            for i in temporary_elements {
            // let mut arbitrary_element = field.random_element();
            let mut arbitrary_element = i;

            // Ensure uniqueness of the elements, and that they are not zero
            while random_elements_set.get(&arbitrary_element).is_some() {
                arbitrary_element = field.random_element();
            }

            x_evaluation_points.push(arbitrary_element.clone());
            // Populate contraint polynomials for target polynomial
            constraint_polynomials.push(Polynomial::new(vec![
                -arbitrary_element.clone(),
                field.element(1),
            ]));
        }

        let x_evaluation_points = vec![field.element(5), field.element(7)];

        let target_polynomial = Self::get_target_polynomial(constraint_polynomials, field.clone());

        // Points to use for interpolation
        let points = Self::determine_points(x_evaluation_points.clone(), field.clone());

        // println!("Got points to interpolate: {:?} \r\n", points);



        // let mut a_polys: Vec<Polynomial<CircuitFieldElement>> = vec![];
        // let mut b_polys: Vec<Polynomial<CircuitFieldElement>> = vec![];
        // let mut c_polys: Vec<Polynomial<CircuitFieldElement>> = vec![];


        // let a_column1 = r1cs.extract_column_from_a(0);
        // let a_column2 = r1cs.extract_column_from_a(1);
        // let a_column3 = r1cs.extract_column_from_a(2);
        // println!("r1cs a complete is: {:?}\r\n", r1cs.a);
        // println!("a column first {:?}", a_column1);
        // println!("a column second {:?}", a_column2);
        // println!("a column third {:?}", a_column3);

        // let maybe_evaluation_points = a_column3.iter().enumerate().map().collect() {
        //     ()
        // }

        // let mut a_3_points_to_interpolate = vec![];
        // for column_value in a_column3.into_iter() {
        //     for x_eval_point in &x_evaluation_points {
        //         println!("value {:?}", (x_eval_point, column_value.clone()));
        //         a_3_points_to_interpolate.push(
        //             (x_eval_point, column_value.clone())
        //         )
        //     }
        // }
        // println!("Points from a column 3 to interpolate are {:?} \r\n", a_3_points_to_interpolate);
        // Next, do the above, but for all columns  - r1cs.extract_all_columns()

        println!("Hello1");

        // let mut points_to_interpolate = vec![];
        // let mut points_to_interpolate = (vec![], vec![], vec![]);
        let (columns_a, columns_b, columns_c) = r1cs.extract_all_columns();

        println!("Hello2");

        // Define a closure to handle the interpolation points for a single matrix
        let handle_matrix = |columns: Vec<Vec<CircuitFieldElement>>, x_evaluation_points: &Vec<CircuitFieldElement>| {
            let mut points = vec![];
            for column in columns {
                let mut column_points = vec![];
                for column_value in column {
                    for x_eval_point in x_evaluation_points {
                        column_points.push((x_eval_point.clone(), column_value.clone()));
                    }
                }
                points.push(column_points);
            }
            points
        };

        // Accumulate points for all matrices
        // points_to_interpolate.0.extend(handle_matrix(columns_a, &x_evaluation_points));
        // points_to_interpolate.1.extend(handle_matrix(columns_b, &x_evaluation_points));
        // points_to_interpolate.2.extend(handle_matrix(columns_c, &x_evaluation_points));


        let interpolation_points = InterpolationPoints {
            a: handle_matrix(columns_a, &x_evaluation_points),
            b: handle_matrix(columns_b, &x_evaluation_points),
            c: handle_matrix(columns_c, &x_evaluation_points)
        };

        let interpolated_polynomials: InterpolatedPolynomials = interpolation_points.interpolate(field);

        let qap = QAP {
            target_polynomial,
            a: interpolated_polynomials.a,
            b: interpolated_polynomials.b,
            c: interpolated_polynomials.c
        };

        Ok(qap)
        // After we interpolate polynomials, re-iterate, placing those polynomials where necessary
        // for i in 0..constraint_count {
            // println!("now checking r1cs a: {:?}", r1cs.a[i]);

            // iterate through a's columns, interpolating any polynomials along the way
            // let a_column: Vec<Polynomial<CircuitFieldElement>> = r1cs.a[i]

            // Transform column into own list
            // let a_column_points: Vec<(CircuitFieldElement, CircuitFieldElement)> = r1cs.a[i]
            //     .iter()
            //     .map(|cell| match cell.is_zero() {
            //         // true => Polynomial::new(vec![field.element(0)]),
            //         true => Polynomial::new(vec![field.element(0)]),
            //         false => {
            //             // println!("Iterating... i field is {:?}. cell is: {:?}", i, cell);
            //             let i_field = field.element(i.try_into().unwrap());
            //             (cell, i_field)
            //             // Polynomial::<CircuitFieldElement>::lagrange_interpolation(
            //             //     vec![(x_evaluation_points[i].clone(), i_field)],
            //             //     field.clone(),
            //             // )
            //         }
            //     })
            //     .collect();

            // let a_column_points: Vec<F> = r1cs.a[i].clone()
            // .into_iter()
            // .map(|cell| cell)
            // .collect();



            // println!("Result column for a {:?}", a_column_points);

            // a_polys.push(a);
            // a_polys.e

            // If any values of a, b or c are 1, we need to assign the interpolated polynomial of the same index we found previously
            // let a = match r1cs.a[i] {
            //     0 => Polynomial::new(vec![field.element(0)]),
            //     // 1 => Polynomial::lagrange_interpolation(points, field)
            //     1 => Polynomial::lagrange_interpolation(vec![(x_evaluation_points[i], i)], field.clone())
            //     _ => panic!("Invalid R1CS used"),
            // };

            // match r1cs.b[i] {}
            // match r1cs.c[i] {}
        // }

        // let testing_five = field.element(5);
        // let testing_seven = field.element(7);

        // let a_2 = Polynomial::<CircuitFieldElement>::lagrange_interpolation(
        //     vec![
        //         (field.element(5), field.element(1)),
        //         (field.element(7), field.element(0)),
        //     ],
        //     field.clone(),
        // );
        // let a_5 = Polynomial::<CircuitFieldElement>::lagrange_interpolation(
        //     vec![
        //         (field.element(5), field.element(0)),
        //         (field.element(7), field.element(1)),
        //     ],
        //     field.clone(),
        // );

        // a_polys.push(a_2);
        // a_polys.push(a_5);

        // let mut interpolated_polynomials = vec![];
        // // populate points we will want to interpolate
        // for i in 0..constraint_count {
        //     x_evaluation_points.iter().for_each(|arb_num| {
        //         interpolated_polynomials.push((arb_num, i));
        //     })
        // }

        // Populate the QAP, interpolating the points where we need to

        // for i in 0..constraint_count {
        //     let a = match r1cs.a[i] {
        //         0 => Polynomial::new(vec![field.element(0)]),
        //         1 => Polynomial::lagrange_interpolation(interpolated_polynomials[i], field),
        //     };
        // }

        // Ok(QAP {
        //     target_polynomial,
        //     a: a_polys,
        //     b: b_polys,
        //     c: c_polys,
        // })
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

    /// Determine the list of points (x, y) which should be used for interpolation for the proposed QAP
    /// given the evaluation points for x values, assuming that the amount of evaluation points is equal to the amount of constraints in the circuit
    fn determine_points(x_evaluation_points: Vec<CircuitFieldElement>, field: CircuitField) -> Vec<(CircuitFieldElement, CircuitFieldElement)> {
        // for (i, x_evaluation_point) in 0..x_evaluation_points.iter().enumerate() {
        let mut points = vec![];
        for x_evaluation_point in x_evaluation_points.iter() {
            // For each evaluation point, we need to add an evaluation from 0..y
            for i in 0..x_evaluation_points.len() {
                let field_element_i = field.element(i.try_into().unwrap());
                points.push((x_evaluation_point.clone(), field_element_i))
            }
        };
        points
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
fn determines_correct_interpolation_points() {
    let field = CircuitField(13);
    let x_evaluation_points = vec![field.element(5), field.element(7)];
    let points = QAP::determine_points(x_evaluation_points, field.clone());
    let expected = vec![
        (field.element(5), field.element(0)),
        (field.element(5), field.element(1)),
        (field.element(7), field.element(0)),
        (field.element(7), field.element(1)),
    ];

    assert_eq!(points, expected)
}


