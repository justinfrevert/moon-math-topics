use crate::{
    circuit::{Circuit, Operation},
    circuit_field::{CircuitFieldElement, FieldElementZero},
    polynomial::Polynomial,
    r1cs::R1CS,
    CircuitField, Node,
};
use std::fmt::Debug;
use std::fmt::Display;
use std::ops::{Add, Mul};

// Definition for quadratic arithmetic program
#[derive(Debug)]
pub struct QAP<F: Display> {
    r1cs: R1CS<F>,
    field: CircuitField,
    pub target_polynomial: Polynomial<CircuitFieldElement>,
    pub a: Vec<Polynomial<CircuitFieldElement>>,
    pub b: Vec<Polynomial<CircuitFieldElement>>,
    pub c: Vec<Polynomial<CircuitFieldElement>>,
}

impl<F: Add + Mul + Clone + Debug + Display + FieldElementZero> QAP<F> {
    pub fn new(r1cs: R1CS<F>, field: CircuitField) -> Result<Self, QAPError> {
        // TODO: remove unwrap
        if field.0 < r1cs.a.len().try_into().unwrap() {
            return Err(QAPError::FieldTooSmall);
        }

        let constraint_count = r1cs.a.len();
        let mut random_elements = vec![];
        let mut constraint_polynomials = vec![];
        // let mut interpolated_polynomials = vec![];

        println!("Constriant count: {}", constraint_count);

        // Populate the polynomials used for finding the t polynomial, and keep their values for later
        for i in 0..constraint_count {
            let arbitrary_element = field.random_element();
            random_elements.push(arbitrary_element.clone());
            // Populate contraint polynomials for target polynomial
            constraint_polynomials.push(Polynomial::new(vec![
                -arbitrary_element.clone(),
                field.element(1),
            ]));

            // populate points we will want to interpolate
            // interpolated_polynomials.push((arbitrary_element, i));
        }

        // get shape of r1cs for placing interpolated polynomials
        // let mut interpolated_polynomials = r1cs.a.clone();

        // let mut points = vec![];

        // for i in 0..constraint_count {
        //     points.push(vec![]);
        //     for (element_idx, arbitrary_element) in random_elements.iter().enumerate() {
        //         // let points = vec![];
        //         points[element_idx].push((arbitrary_element, i))

        //         // interpolated_polynomials[i][element_idx] =
        //         //     Polynomial::lagrange_interpolation(points, field.clone());
        //     }
        // }

        // let temporary_vals = vec![field.element(5), field.element(7)];?

        // for (i, arbitrary_element) in temporary_vals.iter().enumerate() {
        //     // let arbitrary_element = field.random_element();
        //     random_elements.push(arbitrary_element.clone());
        //     // Populate contraint polynomials for target polynomial
        //     constraint_polynomials.push(Polynomial::new(vec![
        //         -arbitrary_element.clone(),
        //         field.element(1),
        //     ]));

        //     // populate points we will want to interpolate
        //     interpolated_polynomials.push((arbitrary_element, i));
        // }

        let target_polynomial = Self::get_target_polynomial(constraint_polynomials, field.clone());

        let mut a_polys: Vec<Polynomial<CircuitFieldElement>> = vec![];
        let mut b_polys: Vec<Polynomial<CircuitFieldElement>> = vec![];
        let mut c_polys: Vec<Polynomial<CircuitFieldElement>> = vec![];

        // After we interpolate polynomials, re-iterate, placing those polynomials where necessary
        for i in 0..constraint_count {
            println!("now checking r1cs a: {:?}", r1cs.a[i]);

            let a: Vec<Polynomial<CircuitFieldElement>> = r1cs.a[i]
                .iter()
                .map(|inner| match inner.is_zero() {
                    true => Polynomial::new(vec![field.element(0)]),
                    false => {
                        let i_field = field.element(i.try_into().unwrap());
                        Polynomial::<CircuitFieldElement>::lagrange_interpolation(
                            vec![(random_elements[i].clone(), i_field)],
                            field.clone(),
                        )
                    }
                })
                .collect();

            // a_polys.push(a);
            // a_polys.e

            // If any values of a, b or c are 1, we need to assign the interpolated polynomial of the same index we found previously
            // let a = match r1cs.a[i] {
            //     0 => Polynomial::new(vec![field.element(0)]),
            //     // 1 => Polynomial::lagrange_interpolation(points, field)
            //     1 => Polynomial::lagrange_interpolation(vec![(random_elements[i], i)], field.clone())
            //     _ => panic!("Invalid R1CS used"),
            // };

            // match r1cs.b[i] {}
            // match r1cs.c[i] {}
        }

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
        //     random_elements.iter().for_each(|arb_num| {
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

        Ok(QAP {
            r1cs,
            field,
            target_polynomial,
            a: a_polys,
            b: b_polys,
            c: c_polys,
        })
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

    let target_polynomial = QAP::<CircuitFieldElement>::get_target_polynomial(polys, field.clone());
    assert_eq!(
        target_polynomial,
        Polynomial::new(vec![field.element(9), field.element(1), field.element(1)])
    );
}
