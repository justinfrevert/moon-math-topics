use crate::{
    circuit::{Circuit, Operation},
    circuit_field::CircuitFieldElement,
    polynomial::Polynomial,
    r1cs::R1CS,
    CircuitField, Node,
};
use std::ops::{Add, Mul};

// Definition for quadratic arithmetic program
pub struct QAP<F> {
    r1cs: R1CS<F>,
    field: CircuitField,
    pub target_polynomial: Polynomial<CircuitFieldElement>,
}

impl<F: Add + Mul + Clone> QAP<F> {
    pub fn new(r1cs: R1CS<F>, field: CircuitField) -> Result<Self, QAPError> {
        // TODO: remove unwrap
        if field.0 < r1cs.a.len().try_into().unwrap() {
            return Err(QAPError::FieldTooSmall);
        }

        // 3 steps:

        // 1:
        // k = amount of constraints.
        // choose k different invertible elements from the field
        let constraint_count = r1cs.a.len();
        let mut arbitrary_numbers = vec![];
        let mut constraint_polynomials = vec![];

        // Populate the polynomials used for finding the t polynomial, and keep their values for later
        for _ in 0..constraint_count {
            let arbitrary_element = field.random_element();
            arbitrary_numbers.push(arbitrary_element.clone());
            constraint_polynomials
                .push(Polynomial::new(vec![-arbitrary_element, field.element(1)]));
        }

        let testing_five = field.element(5);
        let testing_seven = field.element(7);

        // Temporarily use polynomials from exercise for easy follow along
        let temporary_polynomials = vec![
            // Polynomial::new(vec![-field.element(5), field.element(1)]),
            Polynomial::new(vec![-testing_five.clone(), field.element(1)]),
            // Polynomial::new(vec![-field.element(7), field.element(1)]),
            Polynomial::new(vec![-testing_seven, field.element(1)]),
        ];

        // let target_polynomial = Self::get_target_polynomial(constraint_polynomials, field.clone());
        let target_polynomial = Self::get_target_polynomial(temporary_polynomials, field.clone());

        let mut a_polys: Vec<Polynomial<F>> = vec![];
        let mut b_polys: Vec<Polynomial<F>> = vec![];
        let mut c_polys: Vec<Polynomial<F>> = vec![];
        for i in 0..constraint_count {
            // let m_1 =;
            let a_i = r1cs.a[i].clone();
            // a_polys.push(target_polynomial.lagrange_polynomial());
        }

        Ok(QAP {
            r1cs,
            field,
            target_polynomial,
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
