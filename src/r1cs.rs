// TODOS in order to implement R1CS from Alebraic circuits:
// Implement boolean, or equality check constraint
use core::fmt::Debug;

use crate::circuit_field::CircuitFieldElement;

const CONSTANTS_INDEX: usize = 0;
const OUT_INDEX: usize = 1;
const R1CS_OFFSET: usize = 2;

#[derive(Debug)]
pub struct R1CS<F> {
    pub a: Vec<Vec<F>>,
    pub b: Vec<Vec<F>>,
    pub c: Vec<Vec<F>>,
    pub witness: Vec<F>,
}

// impl R1CS<CircuitFieldElement> {
//     pub fn new(a: Vec<Vec<CircuitFieldElement>>, b: Vec<Vec<CircuitFieldElement>>, c: Vec<Vec<CircuitFieldElement>>, witness: Vec<CircuitFieldElement>) -> Self {
//         R1CS { a, b, c, witness }
//     }

//     pub fn extract_column_from_a(&self, column_index: usize) -> Vec<CircuitFieldElement> {
//         self.a.iter().map(|row| row[column_index].clone()).collect()
//     }

//     pub fn extract_column_from_b(&self, column_index: usize) -> Vec<CircuitFieldElement> {
//         self.b.iter().map(|row| row[column_index].clone()).collect()
//     }

//     pub fn extract_column_from_c(&self, column_index: usize) -> Vec<CircuitFieldElement> {
//         self.c.iter().map(|row| row[column_index].clone()).collect()
//     }

//     // Method to extract all columns from matrices 'a', 'b', and 'c'
//     pub fn extract_all_columns(&self) -> (Vec<Vec<CircuitFieldElement>>, Vec<Vec<CircuitFieldElement>>, Vec<Vec<CircuitFieldElement>>) {
//         let num_columns = if !self.a.is_empty() {
//             self.a[0].len()
//         } else {
//             0 // Assuming all matrices have the same number of columns
//         };

//         let mut columns_a = Vec::with_capacity(num_columns);
//         let mut columns_b = Vec::with_capacity(num_columns);
//         let mut columns_c = Vec::with_capacity(num_columns);

//         for column_index in 0..num_columns {
//             columns_a.push(self.extract_column_from_a(column_index));
//             columns_b.push(self.extract_column_from_b(column_index));
//             columns_c.push(self.extract_column_from_c(column_index));
//         }

//         (columns_a, columns_b, columns_c)
//     }
// }

impl<E: Clone> R1CS<E> {
    pub fn new(a: Vec<Vec<E>>, b: Vec<Vec<E>>, c: Vec<Vec<E>>, witness: Vec<E>) -> Self {
        R1CS { a, b, c, witness }
    }

    pub fn extract_column_from_a(&self, column_index: usize) -> Vec<E> {
        self.a.iter().map(|row| row[column_index].clone()).collect()
    }

    pub fn extract_column_from_b(&self, column_index: usize) -> Vec<E> {
        self.b.iter().map(|row| row[column_index].clone()).collect()
    }

    pub fn extract_column_from_c(&self, column_index: usize) -> Vec<E> {
        self.c.iter().map(|row| row[column_index].clone()).collect()
    }

    // Method to extract all columns from matrices 'a', 'b', and 'c'
    pub fn extract_all_columns(&self) -> (Vec<Vec<E>>, Vec<Vec<E>>, Vec<Vec<E>>) {
        let num_columns = if !self.a.is_empty() {
            self.a[0].len()
        } else {
            0 // Assuming all matrices have the same number of columns
        };

        let mut columns_a = Vec::with_capacity(num_columns);
        let mut columns_b = Vec::with_capacity(num_columns);
        let mut columns_c = Vec::with_capacity(num_columns);

        for column_index in 0..num_columns {
            columns_a.push(self.extract_column_from_a(column_index));
            columns_b.push(self.extract_column_from_b(column_index));
            columns_c.push(self.extract_column_from_c(column_index));
        }

        (columns_a, columns_b, columns_c)
    }
}



// #[derive(Debug)]
// pub struct UninitializedR1CS<F> {
//     pub a: Vec<Vec<F>>,
//     pub b: Vec<Vec<F>>,
//     pub c: Vec<Vec<F>>,
//     pub witness: Vec<F>,
//     pub r1cs_length: usize,
//     pub total_gates: usize,
// }

#[derive(Debug)]
pub struct UninitializedR1CS<CircuitFieldElement> {
    pub a: Vec<Vec<CircuitFieldElement>>,
    pub b: Vec<Vec<CircuitFieldElement>>,
    pub c: Vec<Vec<CircuitFieldElement>>,
    pub witness: Vec<CircuitFieldElement>,
    pub r1cs_length: usize,
    pub total_gates: usize,
}

impl UninitializedR1CS<CircuitFieldElement> {
    pub fn new(one: CircuitFieldElement, zero: CircuitFieldElement, r1cs_length: usize, total_gates: usize) -> Self {
        let a: Vec<Vec<CircuitFieldElement>> = vec![];
        let b: Vec<Vec<CircuitFieldElement>> = vec![];
        let c: Vec<Vec<CircuitFieldElement>> = vec![];
        let mut witness: Vec<CircuitFieldElement> = vec![];

        for _ in 0..r1cs_length {
            witness.push(zero.clone());
        }
        witness[CONSTANTS_INDEX] = one;

        UninitializedR1CS {
            a,
            b,
            c,
            witness,
            r1cs_length,
            total_gates,
        }
    }

    // Adds a new row to the current R1CS, which represents an area to store constraints
    pub fn add_constraint(&mut self, zero: CircuitFieldElement) {
        let mut gate = vec![];
        for _ in 0..self.r1cs_length {
            gate.push(zero.clone());
        }
        self.a.push(gate.clone());
        self.b.push(gate.clone());
        self.c.push(gate);
    }

    pub fn add_to_constraint_a(
        &mut self,
        lhs_index: usize,
        lhs: CircuitFieldElement,
        current_constraint: usize,
        one_value: CircuitFieldElement,
    ) {
        // Original index, but bumped right past the const value, as well as the out value(2 places)
        let lhs_r1cs_bumped = lhs_index + R1CS_OFFSET;
        self.a[current_constraint][lhs_r1cs_bumped] = one_value;
        // self.witness.push(lhs);
        self.witness[lhs_r1cs_bumped] = lhs;
    }

    pub fn add_to_constraint_b(
        &mut self,
        rhs_index: usize,
        rhs: CircuitFieldElement, 
        current_constraint: usize,
        one_value: CircuitFieldElement
    ) {
        // Original index, but bumped right past the const value, as well as the out value(2 places)
        let rhs_r1cs_bumped = rhs_index + R1CS_OFFSET;
        self.witness[rhs_r1cs_bumped] = rhs;
        self.b[current_constraint][rhs_r1cs_bumped] = one_value;
    }

    pub fn add_to_constraint_c(
        &mut self,
        node_idx: usize,
        current_constraint: usize,
        calculation:CircuitFieldElement, 
        one: CircuitFieldElement
    ) {
        // Original index, but bumped right past the const value, as well as the out value(2 places)
        // We use the node index because that represents the current node(that's what the output of the gate is...)
        let out_index_r1cs_bumped = node_idx + R1CS_OFFSET;

        // If done, add the final calculated values to a specific location defined by the const. Otherwise, add to the current gate's expected location
        if node_idx == self.total_gates - 1 {
            self.c[current_constraint][OUT_INDEX] = one;
            self.witness[OUT_INDEX] = calculation;
        } else {
            self.c[current_constraint][out_index_r1cs_bumped] = one;
            self.witness[out_index_r1cs_bumped] = calculation;
        }
    }

    pub fn into_r1cs(&self) -> R1CS<CircuitFieldElement> {
        R1CS::new(
            self.a.clone(),
            self.b.clone(),
            self.c.clone(),
            self.witness.clone(),
        )
    }
}
