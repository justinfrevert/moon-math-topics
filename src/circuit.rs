use crate::circuit_field::{CicuitFieldInfo, CircuitFieldElement, NewValue};
use crate::circuit_field::{FieldElementOne, FieldElementZero};
use crate::r1cs::{UninitializedR1CS, R1CS};
use core::fmt::Debug;
use std::{
    collections::HashMap,
    ops::{Add, Mul},
};

#[derive(Clone, Debug)]
pub struct Constant<F>(F);

#[derive(Clone, Debug)]
pub enum Operation {
    Add,
    Multiply,
}

#[derive(Clone, Debug)]
pub struct Node<F: Add<Output = F> + Mul<Output = F> + Clone + Debug + PartialOrd + Ord> {
    value: Option<F>,
    operation: Option<Operation>,
    lhs_index: Option<usize>,
    rhs_index: Option<usize>,
}

impl<F: Add<Output = F> + Mul<Output = F> + Clone + Debug + PartialOrd + Ord> Node<F> {
    pub fn constant(value: F) -> Self {
        Node {
            value: Some(value),
            operation: None,
            lhs_index: None,
            rhs_index: None,
        }
    }
    pub fn operation(operation: Operation, lhs_index: usize, rhs_index: usize) -> Self {
        Node {
            operation: Some(operation),
            lhs_index: Some(lhs_index),
            rhs_index: Some(rhs_index),
            value: None,
        }
    }
}

fn circuit_calculate<F: Add<Output = F> + Mul<Output = F> + Clone + Debug + Ord + PartialOrd>(gates: Vec<Node<F>>) -> Option<F> {
        // A cache for numbers calculated on-the-fly in any Operative nodes. This might support not having to execute multiple
        // loops to order the nodes
        let mut supporting_results: HashMap<usize, F> = HashMap::new();

        for (idx, node) in gates.iter().enumerate() {
            if let Some(operation) = &node.operation {
                let lhs_index = node.lhs_index.unwrap();
                let rhs_index = node.rhs_index.unwrap();

                let lhs = gates[lhs_index]
                    .value
                    .clone()
                    .unwrap_or_else(|| supporting_results.get(&lhs_index).unwrap().clone());
                let rhs = gates[rhs_index]
                    .value
                    .clone()
                    .unwrap_or_else(|| supporting_results.get(&rhs_index).unwrap().clone());

                let calculation = match operation {
                    Operation::Add => lhs + rhs,
                    Operation::Multiply => lhs * rhs,
                };
                supporting_results.insert(idx, calculation);
            }
        }
        supporting_results.get(&(gates.len() - 1)).cloned()
}

fn circuit_calculate_with_trace<F: Add<Output = F> + Mul<Output = F> + Clone + Debug + Ord + PartialOrd   + FieldElementZero + FieldElementOne>(gates: Vec<Node<F>>) -> (Option<F>, R1CS<F>) {
    let mut supporting_results: HashMap<usize, F> = HashMap::new();

    // A value representing zero, contextually
    let zero_value = gates[0].clone().value.unwrap().zero();
    let one_value = gates[0].clone().value.unwrap().one();
    let r1cs_length = gates.len() + 1;

    let mut uninitialized_r1cs = UninitializedR1CS::new(one_value.clone(), zero_value.clone(), r1cs_length, gates.len());

    let mut current_constraint = 0;
    for (node_idx, node) in gates.iter().enumerate() {
        if let Some(operation) = &node.operation {
            // Current node is an operation we'll add a gate
            uninitialized_r1cs.add_constraint(zero_value.clone());

            let lhs_index = node.lhs_index.unwrap();
            let rhs_index = node.rhs_index.unwrap();

            let lhs = gates[lhs_index]
                .value
                .clone()
                .unwrap_or_else(|| supporting_results.get(&lhs_index).unwrap().clone());
            let rhs = gates[rhs_index]
                .value
                .clone()
                .unwrap_or_else(|| supporting_results.get(&rhs_index).unwrap().clone());

            let calculation = match operation {
                // Remove clones...
                Operation::Add => lhs.clone() + rhs.clone(),
                Operation::Multiply => lhs.clone() * rhs.clone(),
            };

            if calculation > lhs.zero() {
                uninitialized_r1cs.add_to_constraint_c(node_idx, current_constraint, calculation.clone(), one_value.clone());
            }

            if lhs > lhs.zero() {
                uninitialized_r1cs.add_to_constraint_a(lhs_index, lhs, current_constraint, one_value.clone());
            }

            if rhs > rhs.zero() {
                uninitialized_r1cs.add_to_constraint_b(rhs_index, rhs, current_constraint, one_value.clone());
            }

            supporting_results.insert(node_idx, calculation);
            current_constraint += 1;
        }
    }
    let r1cs = uninitialized_r1cs.into_r1cs();
    (
        // Final gate result, or output
        supporting_results.get(&(gates.len() - 1)).cloned(),
        r1cs,
    )
}

#[derive(Debug, Clone)]
// Initial circuit implementation, written to more or less fit the description in the text, to get an intuition on it
pub struct Circuit<
    // F: Add<Output = F> + Mul<Output = F> + Clone + Debug + Ord + PartialOrd,
    FI: CicuitFieldInfo,
> {
    // gates: Vec<Node<F>>,
    gates: Vec<Node<CircuitFieldElement>>,
    // Any extra field information
    field: FI,
}

impl<
        FI: CicuitFieldInfo,
    > Circuit<FI>
{
    pub fn new(gates: Vec<Node<CircuitFieldElement>>, field: FI) -> Self {
        Circuit { gates, field }
    }

    pub fn calculate(self) -> Option<CircuitFieldElement> {
        circuit_calculate(self.gates)
    }

    pub fn calculate_with_trace(self) -> (Option<CircuitFieldElement>, R1CS<CircuitFieldElement>) {
        circuit_calculate_with_trace(self.gates)
    }
}


#[derive(Debug, Clone)]
// A circuit that is differentiated from the above circuit that was implemented for learning due to not requiring field info.
pub struct CircuitNoFieldInfo<
    F: Add<Output = F> + Mul<Output = F> + Clone + Debug + Ord + PartialOrd + FieldElementZero + FieldElementOne,
> {
    gates: Vec<Node<F>>
}

impl<F: Add<Output = F> + Mul<Output = F> + Clone + Debug + Ord + PartialOrd + FieldElementZero + FieldElementOne> CircuitNoFieldInfo<F> {
    pub fn new(gates: Vec<Node<F>>) -> Self {
        CircuitNoFieldInfo { gates }
    }

    pub fn calculate_with_trace(&self) -> (Option<F>, R1CS<F>) {
        circuit_calculate_with_trace(self.clone().gates)
    }
}
