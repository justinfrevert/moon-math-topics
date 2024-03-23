use std::collections::HashMap;

#[derive(Clone, Debug)]
pub struct Constant(u32);

#[derive(Clone, Debug)]
pub enum Operation {
    Add,
    Multiply
}

#[derive(Clone, Debug)]
pub struct Node {
    value: Option<u32>,
    operation: Option<Operation>,
    lhs_index: Option<usize>,
    rhs_index: Option<usize>,
}

impl Node {
    pub fn constant(value: u32) -> Self {
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
            value: None
        }
    }
}

#[derive(Debug, Clone)]
pub struct Circuit(Vec<Node>);

impl Circuit {
    pub fn new(nodes: Vec<Node>) -> Self {
        Circuit(nodes.into())
    }

    pub fn calculate(mut self) -> Option<u32> {
        // A cache for numbers calculated on-the-fly in any Operative nodes. This might support not having to execute multiple
        // loops to order the nodes
        let mut supporting_results: HashMap<usize, u32> = HashMap::new();

        for (idx, node) in self.0.iter().enumerate() {
            if let Some(operation) = &node.operation {
                let lhs_index = node.lhs_index.unwrap();
                let rhs_index = node.rhs_index.unwrap();

                let lhs = self.0[lhs_index].value.unwrap_or_else(|| *supporting_results.get(&lhs_index).unwrap());
                let rhs = self.0[rhs_index].value.unwrap_or_else(|| *supporting_results.get(&rhs_index).unwrap());

                let calculation = match operation {
                    Operation::Add => lhs + rhs,
                    Operation::Multiply => lhs * rhs
                };
                supporting_results.insert(idx, calculation);
            }
        };
        supporting_results.get(&(self.0.len() - 1)).cloned()
    }
}
