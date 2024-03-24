pub mod circuit;
pub mod circuit_field;


use circuit::{Circuit, Operation::*};
use crate::{
    circuit::Node,
    circuit_field::CircuitField
};

fn main() {
    // Use the circuit functions
    // Try showing 3 factors of a number 27
    let instructions = vec![
        Node::constant(3),
        Node::constant(3),
        Node::operation(Multiply, 0, 1),

        Node::constant(3),
        Node::operation(Multiply, 2, 3),
    ];

    let mut c = Circuit::new(instructions);
    c.calculate();
}

#[test]
fn three_factors_example_works() {
    let number = 27;
    let x = 3;

    let instructions = vec![
        Node::constant(x),
        Node::constant(x),
        Node::operation(Multiply, 0, 1),
        Node::constant(x),
        Node::operation(Multiply, 2, 3),
    ];

    let mut c = Circuit::new(instructions);
    assert_eq!(c.calculate(), Some(number));
}

#[test]
fn addition() {
    let number = 9;

    let instructions = vec![
        Node::constant(4),
        Node::constant(5),
        Node::operation(Add, 0, 1),
    ];

    let mut c = Circuit::new(instructions);
    assert_eq!(c.calculate(), Some(number));
}

#[test]
fn addition_and_multiplication() {
    let number = 81;

    let instructions = vec![
        Node::constant(4),
        Node::constant(5),
        Node::operation(Add, 0, 1),
        Node::constant(9),
        Node::operation(Multiply, 2, 3)
    ];

    let mut c = Circuit::new(instructions);
    assert_eq!(c.calculate(), Some(number));
}

#[test]
fn with_generic_types() {
    let field = CircuitField(13);
    let x = field.element(3);
    let number = field.element(1);

    let instructions = vec![
        Node::constant(x.clone()),
        Node::constant(x.clone()),
        Node::operation(Multiply, 0, 1),
        Node::constant(x),
        Node::operation(Multiply, 2, 3),
    ];

    let c = Circuit::new(instructions);
    assert_eq!(c.calculate(), Some(number));
}

#[test]
fn addition_and_multiplication_with_modulus() {
    let field = CircuitField(13);
    let number = field.element(81);

    let instructions = vec![ 
        Node::constant(field.element(4)),
        Node::constant(field.element(5)),
        Node::operation(Add, 0, 1),
        Node::constant(field.element(9)),
        Node::operation(Multiply, 2, 3)
    ];

    let mut c = Circuit::new(instructions);
    assert_eq!(c.calculate(), Some(number));
}

#[test]
fn squares() {
    let number = 36;
    let instructions = vec![
        Node::constant(2),
        Node::constant(3),
        Node::operation(Multiply, 0, 0), // 2 ^2
        Node::operation(Multiply, 1, 1), // 3^2
        Node::operation(Multiply, 2, 3), // 2^2 * 3^2
    ];
    let mut c = Circuit::new(instructions);
    assert_eq!(c.calculate(), Some(number));
}

#[test]
fn squares_with_modulus() {
    let field = CircuitField(9);
    let number = field.element(36);

    let instructions = vec![
        Node::constant(field.element(2)),
        Node::constant(field.element(3)),
        Node::operation(Multiply, 0, 0), // 2 ^2
        Node::operation(Multiply, 1, 1), // 3^2
        Node::operation(Multiply, 2, 3), // 2^2 * 3^2
    ];

    let mut c = Circuit::new(instructions);
    assert_eq!(c.calculate(), Some(number));
}

#[test]
fn point_on_curve_circuit() {
    // From text: the tiny jub-jub curve is considered over the prime field \mathbb{F}_{13}
    let field = CircuitField(13);
    // 0 = 1 + 8 · x2 · y2 + 10 · x2 + 12y2 ⇔ 
    let expected_result = field.element(0);

    // Some points to test: (1,2), (1, 11), (4, 0)
    let x = field.element(1);
    let y = field.element(2);

    let instructions = vec![
        Node::constant(field.element(1)), // 0
        Node::constant(field.element(8)), // 1
        Node::constant(x),                // 2
        Node::constant(y),                // 3
        Node::constant(field.element(10)),// 4 
        Node::constant(field.element(12)),// 5

        
        // Initial multiplications
        Node::operation(Multiply, 2, 2), // 1. x * x // idx 6
        Node::operation(Multiply, 3, 3), // 2. y * y
        Node::operation(Multiply, 4, 6), // 3. 10 * x^2
        Node::operation(Multiply, 6, 7), // 4. y^2 * x^2 // idx 9
        Node::operation(Multiply, 5, 7), // 5. y^2 * 12
        Node::operation(Multiply, 1, 9), // 6. 8 * (x^2 * y^2) // 11

        // Additions
        Node::operation(Add, 0, 8), // 7. 1 + (10x^2) // idx 12
        Node::operation(Add, 11, 12), // 8. (8 * x^2 * y^2) + (1 + (10x^2)))
        Node::operation(Add, 10, 13), // 9. ((8 * x^2 * y^2) + (1 + (10x^2))) + (12 *y^2))
    ];

    let c = Circuit::new(instructions);
    assert_eq!(c.calculate(), Some(expected_result));
}