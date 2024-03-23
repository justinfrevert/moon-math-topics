pub mod circuit;
use circuit::{Circuit, Operation::*};
use crate::circuit::Node;

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
fn point_on_curve_circuit() {
    // = 1 + 8 · x2 · y2 + 10 · x2 + 12y2 ⇔  flattened = 
    fn point_on_curve_circuit(x: u32, y: u32) -> Circuit {
        let instructions = vec![
            Node::constant(1), // 0
            Node::constant(8),
            Node::constant(x),
            Node::constant(y),
            Node::constant(10),
            Node::constant(12), // 5

            // Initial multiplications
            Node::operation(Multiply, 2, 2), // 1. x * x // idx 6
            Node::operation(Multiply, 3, 3), // 2. y * y
            Node::operation(Multiply, 4, 6), // 3. 10 * x^2
            Node::operation(Multiply, 6, 7), // 4. y^2 * x^2 // idx 9
            Node::operation(Multiply, 5, 7), // 5.  y^2 * 12
            Node::operation(Multiply, 1, 9), // 6. 8 * (x^2 * y^2) 

            // Additions
            Node::operation(Add, 0, 8), // 7. 1 + (10x^2) // idx 12
            Node::operation(Add, 11, 12), // 8. (8 * x^2 * y^2) + (1 + (10x^2)))
            Node::operation(Add, 5, 13), // 9. ((8 * x^2 * y^2) + (1 + (10x^2))) + (12 *y^2))
        ];
        Circuit::new(instructions)
    }

    let number = 0;
    let x = 1;
    let y = 2;
    let mut c = point_on_curve_circuit(x, y);
    assert_eq!(c.calculate(), Some(number));
}