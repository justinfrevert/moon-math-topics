pub mod circuit;
pub mod circuit_field;
pub mod polynomial;
pub mod qap;
pub mod r1cs;

use crate::circuit_field::CircuitFieldElement;
use crate::{circuit::Node, circuit_field::CircuitField};
use circuit::{Circuit, Operation::*};

use std::ops::Neg;

use polynomial::Polynomial;
use qap::QAP;
fn main() {
    let field = CircuitField(13);
    // let zero = field.element(0);
    // let one = field.element(1);
    let w_0 = field.element(3.clone());
    let w_1 = field.element(4.clone());
    let w_2 = field.element(2.clone());

    let instructions = vec![
        Node::constant(w_0.clone()),
        Node::constant(w_1.clone()),
        Node::constant(w_2.clone()),
        // w_0 * w_1
        Node::operation(Multiply, 0, 1), // v1 // idx 4
        // v_1 * w_2
        Node::operation(Multiply, 2, 3), // v2 // idx 5
    ];

    let circuit = Circuit::new(instructions, field.clone());
    let (_, r1cs) = circuit.calculate_with_trace();

    println!("r1cs {:?}", r1cs.a);

    let qap = QAP::new(r1cs, field.clone()).unwrap();

    // println!("{:?}", r1cs.b);

    // // A = np.array([[0,0,1,0,0,0,0,0],
    // //     [0,0,0,0,1,0,0,0],
    // //     [0,0,0,0,0,0,1,0]])

    // // B = np.array([[0,0,0,1,0,0,0,0],
    // //         [0,0,0,0,0,1,0,0],
    // //         [0,0,0,0,0,0,0,1]])

    // // C = np.array([[0,0,0,0,0,0,1,0],
    // //         [0,0,0,0,0,0,0,1],
    // //         [0,1,0,0,0,0,0,0]])

    // // Witness should appear as:
    // // [1, out, x, y, z, u, v1, v2].

    // let a = vec![
    //     vec![
    //         zero.clone(),
    //         zero.clone(),
    //         one.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //     ],
    //     vec![
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         one.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //     ],
    //     vec![
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         one.clone(),
    //         zero.clone(),
    //     ],
    // ];

    // let b = vec![
    //     vec![
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         one.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //     ],
    //     vec![
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         one.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //     ],
    //     vec![
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         one.clone(),
    //     ],
    // ];
    // let c = vec![
    //     vec![
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         one.clone(),
    //         zero.clone(),
    //     ],
    //     vec![
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         one.clone(),
    //     ],
    //     vec![
    //         zero.clone(),
    //         one.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //         zero.clone(),
    //     ],
    // ];

    // let v1 = field.element(12); // x * y == 12
    // let v2 = field.element(10); // z * u == 10
    // let out = field.element(120); // v1 * v2 == 120

    // let expected_witness = vec![one, out.clone(), x, y, z, u, v1, v2];

    // assert_eq!(result, Some(out));
    // assert_eq!(r1cs.a, a);
    // assert_eq!(r1cs.b, b);
    // assert_eq!(r1cs.c, c);
    // assert_eq!(r1cs.witness, expected_witness);

    // let lhs = Polynomial::new(vec![
    //     field.element(1),
    //     field.element(1),
    //     field.element(5),
    //     field.element(3),
    // ]);
    // let rhs = Polynomial::new(vec![field.element(3), field.element(5), field.element(2)]);
}

#[test]
fn three_factors_example_works() {
    let number = 27;
    let x = 3;
    let unused_field = CircuitField(1);

    let instructions = vec![
        Node::constant(x),
        Node::constant(x),
        Node::operation(Multiply, 0, 1),
        Node::constant(x),
        Node::operation(Multiply, 2, 3),
    ];

    let mut c = Circuit::new(instructions, unused_field);
    assert_eq!(c.calculate(), Some(number));
}

#[test]
fn addition() {
    let number = 9;
    let unused_field = CircuitField(1);
    let instructions = vec![
        Node::constant(4),
        Node::constant(5),
        Node::operation(Add, 0, 1),
    ];

    let mut c = Circuit::new(instructions, unused_field);
    assert_eq!(c.calculate(), Some(number));
}

#[test]
fn addition_and_multiplication() {
    let number = 81;
    let unused_field = CircuitField(1);

    let instructions = vec![
        Node::constant(4),
        Node::constant(5),
        Node::operation(Add, 0, 1),
        Node::constant(9),
        Node::operation(Multiply, 2, 3),
    ];

    let mut c = Circuit::new(instructions, unused_field);
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

    let c = Circuit::new(instructions, field);
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
        Node::operation(Multiply, 2, 3),
    ];

    let mut c = Circuit::new(instructions, field);
    assert_eq!(c.calculate(), Some(number));
}

#[test]
fn squares() {
    let number = 36;
    let unusedfield = CircuitField(13);

    let instructions = vec![
        Node::constant(2),
        Node::constant(3),
        Node::operation(Multiply, 0, 0), // 2 ^2
        Node::operation(Multiply, 1, 1), // 3^2
        Node::operation(Multiply, 2, 3), // 2^2 * 3^2
    ];
    let mut c = Circuit::new(instructions, unusedfield);
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

    let mut c = Circuit::new(instructions, field);
    assert_eq!(c.calculate(), Some(number));
}

#[cfg(test)]
fn point_on_curve_circuit(
    x: CircuitFieldElement,
    y: CircuitFieldElement,
    field: CircuitField,
) -> Circuit<CircuitFieldElement, CircuitField> {
    let instructions: Vec<Node<circuit_field::CircuitFieldElement>> = vec![
        Node::constant(field.element(1)),  // 0
        Node::constant(field.element(8)),  // 1
        Node::constant(x),                 // 2
        Node::constant(y),                 // 3
        Node::constant(field.element(10)), // 4
        Node::constant(field.element(12)), // 5
        // Initial multiplications
        Node::operation(Multiply, 2, 2), // 1. x * x // idx 6
        Node::operation(Multiply, 3, 3), // 2. y * y
        Node::operation(Multiply, 4, 6), // 3. 10 * x^2
        Node::operation(Multiply, 6, 7), // 4. y^2 * x^2 // idx 9
        Node::operation(Multiply, 5, 7), // 5. y^2 * 12
        Node::operation(Multiply, 1, 9), // 6. 8 * (x^2 * y^2) // 11
        // Additions
        Node::operation(Add, 0, 8),   // 7. 1 + (10x^2) // idx 12
        Node::operation(Add, 11, 12), // 8. (8 * x^2 * y^2) + (1 + (10x^2)))
        Node::operation(Add, 10, 13), // 9. ((8 * x^2 * y^2) + (1 + (10x^2))) + (12 *y^2))
    ];
    Circuit::new(instructions, field)
}

#[test]
fn point_on_curve_circuit_works() {
    // From text: the tiny jub-jub curve is considered over the prime field \mathbb{F}_{13}
    let field = CircuitField(13);
    // 0 = 1 + 8 · x^2 · y^2 + 10 · x^2 + 12y^2
    let expected_result = field.element(0);

    // Some points to test: (1,2), (1, 11), (4, 0), (5,2) (5,11), (6,5), (6,8), ...... (12,8)
    let circuit_1_2 = point_on_curve_circuit(field.element(1), field.element(2), field.clone());
    let circuit_1_11 = point_on_curve_circuit(field.element(1), field.element(11), field.clone());
    // Doesn't work:
    let circuit_4_0 = point_on_curve_circuit(field.element(4), field.element(0), field.clone());

    assert_eq!(circuit_1_2.calculate(), Some(expected_result.clone()));
    assert_eq!(circuit_1_11.calculate(), Some(expected_result.clone()));
}

#[test]
// Example from rareskills book
fn simple_circuit() {
    let field = CircuitField(13);

    let x = field.element(2);
    let y = field.element(3);

    // x^2 y
    let instructions = vec![
        Node::constant(x),
        Node::constant(y),
        Node::operation(Multiply, 0, 0), // x^2
        Node::operation(Multiply, 1, 2), // (x^2) * y
    ];
    let c = Circuit::new(instructions, field.clone());

    assert_eq!(c.calculate(), Some(field.element(12)));
}

#[test]
// Example from rareskills book
fn simple_circuit_with_example_witness() {
    let field = CircuitField(13);

    let x = field.element(41);
    let y = field.element(103);

    // x^2 y
    let instructions = vec![
        Node::constant(x),
        Node::constant(y),
        Node::operation(Multiply, 0, 1),
    ];

    let circuit = Circuit::new(instructions, field.clone());
    let zero = field.element(0);
    let one = field.element(1);

    let (result, r1cs) = circuit.calculate_with_trace();

    // [1, out, x, y].
    // A is [0, 0, 1, 0], because x is present, and none of the other variables are.
    // B is [0, 0, 0, 1] because the variables in the right hand side are just y, and
    // C is [0, 1, 0, 0] because we only have the out variable.

    // Witness should appear as:
    // [1, 4223, 41, 103], or [1, out, x, y].
    let a = vec![[zero.clone(), zero.clone(), one.clone(), zero.clone()]];
    let b = vec![[zero.clone(), zero.clone(), zero.clone(), one.clone()]];
    let c = vec![[zero.clone(), one.clone(), zero.clone(), zero]];
    let expected_witness = vec![
        field.element(1),
        field.element(4223),
        field.element(41),
        field.element(103),
    ];

    assert_eq!(r1cs.a, a);
    assert_eq!(r1cs.b, b);
    assert_eq!(r1cs.c, c);
    assert_eq!(r1cs.witness, expected_witness);
    assert_eq!(result, Some(field.element(4223)));
}

// Example from rareskills book x * y * z * u
#[test]
fn example_r1cs_more_terms() {
    let field = CircuitField(42000);
    let zero = field.element(0);
    let one = field.element(1);
    let x = field.element(3.clone());
    let y = field.element(4.clone());
    let z = field.element(2.clone());
    let u = field.element(5.clone());

    let instructions = vec![
        // 4 inputs
        Node::constant(x.clone()), // 0
        Node::constant(y.clone()), // 1
        Node::constant(z.clone()),
        Node::constant(u.clone()),
        // x * y
        Node::operation(Multiply, 0, 1), // v1 // idx 4
        // z * u
        Node::operation(Multiply, 2, 3), // v2 // idx 5
        // v1 * v2 aka ((x * y) * (z * u))
        Node::operation(Multiply, 4, 5),
    ];

    let mut circuit = Circuit::new(instructions, field.clone());
    let (result, r1cs) = circuit.calculate_with_trace();

    // A = np.array([[0,0,1,0,0,0,0,0],
    //     [0,0,0,0,1,0,0,0],
    //     [0,0,0,0,0,0,1,0]])

    // B = np.array([[0,0,0,1,0,0,0,0],
    //         [0,0,0,0,0,1,0,0],
    //         [0,0,0,0,0,0,0,1]])

    // C = np.array([[0,0,0,0,0,0,1,0],
    //         [0,0,0,0,0,0,0,1],
    //         [0,1,0,0,0,0,0,0]])

    // Witness should appear as:
    // [1, out, x, y, z, u, v1, v2].

    let a = vec![
        vec![
            zero.clone(),
            zero.clone(),
            one.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
        ],
        vec![
            zero.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
            one.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
        ],
        vec![
            zero.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
            one.clone(),
            zero.clone(),
        ],
    ];

    let b = vec![
        vec![
            zero.clone(),
            zero.clone(),
            zero.clone(),
            one.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
        ],
        vec![
            zero.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
            one.clone(),
            zero.clone(),
            zero.clone(),
        ],
        vec![
            zero.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
            one.clone(),
        ],
    ];
    let c = vec![
        vec![
            zero.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
            one.clone(),
            zero.clone(),
        ],
        vec![
            zero.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
            one.clone(),
        ],
        vec![
            zero.clone(),
            one.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
            zero.clone(),
        ],
    ];

    let v1 = field.element(12); // x * y == 12
    let v2 = field.element(10); // z * u == 10
    let out = field.element(120); // v1 * v2 == 120

    let expected_witness = vec![one, out, x, y, z, u, v1, v2];

    assert_eq!(r1cs.a, a);
    assert_eq!(r1cs.b, b);
    assert_eq!(r1cs.c, c);

    assert_eq!(r1cs.witness, expected_witness);
}

#[test]
fn qap_works() {
    // Field must be larger than
    let field = CircuitField(420);
}
