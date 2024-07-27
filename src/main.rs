pub mod circuit;
pub mod circuit_field;
pub mod polynomial;
pub mod qap;
pub mod r1cs;
mod groth16;

use blstrs::{pairing, G1Projective, G2Projective, Scalar};
use group::{self, Curve, Group};

use circuit_field::FieldElementZero;
use crypto_bigint::U512;
use itertools::izip;
use groth16::Groth16;
use group::ff::{Field, PrimeField};
use rand::{thread_rng, Rng};

use crate::{circuit::Node, circuit_field::{CircuitField, CircuitFieldElement}};
use circuit::{Circuit, CircuitNoFieldInfo, Operation::*};
use qap::QAP;

fn main() {
    let w_0 = Scalar::from(1);
    let w_1 = Scalar::from(2);

    let instructions = vec![
        Node::constant(w_0.clone()),
        Node::constant(w_1.clone()),
        Node::operation(Multiply, 0, 1), 
    ];

    let circuit: CircuitNoFieldInfo<Scalar> = CircuitNoFieldInfo::new(instructions);
    let (_, r1cs_instance) = circuit.calculate_with_trace();

    let rng = thread_rng();

    let qap_instance = QAP::new_with_scalars(r1cs_instance.clone(), rng).unwrap();

    let num_rows = r1cs_instance.clone().a.len();
    let num_columns = r1cs_instance.a[0].len();

    let groth16_params = Groth16::setup(qap_instance.clone(), num_rows, num_columns);

    let proof = Groth16::prove(qap_instance, groth16_params.clone());
    let public_input = vec![Scalar::from(1)];

    let verification_result = Groth16::verify(public_input, proof, groth16_params);
    println!("Verification was succesful? {:?}", verification_result);
    assert!(verification_result);
    // let (_, r1cs) = circuit.calculate_with_trace();

    // let qap = QAP::new(r1cs, field.clone()).unwrap();
    // // println!("QAP: {:?}", qap);
    // qap.verify(field);

    // let num_rows = 5;
    // let num_columns = 5;

    // let groth16_params = Groth16::setup(qap, num_rows, num_columns);

}

#[test]
// Testing basic circuit which was written for learning only
fn three_factors_example_works() {
    let field = CircuitField(U512::from_u32(29));
    let number = field.element(U512::from_u32(2));
    let x = field.element(U512::from_u32(3));
    let y = field.element(U512::from_u32(4));
    let z = field.element(U512::from_u32(5));

    let instructions = vec![
        Node::constant(x),
        Node::constant(y),
        Node::operation(Multiply, 0, 1),
        Node::constant(z),
        Node::operation(Multiply, 2, 3),
    ];

    let c = Circuit::new(instructions, field);
    assert_eq!(c.calculate(), Some(number));
}

#[test]
fn addition() {
    let field = CircuitField(U512::from_u32(10));
    let number = field.element(U512::from_u32(9));
    let instructions = vec![
        Node::constant(field.element(U512::from_u32(4))),
        Node::constant(field.element(U512::from_u32(5))),
        Node::operation(Add, 0, 1),
    ];

    let mut c = Circuit::new(instructions, field);
    assert_eq!(c.calculate(), Some(number));
}

#[test]
// Testing circuit implementation written to support realistic types - Scalars of the field related to the elliptic curve
fn addition_field_scalar() {
    let number = Scalar::from(9);

    let instructions = vec![
        Node::constant(Scalar::from(4)),
        Node::constant(Scalar::from(5)),
        Node::operation(Add, 0, 1),
    ];

    let mut c = CircuitNoFieldInfo::new(instructions);
    let (result, _) = c.calculate_with_trace();
    assert_eq!(result, Some(number));
}

#[test]
fn multiplication() {
    let field = CircuitField(U512::from_u32(5));
    let answer = field.element(U512::from_u32(0));

    let instructions = vec![
        Node::constant(field.element(U512::from_u32(4))),
        Node::constant(field.element(U512::from_u32(5))),
        Node::operation(Multiply, 0, 1),
    ];

    let mut c = Circuit::new(instructions, field);
    assert_eq!(c.calculate(), Some(answer));
}

#[test]
fn multiplication_field_scalar() {
    let answer = Scalar::from(20);

    let instructions = vec![
        Node::constant(Scalar::from(4)),
        Node::constant(Scalar::from(5)),
        Node::operation(Multiply, 0, 1),
    ];

    let c = CircuitNoFieldInfo::new(instructions);
    let (result, _) = c.calculate_with_trace();
    assert_eq!(result, Some(answer));
}

#[test]
fn multiplication2() {
    let field = CircuitField(U512::from_u32(69));
    let answer = field.element(U512::from_u32(18));

    let instructions = vec![
        Node::constant(field.element(U512::from_u32(400))),
        Node::constant(field.element(U512::from_u32(531))),
        Node::operation(Multiply, 0, 1),
    ];

    let mut c = Circuit::new(instructions, field);
    assert_eq!(c.calculate(), Some(answer));
}

#[test]
fn addition_and_multiplication() {
    let field = CircuitField(U512::from_u32(421));

    let number = field.element(U512::from_u32(176));

    let instructions = vec![
        Node::constant(field.element(U512::from_u32(42))),
        Node::constant(field.element(U512::from_u32(999))),
        Node::operation(Add, 0, 1),
        Node::constant(field.element(U512::from_u32(3))),
        Node::operation(Multiply, 2, 3),
    ];

    let mut c = Circuit::new(instructions, field);
    assert_eq!(c.calculate(), Some(number));
}

#[test]
fn addition_and_multiplication_field_scalar() {
    let number = Scalar::from(3123);

    let instructions = vec![
        Node::constant(Scalar::from(42)),
        Node::constant(Scalar::from(999)),
        Node::operation(Add, 0, 1),
        Node::constant(Scalar::from(3)),
        Node::operation(Multiply, 2, 3),
    ];

    let mut c = CircuitNoFieldInfo::new(instructions);
    let (result, _) = c.calculate_with_trace();
    assert_eq!(result, Some(number));
}

#[test]
fn with_generic_types() {
    let field = CircuitField(U512::from_u32(13));
    let x = field.element(U512::from_u32(3));
    let number = field.element(U512::from_u32(1));

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
fn squares() {
    let field = CircuitField(U512::from_u32(37));
    let number = field.element(U512::from_u32(36));

    let instructions = vec![
        Node::constant(field.element(U512::from_u32(2))),
        Node::constant(field.element(U512::from_u32(3))),
        Node::operation(Multiply, 0, 0), // 2 ^2
        Node::operation(Multiply, 1, 1), // 3^2
        Node::operation(Multiply, 2, 3), // 2^2 * 3^2
    ];
    let c = Circuit::new(instructions, field);
    assert_eq!(c.calculate(), Some(number));
}

#[test]
fn squares_with_modulus() {
    let field = CircuitField(U512::from_u32(9));
    let number = field.element(U512::from_u32(0));

    let instructions = vec![
        Node::constant(field.element(U512::from_u32(2))),
        Node::constant(field.element(U512::from_u32(3))),
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
// ) -> Circuit<CircuitFieldElement, CircuitField> {
) -> Circuit<CircuitField> {
    let instructions: Vec<Node<circuit_field::CircuitFieldElement>> = vec![
        Node::constant(field.element(U512::from_u32(1))),  // 0
        Node::constant(field.element(U512::from_u32(8))),  // 1
        Node::constant(x),                 // 2
        Node::constant(y),                 // 3
        Node::constant(field.element(U512::from_u32(10))), // 4
        Node::constant(field.element(U512::from_u32(12))), // 5
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
    let field = CircuitField(U512::from_u32(13));
    // 0 = 1 + 8 · x^2 · y^2 + 10 · x^2 + 12y^2
    let expected_result = field.element(U512::from_u32(0));

    // Some points to test: (1,2), (1, 11), (4, 0), (5,2) (5,11), (6,5), (6,8), ...... (12,8)
    let circuit_1_2 = point_on_curve_circuit(field.element(U512::from_u32(1)), field.element(U512::from_u32(2)), field.clone());
    let circuit_1_11 = point_on_curve_circuit(field.element(U512::from_u32(1)), field.element(U512::from_u32(11)), field.clone());
    // Doesn't work:
    let circuit_4_0 = point_on_curve_circuit(field.element(U512::from_u32(4)), field.element(U512::from_u32(0)), field.clone());

    assert_eq!(circuit_1_2.calculate(), Some(expected_result.clone()));
    assert_eq!(circuit_1_11.calculate(), Some(expected_result.clone()));
}

#[test]
// Example from rareskills book
fn simple_circuit() {
    let field = CircuitField(U512::from_u32(13));

    let x = field.element(U512::from_u32(2));
    let y = field.element(U512::from_u32(3));

    // x^2 y
    let instructions = vec![
        Node::constant(x),
        Node::constant(y),
        Node::operation(Multiply, 0, 0), // x^2
        Node::operation(Multiply, 1, 2), // (x^2) * y
    ];
    let c = Circuit::new(instructions, field.clone());

    assert_eq!(c.calculate(), Some(field.element(U512::from_u32(12))));
}

#[test]
// Example from rareskills book
fn simple_circuit_with_example_witness() {
    let field = CircuitField(U512::from_u32(13));

    let x = field.element(U512::from_u32(41));
    let y = field.element(U512::from_u32(103));

    // x^2 y
    let instructions = vec![
        Node::constant(x),
        Node::constant(y),
        Node::operation(Multiply, 0, 1),
    ];

    let circuit = Circuit::new(instructions, field.clone());
    let zero = field.element(U512::from_u32(0));
    let one = field.element(U512::from_u32(1));

    let (result, r1cs) = circuit.calculate_with_trace();

    assert_eq!(result, Some(field.element(U512::from_u32(11))));

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
        field.element(U512::from_u32(1)),
        field.element(U512::from_u32(11)),
        field.element(U512::from_u32(41)),
        field.element(U512::from_u32(103)),
    ];

    assert_eq!(r1cs.a, a);
    assert_eq!(r1cs.b, b);
    assert_eq!(r1cs.c, c);
    assert_eq!(r1cs.witness, expected_witness);
    // assert_eq!(result, Some(field.element(U512::from_u32(4223))));
}

#[test]
fn subtract_circuit() {
    assert!(false)
}


#[test]
fn boolean_circuit() {
    let field = CircuitField(U512::from_u32(13));

    let x = field.element(U512::from_u32(41));
    let y = field.element(U512::from_u32(103));

    // x^2 y
    let instructions = vec![
        Node::constant(x),
        Node::constant(y),
        Node::operation(Multiply, 0, 1),
    ];

    let circuit = Circuit::new(instructions, field.clone());
    let zero = field.element(U512::from_u32(0));
    let one = field.element(U512::from_u32(1));

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
        field.element(U512::from_u32(1)),
        field.element(U512::from_u32(4223)),
        field.element(U512::from_u32(41)),
        field.element(U512::from_u32(103)),
    ];

    assert_eq!(r1cs.a, a);
    assert_eq!(r1cs.b, b);
    assert_eq!(r1cs.c, c);
    assert_eq!(r1cs.witness, expected_witness);
    assert_eq!(result, Some(field.element(U512::from_u32(4223))));
}

// Example from rareskills book x * y * z * u
#[test]
fn example_r1cs_more_terms() {
    let field = CircuitField(U512::from_u32(42001));
    let zero = field.element(U512::from_u32(0));
    let one = field.element(U512::from_u32(1));
    let x = field.element(U512::from_u32(3).clone());
    let y = field.element(U512::from_u32(4).clone());
    let z = field.element(U512::from_u32(2).clone());
    let u = field.element(U512::from_u32(5).clone());

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

    let v1 = field.element(U512::from_u32(12)); // x * y == 12
    let v2 = field.element(U512::from_u32(10)); // z * u == 10
    let out = field.element(U512::from_u32(120)); // v1 * v2 == 120

    let expected_witness = vec![one, out, x, y, z, u, v1, v2];

    assert_eq!(r1cs.a, a);
    assert_eq!(r1cs.b, b);
    assert_eq!(r1cs.c, c);

    assert_eq!(r1cs.witness, expected_witness);
}

#[test]
fn qap_works() {
    let field = CircuitField(U512::from_u32(13));
    let w_0 = field.element(U512::from_u32(3).clone());
    let w_1 = field.element(U512::from_u32(4).clone());
    let w_2 = field.element(U512::from_u32(2).clone());

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
    let qap = QAP::new(r1cs, field.clone()).unwrap();
    assert!(qap.verify(field));
}

#[test]
fn qap_does_not_verify_false_proof() {
    let field = CircuitField(U512::from_u32(13));
    let w_0 = field.element(U512::from_u32(3).clone());
    let w_1 = field.element(U512::from_u32(4).clone());
    let w_2 = field.element(U512::from_u32(2).clone());

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
    let (_, mut r1cs) = circuit.calculate_with_trace();

    r1cs.b[0][1] = field.element(U512::from_u32(1));
    r1cs.b[0][2] = field.element(U512::from_u32(1));
    r1cs.b[0][3] = field.element(U512::from_u32(1));

    let qap = QAP::new(r1cs, field.clone()).unwrap();

    assert_eq!(qap.verify(field), false);
}

#[test]
// An intuition for validating a proof which is similar to the QAP's verification, but translated almost directly to the pairing-based equivalent
fn qap_committed_polynomial_evaluation_proof() {
    env_logger::init();

    let w_0 = Scalar::from(1);
    let w_1 = Scalar::from(2);

    let instructions = vec![
        Node::constant(w_0.clone()),
        Node::constant(w_1.clone()),
        Node::operation(Multiply, 0, 1), 
    ];

    let circuit: CircuitNoFieldInfo<Scalar> = CircuitNoFieldInfo::new(instructions);
    let (_, r1cs_instance) = circuit.calculate_with_trace();

    let rng = thread_rng();

    let qap_instance = QAP::new_with_scalars(r1cs_instance.clone(), rng).unwrap();
    let num_rows = r1cs_instance.clone().a.len();
    let num_columns = r1cs_instance.a[0].len();

    let m = num_rows;
    let n = num_columns;

    let crs = {
        let tau_plain = 3;
        let tau = Scalar::from(tau_plain);
        let mut powers_of_tau = vec![];

        // deg(T)−1
        for j in 0..qap_instance.target_polynomial.0.len() - 1 {
            let power_of_tau = Scalar::from(tau_plain.pow(j.try_into().unwrap()));
            let val = G1Projective::generator() * power_of_tau;
            powers_of_tau.push(val);
        }

        // upper right side
        let mut upper_right_side = vec![];
        for j in 0..n {
            let a_tau = qap_instance.a[j].evaluate(tau);
            let b_tau = qap_instance.b[j].evaluate(tau);
            let c_tau = qap_instance.c[j].evaluate(tau);

            let numerator = a_tau + b_tau + c_tau;
            upper_right_side.push(G1Projective::generator() * numerator)
        }

        // lower left side
        let mut lower_left_side = vec![];
        for j in 1..m {
            // Text says j + n... is that right? 
            // Instead, we will just do j for now? 
            let a_tau = qap_instance.a[j].evaluate(tau);
            let b_tau = qap_instance.b[j].evaluate(tau);
            let c_tau = qap_instance.c[j].evaluate(tau);

            let numerator = a_tau + b_tau + c_tau;
            // lower_left_side.push(G1Projective::generator() * numerator);
            lower_left_side.push(G1Projective::generator() * numerator);
        }

        // lower right side
        let mut lower_right_side = vec![];
        for i in 0..qap_instance.target_polynomial.0.len() - 2 {
            let numerator = tau * Scalar::from(i as u64) * qap_instance.target_polynomial.evaluate(tau);
            lower_right_side.push(G1Projective::generator() * numerator);
        }

        // G_2
        let mut g_2_powers_of_tau = vec![];
        // deg(T)−1
        for j in 0..qap_instance.target_polynomial.0.len() - 1 {
            let power_of_tau = Scalar::from(tau_plain.pow(j.try_into().unwrap()));
            let val = G2Projective::generator() * power_of_tau;
            g_2_powers_of_tau.push(val);
        };
        (powers_of_tau, g_2_powers_of_tau)
    };

    assert_eq!(num_rows, 1);
    assert_eq!(num_columns, 4);

    let (a, b, c) = {
        let mut g_1_a = G1Projective::identity();
        let mut g_2_b = G2Projective::identity();
        let mut g_1_w = G1Projective::identity();    
        for (a_poly, b_poly, c_poly) in izip!(qap_instance.a, qap_instance.b, qap_instance.c) {
            g_1_a += a_poly.evaluate_polynomial_commitment(&crs.0);
            g_2_b += b_poly.evaluate_polynomial_commitment_g2(&crs.1);
            g_1_w += c_poly.evaluate_polynomial_commitment(&crs.0);
        }
        (g_1_a, g_2_b, g_1_w)
    };

    let verification_result = pairing(&a.to_affine(), &b.to_affine()) == pairing(&c.to_affine(), &G2Projective::generator().to_affine());
    assert!(verification_result);
}

#[test]
fn qap_committed_polynomial_evaluation_proof_alpha() {
    env_logger::init();

    let w_0 = Scalar::from(1);
    let w_1 = Scalar::from(2);

    let instructions = vec![
        Node::constant(w_0.clone()),
        Node::constant(w_1.clone()),
        Node::operation(Multiply, 0, 1), 
    ];

    let circuit: CircuitNoFieldInfo<Scalar> = CircuitNoFieldInfo::new(instructions);
    let (_, r1cs_instance) = circuit.calculate_with_trace();

    let rng = thread_rng();

    let prover_random_r_term = Scalar::from(6);
    let prover_random_t_term = Scalar::from(7);

    let qap_instance = QAP::new_with_scalars(r1cs_instance.clone(), rng).unwrap();
    let num_rows = r1cs_instance.clone().a.len();
    let num_columns = r1cs_instance.a[0].len();

    let m = num_rows;
    let n = num_columns;
    let alpha = Scalar::from(7);
    let beta = Scalar::from(8);
    let gamma = Scalar::from(9);
    let delta = Scalar::from(1);

    let crs_trap_doors: ((G1Projective, G1Projective, G1Projective), (G2Projective, G2Projective, G2Projective)) = {
        ((
            G1Projective::generator() * alpha,
            G1Projective::generator() * beta,
            G1Projective::generator() * delta,
        ),
        (
            G2Projective::generator() * beta,
            G2Projective::generator() * gamma,
            G2Projective::generator() * delta,
        ))
    };

    let crs: (Vec<G1Projective>, Vec<G2Projective>, Vec<G1Projective>, Vec<G1Projective>) = {
        let tau_plain = 3;
        let tau = Scalar::from(tau_plain);
        let mut powers_of_tau = vec![];

        // deg(T)−1
        for j in 0..qap_instance.target_polynomial.0.len() - 1 {
            let power_of_tau = Scalar::from(tau_plain.pow(j.try_into().unwrap()));
            let val = G1Projective::generator() * power_of_tau;
            powers_of_tau.push(val);
        }

        // upper right side
        let mut upper_right_side = vec![];
        for j in 0..n {
            let a_tau = qap_instance.a[j].evaluate(tau);
            let b_tau = qap_instance.b[j].evaluate(tau);
            let c_tau = qap_instance.c[j].evaluate(tau);

            let numerator = a_tau + b_tau + c_tau;
            // let result = numerator * gamma.invert().unwrap();
            upper_right_side.push(G1Projective::generator() * numerator)
        }

        // lower left side
        let mut lower_left_side = vec![];
        for j in 1..m {
            let a_tau = qap_instance.a[j].evaluate(tau);
            let b_tau = qap_instance.b[j].evaluate(tau);
            let c_tau = qap_instance.c[j].evaluate(tau);

            let numerator = a_tau + b_tau + c_tau;
            // let result = numerator * delta.invert().unwrap();
            lower_left_side.push(G1Projective::generator() * numerator);
        }

        // lower right side
        let mut lower_right_side = vec![];
        for i in 0..qap_instance.target_polynomial.0.len() - 1 {
            let numerator = tau * Scalar::from(i as u64) * qap_instance.target_polynomial.evaluate(tau);
            lower_right_side.push(G1Projective::generator() * numerator);
        }

        // G_2
        let mut g_2_powers_of_tau = vec![];
        // deg(T)−1
        for j in 0..qap_instance.target_polynomial.0.len() - 1 {
            let power_of_tau = Scalar::from(tau_plain.pow(j.try_into().unwrap()));
            let val = G2Projective::generator() * power_of_tau;
            g_2_powers_of_tau.push(val);
        };
        (powers_of_tau, g_2_powers_of_tau, lower_right_side, lower_left_side)
    };

    assert_eq!(num_rows, 1);
    assert_eq!(num_columns, 4);

    let (is_verified, h) = qap_instance.verify_with_data();
    assert!(is_verified);

    let (crs_g1, crs_g2) = crs_trap_doors;
    let (_, _, g1_delta) = crs_g1;
    let (_, _, g2_delta) = crs_g2;

    let (a, b_g_1, b, calculated_w) = {
        let mut g_1_a = G1Projective::identity();
        let mut g_1_b = G1Projective::identity();
        let mut g_2_b = G2Projective::identity();
        let mut g_1_w = G1Projective::identity();    

        for (a_poly, b_poly, c_poly) in izip!(qap_instance.a, qap_instance.b, qap_instance.c) {
            g_1_a += a_poly.evaluate_polynomial_commitment(&crs.0);
            g_1_b += b_poly.evaluate_polynomial_commitment(&crs.0);
            g_2_b += b_poly.evaluate_polynomial_commitment_g2(&crs.1);
            g_1_w += c_poly.evaluate_polynomial_commitment(&crs.0);
        }

        // Scaling by prover-side random terms
        // g_1_a += g1_delta * prover_random_r_term;
        // g_1_b += g1_delta * prover_random_t_term;
        // g_2_b += g2_delta * prover_random_t_term;

        (g_1_a, g_1_b, g_2_b, g_1_w)
    };

    // Collect collection of w group elements into one 
    // let w: G1Projective = crs.3.iter().sum();
    // assert_eq!(w, calculated_w);

    // Below works for now, but will need to incorporate A and B
    let c = h.evaluate_polynomial_commitment(&crs.2) + calculated_w;
    // let c = h.evaluate_polynomial_commitment(&crs.2) + w *- prover_random_r_term;

    // let c = h.evaluate_polynomial_commitment(&crs.2) + w + (g1_delta * -(prover_random_r_term * prover_random_t_term));
    // let c: G1Projective = h.evaluate_polynomial_commitment(&crs.2) + calculated_w + a + b_g_1 + (g1_delta * (-prover_random_r_term * prover_random_t_term));

    let verification_result = 
    pairing(&a.to_affine(), &b.to_affine()) == 
    pairing(&c.to_affine(), &G2Projective::generator().to_affine());

    // let verification_result = 
    // pairing(&a.to_affine(), &b.to_affine()) == 
    // pairing(&c.to_affine(), &g2_delta.to_affine());

    assert!(verification_result);
}