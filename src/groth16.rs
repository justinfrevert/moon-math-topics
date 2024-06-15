use crate::{QAP, r1cs::R1CS, circuit_field::CircuitFieldElement};
use blstrs::{G1Projective, G2Projective, Scalar };
use group::{self, Group};
use rand::rngs::OsRng;

pub struct SimulationTrapdoor<FE>(FE, FE, FE, FE, FE);

// pub struct CRS(
//     ((G1Projective, G1Projective, G1Projective, Vec<G1Projective>), G1Projective, G1Projective, G1Projective),
//     (G2Projective, G2Projective, G2Projective, Vec<G2Projective>)
// );

// pub struct CRS_G1 {
//     trapdoor_values: G1Projective
// }

pub struct TrapDoorValues {
    alpha: G1Projective,
    beta: G1Projective, 
    delta: G1Projective,
    taus: Vec<G1Projective>
}

impl TrapDoorValues {
    fn new(
        alpha: G1Projective,
        beta: G1Projective, 
        delta: G1Projective,
        taus: Vec<G1Projective>
    ) -> Self {
        TrapDoorValues {
            alpha,
            beta, 
            delta,
            taus
        }
    }
}

pub struct CRS {
    g_1_trapdoor_values: TrapDoorValues,
    /// Accumulated calculations over each column(the upper right side of the definition)
    g_1_columns_values: Vec<G1Projective>,
    /// Accumulated calculations over each row(the lower left side of the definition)
    g_1_rows_values: Vec<G1Projective>,
    /// Accumulated calculations against the target polynomial(lower right side of the definition)
    g_1_target_polynomial_values: Vec<G1Projective>,
    g_2_values: (G2Projective, G2Projective, G2Projective, Vec<G2Projective>)
}

impl CRS {
    fn new(
        g_1_trapdoor_values: TrapDoorValues,
        g_1_columns_values: Vec<G1Projective>,
        g_1_rows_values: Vec<G1Projective>,
        g_1_target_polynomial_values: Vec<G1Projective>,
        g_2_values: (G2Projective,G2Projective,G2Projective,Vec<G2Projective>)
    ) -> Self {
        CRS {
            g_1_trapdoor_values, g_1_columns_values, g_1_rows_values, g_1_target_polynomial_values, g_2_values
        }
    }
}

// pub struct Groth16<G>(PhantomData<G>);
pub struct Groth16();

// impl<FE, F: Field<FE>> Groth16<FE, F> {
// impl<G: Group + std::ops::Mul<Scalar, Output = G>> Groth16<G> {
impl Groth16 {
    #[cfg(feature = "proving")]
    // note: n = num rows, m = num columns
    pub fn setup(qap: QAP<Scalar>, n: usize, m: usize) -> CRS {
        use blstrs::G2Projective;
        use group::ff::Field;
        use rand::RngCore;

        let rng = OsRng;

        // let alpha = Scalar::from(rng.next_u64());
        // let beta = Scalar::from(rng.next_u64());
        // let gamma = Scalar::from(rng.next_u64());
        // let delta = Scalar::from(rng.next_u64());
        // let tau = Scalar::from(rng.next_u64());
        // Chosen values for comparison with text
        let alpha = Scalar::from(6);
        let beta = Scalar::from(5);
        let gamma = Scalar::from(4);
        let delta = Scalar::from(3);
        let tau_plain = 2;
        let tau = Scalar::from(tau_plain);

        let simulation_trapdoor = SimulationTrapdoor(alpha, beta, gamma, delta, tau);

        // Upper left side
        let g_alpha = G1Projective::generator() * alpha;
        let g_beta = G1Projective::generator() * beta;
        let g_delta = G1Projective::generator() * delta;

        let mut powers_of_tau = vec![];

        // deg(T )−1
        for j in 0..qap.target_polynomial.0.len() - 1 {
            let power_of_tau = Scalar::from(tau_plain.pow(j.try_into().unwrap()));
            let val = G1Projective::generator() * power_of_tau;
            powers_of_tau.push(val);
        }

        // let upper_left_side = (g_alpha, g_beta, g_delta, powers_of_tau);
        let upper_left_side = TrapDoorValues::new(g_alpha, g_beta, g_delta, powers_of_tau);

        // upper right side
        // let mut upper_right_side = Scalar::ZERO;
        let mut upper_right_side = vec![];
        for j in 0..n {
            let a_tau = qap.a[j].evaluate(tau);
            let b_tau = qap.b[j].evaluate(tau);
            let c_tau = qap.c[j].evaluate(tau);

            let numerator = (beta * a_tau) + (alpha * b_tau) + c_tau;

            let result = numerator * gamma.invert().unwrap();
            let result = G1Projective::generator() * result;
            // upper_right_side = upper_right_side + ans;
            upper_right_side.push(result)
        }
        // let upper_right_side = G1Projective::generator() * upper_right_side;

        // lower left side
        // let mut lower_left_side = Scalar::ZERO;
        let mut lower_left_side = vec![];
        for j in 1..m {
            let a_tau = qap.a[j + n].evaluate(tau);
            let b_tau = qap.b[j + n].evaluate(tau);
            let c_tau = qap.b[j + n].evaluate(tau);

            let numerator = (beta * a_tau) + (alpha * b_tau) + c_tau;

            let result = numerator * delta.invert().unwrap();
            lower_left_side.push(G1Projective::generator() * result);
            // let lower_left_side = G1Projective::generator() * lower_left_side;
        }
        // let lower_left_side = G1Projective::generator() * lower_left_side;

        // lower right side
        let mut lower_right_side = vec![];
        for i in 0..qap.target_polynomial.0.len() - 2 {
            let numerator = tau * Scalar::from(i as u64) * qap.target_polynomial.evaluate(tau);
            let result  = numerator * delta.invert().unwrap();
            lower_right_side.push(G1Projective::generator() * result);
        }

        // let lower_right_side = G1Projective::generator() * lower_right_side;

        // let crs_g1 = (upper_left_side, upper_right_side, lower_left_side, lower_right_side);

        // G_2
        let g_2_beta = G2Projective::generator() * beta;
        let g_2_gamma = G2Projective::generator() * gamma;
        let g_2_delta = G2Projective::generator() * delta;
        let mut g_2_powers_of_tau = vec![];
        // deg(T )−1
        for j in 0..qap.target_polynomial.0.len() - 1 {
            let power_of_tau = Scalar::from(tau_plain.pow(j.try_into().unwrap()));
            let val = G2Projective::generator() * power_of_tau;
            g_2_powers_of_tau.push(val);
        }

        let crs_g2 = (g_2_beta, g_2_gamma, g_2_delta, g_2_powers_of_tau);

        // CRS::new(crs_g1, crs_g2)
        CRS::new(
            upper_left_side, 
            upper_right_side, 
            lower_left_side, 
            lower_right_side, 
            crs_g2
        )
        // (crs_g1, crs_g2)
    }

    #[cfg(feature = "proving")]
    fn prove(r1cs_instance: R1CS<Scalar>, qap_instance: QAP<Scalar>, crs: CRS) {
        // use crate::r1cs::R1CS;

        let mut g_1_w = G1Projective::identity();
        for (idx, w) in r1cs_instance.witness.into_iter().enumerate() {
            // TODO: check if its the correct idx according to the spec
            g_1_w += crs.g_1_rows_values[idx] * w;
        }

        // let mut g_1_a = G1Projective::identity();

        for (idx, row) in crs.g_1_rows_values.iter().enumerate() {

            let tau_at_idx = crs.g_1_trapdoor_values.taus[idx];
            // qap_instance.a[idx].evaluate(tau_at_idx);
            qap_instance.a[idx].evaluate_commitment(tau_at_idx);

        }


    //     let g_w = ()
    }
    fn verify() -> bool {
        false
    }
}

#[test]
fn cancel_test() {

}