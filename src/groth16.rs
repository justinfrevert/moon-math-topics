use crate::QAP;
use blstrs::{pairing, G1Projective, G2Projective, Scalar };
use group::{self, Curve, Group};
use itertools::izip;

#[derive(Clone)]
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

#[derive(Clone)]
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

pub struct Groth16();

impl Groth16 {
    #[cfg(feature = "proving")]
    pub fn setup(qap: QAP<Scalar>, n: usize, m: usize) -> CRS {
        use group::ff::Field;

        // Chosen values for comparison with text. This is obviously a temporary measure 
        let alpha = Scalar::from(6);
        let beta = Scalar::from(5);
        let gamma = Scalar::from(4);
        let delta = Scalar::from(3);
        let tau_plain = 2;
        let tau = Scalar::from(tau_plain);

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

        let upper_left_side = TrapDoorValues::new(g_alpha, g_beta, g_delta, powers_of_tau);

        // upper right side
        let mut upper_right_side = vec![];
        for j in 0..n {
            let a_tau = qap.a[j].evaluate(tau);
            let b_tau = qap.b[j].evaluate(tau);
            let c_tau = qap.c[j].evaluate(tau);

            let numerator = (beta * a_tau) + (alpha * b_tau) + c_tau;

            let result = numerator * gamma.invert().unwrap();
            upper_right_side.push(G1Projective::generator() * result)
        }

        // lower left side
        // let mut lower_left_side = Scalar::ZERO;
        let mut lower_left_side = vec![];
        for j in 1..m {
            // Text says j + n... is that right? 
            // let a_tau = qap.a[j + n].evaluate(tau);
            // let b_tau = qap.b[j + n].evaluate(tau);
            // let c_tau = qap.c[j + n].evaluate(tau);
            // Instead, we will just do j for now? 
            let a_tau = qap.a[j].evaluate(tau);
            let b_tau = qap.b[j].evaluate(tau);
            let c_tau = qap.c[j].evaluate(tau);

            let numerator = (beta * a_tau) + (alpha * b_tau) + c_tau;

            let result = numerator * delta.invert().unwrap();
            lower_left_side.push(G1Projective::generator() * result);
        }

        // lower right side
        let mut lower_right_side = vec![];
        for i in 0..qap.target_polynomial.0.len() - 2 {
            let numerator = tau * Scalar::from(i as u64) * qap.target_polynomial.evaluate(tau);
            let result  = numerator * delta.invert().unwrap();
            lower_right_side.push(G1Projective::generator() * result);
        }

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
        CRS::new(
            upper_left_side, 
            upper_right_side, 
            lower_left_side, 
            lower_right_side, 
            crs_g2
        )
    }

    #[cfg(feature = "proving")]
    pub fn prove(qap_instance: QAP<Scalar>, crs: CRS) -> (G1Projective, G2Projective, G1Projective) {
        let mut g_1_a = G1Projective::identity();
        let mut g_2_b = G2Projective::identity();
        let mut g_1_w = G1Projective::identity();

        for (a_poly, b_poly, c_poly) in izip!(qap_instance.a, qap_instance.b, qap_instance.c) {
            g_1_a += a_poly.evaluate_polynomial_commitment(&crs.g_1_trapdoor_values.taus);
            g_2_b += b_poly.evaluate_polynomial_commitment_g2(&crs.g_2_values.3);
            g_1_w += c_poly.evaluate_polynomial_commitment(&crs.g_1_trapdoor_values.taus);
        }

        g_1_w += qap_instance.target_polynomial.evaluate_polynomial_commitment(&crs.g_1_trapdoor_values.taus);
        (g_1_a, g_2_b, g_1_w)
    }

    pub fn verify(public_input: Vec<Scalar>, proof: (G1Projective, G2Projective, G1Projective), crs: CRS) -> bool {
        false
    }
}

#[test]
fn cancel_test() {

}