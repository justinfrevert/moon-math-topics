use crate::QAP;
use blstrs::{G1Projective, Scalar};
use group::{self, Group};
use rand::rngs::OsRng;

pub struct SimulationTrapdoor<FE>(FE, FE, FE, FE, FE);

pub struct CRS();

// pub struct Groth16<G>(PhantomData<G>);
pub struct Groth16();

// impl<FE, F: Field<FE>> Groth16<FE, F> {
// impl<G: Group + std::ops::Mul<Scalar, Output = G>> Groth16<G> {
impl Groth16 {
    #[cfg(feature = "proving")]
    // note: n = num rows, m = num columns
    pub fn setup(qap: QAP, n: usize, m: usize) -> CRS {
        // let alpha = G::random(OsRng);
        // let beta = G::random(OsRng);
        // let gamma = G::random(OsRng);
        // let delta = G::random(OsRng);
        // let tau = G::random(OsRng);

        use rand::RngCore;

        let rng = OsRng;

        // let alpha = Scalar::from(rng.next_u64());
        // let beta = Scalar::from(rng.next_u64());
        // let gamma = Scalar::from(rng.next_u64());
        // let delta = Scalar::from(rng.next_u64());
        // let tau = Scalar::from(rng.next_u64());
        // Chosen values for comparison with text
        let alpha = Scalar::from(5);
        let beta = Scalar::from(6);
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

        let upper_left_side = (g_alpha, g_beta, g_delta, powers_of_tau);

        // Get points on the right
        // for n_j in 0..n {
        //     let first = beta * qap.a[n_j].evaluate(tau);
        // }

        // upper right side
        // let mut upper_right_elements = vec![];

        // for j in 0..num_rows {
        //     let numerator = beta * qap.a[j] + (alpha * qap.b[j] * tau) + (qap.c[j] * tau);
        //     let exponent = numerator / delta;

        //     upper_right_elements.push(G::generator() * exponent);
        // }

        // let numerator = ;
        // let denominator = ;

        CRS()
    }

    #[cfg(feature = "proving")]
    fn prove() {}
    fn verify() -> bool {
        false
    }
}