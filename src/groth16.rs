use std::marker::PhantomData;

use crate::circuit_field::Field;

pub struct SimulationTrapdoor<FE>(FE, FE, FE, FE, FE);

pub struct Groth16<FE, F>(PhantomData<(FE, F)>);

impl <FE, F: Field<FE>>Groth16<FE, F> {
    #[cfg(feature = "proving")]
    pub fn setup(field: F) -> SimulationTrapdoor<FE> {
        // α , β ,γ , δ and τ 
        let a: FE = field.random_element();
        let b: FE = field.random_element();
        let c: FE = field.random_element();
        let d: FE = field.random_element();
        let e: FE = field.random_element();

        SimulationTrapdoor(a, b, c, d, e)
    }

    #[cfg(feature = "proving")]
    fn prove() {}
    fn verify() -> bool {
        false
    }
}