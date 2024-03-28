use crate::{CircuitField, r1cs::R1CS};

// Definition for quadratic arithmetic program
pub struct QAP<F>{
    r1cs: R1CS<F>,
    field: CircuitField
}

impl<F> QAP<F> {
    fn new(r1cs: R1CS<F>, field: CircuitField) -> Result<Self, QAPError> {
        // TODO: remove unwrap
        if field.0 < r1cs.a.len().try_into().unwrap() {
            return Err(QAPError::FieldTooSmall);
        }
        Ok(QAP {
            r1cs,
            field
        })
    }
}

pub enum QAPError {
    FieldTooSmall
}