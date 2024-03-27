

// TODOS in order to implement R1CS from Alebraic circuits:
// Implement boolean, or equality check constraint

#[derive(Debug)]
pub struct R1CS<F> {
    pub a: Vec<Vec<F>>,
    pub b: Vec<Vec<F>>,
    pub c: Vec<Vec<F>>,
    pub witness: Vec<F>
}

impl <F> R1CS<F> {
    pub fn new(
        a: Vec<Vec<F>>,
        b: Vec<Vec<F>>,
        c: Vec<Vec<F>>,
        witness: Vec<F>
    ) -> Self {
        R1CS { a, b, c, witness }
    }
}