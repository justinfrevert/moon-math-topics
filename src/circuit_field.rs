// use std::ops::{Add, Mul};
// pub struct CircuitField(u64);

// impl Add for CircuitField {
//     fn add(self, rhs: Self) -> Self::Output {
//         self
//     }
// }

use std::ops::AddAssign;
use std::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CircuitField(pub u64);

impl CircuitField {
    pub fn element(&self, value: u64) -> CircuitFieldElement {
        CircuitFieldElement {
            field: self.clone(),
            value: value % self.0
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CircuitFieldElement {
    pub value: u64,
    field: CircuitField,
}

impl CircuitFieldElement {
    pub fn new(value: u64, field: CircuitField) -> Self {
        CircuitFieldElement { value, field }
    }
}

impl Add for CircuitFieldElement {
    type Output = CircuitFieldElement;
    fn add(self, rhs: Self) -> Self::Output {
        let ans = (self.value + rhs.value) % self.field.clone().0;
        CircuitFieldElement::new(ans, self.field)
    }
}

impl AddAssign for CircuitFieldElement {
    fn add_assign(&mut self, rhs: Self) {
        // let result = (self.value + rhs.value) % self.field.0;
        self.value = (self.clone().value + rhs.value) % self.clone().field.0;
    }
}

impl Neg for CircuitFieldElement {
    type Output = CircuitFieldElement;

    fn neg(self) -> Self::Output {
        let result = (self.clone().field.0 - self.clone().value) % self.clone().field.clone().0;
        // let result = (self.field.0 - self.value) % self.field.0;
        CircuitFieldElement::new(result, self.clone().field)
    }
}

impl Mul for CircuitFieldElement {
    type Output = CircuitFieldElement;

    fn mul(self, rhs: Self) -> Self::Output {
        let result = (self.clone().value * rhs.value) % self.clone().field.0;
        CircuitFieldElement::new(result, self.field)
    }
}


#[test]
fn adds() {
    let field = CircuitField(u64::from(101_u32));

    let field_element_lower = CircuitFieldElement {
        value: u64::from(100_u32),
        field: field.clone(),
    };

    let field_element_higher = CircuitFieldElement {
        value: u64::from(2_u32),
        field: field.clone(),
    };

    assert_eq!(
        field_element_lower + field_element_higher,
        CircuitFieldElement::new(u64::from(1_u32), field)
    );
}

#[test]
fn example_field() {
    let field = CircuitField(u64::from(41_u32));

    let field_element_lower = CircuitFieldElement {
        value: u64::from(1_u32),
        field: field.clone(),
    };

    let field_element_higher = CircuitFieldElement {
        value: u64::from(40_u32),
        field: field.clone(),
    };

    assert_eq!(
        field_element_lower + field_element_higher,
        CircuitFieldElement::new(u64::from(0_u32), field)
    );
}

#[test]
fn low_field() {
    let field = CircuitField(u64::from(13_u32));

    let field_element_lower = CircuitFieldElement {
        value: u64::from(3_u32),
        field: field.clone(),
    };

    let field_element_higher = CircuitFieldElement {
        value: u64::from(10_u32),
        field: field.clone(),
    };

    assert_eq!(
        field_element_lower + field_element_higher,
        CircuitFieldElement::new(u64::from(0_u32), field)
    );
}
