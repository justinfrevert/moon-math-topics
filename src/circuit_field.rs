use std::cmp::Ordering;
use std::ops::AddAssign;
use std::ops::{Add, Mul, Neg};

#[derive(Clone, Debug, PartialEq, Eq)]
// Field with sole value of modulus
pub struct CircuitField(pub u64);

impl CircuitField {
    pub fn element(&self, value: u64) -> CircuitFieldElement {
        CircuitFieldElement {
            field: self.clone(),
            value: value % self.0
        }
    }
}

pub trait CicuitFieldInfo {
    fn modulus(&self) -> u64;
 }

 impl CicuitFieldInfo for CircuitField {
    fn modulus(&self) -> u64 {
        self.0
    }
 }

 pub trait FieldElementZero {
    fn zero(&self) -> Self;
 }

 impl FieldElementZero for CircuitFieldElement {
    fn zero(&self) -> Self {
        CircuitFieldElement { value: 0, field: self.clone().field }
    }
 }

 pub trait FieldElementOne {
    fn one(&self) -> Self;
 }

 impl FieldElementOne for CircuitFieldElement {
    fn one(&self) -> Self {
        CircuitFieldElement { value: 1, field: self.clone().field }
    }
 }

 // For testing only:
 impl FieldElementZero for u32 {
    fn zero(&self) -> Self {
        0_u32
    }
 }

 impl FieldElementOne for u32 {
    fn one(&self) -> Self {
        1_u32
    }
 }

 pub trait NewValue<T> {
    fn new(value: T, supporting_value: CircuitField) -> Self;
 }

 impl <T: Into<u64>> NewValue<T> for CircuitFieldElement {
    fn new(value: T, supporting_value: CircuitField) -> Self {
        CircuitFieldElement {
            value:  value.into(),
            field: supporting_value
        }
    }
 }

 impl <T: Into<u32>> NewValue<T> for u32 {
    fn new(value: T, _supporting_value: CircuitField) -> Self {
        value.into()
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

    pub fn circuit_field_zero(&self) -> Self {
        self.field.element(0)
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

impl Ord for CircuitFieldElement {
    fn cmp(&self, other: &Self) -> Ordering {
        self.value.cmp(&other.value)
    }
}

impl PartialOrd for CircuitFieldElement {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
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
