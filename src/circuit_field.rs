#[cfg(feature = "proving")]
use rand::Rng;
use std::cmp::Ordering;
use std::fmt::Display;
use std::ops::{AddAssign, Sub, SubAssign};
use std::ops::{Add, Mul, Neg};

#[derive(Clone, Debug, PartialEq, Eq)]
// Field with sole value of modulus
pub struct CircuitField(pub u64);

impl CircuitField {
    pub fn element(&self, value: u64) -> CircuitFieldElement {
        CircuitFieldElement {
            field: self.clone(),
            value: value % self.0,
        }
    }
    #[cfg(feature = "proving")]
    pub fn random_element(&self) -> CircuitFieldElement {
        let mut rng = rand::thread_rng();
        let value: u64 = rng.gen_range(0..self.0);
        self.element(value)
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
    fn is_zero(&self) -> bool;
}

impl FieldElementZero for CircuitFieldElement {
    fn zero(&self) -> Self {
        CircuitFieldElement {
            value: 0,
            field: self.clone().field,
        }
    }
    
    fn is_zero(&self) -> bool {
        self.value == 0_u64
    }
}

pub trait FieldElementOne {
    fn one(&self) -> Self;
}

impl FieldElementOne for CircuitFieldElement {
    fn one(&self) -> Self {
        CircuitFieldElement {
            value: 1,
            field: self.clone().field,
        }
    }
}

// For testing only:
impl FieldElementZero for u32 {
    fn zero(&self) -> Self {
        0_u32
    }
    fn is_zero(&self) -> bool {
        self == &0_u32
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

impl<T: Into<u64>> NewValue<T> for CircuitFieldElement {
    fn new(value: T, supporting_value: CircuitField) -> Self {
        CircuitFieldElement {
            value: value.into(),
            field: supporting_value,
        }
    }
}

impl<T: Into<u32>> NewValue<T> for u32 {
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

    pub fn invert(&self) -> Self {
        // TODO: Replace field modulus type with something which would allow larger exponentiations
        let small_modulus: u32 = self.field.0.try_into().unwrap();
        let value = self.value.pow(small_modulus - 2_u32) % self.field.0;

        CircuitFieldElement {
            field: self.field.clone(),
            value
        }
    }
}

impl Add for CircuitFieldElement {
    type Output = CircuitFieldElement;
    fn add(self, rhs: Self) -> Self::Output {
        let ans = (self.value + rhs.value) % self.field.clone().0;
        CircuitFieldElement::new(ans, self.field)
    }
}

// impl AddAssign for CircuitFieldElement {
//     fn add_assign(&mut self, rhs: Self) {
//         let ans = (self.value + rhs.value) % self.field.clone().0;
//         self = ans;
//     }
// }

impl Sub for CircuitFieldElement {
    type Output = CircuitFieldElement;
    fn sub(self, rhs: Self) -> Self::Output {
        self + -rhs
    }
}

impl SubAssign<&CircuitFieldElement> for CircuitFieldElement {
    fn sub_assign(&mut self, rhs: &CircuitFieldElement) {
        // self = self + -rhs

        if self.value < rhs.value {
            self.value = (self.field.0 - (rhs.value - self.value)) % self.field.0;
        } else {
            self.value = (self.value - rhs.value) % self.field.0;
        }
    }
}

// impl SubAssign for CircuitFieldElement {
//     fn sub_assign(&mut self, rhs: Self) {
//         // self = self + -rhs;

//         let mut subtraction = self.clone() + -rhs;
//         self = &subtraction;
//     }
// }

impl<'a, 'b> Add<&'b CircuitFieldElement> for &'a CircuitFieldElement {
    type Output = CircuitFieldElement;

    fn add(self, rhs: &'b CircuitFieldElement) -> Self::Output {
        assert_eq!(
            self.field.0, rhs.field.0,
            "Fields must be the same for addition"
        );
        let ans = (self.value + rhs.value) % self.field.0;
        CircuitFieldElement::new(ans, self.field.clone())
    }
}

impl AddAssign for CircuitFieldElement {
    fn add_assign(&mut self, rhs: Self) {
        self.value = (self.clone().value + rhs.value) % self.clone().field.0;
    }
}

impl Neg for CircuitFieldElement {
    type Output = CircuitFieldElement;

    fn neg(self) -> Self::Output {
        let result = (self.clone().field.0 - self.clone().value) % self.clone().field.clone().0;
        CircuitFieldElement::new(result, self.clone().field)
    }
}

impl<'a> Neg for &'a CircuitFieldElement {
    type Output = CircuitFieldElement;

    fn neg(self) -> Self::Output {
        let result = (self.field.0 - self.value) % self.field.0;
        CircuitFieldElement::new(result, self.field.clone())
    }
}

impl Mul for CircuitFieldElement {
    type Output = CircuitFieldElement;

    fn mul(self, rhs: Self) -> Self::Output {
        let result = (self.clone().value * rhs.value) % self.clone().field.0;
        CircuitFieldElement::new(result, self.field)
    }
}

impl<'a, 'b> Mul<&'b CircuitFieldElement> for &'a CircuitFieldElement {
    type Output = CircuitFieldElement;

    fn mul(self, rhs: &'b CircuitFieldElement) -> Self::Output {
        assert_eq!(
            self.field.0, rhs.field.0,
            "Fields must be the same for multiplication"
        );
        let result = (self.value * rhs.value) % self.field.0;
        CircuitFieldElement::new(result, self.field.clone())
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

impl Display for CircuitFieldElement {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.value)
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
fn subtracts_including_negative_case() {
    let field = CircuitField(8);
    let a = field.element(4);
    let b = field.element(5);
    assert_eq!(a - b, field.element(7));
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
