use blstrs::Scalar;
#[cfg(feature = "proving")]
use rand::Rng;
use std::cmp::Ordering;
use std::fmt::Display;
use std::hash::Hash;
use std::ops::{AddAssign, Sub, SubAssign};
use std::ops::{Add, Mul, Neg};
use group::ff::Field as FieldT;

use crypto_bigint::{AddMod, ConstZero, Constants, Integer, MulMod, NegMod, NonZero, Random, SubMod, Zero, U512};

#[derive(Clone, Debug, PartialEq, Eq)]
// Field with sole value of modulus
pub struct CircuitField(pub U512);

pub trait Field<F> {
    fn element(&self, value: U512) -> F;
    #[cfg(feature = "proving")]
    fn random_element(&self) -> F;
}

impl CircuitField {
    pub fn element(&self, value: U512) -> CircuitFieldElement {
        CircuitFieldElement {
            field: self.clone(),
            // value: value % self.0,
            value: value.add_mod(&U512::ZERO, &self.0)
        }
    }

    #[cfg(feature = "proving")]
    pub fn random_element(&self) -> CircuitFieldElement {
        use crypto_bigint::{NonZero, RandomMod};

        let mut rng = rand::thread_rng();
        // let value: U512 = rng.gen_range(0..self.0);
        // let value: U512 = rand_mod(0..self.0);
        let value: U512 = U512::random_mod(&mut rng, &NonZero::new(self.0).unwrap());
        self.element(value)
    }
}

impl Hash for CircuitField {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.0.hash(state);
    }
}


pub trait CicuitFieldInfo {
    fn modulus(&self) -> U512;
}

impl CicuitFieldInfo for CircuitField {
    fn modulus(&self) -> U512 {
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
            value: U512::ZERO,
            field: self.clone().field,
        }
    }
    
    fn is_zero(&self) -> bool {
        self.value == U512::ZERO
    }
}

pub trait FieldElementOne {
    fn one(&self) -> Self;
}

impl FieldElementOne for CircuitFieldElement {
    fn one(&self) -> Self {
        CircuitFieldElement {
            value: U512::ONE,
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

impl FieldElementOne for Scalar {
    fn one(&self) -> Self {
        Scalar::ONE
    }
}

impl FieldElementZero for Scalar {
    fn is_zero(&self) -> bool {
        *self == Scalar::ZERO
    }

    fn zero(&self) -> Self {
        Scalar::ZERO
    }
}

pub trait NewValue<T> {
    fn new(value: T, supporting_value: CircuitField) -> Self;
}

impl<T: Into<U512>> NewValue<T> for CircuitFieldElement {
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
    pub value: U512,
    // pub value: DynResidue<{U512::LIMBS}>,
    pub field: CircuitField,
}

impl CircuitFieldElement {
    pub fn new(value: U512, field: CircuitField) -> Self {
        // let dyn_params = DynResidueParams::new(&field.modulus());
        // let field_value: DynResidue<{U512::LIMBS}> = DynResidue::new(&value, dyn_params);
        let is_modded_val = value.add_mod(&U512::ZERO, &field.0);
        CircuitFieldElement { value: is_modded_val, field }
    }

    pub fn circuit_field_zero(&self) -> Self {
        self.field.element(U512::ZERO)
    }

    pub fn invert(&self) -> Self {
        // TODO: Replace field modulus type with something which would allow larger exponentiations
        // let small_modulus = self.field.0;

        let subtrahend = U512::from_u32(2);
        let exponent: CircuitFieldElement = CircuitFieldElement::new(self.field.0 - subtrahend, self.clone().field);

        // let value = self.value.pow(small_modulus - 2_u32) % self.field.0;

        // let (_, value) = self.pow(exponent) % self.field;
        let (_, rem) = self.pow(exponent).div_rem(&NonZero::new(self.field.0).unwrap());

        CircuitFieldElement {
            field: self.field.clone(),
            value: rem
        }
    }

    // pub fn pow(&self, rhs: Self) -> Self {
    //     // self.value.mul_mod(&rhs.value, &self.field.0);


    // }

    // pub fn pow(&self, n: Self) -> Self {

    //     // CircuitFieldElement::one(&self)

    //     let two = NonZero::new(U512::from_u32(2)).unwrap();
    //     match n.value {
    //         U512::ZERO => CircuitFieldElement::one(&self),
    //         U512::ONE => *self,


    //         // i if i % two == 0 => exp(x * x, n / 2),
    //         i if i.div_rem(&two).1 == U512::ZERO => {
    //             let divided = CircuitFieldElement::new(i / two, self.field);
    //             (self *  self).pow(divided)
    //             // pow(x * x, n / 2)
    //         },

    //         // _ => self * (self * self).pow(n - 1) / 2,
    //         _ => {
    //             let initial_pow = self * self;
    //             self * initial_pow.pow(n - 1) / 2
    //         },
    //     }     
    // }

    pub fn pow(&self, rhs: Self) -> U512 {
        // // self.value.as_
        // let x = self.value;
        // let n = rhs.value;

        // let p = U512::from_u32(3);
        // let two = NonZero::new(U512::from_u16(2)).unwrap();
        // if n == U512::ZERO {
        //     U512::ONE
        // } else if n == U512::ONE {
        //     x
        // } else if n % two == U512::ZERO {
        //     pow(x.mul_mod(&x, &p), n / two)
        // } else {
        //     x.mul_mod(&pow(x.mul_mod(&x, &p), (n - U512::ONE) / two), &p) 
        // } 
        pow(self.value, rhs.value, self.field.0)
    }

}


fn pow(x: U512, n: U512, p: U512) -> U512 {
    // let p = U512::from_u32(3);
    let two = NonZero::new(U512::from_u16(2)).unwrap();
    if n == U512::ZERO {
        U512::ONE
    } else if n == U512::ONE {
        x
    } else if n % two == U512::ZERO {
        pow(x.mul_mod(&x, &p), n / two, p)
    } else {
        x.mul_mod(&pow(x.mul_mod(&x, &p), (n - U512::ONE) / two, p), &p) 
    } 
}

impl Hash for CircuitFieldElement {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.field.hash(state);
        self.value.hash(state);
    }
}

impl Add for CircuitFieldElement {
    type Output = CircuitFieldElement;
    fn add(self, rhs: Self) -> Self::Output {
        let ans = self.value.add_mod(&rhs.value, &self.field.0);
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

        // if self.value < rhs.value {
        //     self.value = (self.field.0 - (rhs.value - self.value)) % self.field.0;
        // } else {
        //     self.value = (self.value - rhs.value) % self.field.0;
        // }

        self.value = self.value.sub_mod(&rhs.value, &self.field.0)
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
        // let ans = (self.value + rhs.value) % self.field.0;
        let ans = self.value.add_mod(&rhs.value, &self.field.0);
        CircuitFieldElement::new(ans, self.field.clone())
    }
}

impl AddAssign for CircuitFieldElement {
    fn add_assign(&mut self, rhs: Self) {
        // self.value = (self.clone().value + rhs.value) % self.clone().field.0;
        self.value = self.value.add_mod(&rhs.value, &self.field.0);
    }
}

impl Neg for CircuitFieldElement {
    type Output = CircuitFieldElement;

    fn neg(self) -> Self::Output {
        // let result = (self.clone().field.0 - self.clone().value) % self.clone().field.clone().0;
        let result  = self.value.neg_mod(&self.field.0);
        CircuitFieldElement::new(result, self.clone().field)
    }
}

impl<'a> Neg for &'a CircuitFieldElement {
    type Output = CircuitFieldElement;

    fn neg(self) -> Self::Output {
        // let result = (self.field.0 - self.value) % self.field.0;
        let result = self.value.neg_mod(&self.field.0);
        CircuitFieldElement::new(result, self.field.clone())
    }
}

impl Mul for CircuitFieldElement {
    type Output = CircuitFieldElement;

    fn mul(self, rhs: Self) -> Self::Output {
        // let result = (self.clone().value * rhs.value) % self.clone().field.0;
        let result = self.value.mul_mod(&rhs.value, &self.field.0);
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
        // let result = (self.value * rhs.value) % self.field.0;
        let result = self.value.mul_mod(&rhs.value, &self.field.0);
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
    let field = CircuitField(U512::from(101_u32));

    let field_element_lower = CircuitFieldElement {
        value: U512::from(100_u32),
        field: field.clone(),
    };

    let field_element_higher = CircuitFieldElement {
        value: U512::from(2_u32),
        field: field.clone(),
    };

    assert_eq!(
        field_element_lower + field_element_higher,
        CircuitFieldElement::new(U512::from(1_u32), field)
    );
}

#[test]
fn subtracts_including_negative_case() {
    let field = CircuitField(U512::from_u32(8));
    let a = field.element(U512::from_u32(4));
    let b = field.element(U512::from_u32(5));
    assert_eq!(a - b, field.element(U512::from_u32(7)));
}

#[test]
fn example_field() {
    let field = CircuitField(U512::from(41_u32));

    let field_element_lower = CircuitFieldElement {
        value: U512::from(1_u32),
        field: field.clone(),
    };

    let field_element_higher = CircuitFieldElement {
        value: U512::from(40_u32),
        field: field.clone(),
    };

    assert_eq!(
        field_element_lower + field_element_higher,
        CircuitFieldElement::new(U512::from(0_u32), field)
    );
}

#[test]
fn low_field() {
    let field = CircuitField(U512::from(13_u32));

    let field_element_lower = CircuitFieldElement {
        value: U512::from(3_u32),
        field: field.clone(),
    };

    let field_element_higher = CircuitFieldElement {
        value: U512::from(10_u32),
        field: field.clone(),
    };

    assert_eq!(
        field_element_lower + field_element_higher,
        CircuitFieldElement::new(U512::from(0_u32), field)
    );
}

#[test]
fn pow_works() {
    let field = CircuitField(U512::from_u128(997003001));
    let element = field.element(U512::from_u32(999));
    let element_right = field.element(U512::from_u32(3));
    assert_eq!(element.pow(element_right), U512::from_u128(997002999));
}

#[test]
fn pow_with_modulus_works() {
    let field = CircuitField(U512::from_u128(5));
    let element = field.element(U512::from_u32(999));
    let element_right = field.element(U512::from_u32(3));
    assert_eq!(element.pow(element_right), U512::from_u128(4));
}
