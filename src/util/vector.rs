#![allow(dead_code)]

use rug::{Integer, Complete};
use rug::rand::RandState;
use crate::util::group::Group;
use rayon::prelude::*;

pub fn gen_random_vector(dim: usize, bound: &Integer, rand: &mut RandState<'_>) -> Vec<Integer> {
    let mut vec = vec![Integer::from(0); dim];
    for i in 0..dim {
        vec[i] = bound.clone().random_below(rand);
    }
    vec
}

pub fn gen_random_vector_signed(dim: usize, bound: &Integer, rand: &mut RandState<'_>) -> Vec<Integer> {
    let mut vec = vec![Integer::from(0); dim];
    let bound2 = bound.clone() * Integer::from(2); 
    for i in 0..dim {
        vec[i] = bound2.clone().random_below(rand) - bound;
    }
    vec
}

pub fn vec_mul_scalar(vec: &Vec<Integer>, scalar: &Integer) -> Vec<Integer> {
    let mut res = vec![Integer::from(0); vec.len()];
    for i in 0..vec.len() {
        res[i] = vec[i].clone() * scalar.clone();
    }
    res
}

pub fn vec_add(vec1: &Vec<Integer>, vec2: &Vec<Integer>) -> Vec<Integer> {
    assert!(vec1.len() == vec2.len());
    let mut res = vec![Integer::from(0); vec1.len()];
    for i in 0..vec1.len() {
        res[i] = vec1[i].clone() + vec2[i].clone();
    }
    res
}

pub fn vec_mod(vec: &mut Vec<Integer>, modulus: &Integer) {
    for i in 0..vec.len() {
        // vec[i] = vec[i].clone().div_rem_euc(modulus.clone()).1;
        vec[i] = int_mod(&vec[i], modulus);
    }
}

pub fn vec_pow(vec: &[Integer], exp: &Integer, modulus: &Integer) -> Vec<Integer> {
    vec.par_iter()
        .map(|val| val.clone().pow_mod(exp, modulus).unwrap())
        .collect()
}

pub fn vec_pow_by_vec(vec: &[Integer], exp: &[Integer], modulus: &Integer) -> Vec<Integer> {
    assert_eq!(vec.len(), exp.len());   
    vec.par_iter()
        .zip(exp.par_iter())
        .map(|(val, exp)| val.clone().pow_mod(exp, modulus).unwrap())
        .collect()
}

pub fn vec_inner_pow(v_base: &Vec<Integer>, v_exp: &Vec<Integer>, grp: &Group) -> Integer {
    assert_eq!(v_base.len(), v_exp.len());
    
    let modulo = &grp.n_sq;
    let val = v_base.par_iter().zip(v_exp.par_iter()).map(|(base, exp)| {
        base.clone().pow_mod(exp, modulo).unwrap()
    }).collect::<Vec<Integer>>();

    let mut out = Integer::from(1);
    for i in 0..val.len() {
        out *= val[i].clone();
        out = out % modulo;
    }
    out
}

pub fn vec_exp_with_base(base: &Integer, v_exp: &Vec<Integer>, modulo: &Integer) -> Vec<Integer> {
    v_exp.par_iter()
        .map(|exp| base.clone().pow_mod(exp, modulo).unwrap())
        .collect()
}


pub fn tensor_product_vecs(vec1: &Vec<Integer>, vec2: &Vec<Integer>, modulo: &Integer) -> Vec<Integer> {
    let mut res = vec![Integer::from(0); vec1.len() * vec2.len()];
    for i in 0..vec1.len() {
        for j in 0..vec2.len() {
            let mul = vec1[i].clone() * vec2[j].clone();
            res[i * vec2.len() + j] = int_mod(&mul, &modulo);
            // res[i * vec2.len() + j] = vec1[i].clone() * vec2[j].clone() % modulo;
        }
    }
    res
}

pub fn eval_quadratic(
    x: &Vec<Integer>, 
    y: &Vec<Integer>, 
    f: &Vec<Integer>, 
) -> Integer {
    assert_eq!(f.len(), x.len() * y.len());
    let mut out = Integer::from(0);
    for i in 0..x.len() {
        for j in 0..y.len() {
            out += f[i * y.len() + j].clone() * x[i].clone() * y[j].clone();
        }
    }
    out
}

pub fn int_mod(input: &Integer, modulo: &Integer) -> Integer {
    let mut val = input.div_rem_ref(modulo).complete().1;
    
    let val2 = val.clone() * Integer::from(2);
    if val2.ge(modulo){
        val = val.clone() - modulo;
    } else if val2.le(&(modulo * Integer::from(-1))) {
        val = val.clone() + modulo;
    }
    val
}