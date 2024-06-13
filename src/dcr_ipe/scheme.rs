// dcr_ipe.rs
extern crate rug;
use rug::Integer;
use rug::rand::RandState;
use std::ops::MulAssign;

use crate::util::group::Group;
use crate::util::vector::{vec_exp_with_base, vec_pow, int_mod, vec_mod};

// takes input of dimension, secret key bound, and group
// returns a key pair (sk, pk)
pub fn dcr_setup(
    dim: usize, 
    sk_bound: &Integer, 
    grp: &Group,
    rng: &mut RandState<'_>,
) -> (Vec<Integer>, Vec<Integer>) {
    let mut sk = vec![Integer::from(0); dim];
    let mut bound2 = sk_bound.clone();
    bound2.mul_assign(2);
    for i in 0..dim {
        sk[i] = bound2.clone().random_below(rng);
        sk[i] -= sk_bound;
    }
    let pk = vec_exp_with_base(&grp.g, &sk, &grp.n_sq);

    (sk, pk)
}

pub fn dcr_keygen(sk: &Vec<Integer>, y: &Vec<Integer>) -> Integer {
    // sk_y = inner product between sk and y
    assert_eq!(sk.len(), y.len(), "sk and y must have the same length");
    let mut sk_y = Integer::from(0);
    for i in 0..sk.len() {
        sk_y += sk[i].clone() * y[i].clone();
    }
    sk_y
}

pub fn dcr_enc(pk: &Vec<Integer>, x: &Vec<Integer>, grp: &Group, rand: &mut RandState<'_>) -> Vec<Integer> {
    assert!(pk.len() == x.len());

    let r: Integer = grp.n_sq.clone().random_below(rand);
    let ct_x_last = grp.g.clone().pow_mod(&r, &grp.n_sq).unwrap();
    let mut x_encode = vec![Integer::from(0); x.len()];
    for i in 0..x.len() {
        x_encode[i] = x[i].clone() * &grp.n + Integer::from(1);
    }
    let mut ct_x = vec_pow(&pk, &r, &grp.n_sq);
    ct_x.iter_mut()
    .zip(x_encode.iter_mut())
    .for_each(|(val, x_encode)| val.mul_assign(x_encode.clone()));
    ct_x.push(ct_x_last);
    vec_mod(&mut ct_x, &grp.n_sq);
    ct_x
}

pub fn dcr_dec(ct_x: &Vec<Integer>, y: &Vec<Integer>, sk_y: &Integer, grp: &Group) -> Integer {
    assert!(ct_x.len() == y.len() + 1);

    let dim = y.len();
    let mut out = Integer::from(1);
    
    for i in 0..ct_x.len()-1 {
        out *= ct_x[i].clone().pow_mod(&y[i], &grp.n_sq).unwrap();
        out = int_mod(&out, &grp.n_sq);
    }
    let num = 
            ct_x[dim].clone().pow_mod(&sk_y, &grp.n_sq).unwrap()
            .invert(&grp.n_sq).unwrap();
    out *= num;
    // out = out.div_rem_euc(grp.n_sq.clone()).1;
    out = int_mod(&out, &grp.n_sq);
    // out = (out - 1) / &grp.n;
    out = out - 1;
    let mut quo = out.clone();
    let mut rem = grp.n.clone();
    quo.div_rem_mut(&mut rem);
    // let (quo, rem) = out.div_rem(grp.n.clone());
    assert!(rem == Integer::from(0), "Remainder is not 0");
    quo
}