// ipe.rs
extern crate rug;
use rug::Integer;
use rug::rand::RandState;

use crate::util::group::Group;
use crate::util::group::discrete_logarithm;
use crate::util::matrix::*;
use crate::util::vector::*;
use super::keys::IpeSk;
use std::time::{Duration, SystemTime};

pub fn ipe_setup(group: &Group, dim: usize, b: usize, rng: &mut RandState<'_>) -> IpeSk {
    IpeSk::new(dim, b, group, rng)
}


pub fn ipe_keygen(sk: &IpeSk, y: &Vec<Integer>, grp: &Group) -> Vec<Integer> {
    let mut val;
    let mod_val = grp.delta.clone();
    
    let mut u_t_y_left = Vec::with_capacity(sk.u_t.rows);
    let mut u_t_y = Vec::with_capacity(u_t_y_left.len() + y.len());
    // U_t_y = (-1 * sk->U_t * y) || y
    for i in 0..sk.u_t.rows {
        val = Integer::from(0);
        for j in 0..sk.u_t.cols {
            val -= sk.u_t.get(i, j).clone() * y[j].clone();
            // val = val.clone().div_rem_euc(mod_val.clone()).1;
            val = int_mod(&val, &mod_val);
        }
        u_t_y_left.push(val);
    }

    u_t_y.extend_from_slice(&u_t_y_left);
    u_t_y.extend_from_slice(y);

    // sk_f_mat = (sk->D_inv_left + sk->D_inv_right * sk->U) * u_t_y
    let mut left = Matrix::new(1, u_t_y.len());
    for i in 0..u_t_y.len() {
        let val = &u_t_y[i];
        left.set(0, i, val.clone());
    }

    let sk_f_mat = (left * &sk.d) % &mod_val;

    let sk_f_col = sk_f_mat.get_row(0);

    let start = SystemTime::now();
    let sk_f: Vec<Integer> = sk_f_col
    .iter()
    .map(|val| {
        grp.g
            .clone()
            .pow_mod(val, &grp.n_sq)
            .unwrap_or_else(|e| panic!("Error in ipe_keygen: {}", e))
    })
    .collect();
    let end = start.elapsed();
    let num_pow_mod = sk_f.len();
    // println!("Time elapsed in ipe_keygen::pow_mod is: {:?} for {} pow_mod", end, num_pow_mod);
    
    sk_f
}

pub fn ipe_enc(
    sk: &IpeSk,
    x: &Vec<Integer>,
    grp: &Group,
    mult_mu: bool,
    rng: &mut RandState<'_>,
) -> Vec<Integer> {
    let mod_val = grp.delta.clone();
    let r = mod_val.clone().random_below(rng);

    let rand = gen_random_vector(sk.d_perp.cols, &mod_val, rng);
    let mut d_perp_rand = sk.d_perp.mul_vec(&rand);
    vec_mod(&mut d_perp_rand, &mod_val);

    let mut x_mu = if mult_mu {
        vec_mul_scalar(x, &grp.mu)
    } else {
        x.to_vec()
    };
    vec_mod(&mut x_mu, &mod_val);

    let right = vec_mul_scalar(&sk.a, &r);
    let mut right2 = sk.u.mul_vec(&sk.a);
    right2 = vec_mul_scalar(&right2, &r);
    right2 = vec_add(&right2, &x_mu);

    let mut right_joined = Vec::new();
    right_joined.extend(right);
    right_joined.extend(right2);

    let mut ctxt = sk.d_inv.mul_vec(&right_joined);
    ctxt = vec_add(&ctxt, &d_perp_rand);
    vec_mod(&mut ctxt, &mod_val);

    ctxt
}

pub fn ipe_enc_matrix_expression(sk: &IpeSk, grp: &Group, mult_mu: bool, rand: &mut RandState<'_>) -> Matrix {
    let mod_val = grp.delta.clone();
    let r = mod_val.clone().random_below(rand);

    // enc(x) = sk_enc * mu * x + (D_perp * rvec + sk1 * r)
    // matrix = (sk_enc * mu || (D_perp * rvec + sk1 * r))

    // Compute D_perp_rand = D_null_space_basis * rvec
    let rvec = gen_random_vector(sk.d_perp.cols, &mod_val, rand);
    let d_perp_rand = sk.d_perp.clone() * rvec;

    // tmp = sk1 * r
    let tmp = vec_mul_scalar(&sk.sk1, &r);

    // x_b = D_perp_rand + sk1 * r
    let mut enc_x_b = vec_add(&d_perp_rand, &tmp);
    vec_mod(&mut enc_x_b, &mod_val);

    let sk_enc_mu = if mult_mu {
        let mut sk_enc_mu = sk.sk_enc.clone();
        sk_enc_mu.mul_scalar_inplace(&grp.mu);
        sk_enc_mu.mod_inplace(&mod_val);
        sk_enc_mu
    } else {
        sk.sk_enc.clone()
    };

    let mut ipe_enc_mat = concatenate_vec_col(&sk_enc_mu, &enc_x_b);
    ipe_enc_mat.mod_inplace(&mod_val);

    ipe_enc_mat 
}

pub fn ipe_dec(sk_f: &Vec<Integer>, ctxt: &Vec<Integer>, grp: &Group, solve_dl: bool) -> Integer {
    let start = SystemTime::now();
    let mut out = vec_inner_pow(&sk_f, &ctxt, &grp);
    let end = start.elapsed();
    // println!("Time elapsed in ipe_dec::vec_inner_pow is: {:?} for {} pow_mod", end, sk_f.len());
    if solve_dl {
        out = discrete_logarithm(out.clone(), &grp);
    }
    out
}


