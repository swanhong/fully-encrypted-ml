// ipe.rs
extern crate rug;
use rug::Integer;
use rug::rand::RandState;

use std::time::SystemTime;

use crate::util::group::Group;
use crate::util::group::discrete_logarithm;
use crate::util::matrix::*;
use crate::util::vector::*;
use super::keys::IPE_sk;

pub fn ipe_setup(group: &Group, dim: usize, b: usize) -> IPE_sk {
    IPE_sk::new(dim, b, group)
}


pub fn ipe_keygen(sk: &IPE_sk, y: &Vec<Integer>, grp: &Group) -> Vec<Integer> {
    let mut val = Integer::new();
    let mut val2 = Integer::new();
    let mod_val = grp.delta.clone();
    
    let mut u_t_y_left = Vec::with_capacity(sk.u_t.rows);
    let mut u_t_y = Vec::with_capacity(u_t_y_left.len() + y.len());
    // let mut sk_f_mat = Matrix::new(1, sk.d.cols);
    // U_t_y = (-1 * sk->U_t * y) || y
    for i in 0..sk.u_t.rows {
        val = Integer::from(0);
        for j in 0..sk.u_t.cols {
            val -= sk.u_t.get(i, j).clone() * y[j].clone();
            val = val.clone().div_rem_euc(mod_val.clone()).1;
        }
        u_t_y_left.push(val);
    }

    u_t_y.extend_from_slice(&u_t_y_left);
    u_t_y.extend_from_slice(y);

    println!("sk->u_t = {}", sk.u_t);
    println!("y = {:?}", y);
    println!("u_t_y_left = {:?}", u_t_y_left);
    println!("u_t_y = {:?}", u_t_y);

    // sk_f_mat = (sk->D_inv_left + sk->D_inv_right * sk->U) * u_t_y
    let mut left = Matrix::new(1, u_t_y.len());
    for i in 0..u_t_y.len() {
        let val = &u_t_y[i];
        left.set(0, i, val.clone());
    }

    let mut sk_f_mat = (left * &sk.d) % &mod_val;

    println!("sk_f_mat = {}", sk_f_mat);

    
    let sk_f_col = sk_f_mat.get_row(0);

    // let mut sk_f: Vec<Integer> = Vec::new();
    // for i in 0..sk_f_col.len() {
    //     let power = match grp.g.clone().pow_mod(&sk_f_col[i], &grp.delta) {
    //         Ok(x) => x,
    //         Err(e) => panic!("Error in ipe_keygen: {}", e),
    //     };
    //     sk_f.push(power);
    // }
    let sk_f: Vec<Integer> = sk_f_col
    .iter()
    .map(|val| {
        grp.g
            .clone()
            .pow_mod(val, &grp.n_sq)
            .unwrap_or_else(|e| panic!("Error in ipe_keygen: {}", e))
    })
    .collect();
    
    sk_f
}

pub fn ipe_enc_matrix_expression(sk: &IPE_sk, grp: &Group, mult_mu: bool, rand: &mut RandState<'_>) -> Matrix {
    let mod_val = grp.delta.clone();
    let r = mod_val.clone().random_below(rand);

    // enc(x) = sk_enc * mu * x + (D_perp * rvec + sk1 * r)
    // matrix = (sk_enc * mu || (D_perp * rvec + sk1 * r))

    // Compute D_perp_rand = D_null_space_basis * rvec
    let rvec = gen_random_vector(sk.d_perp.cols, &mod_val, rand);
    let d_perp_rand = sk.d_perp.clone() * rvec;

    // tmp = sk1 * r
    let tmp = mul_vec_scalar(&sk.sk1, &r);

    // x_b = D_perp_rand + sk1 * r
    let mut enc_x_b = add_vec(&d_perp_rand, &tmp);
    mod_vec(&mut enc_x_b, &mod_val);

    let sk_enc_mu = if mult_mu {
        let mut sk_enc_mu = sk.sk_enc.clone();
        sk_enc_mu.mul_scalar_inplace(&grp.mu);
        sk_enc_mu.mod_inplace(&mod_val);
        sk_enc_mu
    } else {
        sk.sk_enc.clone()
    };

    println!("sk_enc_mu size = {} x {}", sk_enc_mu.rows, sk_enc_mu.cols);
    println!("enc_x_b len = {}", enc_x_b.len());
    let mut ipe_enc_mat = concatenate_vec_col(&sk_enc_mu, &enc_x_b);
    ipe_enc_mat.mod_inplace(&mod_val);

    ipe_enc_mat 
}

pub fn ipe_dec(sk_f: &Vec<Integer>, ctxt: &Vec<Integer>, grp: &Group, solve_dl: bool) -> Integer {
    let mut out = vec_inner_pow(&sk_f, &ctxt, &grp);
    println!("out before dl = {}", out);
    if solve_dl {
        out = discrete_logarithm(out.clone(), &grp);
    }
    out
}


