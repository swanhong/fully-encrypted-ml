#![allow(dead_code)]

extern crate rug;
use rug::Integer;
use rug::rand::RandState;

use crate::util::group::Group;
use crate::util::group::discrete_logarithm;
use crate::util::matrix::*;
use crate::util::vector::*;
use crate::ipfe::scheme::{ipfe_keygen, ipfe_enc, ipfe_dec};
use super::keys::QfeSk;

pub fn qfe_setup(grp: &Group, dim: usize, q: usize, rng: &mut RandState<'_>) -> QfeSk {
    QfeSk::new(grp, dim, q, rng)
}

pub fn qfe_keygen(
    sk: &QfeSk,
    f: &Vec<Integer>,
    grp: &Group,
) -> Vec<Integer> {
    // f.len() = dim^2
    let modulo = grp.delta.clone();

    // sk_f = sk->M_f * f
    // sk.m : 4dim * dim^2
    // f : dim^2
    // mf: 4 * dim
    let mut mf = sk.m.mul_vec(&f);
    vec_mod(&mut mf, &modulo);

    // for inverse decomposition, add power of base to sk_f
    // sk_f.len() = 4*dim + 2
    let sk_f = ipfe_keygen(&sk.ipfe_sk, &mf, grp);
    let dx_dy = Matrix::tensor_product(&sk.d_x, &sk.d_y, &modulo);
    let dx_dy_t = dx_dy.transpose();

    // f_dx_dy.len() = (dim + q)^2
    let mut f_dx_dy = dx_dy_t.mul_vec(&f);
    vec_mod(&mut f_dx_dy, &modulo);

    let sk_red = vec_exp_with_base(&grp.g, &f_dx_dy, &grp.n_sq);

    // sf_f.len() = 4 * dim + q + 2
    // sk_red.len() = (dim + q)^2
    let mut res = Vec::with_capacity(sk_f.len() + sk_red.len());
    res.extend_from_slice(&sk_f);
    res.extend_from_slice(&sk_red);
    // final size = 4 * dim + 2 + (dim + q)^2
    res
}

pub fn get_funcional_key_len(dim: usize, q: usize) -> usize {
    4 * dim + q + 2 + (dim + q) * (dim + q)
}

pub fn divide_vector_for_functional_key(
    sk: &Vec<Integer>,
    dim: usize,
    q: usize,
) -> (Vec<Integer>, Vec<Integer>) {
    assert_eq!(
        sk.len(), 
        4 * dim + q + 2 + (dim + q) * (dim + q),
        "divide_fk... sk.len() = {} != 4 * dim + q + 2 + (dim + q)^2 = {}",
        sk.len(),
        4 * dim + q + 2 + (dim + q) * (dim + q)
    );
    let sk_f = &sk[0..(4 * dim + q + 2)];
    let sk_red = &sk[(4 * dim + q + 2)..sk.len()];
    (sk_f.to_vec(), sk_red.to_vec())
}

pub fn qfe_enc_for_x(
    v: &Matrix,
    d_null: &Matrix,
    d_inv: &Matrix,
    x: &Vec<Integer>,
    r_x: &Vec<Integer>,
    grp: &Group,
    rng: &mut RandState<'_>,
    mult_mu: bool,
    use_sk_v: bool,
) -> Vec<Integer> {
    // right_x = x * mu + (sk->V * r_x)
    // ctxt_x = D_x_null * rand + D_inv_x * right_x
    let modulo = grp.delta.clone();

    // Generate random vectors
    let rand = gen_random_vector(d_null.cols, &modulo, rng);

    let mut x_mu = if mult_mu {
        vec_mul_scalar(&x, &grp.mu)
    } else {
        x.clone()
    };
    vec_mod(&mut x_mu, &modulo);

    // Compute right_x
    // right_x = x_mu + V * r_x, if use_sk_v = false, then remove V * r_x
    let mut right_x = if use_sk_v {
        vec_add(&v.mul_vec(&r_x), &x_mu)
    } else {
        x_mu
    };
    vec_mod(&mut right_x, &modulo);

    // Compute ctxt_x
    let mut ctxt_x = vec_add(
        &d_inv.mul_vec(&right_x), &d_null.mul_vec(&rand));
    vec_mod(&mut ctxt_x, &modulo);
    
    ctxt_x
}

pub fn qfe_enc_for_h(
    sk: &QfeSk, 
    x: &Vec<Integer>, 
    r_x: &Vec<Integer>, 
    r_y: &Vec<Integer>, 
    grp: &Group, 
    rng: &mut RandState,
    use_sk_v: bool,
) -> Vec<Integer> {
    let modulo = grp.delta.clone();
    // tmp = x_mu (+ V * r_x)
    let x_mu = vec_mul_scalar(&x, &grp.mu);
    let mut tmp = if use_sk_v {
        let tmp2 = sk.v.mul_vec(&r_x);
        vec_add(&x_mu, &tmp2)
    } else {
        x_mu
    };
    vec_mod(&mut tmp, &modulo);

    let h_up = tensor_product_vecs(&r_x, &x, &modulo);
    let h_down = tensor_product_vecs(&tmp, &r_y, &modulo);
    let mut h_join = Vec::with_capacity(h_up.len() + h_down.len());
    h_join.extend_from_slice(&h_up);
    h_join.extend_from_slice(&h_down);

    let ctxt_ipfe = ipfe_enc(&sk.ipfe_sk, &h_join, grp, false, rng);
    ctxt_ipfe
}

pub fn qfe_enc(
    sk: &QfeSk,
    x: &Vec<Integer>,
    grp: &Group,
    rng: &mut RandState<'_>,
) -> Vec<Integer> {
    assert_eq!(x.len(), sk.dim, "x.len() = {} != sk.dim = {}", x.len(), sk.dim);
    let r_x = gen_random_vector(2, &grp.delta, rng);
    let r_y = gen_random_vector(2, &grp.delta, rng);
    let ctxt_x = qfe_enc_for_x(
            &sk.v, 
            &sk.d_x_null, 
            &sk.d_x_inv, 
            &x, 
            &r_x,
            &grp, 
            rng, 
            true,
            true,
        );
    let ctxt_y = qfe_enc_for_x(
            &sk.w, 
            &sk.d_y_null, 
            &sk.d_y_inv, 
            &x, 
            &r_y,
            &grp, 
            rng, 
            false,
            true,
        );

    let ctxt_ipfe = qfe_enc_for_h(
        sk, 
        x, 
        &r_x, 
        &r_y, 
        grp, 
        rng,
        true,
    );
    
    assert_eq!(ctxt_x.len(), sk.dim + sk.q, "ctxt_x.len() = {} != sk.dim + sk.q = {}", ctxt_x.len(), sk.dim + sk.q);
    assert_eq!(ctxt_y.len(), sk.dim + sk.q, "ctxt_y.len() = {} != sk.dim + sk.q = {}", ctxt_y.len(), sk.dim + sk.q);
    assert_eq!(ctxt_ipfe.len(), 4 * sk.dim + sk.q + 2, "ctxt_ipfe.len() = {} != 4 * sk.dim + sk.q + 2 = {}", ctxt_ipfe.len(), 4 * sk.dim + sk.q + 2);
    
    let mut res = Vec::with_capacity(ctxt_x.len() + ctxt_y.len() + ctxt_ipfe.len());
    res.extend_from_slice(&ctxt_x);
    res.extend_from_slice(&ctxt_y);
    res.extend_from_slice(&ctxt_ipfe);
    res
}

pub fn qfe_cenc(
    sk: &QfeSk,
    m_f: &Matrix,
    grp: &Group,
    rng: &mut RandState<'_>,
) -> Matrix {
    let r_x = gen_random_vector(2, &grp.delta, rng);
    let r_y = gen_random_vector(2, &grp.delta, rng);
    // rmove random
    // let r_x = vec![Integer::from(0); 2];
    // let r_y = vec![Integer::from(0); 2];

    let mut ct_mat = Matrix::new(m_f.cols, get_ctxt_len(sk.dim, sk.q));
    for i in 0..m_f.cols {
        let f_col = m_f.get_col(i);
        let ct_0 = qfe_enc_for_x(
            &sk.v,
            &sk.d_x_null,
            &sk.d_x_inv,
            &f_col,
            &r_x,
            &grp,
            rng,
            true,
            i == m_f.cols - 1,
        );
        let ct_1 = qfe_enc_for_x(
            &sk.w,
            &sk.d_y_null,
            &sk.d_y_inv,
            &f_col,
            &r_y,
            &grp,
            rng,
            false,
            i == m_f.cols - 1,
        );
        let ct_h = qfe_enc_for_h(
            sk, &f_col, &r_x, &r_y, &grp, rng, i == m_f.cols - 1,
        );

        let mut row = Vec::with_capacity(ct_0.len() + ct_1.len() + ct_h.len());
        row.extend_from_slice(&ct_0);
        row.extend_from_slice(&ct_1);
        row.extend_from_slice(&ct_h);
        ct_mat.set_row(i, &row);
    }
    ct_mat.transpose()
}

pub fn get_ctxt_len(dim: usize, q: usize) -> usize {
    2 * (dim + q) + 4 * dim + q + 2
}

pub fn divide_vector_for_qfe_ctxt(
    ctxt: &Vec<Integer>,
    dim: usize,
    q: usize,
) -> (Vec<Integer>, Vec<Integer>, Vec<Integer>) {
    assert_eq!(
        ctxt.len(),
        2 * (dim + q) + 4 * dim + q + 2,
        "dividing ctxt.. ctxt.len() = {} != 2 * (dim + q) + 4 * dim + q + 2 = {}",
        ctxt.len(),
        2 * (dim + q) + 4 * dim + q + 2
    );
    let ctxt_x = &ctxt[0..(dim + q)];
    let ctxt_y = &ctxt[(dim + q)..(2 * (dim + q))];
    let ctxt_h = &ctxt[(2 * (dim + q))..ctxt.len()];
    (ctxt_x.to_vec(), ctxt_y.to_vec(), ctxt_h.to_vec())
}

pub fn qfe_dec(
    fk: &Vec<Integer>,
    ctxt: &Vec<Integer>,
    dim: usize,
    q: usize,
    grp: &Group,
) -> Integer {
    // println!("qfe_dec, input dim, q = {}, {}", dim, q);
    let (enc_x, enc_y, enc_h) = divide_vector_for_qfe_ctxt(ctxt, dim, q);
    let (sk_f, sk_red) = divide_vector_for_functional_key(fk, dim, q);
    let ctxt_tensor = tensor_product_vecs(&enc_x, &enc_y, &grp.delta);
    let val_mult = vec_inner_pow(&sk_red, &ctxt_tensor, &grp);



    let out_ipfe = ipfe_dec(&sk_f, &enc_h, &grp, false);
    let out_ipfe_inv = out_ipfe.clone().invert(&grp.n_sq).unwrap();

    let val = val_mult * out_ipfe_inv.clone();
    let val_mult = int_mod(&val, &grp.n_sq);
    discrete_logarithm(val_mult, &grp)
}