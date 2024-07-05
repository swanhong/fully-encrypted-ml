#![allow(dead_code)]

extern crate rug;
use rug::Integer;
use rug::rand::RandState;

use crate::util::group::Group;
use crate::util::group::discrete_logarithm;
use crate::util::matrix::*;
use crate::util::vector::*;
use crate::ipfe::scheme::{ipfe_keygen, ipfe_enc, ipfe_enc_matrix_expression, ipfe_dec};
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

pub fn qfe_enc(
    sk: &QfeSk,
    x: &Vec<Integer>,
    grp: &Group,
    rng: &mut RandState<'_>,
) -> Vec<Integer> {
    let mod_val = grp.delta.clone();

    fn qfe_enc_for_x(
        v: &Matrix,
        d_null: &Matrix,
        d_inv: &Matrix,
        x: &Vec<Integer>,
        grp: &Group,
        rng: &mut RandState<'_>,
        mult_mu: bool,
    ) -> (Vec<Integer>, Vec<Integer>) {
        // right_x = x * mu + (sk->V * r_x)
        // ctxt_x = D_x_null * rand + D_inv_x * right_x
        let modulo = grp.delta.clone();
        let mu = grp.mu.clone();

        let x_mu = if mult_mu {
            let mut x_temp = vec_mul_scalar(&x, &mu);
            vec_mod(&mut x_temp, &modulo);
            x_temp
        } else {
            x.clone()
        };
        
        // Generate random vectors
        let r_x = gen_random_vector(2, &modulo, rng);
        let r_x = vec_mul_scalar(&r_x, &(Integer::from(2) * &grp.n));
        let rand = gen_random_vector(d_null.cols, &modulo, rng);

        // Compute right_x
        let mut right_x = vec_add(&v.mul_vec(&r_x), &x_mu);
        vec_mod(&mut right_x, &modulo);

        // Compute ctxt_x
        let mut ctxt_x = vec_add(&d_inv.mul_vec(&right_x), &d_null.mul_vec(&rand));
        vec_mod(&mut ctxt_x, &modulo);
        
        (ctxt_x, r_x)
    }

    let (ctxt_x, r_x) 
        = qfe_enc_for_x(
            &sk.v, 
            &sk.d_x_null, 
            &sk.d_x_inv, 
            &x, 
            &grp, 
            rng, 
            true
        );
    let (ctxt_y, r_y) 
        = qfe_enc_for_x(
            &sk.w, 
            &sk.d_y_null, 
            &sk.d_y_inv, 
            &x, 
            &grp, 
            rng, 
            false
        );

    // tmp = x_mu + V * r_x
    let tmp1 = vec_mul_scalar(&x, &grp.mu);
    let tmp2 = sk.v.mul_vec(&r_x);
    let mut tmp = vec_add(&tmp1, &tmp2);
    vec_mod(&mut tmp, &mod_val);

    let h_up = tensor_product_vecs(&r_x, &x, &mod_val);
    let h_down = tensor_product_vecs(&tmp, &r_y, &mod_val);
    let mut h_join = Vec::with_capacity(h_up.len() + h_down.len());
    h_join.extend_from_slice(&h_up);
    h_join.extend_from_slice(&h_down);

    let ctxt_ipfe = ipfe_enc(&sk.ipfe_sk, &h_join, grp, false, rng);
    
    // size of ctxt_x = dim + q
    // size of ctxt_y = dim + q
    // size of ctxt_ipfe = 4 * dim + q + 2
    // println!("end of qfe_enc");
    println!("ctxt_x size = {} = {}", ctxt_x.len(), sk.dim + sk.q);
    println!("ctxt_y size = {} = {}", ctxt_y.len(), sk.dim + sk.q);
    println!("ctxt_ipfe size = {} = {}", ctxt_ipfe.len(), 4 * sk.dim + sk.q + 2);
    
    let mut res = Vec::with_capacity(ctxt_x.len() + ctxt_y.len() + ctxt_ipfe.len());
    res.extend_from_slice(&ctxt_x);
    res.extend_from_slice(&ctxt_y);
    res.extend_from_slice(&ctxt_ipfe);
    res
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

fn gen_enc_matrix_for_x(
    d_null: &Matrix,
    d_inv: &Matrix,
    v: &Matrix,
    grp: &Group,
    rng: &mut RandState<'_>,
    mult_mu: bool,
) -> (Matrix, Vec<Integer>) {
    let modulo = grp.delta.clone();
    let mu = grp.mu.clone();
    
    // ctxt_x = d_inv * x + (d_null * r_null + d_inv * v * r_v)
    let r_v = gen_random_vector(v.cols, &modulo, rng);
    let r_null = gen_random_vector(d_null.cols, &modulo, rng);
    let mut right = v.mul_vec(&r_v);
    vec_mod(&mut right, &modulo);
    right = d_inv.mul_vec(&right);
    vec_mod(&mut right, &modulo);

    let mut mat_const_term = d_null.mul_vec(&r_null);
    mat_const_term = vec_add(&mat_const_term, &right);
    vec_mod(&mut mat_const_term, &modulo);
    // println!("enc_b = {:?}", mat_const_term);
    // println!("sk.d_inv = {}", d_inv);

    let mut mat_left = d_inv.clone();
    if mult_mu {
        mat_left.mul_scalar_inplace(&mu);
        mat_left.mod_inplace(&modulo);
    }
    let enc_x = concatenate_vec_col(&mat_left, &mat_const_term);

    (enc_x, r_v)
}


fn compute_m_h_b_1(
    qfe_sk: &QfeSk,
    n_x: usize,
    n_y: usize,
    r_x: &Vec<Integer>,
    r_y: &Vec<Integer>,
    grp: &Group,
) -> Matrix {
    let modulo = grp.delta.clone();

    // M_h = (r_x \otimes I_{n_y}      0)
    //       (0      I_{n_x} \otimes r_y)

    // Initialize M_h_ul, M_h_ur, M_h_ll, M_h_lr matrices
    let m_h_ur = Matrix::new(2 * n_y, n_x);
    let m_h_ll = Matrix::new(2 * n_x, n_y);

    let i_n_x = Matrix::get_identity(n_x);
    let i_n_y = Matrix::get_identity(n_y);
    
    // Compute M_h_ul
    let m_h_ul = Matrix::tensor_product_vec_left(&r_x, &i_n_y, &modulo);

    // Compute M_h_lr
    let m_h_lr = Matrix::tensor_product_vec_right(&i_n_x, &r_y, &modulo);

    // Concatenate M_h_ul and M_h_ur
    let m_h_u = concatenate_col(&m_h_ul, &m_h_ur);

    // Concatenate M_h_ll and M_h_lr
    let m_h_l = concatenate_col(&m_h_ll, &m_h_lr);

    // Concatenate M_h_u and M_h_l
    let m_h = concatenate_row(&m_h_u, &m_h_l);

    // Create rx_ry vector and compute V_rx_ry
    let rx_ry = tensor_product_vecs(r_x, r_y, &modulo);

    let i_2 = Matrix::get_identity(2);
    let v_i_2 = Matrix::tensor_product(&qfe_sk.v, &i_2, &modulo);

    // Compute V_rx_ry
    let mut v_rx_ry = v_i_2.mul_vec(&rx_ry);
    vec_mod(&mut v_rx_ry, &modulo);

    // Extract values from V_rx_ry and set in h_b
    // h_b = (0, ... 0) || v_rx_ry (0: 2*n_x times)
    let mut h_b = vec![Integer::from(0); m_h.rows];
    for i in 0..v_rx_ry.len() {
        let val = v_rx_ry[i].clone();
        h_b[2 * n_x + i] = val.clone();
    }

    // Concatenate M_h and h_b to create M_h_b
    let m_h_b = concatenate_vec_col(&m_h, &h_b);

    // Create M_h_b_1 as a diagonal matrix with M_h_b & 1
    let mut m_h_b_1 = Matrix::new(m_h_b.rows + 1, m_h_b.cols + 1);
    m_h_b_1.set(m_h_b.rows, m_h_b.cols, Integer::from(1));

    for i in 0..m_h_b.rows {
        for j in 0..m_h_b.cols {
            m_h_b_1.set(i, j, m_h_b.get(i, j));
        }
    }
    m_h_b_1
}

pub fn qfe_enc_matrix_expression(
    qfe_sk: &QfeSk,
    n_x: usize,
    n_y: usize,
    grp: &Group,
    rng: &mut RandState<'_>,
) -> (Matrix, Matrix, Matrix) {
    let (enc_x, r_x) = gen_enc_matrix_for_x(
        &qfe_sk.d_x_null,
        &qfe_sk.d_x_inv,
        &qfe_sk.v,
        grp,
        rng,
        true,
    );

    let (enc_y, r_y) = gen_enc_matrix_for_x(
        &qfe_sk.d_y_null,
        &qfe_sk.d_y_inv,
        &qfe_sk.w,
        grp,
        rng,
        false,
    );

    let m_h_b_1 = compute_m_h_b_1(qfe_sk, n_x, n_y, &r_x, &r_y, grp);

    let ipfe_enc_mat = ipfe_enc_matrix_expression(&qfe_sk.ipfe_sk, grp, false, rng);

    let mut enc_h = ipfe_enc_mat * m_h_b_1;
    enc_h.mod_inplace(&grp.delta);

    // enc_h = decomp.matrix_col(&enc_h);
    (enc_x, enc_y, enc_h)
}

pub fn qfe_enc_matrix_same_xy(
    qfe_sk: &QfeSk,
    n_x: usize,
    grp: &Group,
    rng: &mut RandState<'_>,
) -> (Matrix, Matrix, Matrix) {
    let modulo = grp.delta.clone();


    let (enc_x, r_x) = gen_enc_matrix_for_x(
        &qfe_sk.d_x_null,
        &qfe_sk.d_x_inv,
        &qfe_sk.v,
        grp,
        rng,
        false,
    );
    let (enc_y, r_y) = gen_enc_matrix_for_x(
        &qfe_sk.d_y_null,
        &qfe_sk.d_y_inv,
        &qfe_sk.w,
        grp,
        rng,
        true,
    );

    let ipfe_enc_mat = ipfe_enc_matrix_expression(&qfe_sk.ipfe_sk, grp, false, rng);
    let m_h_b_1 = compute_m_h_b_1(qfe_sk, n_x, n_x, &r_x, &r_y, grp);
    let mut qfe_enc_h_origin = ipfe_enc_mat * m_h_b_1;
    qfe_enc_h_origin.mod_inplace(&modulo);

    let mut qfe_enc_h = Matrix::new(qfe_enc_h_origin.rows, n_x + 1);
    for i in 0..qfe_enc_h_origin.rows {
        for j in 0..n_x {
            let mut val1 = qfe_enc_h_origin.get(i, j);
            let val2 = qfe_enc_h_origin.get(i, j + n_x);
            let tmp = (val1 * &grp.mu) + val2;
            val1 = int_mod(&tmp, &modulo);
            // val1 = ((val1 * &grp.mu) + val2) % &modulo;
            qfe_enc_h.set(i, j, val1);
        }
        let mut val1 = qfe_enc_h_origin.get(i, 2 * n_x);
        let val2 = qfe_enc_h_origin.get(i, 2 * n_x + 1);
        let tmp = val1 + val2;
        val1 = int_mod(&tmp, &modulo);
        // val1 = (val1 + val2) % &modulo;
        qfe_enc_h.set(i, n_x, val1);
    }

    (enc_x, enc_y, qfe_enc_h)
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
    // println!("qfe_dec");
    // println!("enc_x size = {}", enc_x.len());
    // println!("enc_y size = {}", enc_y.len());
    // println!("enc_h size = {}", enc_h.len());
    // println!("sk_f size = {}", sk_f.len());
    // println!("sk_red size = {}", sk_red.len());
    // println!("enc_x.len * enc_y.len = sk_red len");
    // println!("enc_h len = sk_f len");
    let ctxt_tensor = tensor_product_vecs(&enc_x, &enc_y, &grp.delta);
    let val_mult = vec_inner_pow(&sk_red, &ctxt_tensor, &grp);

    let out_ipfe = ipfe_dec(&sk_f, &enc_h, &grp, false);
    let out_ipfe_inv = out_ipfe.clone().invert(&grp.n_sq).unwrap();

    let val = val_mult * out_ipfe_inv.clone();
    let val_mult = int_mod(&val, &grp.n_sq);
    discrete_logarithm(val_mult, &grp)
}