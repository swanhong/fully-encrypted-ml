// ipe.rs
extern crate rug;
use rug::Integer;
use rug::rand::RandState;

use std::time::SystemTime;

use crate::util::group::Group;
use crate::util::group::discrete_logarithm;
use crate::util::matrix::*;
use crate::util::vector::*;
use crate::util::decomp::Decomp;
use crate::ipe::keys::IpeSk;
use crate::ipe::scheme::{ipe_keygen, ipe_enc, ipe_enc_matrix_expression, ipe_dec};
use super::keys::QeSk;

pub fn qe_setup(grp: &Group, n_x: usize, n_y: usize, B: usize, rng: &mut RandState<'_>) -> QeSk {
    QeSk::new(grp, n_x, n_y, B, rng)
}

pub fn qe_keygen(
    sk: &QeSk,
    f: &Vec<Integer>,
    grp: &Group,
    decomp: &Decomp,
) -> (Vec<Integer>, Vec<Integer>) {
    let modulo = grp.delta.clone();

    let mut mf = sk.m.mul_vec(&f);
    vec_mod(&mut mf, &modulo);

    // for inverse decomposition, add power of base to sk_f
    let sk_f_nonpower = ipe_keygen(&sk.ipe_sk, &mf, grp);
    let sk_f = decomp.vector_pow_exp(&sk_f_nonpower);
    
    let dx_power = decomp.matrix_pow_col(&sk.d_x);
    let dy_power = decomp.matrix_pow_col(&sk.d_y);
    let dx_dy: Matrix = Matrix::tensor_product(&dx_power, &dy_power, &modulo);
    let dx_dy_t = dx_dy.transpose();

    // f_dx_dy = dx_dy_t * f
    let mut f_dx_dy = dx_dy_t.mul_vec(&f);
    vec_mod(&mut f_dx_dy, &modulo);

    let sk_red = vec_exp_with_base(&grp.g, &f_dx_dy, &grp.n_sq);

    (sk_f, sk_red)
}

pub fn qe_enc(
    sk: &QeSk,
    x: &Vec<Integer>,
    y: &Vec<Integer>,
    grp: &Group,
    decomp: &Decomp,
    rng: &mut RandState<'_>,
) -> (Vec<Integer>, Vec<Integer>, Vec<Integer>) {
    let mod_val = grp.delta.clone();

    fn qe_enc_for_x(
        v: &Matrix,
        d_null: &Matrix,
        d_inv: &Matrix,
        x: &Vec<Integer>,
        grp: &Group,
        decomp: &Decomp,
        rng: &mut RandState<'_>,
        mult_mu: bool,
    ) -> (Vec<Integer>, Vec<Integer>) {
        // right_x = x * mu + (sk->V * r_x)
        // ctxt_x = D_x_null_space_basis * rand + D_inv_x * right_x
        let modulo = grp.delta.clone();
        let mu = grp.mu.clone();

        let mut x_mu = x.clone();
        if mult_mu {
            x_mu = vec_mul_scalar(&x_mu, &mu);
            vec_mod(&mut x_mu, &modulo);
        }

        // right_x = x * mu + (sk->V * r_x)
        let r_x = gen_random_vector(2, &modulo, rng);
        let mut right_x = v.mul_vec(&r_x);
        right_x = vec_add(&mut right_x, &x_mu);
        vec_mod(&mut right_x, &modulo);

        // ctxt_x = D_x_null_space_basis * rand + D_inv_x * right_x
        let rand = gen_random_vector(d_null.cols, &modulo, rng);
        let mut d_perp_rand = d_null.mul_vec(&rand);
        vec_mod(&mut d_perp_rand, &modulo);

        let mut ctxt_x = d_inv.mul_vec(&right_x);
        ctxt_x = vec_add(&ctxt_x, &d_perp_rand);
        vec_mod(&mut ctxt_x, &modulo);

        (decomp.vector(&ctxt_x), r_x)
    }

    let (ctxt_x, r_x) 
        = qe_enc_for_x(
            &sk.v, 
            &sk.d_x_null, 
            &sk.d_x_inv, 
            &x, 
            &grp, 
            decomp, 
            rng, 
            true
        );
    let (ctxt_y, r_y) 
        = qe_enc_for_x(
            &sk.w, 
            &sk.d_y_null, 
            &sk.d_y_inv, 
            &y, 
            &grp, 
            decomp, 
            rng, 
            false
        );

    let h_up = tensor_product_vecs(&r_x, &y, &mod_val);
    let h_down = tensor_product_vecs(&x, &r_y, &mod_val);

    let mut h_join = vec![Integer::from(0); h_up.len() + h_down.len()];
    for i in 0..h_up.len() {
        h_join[i] = h_up[i].clone();
    }
    for i in 0..h_down.len() {
        h_join[h_up.len() + i] = h_down[i].clone();
    }

    let ctxt_ipe_nondecomp = ipe_enc(&sk.ipe_sk, &h_join, grp, false, rng);
    let ctxt_ipe = decomp.vector(&ctxt_ipe_nondecomp);
    (ctxt_x, ctxt_y, ctxt_ipe)
}


fn gen_enc_matrix_for_x(
    d_null: &Matrix,
    d_inv: &Matrix,
    v: &Matrix,
    grp: &Group,
    decomp: &Decomp,
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
    let enc_x_nondecomp = concatenate_vec_col(&mat_left, &mat_const_term);

    // println!("enc_x_non_decomp = {}", enc_x_nondecomp);

    (decomp.matrix_col(&enc_x_nondecomp), r_v)
}


fn compute_m_h_b_1(
    qe_sk: &QeSk,
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

    // println!("i_n_x = {}", i_n_x);
    // println!("r_y = {:?}", r_y);
    // println!("m_h_lr = {}", m_h_lr);

    // Concatenate M_h_ul and M_h_ur
    let m_h_u = concatenate_col(&m_h_ul, &m_h_ur);

    // Concatenate M_h_ll and M_h_lr
    let m_h_l = concatenate_col(&m_h_ll, &m_h_lr);

    // Concatenate M_h_u and M_h_l
    let m_h = concatenate_row(&m_h_u, &m_h_l);
    // println!("m_h_l = {}", m_h_l);
    // println!("m_h = {}", m_h);
    // Create rx_ry vector and compute V_rx_ry
    let rx_ry = tensor_product_vecs(r_x, r_y, &modulo);

    let i_2 = Matrix::get_identity(2);
    let v_i_2 = Matrix::tensor_product(&qe_sk.v, &i_2, &modulo);

    // Compute V_rx_ry
    let mut v_rx_ry = v_i_2.mul_vec(&rx_ry);
    vec_mod(&mut v_rx_ry, &modulo);

    // Extract values from V_rx_ry and set in h_b
    // h_b = (0, ... 0) || v_rx_ry (0: 2*n_x times)
    let mut h_b = vec![Integer::from(0); m_h.rows];
    for i in 0..(2 * n_x) {
        let val = v_rx_ry.get(i).unwrap();
        h_b[2 * n_y + i] = val.clone();
    }

    // Concatenate M_h and h_b to create M_h_b
    let m_h_b = concatenate_vec_col(&m_h, &h_b);

    // println!("m_h = {}", m_h);
    // println!("h_b = {:?}", h_b);
    // println!("m_h_b = {}", m_h_b);

    // Create M_h_b_1 as a diagonal matrix with M_h_b & 1
    let mut m_h_b_1 = Matrix::new(m_h_b.rows + 1, m_h_b.cols + 1);
    m_h_b_1.set(m_h_b.rows, m_h_b.cols, Integer::from(1));

    for i in 0..m_h_b.rows {
        for j in 0..m_h_b.cols {
            m_h_b_1.set(i, j, m_h_b.get(i, j));
        }
    }
    // println!("m_h_b_1 = {}", m_h_b_1);

    m_h_b_1
}

pub fn qe_enc_matrix_expression(
    qe_sk: &QeSk,
    n_x: usize,
    n_y: usize,
    grp: &Group,
    decomp: &Decomp,
    rng: &mut RandState<'_>,
) -> (Matrix, Matrix, Matrix) {
    let (enc_x, r_x) = gen_enc_matrix_for_x(
        &qe_sk.d_x_null,
        &qe_sk.d_x_inv,
        &qe_sk.v,
        grp,
        decomp,
        rng,
        false,
    );

    let (enc_y, r_y) = gen_enc_matrix_for_x(
        &qe_sk.d_y_null,
        &qe_sk.d_y_inv,
        &qe_sk.w,
        grp,
        decomp,
        rng,
        false,
    );

    // println!("enc_x size = {} x {}", enc_x.rows, enc_x.cols);
    // println!("enc_y size = {} x {}", enc_y.rows, enc_y.cols);

    let m_h_b_1 = compute_m_h_b_1(qe_sk, n_x, n_y, &r_x, &r_y, grp);

    // println!("m_h_b_1 size = {} x {}", m_h_b_1.rows, m_h_b_1.cols);

    let ipe_enc_mat = ipe_enc_matrix_expression(&qe_sk.ipe_sk, grp, false, rng);

    // println!("ipe_enc_mat size = {} x {}", ipe_enc_mat.rows, ipe_enc_mat.cols);

    let mut enc_h = ipe_enc_mat * m_h_b_1;
    enc_h.mod_inplace(&grp.delta);

    // println!("enc_h size = {} x {}", enc_h.rows, enc_h.cols);

    (enc_x, enc_y, enc_h)
}

pub fn qe_enc_matrix_same_xy(
    qe_sk: &QeSk,
    n_x: usize,
    grp: &Group,
    decomp: &Decomp,
    rng: &mut RandState<'_>,
) -> (Matrix, Matrix, Matrix) {
    let modulo = grp.delta.clone();


    let (enc_x, r_x) = gen_enc_matrix_for_x(
        &qe_sk.d_x_null,
        &qe_sk.d_x_inv,
        &qe_sk.v,
        grp,
        decomp,
        rng,
        false,
    );
    let (enc_y, r_y) = gen_enc_matrix_for_x(
        &qe_sk.d_y_null,
        &qe_sk.d_y_inv,
        &qe_sk.w,
        grp,
        decomp,
        rng,
        false,
    );



    let ipe_enc_mat = ipe_enc_matrix_expression(&qe_sk.ipe_sk, grp, false, rng);
    let m_h_b_1 = compute_m_h_b_1(qe_sk, n_x, n_x, &r_x, &r_y, grp);

    let mut qe_enc_h_origin = ipe_enc_mat * m_h_b_1;
    qe_enc_h_origin.mod_inplace(&modulo);

    let mut qe_enc_h_nodecomp = Matrix::new(qe_enc_h_origin.rows, n_x + 1);
    for i in 0..qe_enc_h_origin.rows {
        for j in 0..n_x {
            let mut val1 = qe_enc_h_origin.get(i, j);
            let val2 = qe_enc_h_origin.get(i, j + n_x);
            val1 = ((val1 * &grp.mu) + val2) % &modulo;
            qe_enc_h_nodecomp.set(i, j, val1);
        }
        let mut val1 = qe_enc_h_origin.get(i, 2 * n_x);
        let val2 = qe_enc_h_origin.get(i, 2 * n_x + 1);
        val1 = (val1 + val2) % &modulo;
        qe_enc_h_nodecomp.set(i, n_x, val1);
    }
    println!("enc_h_nondecomp size = {} x {}", qe_enc_h_nodecomp.rows, qe_enc_h_nodecomp.cols);
    let enc_h = decomp.matrix_col(&qe_enc_h_nodecomp);

    (enc_x, enc_y, enc_h)
}

pub fn qe_dec(
    (sk_f, sk_red): (&Vec<Integer>, &Vec<Integer>),
    (enc_x, enc_y, enc_h): (&Vec<Integer>, &Vec<Integer>, &Vec<Integer>),
    grp: &Group,
) -> Integer {
    let ctxt_tensor = tensor_product_vecs(&enc_x, &enc_y, &grp.delta);
    let val_mult = vec_inner_pow(&sk_red, &ctxt_tensor, &grp);

    let out_ipe = ipe_dec(&sk_f, &enc_h, &grp, false);
    let out_ipe_inv = out_ipe.clone().invert(&grp.n_sq).unwrap();
    // println!("out_ipe = {}", out_ipe);
    // println!("out_ipe_inv = {}", out_ipe_inv.clone());


    
    let val_mult = (val_mult * out_ipe_inv) % &grp.n_sq;
    // println!("val_mult = {}", val_mult.clone());
    discrete_logarithm(val_mult, &grp)
}