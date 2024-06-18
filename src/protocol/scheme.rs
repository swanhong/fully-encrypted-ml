extern crate rug;
use rug::Integer;
use rug::rand::RandState;

use crate::util::group::Group;
use crate::util::matrix::{Matrix, remove_diag_one};
use crate::util::vector::{vec_add, vec_mod};
use crate::util::decomp::Decomp;
use crate::dcr_ipe::scheme::{dcr_setup, dcr_enc, dcr_keygen, dcr_dec};
use crate::qe::keys::QeSk;
use crate::qe::scheme::{qe_setup, qe_keygen, qe_enc_matrix_same_xy, qe_dec};
use std::time::SystemTime;

pub fn protocol_setup(
    dim_vec: Vec<usize>,
    f_num: usize,
    k: usize,
    sk_bound: &Integer,
    grp: &Group,
    rng: &mut RandState<'_>,
) -> ((Vec<Integer>, Vec<Integer>), QeSk, Vec<QeSk>, QeSk) {
    let dim = dim_vec[0];
    let (dcr_sk, dcr_pk) = dcr_setup(dim + dim + 1, sk_bound, grp, rng);
    let qe_sk_init = qe_setup(grp, dim + k, dim + k, 2 * (dim + k) + 1, rng);
    let mut qe_sk_fcn = Vec::with_capacity(f_num);
    for i in 1..dim_vec.len()-1 {
        let dim = dim_vec[i];
        qe_sk_fcn.push(qe_setup(grp, dim + k, dim + k, 2 * (dim + k) + 1, rng));
    }
    let dim = dim_vec[dim_vec.len()-1];
    let qe_sk_end = qe_setup(grp, dim + k, dim + k, 2 * (dim + k) + 1, rng);

    ((dcr_sk, dcr_pk), qe_sk_init, qe_sk_fcn, qe_sk_end)
}

pub fn protocol_enc_init(
    dcr_pk: &Vec<Integer>,
    gamma_right: &Matrix,
    x: &Vec<Integer>,
    grp: &Group,
    rng: &mut RandState<'_>,
) -> Vec<Integer> {
    let mut x1 = x.clone();
    x1.push(Integer::from(1));

    let mut gamma_right_x = gamma_right.mul_vec(&x1);
    vec_mod(&mut gamma_right_x, &grp.delta);
    dcr_enc(dcr_pk, &gamma_right_x, &grp, rng)
}

pub fn protocol_keygen_switch(
    qe_sk: &QeSk,
    dcr_sk: &Vec<Integer>,
    h_right: &Matrix,
    gamma_left: &Matrix,
    dim: usize,
    k: usize,
    decomp: &Decomp,
    grp: &Group,
    rng: &mut RandState<'_>,
) -> (
    (Matrix, Matrix, Matrix),
    (Vec<Integer>, Vec<Integer>, Vec<Integer>),
) {
    let modulo = grp.delta.clone();
    let (qe_enc_mat_x, qe_enc_mat_y, qe_enc_mat_h)
        = qe_enc_matrix_same_xy(&qe_sk, dim + k, &grp, &decomp, rng, true);
    
    fn gen_switch_key_for_enc_and_dcr(
        dcr_sk: &Vec<Integer>,
        qe_enc_mat: &Matrix,
        h_right: &Matrix,
        gamma_left: &Matrix,
        decomp: &Decomp,
        modulo: &Integer,
    ) -> (Matrix, Vec<Integer>) {
        let mut xh = qe_enc_mat * h_right;
        xh.mod_inplace(modulo);
        let xh_decomp = decomp.matrix_col(&xh);
        let mut switch_key_x = xh_decomp * gamma_left;

        switch_key_x.mod_inplace(modulo);

        let mut switch_key_dcr_x = vec![Integer::from(0); switch_key_x.rows];
        for i in 0..switch_key_x.rows {
            let row = switch_key_x.get_row(i);
            switch_key_dcr_x[i] = dcr_keygen(&dcr_sk, &row);
        }

        (switch_key_x, switch_key_dcr_x)                
    }
    let (switch_key_x, switch_key_dcr_x) = gen_switch_key_for_enc_and_dcr(
        &dcr_sk,
        &qe_enc_mat_x,
        &h_right,
        &gamma_left,
        &decomp,
        &modulo,
    );

    let (switch_key_y, switch_key_dcr_y) = gen_switch_key_for_enc_and_dcr(
        &dcr_sk,
        &qe_enc_mat_y,
        &h_right,
        &gamma_left,
        &decomp,
        &modulo,
    );

    let (switch_key_h, switch_key_dcr_h) = gen_switch_key_for_enc_and_dcr(
        &dcr_sk,
        &qe_enc_mat_h,
        &h_right,
        &gamma_left,
        &decomp,
        &modulo,
    );

    (
        (switch_key_x, switch_key_y, switch_key_h),
        (switch_key_dcr_x, switch_key_dcr_y, switch_key_dcr_h),
    )
}

pub fn protocol_keyswitch(
    ct_in: &Vec<Integer>,
    (switch_key_x, switch_key_y, switch_key_h): (&Matrix, &Matrix, &Matrix),
    (switch_key_dcr_x, switch_key_dcr_y, switch_key_dcr_h): (&Vec<Integer>, &Vec<Integer>, &Vec<Integer>),
    decomp: &Decomp,
    grp: &Group,
) -> (Vec<Integer>, Vec<Integer>, Vec<Integer>) {
    
    pub fn dcr_dec_multi(
        ct_in: &Vec<Integer>,
        switch_key: &Matrix,
        switch_key_dcr: &Vec<Integer>,
        decomp: &Decomp,
        grp: &Group,
    ) -> Vec<Integer> {
        assert_eq!(switch_key.rows, switch_key_dcr.len(), "error in dcr_dec_multi inputs");
        assert_eq!(ct_in.len(), switch_key.cols + 1, "error in dcr_dec_multi inputs");
        let mut ct_out = vec![Integer::from(0); switch_key.rows];
        for i in 0..switch_key.rows {
            let row = switch_key.get_row(i);
            ct_out[i] = dcr_dec(&ct_in, &row, &switch_key_dcr[i], grp);
        }
        vec_mod(&mut ct_out, &grp.n);
        decomp.vector_inv(&ct_out)
    }

    println!("do dcr_dec_multi for {} x {} times", switch_key_x.rows, switch_key_x.cols);
    println!("do dcr_dec_multi for {} y {} times", switch_key_y.rows, switch_key_y.cols);
    println!("do dcr_dec_multi for {} h {} times", switch_key_h.rows, switch_key_h.cols);
    let ct_out_x = dcr_dec_multi(&ct_in, &switch_key_x, &switch_key_dcr_x, &decomp, &grp);
    let ct_out_y = dcr_dec_multi(&ct_in, &switch_key_y, &switch_key_dcr_y, &decomp, &grp);
    let ct_out_h = dcr_dec_multi(&ct_in, &switch_key_h, &switch_key_dcr_h, &decomp, &grp);
    (ct_out_x, ct_out_y, ct_out_h)
}

pub fn protocol_keygen_i(
    qe_sk_enc: &QeSk,
    qe_sk_keygen: &QeSk,
    h_right: &Matrix,
    hm_left: &Matrix,
    dim: usize,
    k: usize,
    f: &Matrix,
    decomp: &Decomp,
    grp: &Group,
    rng: &mut RandState<'_>,
) -> (
    (Vec<Integer>, Vec<Integer>, Vec<Integer>), 
    (Matrix, Matrix, Matrix), 
    (Matrix, Matrix, Matrix),
) {
    let (qe_enc_mat_x, qe_enc_mat_y, qe_enc_mat_h) = qe_enc_matrix_same_xy(&qe_sk_enc, dim + k, &grp, &decomp, rng, false);
    
    // divide enc_mats into M * b
    pub fn divide_mat_into_m_b(
        mat: &Matrix,
    ) -> (Matrix, Vec<Integer>) {
        let mut mat_left = Matrix::new(mat.rows, mat.cols - 1);
        let mut vec_right = vec![Integer::from(0); mat.rows];

        for i in 0..mat.rows {
            for j in 0..mat.cols - 1 {
                let val = mat.get(i, j);
                mat_left.set(i, j, val);
            }
            let val = mat.get(i, mat.cols - 1);
            vec_right[i] = val;
        }
        (mat_left, vec_right)
    }
    let (qe_enc_mat_x, qe_b_x) = divide_mat_into_m_b(&qe_enc_mat_x);
    let (qe_enc_mat_y, qe_b_y) = divide_mat_into_m_b(&qe_enc_mat_y);
    let (qe_enc_mat_h, qe_b_h) = divide_mat_into_m_b(&qe_enc_mat_h);

    let qe_b_x = decomp.vector(&qe_b_x);
    let qe_b_y = decomp.vector(&qe_b_y);
    let qe_b_h = decomp.vector(&qe_b_h);

    let h_right_origin = remove_diag_one(&h_right);
    let hm_left_origin = remove_diag_one(&hm_left);
    let hmhm = Matrix::tensor_product(&hm_left_origin, &hm_left_origin, &grp.delta);

    fn mat_mul_4(
        a: &Matrix,
        b: &Matrix,
        c: &Matrix,
        d: &Matrix,
        decomp: &Decomp,
        grp: &Group,
    ) -> Matrix {
        // A = Enc, B = H, C = F, D = (H' tensor H')
        // output decomp(A * B) * (C * D)
        let modulo = grp.delta.clone();
        let mut ab = a * b;
        ab.mod_inplace(&modulo);
        ab = decomp.matrix_col(&ab);
        let mut cd = c * d;
        cd.mod_inplace(&modulo);
        let mut out = ab * cd;
        out.mod_inplace(&modulo);
        out
    }

    let total_mat_x = mat_mul_4(
        &qe_enc_mat_x, &h_right_origin, &f, &hmhm, &decomp, &grp);
    let total_mat_y = mat_mul_4(
        &qe_enc_mat_y, &h_right_origin, &f, &hmhm, &decomp, &grp);
    let total_mat_h = mat_mul_4(
        &qe_enc_mat_h, &h_right_origin, &f, &hmhm, &decomp, &grp);
    
    fn gen_f_and_red(
        qe_sk: &QeSk,
        total_mat: Matrix,
        grp: &Group,
        decomp: &Decomp,
    ) -> (Matrix, Matrix) {
        let mut sk_f_mat = Matrix::new(1, 1);
        let mut sk_red_mat = Matrix::new(1, 1);
        println!("do qe_keygen of dim {} for {} times", total_mat.cols, total_mat.rows);
        let start = SystemTime::now();
        for i in 0..total_mat.rows {
            let row = total_mat.get_row(i);
            let (sk_f, sk_red) = qe_keygen(&qe_sk, &row, grp, decomp);
            if i == 0 {
                sk_f_mat = Matrix::new(total_mat.rows, sk_f.len());
                sk_red_mat = Matrix::new(total_mat.rows, sk_red.len());
            }
            sk_f_mat.set_row(i, &sk_f);
            sk_red_mat.set_row(i, &sk_red);
        }
        let end = SystemTime::now();
        let elapsed = end.duration_since(start).unwrap();
        println!("qe_keygen for {} times, time: {:?}", total_mat.rows, elapsed);
        (sk_f_mat, sk_red_mat)
    }

    let (sk_f_mat_x, sk_red_mat_x) = gen_f_and_red(&qe_sk_keygen, total_mat_x, &grp, &decomp);
    let (sk_f_mat_y, sk_red_mat_y) = gen_f_and_red(&qe_sk_keygen, total_mat_y, &grp, &decomp);
    let (sk_f_mat_h, sk_red_mat_h) = gen_f_and_red(&qe_sk_keygen, total_mat_h, &grp, &decomp);
    (
        (qe_b_x, qe_b_y, qe_b_h), 
        (sk_f_mat_x, sk_f_mat_y, sk_f_mat_h), 
        (sk_red_mat_x, sk_red_mat_y, sk_red_mat_h),
    )
}

pub fn protocol_dec_i(
    ctxt_triple: (&Vec<Integer>, &Vec<Integer>, &Vec<Integer>),
    (qe_b_x, qe_b_y, qe_b_h): (&Vec<Integer>, &Vec<Integer>, &Vec<Integer>),
    (sk_f_mat_x, sk_f_mat_y, sk_f_mat_h): (&Matrix, &Matrix, &Matrix),
    (sk_red_mat_x, sk_red_mat_y, sk_red_mat_h): (&Matrix, &Matrix, &Matrix),
    decomp: &Decomp,
    grp: &Group,
) -> (Vec<Integer>, Vec<Integer>, Vec<Integer>) {
    fn compute_f_red_out (
        sk_f_mat: &Matrix,
        sk_red_mat: &Matrix,
        ctxt_triple: (&Vec<Integer>, &Vec<Integer>, &Vec<Integer>),
        qe_b: &Vec<Integer>,
        grp: &Group,
        decomp: &Decomp,
    ) -> Vec<Integer> {
        let mut ct_out = vec![Integer::from(0); sk_f_mat.rows];
        println!("run qe_dec of dim {} for {} times", sk_f_mat.cols, sk_f_mat.rows);
        for i in 0..sk_f_mat.rows {
            let sk_f = sk_f_mat.get_row(i);
            let sk_red = sk_red_mat.get_row(i);
            ct_out[i] =qe_dec((&sk_f, &sk_red), ctxt_triple, grp, decomp);
        }
        ct_out = vec_add(&ct_out, &qe_b);
        vec_mod(&mut ct_out, &grp.n);
        ct_out
    }

    let ct_out_x = compute_f_red_out(
        &sk_f_mat_x, &sk_red_mat_x,
        ctxt_triple,
        qe_b_x,
        grp,
        decomp,
    );
    let ct_out_y = compute_f_red_out(
        &sk_f_mat_y, &sk_red_mat_y,
        ctxt_triple,
        qe_b_y,
        grp,
        decomp,
    );
    let ct_out_h = compute_f_red_out(
        &sk_f_mat_h, &sk_red_mat_h,
        ctxt_triple,
        qe_b_h,
        grp,
        decomp,
    );
    (ct_out_x, ct_out_y, ct_out_h)
}

pub fn protocol_keygen_end(
    qe_sk: &QeSk,
    hm_left: &Matrix,
    f: &Matrix,
    decomp: &Decomp,
    grp: &Group,
) -> (Matrix, Matrix) {
    let modulo = grp.delta.clone();

    let hm_origin = remove_diag_one(&hm_left);
    let hmhm = Matrix::tensor_product(&hm_origin, &hm_origin, &grp.delta);
    
    let mut fhmhm = f * &hmhm;
    fhmhm.mod_inplace(&modulo);

    let mut sk_f_mat = Matrix::new(1, 1);
    let mut sk_red_mat = Matrix::new(1, 1);
    println!("do qe_keygen of dim {} for {} times", fhmhm.cols, fhmhm.rows);
    for i in 0..fhmhm.rows {
        let row = fhmhm.get_row(i);
        let (sk_f, sk_red) = qe_keygen(&qe_sk, &row, grp, decomp);
        if i == 0 {
            sk_f_mat = Matrix::new(fhmhm.rows, sk_f.len());
            sk_red_mat = Matrix::new(fhmhm.rows, sk_red.len());
        }
        sk_f_mat.set_row(i, &sk_f);
        sk_red_mat.set_row(i, &sk_red);
    }
    (sk_f_mat, sk_red_mat)
}

pub fn protocol_dec_end(
    ctxt_triple: (&Vec<Integer>, &Vec<Integer>, &Vec<Integer>),
    (sk_f_mat, sk_red_mat): (&Matrix, &Matrix),
    decomp: &Decomp,
    grp: &Group,
) -> Vec<Integer> {
    let mut ct_out = vec![Integer::from(0); sk_f_mat.rows];
    println!("run qe_dec of dim {} for {} times", sk_f_mat.cols, sk_f_mat.rows);
    for i in 0..sk_f_mat.rows {
        let sk_f = sk_f_mat.get_row(i);
        let sk_red = sk_red_mat.get_row(i);
        ct_out[i] =qe_dec((&sk_f, &sk_red), ctxt_triple, grp, decomp);
    }
    ct_out
}
