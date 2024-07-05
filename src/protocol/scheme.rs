extern crate rug;
use rug::Integer;
use rug::rand::RandState;

use crate::qfe;
use crate::util::group::Group;
use crate::util::matrix::{concatenate_row, concatenate_vec_row, remove_diag_one, Matrix};
use crate::util::vector;
use crate::util::vector::{gen_random_vector, vec_add, vec_mod, vec_mul_scalar, int_mod};
use crate::util::decomp::Decomp;
use crate::dcr::scheme::{dcr_setup, dcr_enc, dcr_keygen, dcr_dec};
use crate::qfe::keys::QfeSk;
use crate::qfe::scheme::{divide_vector_for_functional_key, get_funcional_key_len, qfe_dec, qfe_enc_matrix_same_xy, qfe_keygen, qfe_setup, get_ctxt_len};
use crate::ipfe::scheme::{ipfe_enc, ipfe_keygen};
use std::time::SystemTime;

pub fn protocol_setup(
    dim_vec: Vec<usize>,
    f_num: usize,
    k: usize,
    sk_bound: &Integer,
    grp: &Group,
    rng: &mut RandState<'_>,
) -> ((Vec<Integer>, Vec<Integer>), QfeSk, Vec<QfeSk>, QfeSk) {
    let dim = dim_vec[0];
    let (dcr_sk, dcr_pk) = dcr_setup(2 * dim + 1, sk_bound, grp, rng);
    let qfe_sk_init = qfe_setup(grp, dim + k + 1, 2 * (dim + k + 1) + 1, rng);
    let mut qfe_sk_fcn = Vec::with_capacity(f_num);
    for i in 1..dim_vec.len()-1 {
        let dim = dim_vec[i];
        qfe_sk_fcn.push(qfe_setup(grp, dim + k, 2 * (dim + k) + 1, rng));
    }
    let dim = dim_vec[dim_vec.len()-1];
    let qfe_sk_end = qfe_setup(grp, dim + k, 2 * (dim + k) + 1, rng);

    ((dcr_sk, dcr_pk), qfe_sk_init, qfe_sk_fcn, qfe_sk_end)
}

pub fn protocol_enc_init(
    dcr_pk: &Vec<Integer>,
    gamma_right: &Matrix,
    x: &Vec<Integer>,
    grp: &Group,
    rng: &mut RandState<'_>,
) -> Vec<Integer> {
    let mut x1 = x.clone();
    vec_mod(&mut x1, &grp.delta);
    x1.push(Integer::from(1));

    let mut gamma_right_x = gamma_right.mul_vec(&x1);
    vec_mod(&mut gamma_right_x, &grp.delta);
    dcr_enc(dcr_pk, &gamma_right_x, &grp, rng)
}

pub fn protocol_keygen_switch(
    qfe_sk: &QfeSk,
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
    let (qfe_enc_mat_x, qfe_enc_mat_y, qfe_enc_mat_h)
        = qfe_enc_matrix_same_xy(&qfe_sk, dim + k, &grp, rng);
    
    fn gen_switch_key_for_enc_and_dcr(
        dcr_sk: &Vec<Integer>,
        qfe_enc_mat: &Matrix,
        h_right: &Matrix,
        gamma_left: &Matrix,
        decomp: &Decomp,
        modulo: &Integer,
    ) -> (Matrix, Vec<Integer>) {
        let mut xh = qfe_enc_mat * h_right;
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
        &qfe_enc_mat_x,
        &h_right,
        &gamma_left,
        &decomp,
        &modulo,
    );

    let (switch_key_y, switch_key_dcr_y) = gen_switch_key_for_enc_and_dcr(
        &dcr_sk,
        &qfe_enc_mat_y,
        &h_right,
        &gamma_left,
        &decomp,
        &modulo,
    );

    let (switch_key_h, switch_key_dcr_h) = gen_switch_key_for_enc_and_dcr(
        &dcr_sk,
        &qfe_enc_mat_h,
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
    qfe_sk_enc: &QfeSk,
    qfe_sk_keygen: &QfeSk,
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
    let (qfe_enc_mat_x, qfe_enc_mat_y, qfe_enc_mat_h) = qfe_enc_matrix_same_xy(&qfe_sk_enc, dim + k, &grp, rng);
    
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
    let (qfe_enc_mat_x, qfe_b_x) = divide_mat_into_m_b(&qfe_enc_mat_x);
    let (qfe_enc_mat_y, qfe_b_y) = divide_mat_into_m_b(&qfe_enc_mat_y);
    let (qfe_enc_mat_h, qfe_b_h) = divide_mat_into_m_b(&qfe_enc_mat_h);

    let qfe_b_x = decomp.vector(&qfe_b_x);
    let qfe_b_y = decomp.vector(&qfe_b_y);
    let qfe_b_h = decomp.vector(&qfe_b_h);

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
        &qfe_enc_mat_x, &h_right_origin, &f, &hmhm, &decomp, &grp);
    let total_mat_y = mat_mul_4(
        &qfe_enc_mat_y, &h_right_origin, &f, &hmhm, &decomp, &grp);
    let total_mat_h = mat_mul_4(
        &qfe_enc_mat_h, &h_right_origin, &f, &hmhm, &decomp, &grp);
    
    fn gen_f_and_red(
        qfe_sk: &QfeSk,
        total_mat: Matrix,
        grp: &Group,
    ) -> (Matrix, Matrix) {
        let mut sk_f_mat = Matrix::new(1, 1);
        let mut sk_red_mat = Matrix::new(1, 1);
        println!("do qfe_keygen of dim {} for {} times", total_mat.cols, total_mat.rows);
        let start = SystemTime::now();
        for i in 0..total_mat.rows {
            let row = total_mat.get_row(i);
            let fk = qfe_keygen(&qfe_sk, &row, grp);
            let (sk_f, sk_red) = divide_vector_for_functional_key(&fk, qfe_sk.dim, qfe_sk.q);
            if i == 0 {
                sk_f_mat = Matrix::new(total_mat.rows, sk_f.len());
                sk_red_mat = Matrix::new(total_mat.rows, sk_red.len());
            }
            sk_f_mat.set_row(i, &sk_f);
            sk_red_mat.set_row(i, &sk_red);
        }
        let end = SystemTime::now();
        let elapsed = end.duration_since(start).unwrap();
        println!("qfe_keygen for {} times, time: {:?}", total_mat.rows, elapsed);
        (sk_f_mat, sk_red_mat)
    }

    let (sk_f_mat_x, sk_red_mat_x) = gen_f_and_red(&qfe_sk_keygen, total_mat_x, &grp);
    let (sk_f_mat_y, sk_red_mat_y) = gen_f_and_red(&qfe_sk_keygen, total_mat_y, &grp);
    let (sk_f_mat_h, sk_red_mat_h) = gen_f_and_red(&qfe_sk_keygen, total_mat_h, &grp);
    (
        (qfe_b_x, qfe_b_y, qfe_b_h), 
        (sk_f_mat_x, sk_f_mat_y, sk_f_mat_h), 
        (sk_red_mat_x, sk_red_mat_y, sk_red_mat_h),
    )
}

pub fn protocol_dec_i(
    ctxt_triple: (&Vec<Integer>, &Vec<Integer>, &Vec<Integer>),
    (qfe_b_x, qfe_b_y, qfe_b_h): (&Vec<Integer>, &Vec<Integer>, &Vec<Integer>),
    (sk_f_mat_x, sk_f_mat_y, sk_f_mat_h): (&Matrix, &Matrix, &Matrix),
    (sk_red_mat_x, sk_red_mat_y, sk_red_mat_h): (&Matrix, &Matrix, &Matrix),
    dim: usize,
    q: usize,
    decomp: &Decomp,
    grp: &Group,
) -> (Vec<Integer>, Vec<Integer>, Vec<Integer>) {
    fn compute_f_red_out(
        sk_f_mat: &Matrix,
        sk_red_mat: &Matrix,
        ctxt_triple: (&Vec<Integer>, &Vec<Integer>, &Vec<Integer>),
        qfe_b: &Vec<Integer>,
        dim: usize,
        q: usize,
        grp: &Group,
    ) -> Vec<Integer> {
        let mut ct_out = vec![Integer::from(0); sk_f_mat.rows];
        // println!("run qfe_dec of dim {} for {} times", sk_f_mat.cols, sk_f_mat.rows);
        for i in 0..sk_f_mat.rows {
            let sk_f = sk_f_mat.get_row(i);
            // let sk_f = decomp.vector_pow_exp(&sk_f);
            let sk_red = sk_red_mat.get_row(i);
            // fk = sk_f || sk_red (concatenation)
            let mut fk = Vec::with_capacity(sk_f.len() + sk_red.len());
            fk.extend_from_slice(&sk_f);
            fk.extend_from_slice(&sk_red);
            // println!("fk.len: {}", fk.len());
            // println!(" = sk_f.len: {} + sk_red.len: {}", sk_f.len(), sk_red.len());

            let enc_x = ctxt_triple.0.clone();
            let enc_y = ctxt_triple.1.clone();
            let enc_h = ctxt_triple.2.clone();
            // ctxt = (enc_x, enc_y, enc_h)
            let mut ctxt = Vec::with_capacity(enc_x.len() + enc_y.len() + enc_h.len());
            ctxt.extend_from_slice(&enc_x);
            ctxt.extend_from_slice(&enc_y);
            ctxt.extend_from_slice(&enc_h);
            ct_out[i] =qfe_dec(
                &fk, &ctxt, dim + 1, 2 * (dim + 1) + 1, grp,
            );
        }
        ct_out = vec_add(&ct_out, &qfe_b);
        vec_mod(&mut ct_out, &grp.n);
        ct_out
    }

    let ct_out_x = compute_f_red_out(
        &sk_f_mat_x, &sk_red_mat_x,
        ctxt_triple,
        qfe_b_x,
        dim,
        q,
        grp,
    );
    let ct_out_y = compute_f_red_out(
        &sk_f_mat_y, &sk_red_mat_y,
        ctxt_triple,
        qfe_b_y,
        dim,
        q,
        grp,
    );
    let ct_out_h = compute_f_red_out(
        &sk_f_mat_h, &sk_red_mat_h,
        ctxt_triple,
        qfe_b_h,
        dim,
        q,
        grp,
    );
    (ct_out_x, ct_out_y, ct_out_h)
}

pub fn protocol_keygen_end(
    qfe_sk: &QfeSk,
    hm_left: &Matrix,
    f: &Matrix,
    grp: &Group,
) -> (Matrix, Matrix) {
    let modulo = grp.delta.clone();

    let hm_origin = remove_diag_one(&hm_left);
    let hmhm = Matrix::tensor_product(&hm_origin, &hm_origin, &grp.delta);
    
    let mut fhmhm = f * &hmhm;
    fhmhm.mod_inplace(&modulo);

    let mut sk_f_mat = Matrix::new(1, 1);
    let mut sk_red_mat = Matrix::new(1, 1);
    println!("do qfe_keygen of dim {} for {} times", fhmhm.cols, fhmhm.rows);
    for i in 0..fhmhm.rows {
        let row = fhmhm.get_row(i);
        let fk = qfe_keygen(&qfe_sk, &row, grp);
        let (sk_f, sk_red) = divide_vector_for_functional_key(&fk, qfe_sk.dim, qfe_sk.q);
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
    dim: usize,
    q: usize,
    grp: &Group,
) -> Vec<Integer> {
    let mut ct_out = vec![Integer::from(0); sk_f_mat.rows];
    println!("run qfe_dec of dim {} for {} times", sk_f_mat.cols, sk_f_mat.rows);
    for i in 0..sk_f_mat.rows {
        let sk_f = sk_f_mat.get_row(i);
        // let sk_f = decomp.vector_pow_exp(&sk_f);
        let sk_red = sk_red_mat.get_row(i);

        let enc_x = ctxt_triple.0.clone();
        let enc_y = ctxt_triple.1.clone();
        let enc_h = ctxt_triple.2.clone();
        let enc_x = decomp.vector_inv(&enc_x);
        let enc_y = decomp.vector_inv(&enc_y);
        let enc_h = decomp.vector_inv(&enc_h);

        let mut fk = Vec::with_capacity(sk_f.len() + sk_red.len());
        fk.extend_from_slice(&sk_f);
        fk.extend_from_slice(&sk_red);

        let mut ctxt = Vec::with_capacity(enc_x.len() + enc_y.len() + enc_h.len());
        ctxt.extend_from_slice(&enc_x);
        ctxt.extend_from_slice(&enc_y);
        ctxt.extend_from_slice(&enc_h);
        println!("ctxt size = {}", ctxt.len());
        println!(" = enc_x.len: {} + enc_y.len: {} + enc_h.len: {}", enc_x.len(), enc_y.len(), enc_h.len());
        ct_out[i] = qfe_dec(
            &fk, &ctxt, dim + 1, 2 * (dim + 1) + 1, &grp,
        );
    }
    ct_out
}

pub fn protocol_setup_new(
    dim_vec: Vec<usize>,
    f_num: usize,
    k: usize,
    sk_bound: &Integer,
    grp: &Group,
    rng: &mut RandState<'_>,
) -> ((Vec<Integer>, Vec<Integer>), QfeSk, Vec<QfeSk>, QfeSk) {
    let dim = dim_vec[0];
    let (dcr_sk, dcr_pk) = dcr_setup(2 * dim + 1, sk_bound, grp, rng);
    let qfe_sk_dcr_to_qe = qfe_setup(grp, dim + k + 1, 2 * (dim + k + 1) + 1, rng);
    let mut qfe_sk_qe_to_qe = Vec::new();
    for i in 0..f_num {
        let dim = dim_vec[0];
        qfe_sk_qe_to_qe.push(qfe_setup(grp, dim + k + 1, 2 * (dim + k + 1) + 1, rng));
    }
    let dim: usize = dim_vec[0];
    let qfe_sk_end = qfe_setup(grp, dim + k + 1, 2 * (dim + k + 1) + 1, rng);

    ((dcr_sk, dcr_pk), qfe_sk_dcr_to_qe, qfe_sk_qe_to_qe, qfe_sk_end)
}

pub fn composite_enc_and_f(
    qfe_sk: &QfeSk,
    f: &Matrix,
    grp: &Group,
    dcp: &Decomp,
    rng: &mut RandState<'_>,
) -> Matrix {
    println!("start composite_enc_and_f");
    let dim = qfe_sk.dim;
    let q = qfe_sk.q;
    println!("input f size = {} x {}", f.rows, f.cols);
    println!("decomposed f size = {} x {}", f.rows, f.cols);

    let modulo = grp.delta.clone();

    let mut mat_ctxts = Matrix::new(f.cols, get_ctxt_len(dim, q));
    for i in 0..f.cols {
        let f_col = f.get_col(i);
        let mu_f_col = vec_mul_scalar(&f_col, &grp.mu);
        let rand_x = gen_random_vector(q, &modulo, rng);
        let rand_y = gen_random_vector(q, &modulo, rng);
        let r_x = gen_random_vector(2, &modulo, rng);
        let r_y = gen_random_vector(2, &modulo, rng);

        // remove random for testing
        // println!("remove random for testing");
        // let rand_x = vec![Integer::from(0); q];
        // let rand_y = vec![Integer::from(0); q];
        // let r_x = vec![Integer::from(0); 2];
        // let r_y = vec![Integer::from(0); 2];
        
        let mut r_x = vec_mul_scalar(&r_x, &(Integer::from(2) * &grp.n));
        let mut r_y = vec_mul_scalar(&r_y, &(Integer::from(2) * &grp.n));
        vec_mod(&mut r_x, &modulo);
        vec_mod(&mut r_y, &modulo);
        
        // ct0 = d_x_null * rand_x + d_x_inv * mu * f_col
        // ct1 = d_y_null * rand_y + d_y_inv * f_col
        // if i = f.cols:
        //  ct0 = d_x_null * rand_x + d_x_inv * (mu * f_col + sk->V * r_x)
        //  ct1 = d_y_null * rand_y + d_y_inv * (f_col + sk->W * r_y)
        let ct0_left = qfe_sk.d_x_null.mul_vec(&rand_x);
        let ct0_right = if i < f.cols - 1 {
            qfe_sk.d_x_inv.mul_vec(&mu_f_col)
        } else {
            qfe_sk.d_x_inv.mul_vec(
                &vec_add(&mu_f_col, &qfe_sk.v.mul_vec(&r_x))
            )
        };
        let mut ct0: Vec<Integer> = vec_add(&ct0_left, &ct0_right);
        vec_mod(&mut ct0, &modulo);

        let ct1_left = qfe_sk.d_y_null.mul_vec(&rand_y);
        let ct1_right = if i < f.cols - 1 {
            qfe_sk.d_y_inv.mul_vec(&f_col)
        } else {
            qfe_sk.d_y_inv.mul_vec(
                &vec_add(&f_col, &qfe_sk.w.mul_vec(&r_y))
            )
        };
        let mut ct1 = vec_add(&ct1_left, &ct1_right);
        vec_mod(&mut ct1, &modulo);

        // h = (r_x tensor f_col) || (mu_f_col tensor r_y)
        // if i == f.cols:
        //   h = (r_x tensor f_col) || ((mu_f_col + sk->v * r_x) tensor r_y)
        let h_left = vector::tensor_product_vecs(&r_x, &f_col, &modulo);
        let mut tmp = if i < f.cols - 1 {
            mu_f_col.clone()
        } else {
            vec_add(&mu_f_col, &qfe_sk.v.mul_vec(&r_x))
        };
        vec_mod(&mut tmp, &modulo);

        let h_right = vector::tensor_product_vecs(&tmp, &r_y, &modulo);

        let mut h: Vec<Integer> = Vec::with_capacity(h_left.len() + h_right.len());
        h.extend_from_slice(&h_left);
        h.extend_from_slice(&h_right);
        // println!("h_left = {:?}", h_left);
        // println!("h = {:?}", h);
        // println!("=============");
        let ctxt_ipfe = ipfe_enc(&qfe_sk.ipfe_sk, &h, grp, false, rng);

        // ith row of mat_ctxts = (ct0, ct1, ctxt_ipfe)
        let mut row = Vec::with_capacity(ct0.len() + ct1.len() + ctxt_ipfe.len());
        row.extend_from_slice(&ct0);
        row.extend_from_slice(&ct1);
        row.extend_from_slice(&ctxt_ipfe);
        mat_ctxts.set_row(i, &row);
    }

    let mat_ctxts = dcp.matrix_row(&mat_ctxts);

    println!("end composite_enc_and_f");
    println!("output mat size = {} x {}", mat_ctxts.rows, mat_ctxts.cols);
    // output size = L * (6 * dim + 3 * q + 2) x ((dim + 1)^2 + 1)
    mat_ctxts.transpose()
}

pub fn protocol_keygen_dcr_to_qfe(
    dcr_sk: &Vec<Integer>,
    qfe_sk: &QfeSk,
    h_right: &Matrix,
    gamma_left: &Matrix,
    decomp: &Decomp,
    grp: &Group,
    rng: &mut RandState<'_>,
) -> (Matrix, Vec<Integer>) {
    let mut total_mat = h_right * gamma_left;
    total_mat.mod_inplace(&grp.n);
    
    // composite enc and f
    let fk_mat = composite_enc_and_f(
        qfe_sk, &total_mat, grp, decomp, rng);

    // keygen_dcr for each row of mat_ctxts
    let mut fk = vec![Integer::from(0); fk_mat.rows];
    for i in 0..fk_mat.rows {
        let row = fk_mat.get_row(i);
        fk[i] = dcr_keygen(dcr_sk, &row);
    }
    println!("done keygen_dcr_to_qfe");

    (fk_mat, fk)
}

pub fn protocol_keygen_qfe_to_qfe(
    qfe_sk_enc: &QfeSk,
    qfe_sk_keygen: &QfeSk,
    h_right: &Matrix,
    hm_left: &Matrix,
    f: &Matrix,
    decomp: &Decomp,
    grp: &Group,
    rng: &mut RandState<'_>,
) -> Matrix {
    // function: h_right * f * (hm_left tensor hm-left)
    // let hm_origin = remove_diag_one(&hm_left);
    // let hmhm = Matrix::tensor_product(&hm_origin, &hm_origin, &grp.delta);
    let hmhm = Matrix::tensor_product(&hm_left, &hm_left, &grp.delta);
    println!("h_right dim = {} x {}", h_right.rows, h_right.cols);
    println!("f dim = {} x {}", f.rows, f.cols);
    println!("hmhm dim = {} x {}", hmhm.rows, hmhm.cols);
    let total_mat = h_right * f * &hmhm;

    // composite enc and f
    let mat_ctxts = composite_enc_and_f(
        qfe_sk_enc, &total_mat, grp, decomp, rng);   
    // keygen_qfe for each row of mat_ctxts
    let mut fk_mat = Matrix::new(
        mat_ctxts.rows,
        get_funcional_key_len(qfe_sk_keygen.dim, qfe_sk_keygen.q)
    );
    println!("run qfe_keygen of dim {} for {} times", mat_ctxts.cols, mat_ctxts.rows);
    for i in 0..mat_ctxts.rows {
        let row = mat_ctxts.get_row(i);
        let fk = qfe_keygen(qfe_sk_keygen, &row, grp);
        fk_mat.set_row(i, &fk);
    }
    
    fk_mat
}

pub fn protocol_keygen_qfe_to_plain(
    qfe_sk: &QfeSk,
    hm_left: &Matrix,
    f: &Matrix,
    decomp: &Decomp,
    grp: &Group,
) -> Matrix {
    let hmhm = Matrix::tensor_product(&hm_left, &hm_left, &grp.delta);
    println!("f size = {} x {}", f.rows, f.cols);
    println!("hmhm size = {} x {}", hmhm.rows, hmhm.cols);
    let mut total_mat = f * &hmhm;
    total_mat.mod_inplace(&grp.delta);
    // let total_mat = decomp.matrix_col(&total_mat);

    let mut fk_mat = Matrix::new(
        total_mat.rows,
        get_funcional_key_len(qfe_sk.dim, qfe_sk.q)
    );
    for i in 0..total_mat.rows {
        let row = total_mat.get_row(i);
        let fk = qfe_keygen(qfe_sk, &row, grp);
        fk_mat.set_row(i, &fk);
    }
    fk_mat
}


pub fn protocol_dec_dcr_to_qfe(
    ctxt: &Vec<Integer>,
    fk_mat: &Matrix,
    fk_vec: &Vec<Integer>,
    decomp: &Decomp,
    grp: &Group,
) -> Vec<Integer> {
    let mut ct_out = vec![Integer::from(0); fk_mat.rows];
    println!("run dcr_dec of dim {} for {} times", fk_mat.cols, fk_mat.rows);
    for i in 0..fk_mat.rows {
        ct_out[i] = dcr_dec(ctxt, &fk_mat.get_row(i), &fk_vec[i], grp);
    }
    vec_mod(&mut ct_out, &grp.n);
    decomp.vector_inv(&ct_out)
}

pub fn protocol_dec_qfe(
    ctxt: &Vec<Integer>,
    fk_mat: &Matrix,
    dim: usize,
    q: usize,
    decomp: &Decomp,
    grp: &Group,
    is_output_decomposed: bool,
) -> Vec<Integer> {
    println!("fk_mat size in dec_qfe = {} x {}", fk_mat.rows, fk_mat.cols);
    let mut ct_out = vec![Integer::from(0); fk_mat.rows];
    println!("run qfe_dec of dim {} for {} times", fk_mat.cols, fk_mat.rows);
    for i in 0..fk_mat.rows {
        let fk = fk_mat.get_row(i);
        ct_out[i] = qfe_dec(&fk, &ctxt, dim + 1, 2 * (dim + 1) + 1, grp);
    }
    vec_mod(&mut ct_out, &grp.n);
    if is_output_decomposed {
        decomp.vector_inv(&ct_out)
    } else {
        ct_out
    }
    
}

// pub fn protocol_dec_qfe_to_plain(
//     ctxt: &Vec<Integer>,
//     fk_mat: &Matrix,
//     dim: usize,
//     q: usize,
//     grp: &Group,
// ) -> Vec<Integer> {
//     let mut res = vec![Integer::from(0); fk_mat.rows];
//     println!("run qfe_dec of dim {} for {} times", fk_mat.cols, fk_mat.rows);
//     for i in 0..fk_mat.rows {
//         let fk = fk_mat.get_row(i);
//         res[i] = qfe_dec(&fk, &ctxt, dim + 1, 2 * (dim + 1) + 1, grp);
//     }
//     vec_mod(&mut res, &grp.n);
//     res
// }