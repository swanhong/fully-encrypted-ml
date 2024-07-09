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
use crate::qfe::scheme::{divide_vector_for_functional_key, get_ctxt_len, get_funcional_key_len, qfe_dec, qfe_enc, qfe_cenc, qfe_enc_for_x, qfe_enc_for_h, qfe_enc_matrix_same_xy, qfe_keygen, qfe_setup};
use crate::ipfe::scheme::{ipfe_enc, ipfe_keygen};
use std::time::SystemTime;

pub fn protocol_setup(
    dim_vec: &Vec<usize>,
    k: usize,
    sk_bound: &Integer,
    grp: &Group,
    rng: &mut RandState<'_>,
) -> (
    (Vec<Integer>, Vec<Integer>),
    Vec<QfeSk>,
) {
    let (dcr_sk, dcr_pk) = dcr_setup(2 * dim_vec[0] + 1, sk_bound, grp, rng);
    let mut qfe_sk = Vec::new();
    for i in 0..dim_vec.len() - 1  {
        qfe_sk.push(qfe_setup(grp, dim_vec[i] + k + 1, 2 * (dim_vec[i] + k + 1) + 1, rng));
    }
    (
        (dcr_sk, dcr_pk),
        qfe_sk
    )
}
pub fn protocol_enc_init(
    dcr_pk: &Vec<Integer>,
    gamma_right: &Matrix,
    x: &Vec<Integer>,
    grp: &Group,
    rng: &mut RandState<'_>,
) -> Vec<Integer> {
    let mut x1 = x.clone();
    vec_mod(&mut x1, &grp.n);
    x1.push(Integer::from(1));
    let gamma_right_x = gamma_right.mul_vec(&x1);
    dcr_enc(dcr_pk, &gamma_right_x, &grp, rng)
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
    let fk_mat = qfe_cenc(&qfe_sk, &total_mat, grp, rng);
    let fk_mat = decomp.matrix_col(&fk_mat);

    // keygen_dcr for each row of mat_ctxts
    let mut fk = vec![Integer::from(0); fk_mat.rows];
    for i in 0..fk_mat.rows {
        let row = fk_mat.get_row(i);
        fk[i] = dcr_keygen(dcr_sk, &row);
    }
    (fk_mat, fk)
}

pub fn protocol_keygen_qfe_to_qfe(
    qfe_sk_keygen: &QfeSk,
    qfe_sk_enc: &QfeSk,
    h_right: &Matrix,
    hm_left: &Matrix,
    f: &Matrix,
    decomp: &Decomp,
    grp: &Group,
    rng: &mut RandState<'_>,
) -> Matrix {
    // function: h_right * f * (hm_left tensor hm_left)
    let hmhm = Matrix::tensor_product(&hm_left, &hm_left, &grp.delta);
    println!("h_right = {} x {}", h_right.rows, h_right.cols);
    println!("f = {} x {}", f.rows, f.cols);
    println!("hmhm = {} x {}", hmhm.rows, hmhm.cols);
    let mut total_mat: Matrix = h_right * &(f * &hmhm);
    total_mat.mod_inplace(&grp.n);

    // composite enc and f
    let mat_ctxts = qfe_cenc(&qfe_sk_enc, &total_mat, grp, rng);
    let mat_ctxts = decomp.matrix_col(&mat_ctxts);
    
    // keygen_qfe for each row of mat_ctxts
    println!("run qfe_keygen of dim {} for {} times", mat_ctxts.cols, mat_ctxts.rows);
    let mut fk_mat = Matrix::new(
        mat_ctxts.rows,
        get_funcional_key_len(qfe_sk_keygen.dim, qfe_sk_keygen.q)
    );
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
    grp: &Group,
) -> Matrix {
    let hmhm = Matrix::tensor_product(&hm_left, &hm_left, &grp.delta);
    println!("hmhm size = {} x {}", hmhm.rows, hmhm.cols);
    println!("f size = {} x {}", f.rows, f.cols);
    let mut total_mat = f * &hmhm;
    total_mat.mod_inplace(&grp.n);

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
        ct_out[i] = qfe_dec(&fk, &ctxt, dim, q, grp);
    }
    vec_mod(&mut ct_out, &grp.n);
    if is_output_decomposed {
        decomp.vector_inv(&ct_out)
    } else {
        ct_out
    }
    
}