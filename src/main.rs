extern crate rug;

mod util;
mod parse;
mod dcr_ipe;
mod ipe;
mod qe;
mod protocol;
mod test;

use rug::{Integer, Float};
use rug::ops::Pow;
use rug::rand::RandState;
use std::time::SystemTime;
use crate::protocol::scheme::{protocol_setup, protocol_keygen_switch, protocol_enc_init, protocol_keyswitch, protocol_keygen_i, protocol_dec_i, protocol_keygen_end, protocol_dec_end};
use crate::util::group::Group;
use crate::util::matrix::{Matrix, eval_quadratic_multivariate, concatenate_col, concatenate_vec_row};
use crate::util::vector::{gen_random_vector, vec_mod, vec_mul_scalar};
use crate::util::decomp::Decomp;
use crate::protocol::algorithm::{sample_h, sample_gamma};
use crate::parse::file_reader::read_matrix;

fn _run_protocol_start_to_end() {
    println!("run test_ipe_start_to_end");

    let mut rng = RandState::new(); // Create a single RandState object
    let d = SystemTime::now()
    .duration_since(SystemTime::UNIX_EPOCH)
    .expect("Duration since UNIX_EPOCH failed");
    rng.seed(&Integer::from(d.as_secs()));
    
    let grp = Group::new(10); // Initialize the group        
    // println!("{}", grp);
    
    let dim = 2;
    let dim2 = 2;
    let dim3 = 2;
    let k = 1;

    // let b = 1;
    let bound = 3;
    let sk_bound = Integer::from(10);

    let decomp = Decomp::new(4, &grp);

    // x = dim1
    let x = gen_random_vector(dim, &Integer::from(bound), &mut rng);
    
    // f: dim1 -> dim2
    // let mut f1 = Matrix::new(dim2, dim * dim); 
    // for i in 0..f1.rows {
    //     let f = gen_random_vector(dim, &Integer::from(bound), &mut rng);
    //     let f_tensor = tensor_product_vecs(&f, &f, &grp.delta);
    //     f1.set_row(i, &f_tensor);
    // }
    let f1 = Matrix::random_quadratic_tensored(dim, dim2, &Integer::from(bound), &grp.delta, &mut rng);
    let f2 = Matrix::random_quadratic_tensored(dim2, dim3, &Integer::from(bound), &grp.delta, &mut rng);

    println!("x = {:?}", x);
    println!("f1 = {:?}", f1);
    println!("f2 = {:?}", f2);

    // Perform setup
    println!("start protocol_setup");
    let start = SystemTime::now();
    let ((dcr_sk, dcr_pk), 
        qe_sk_init, 
        qe_sk_vec,
        qe_sk_end,
    ) = protocol_setup(
        vec![dim, dim2, dim3], 
        0, 
        k, 
        &sk_bound, 
        &grp, 
        &mut rng);

    let (h0_left, h0_right) = sample_h(dim, k, &grp.delta, &mut rng);
    let (h1_left, h1_right) = sample_h(dim2, k, &grp.delta, &mut rng);
    println!("h1_left shape = {} {} ", h1_left.rows, h1_left.cols);
    let (gamma_left, gamma_right) = sample_gamma(dim, dim, &grp.delta, &mut rng);
    let end = start.elapsed();
    println!("Time elapsed in protocol_setup is: {:?}", end);
    
    println!("start protocol_keygen_switch");
    let start = SystemTime::now();
    let ((switch_key_x, switch_key_y, switch_key_h), (switch_key_dcr_x, switch_key_dcr_y, switch_key_dcr_h))
    = protocol_keygen_switch(
        &qe_sk_init,
        &dcr_sk, 
        &h0_right, 
        &gamma_left, 
        dim, 
        k, 
        &decomp, 
        &grp, 
        &mut rng);
    let end = start.elapsed();
    println!("Time elapsed in protocol_keygen_switch is: {:?}", end);

    println!("start protocol_enc_init");
    let start = SystemTime::now();
    let ctxt_x = protocol_enc_init(&dcr_pk, &gamma_right, &x, &grp);
    let end = start.elapsed();
    println!("Time elapsed in protocol_enc_init is: {:?}", end);

    println!("start protocol_keyswitch");
    let start = SystemTime::now();
    //  run keyswitch
    let (ct_out_x, ct_out_y, ct_out_h) = protocol_keyswitch(
        &ctxt_x, 
        (&switch_key_x, &switch_key_y, &switch_key_h),
        (&switch_key_dcr_x, &switch_key_dcr_y, &switch_key_dcr_h),
        &grp
    );
    let end = start.elapsed();
    println!("Time elapsed in protocol_keyswitch is: {:?}", end);

    println!("start protocol_keygen_i");
    let start = SystemTime::now();
    let ((qe_b_x, qe_b_y, qe_b_h), 
         (sk_f_mat_x, sk_f_mat_y, sk_f_mat_h), 
         (sk_red_mat_x, sk_red_mat_y, sk_red_mat_h)
        ) = protocol_keygen_i(
            &qe_sk_vec[0], 
            &qe_sk_init, 
            &h1_right, 
            &h0_left, 
            dim2, 
            k, 
            &f1, 
            &decomp, 
            &grp, 
            &mut rng);
    let end = start.elapsed();
    println!("Time elapsed in protocol_keygen_i is: {:?}", end);

    println!("start protocol_dec_i");
    let start = SystemTime::now();
    let (ct_out_x2, ct_out_y2, ct_out_h2) = protocol_dec_i(
        (&ct_out_x, &ct_out_y, &ct_out_h),
        (&qe_b_x, &qe_b_y, &qe_b_h),
        (&sk_f_mat_x, &sk_f_mat_y, &sk_f_mat_h),
        (&sk_red_mat_x, &sk_red_mat_y, &sk_red_mat_h),
        &decomp,
        &grp
    );
    let end = start.elapsed();
    println!("Time elapsed in protocol_dec_i is: {:?}", end);

    println!("start protocol_keygen_end");
    let start = SystemTime::now();
    let (sk_f_mat, sk_f_red) = protocol_keygen_end(&qe_sk_vec[0], &h1_left, &f2, &decomp, &grp);
    let end = start.elapsed();
    println!("Time elapsed in protocol_keygen_end is: {:?}", end);

    println!("start protocol_dec_end");
    let start = SystemTime::now();
    let val_end = protocol_dec_end(
        (&ct_out_x2, &ct_out_y2, &ct_out_h2),
        (&sk_f_mat, &sk_f_red),
        &decomp,
        &grp
    );
    let end = start.elapsed();
    println!("Time elapsed in protocol_dec_end is: {:?}", end);

    println!("reprint inputs");
    println!("x: {:?}", x);
    println!("f1: {:?}", f1);
    println!("f2: {:?}", f2);
    let fx = eval_quadratic_multivariate(&x, &x, &f1);
    println!("f1(x) = {:?}", fx);

    let mut ffx = eval_quadratic_multivariate(&fx, &fx, &f2);
    println!("f2(f1(x)) = {:?}", ffx);
    vec_mod(&mut ffx, &grp.n);
    println!("real mod n = {:?}", ffx);

    println!("eval result: {:?}", val_end);
}

fn run_protocol_with_ml_data() {
    println!("run test_ipe_start_to_end");

    let mut rng = RandState::new(); // Create a single RandState object
    let d = SystemTime::now()
    .duration_since(SystemTime::UNIX_EPOCH)
    .expect("Duration since UNIX_EPOCH failed");
    rng.seed(&Integer::from(d.as_secs()));
    
    let grp = Group::new(100); // Initialize the group        
    // println!("{}", grp);

    let scale = 1048576; // 2^20
    // let scale = 100; // 2^20
    let sk_bound = Integer::from(10);
    let decomp = Decomp::new(4, &grp);

    let dataset = "breast";
    let file_x_test = format!("plain_model/{}/x_test.csv", dataset);
    let file_layer1_weight = format!("plain_model/{}/input_layer.weight.csv", dataset);
    let file_layer1_bias = format!("plain_model/{}/input_layer.bias.csv", dataset);
    let file_layer2_weight = format!("plain_model/{}/output_layer.weight.csv", dataset);
    let file_layer2_bias = format!("plain_model/{}/output_layer.bias.csv", dataset);

    let x_test = read_matrix(&file_x_test, scale);
    let layer1_weight = read_matrix(&file_layer1_weight, scale);
    let layer2_weight = read_matrix(&file_layer2_weight, scale);
    let mut layer1_bias = read_matrix(&file_layer1_bias, scale);
    let mut layer2_bias = read_matrix(&file_layer2_bias, scale);

    let scale_int = Integer::from(scale);
    layer1_bias.mul_scalar_inplace(&scale_int);
    let scale_int_4 = scale_int.pow(4);
    layer2_bias.mul_scalar_inplace(&scale_int_4);

    let mut x = x_test.get_row(0);
    x.push(Integer::from(1));

    let f1_notensor = concatenate_col(&layer1_weight, &layer1_bias);
    let f1_no_one = f1_notensor.gen_tensored_matrix(&grp.delta);
    let mut vec_one = vec![Integer::from(0); f1_no_one.cols];
    vec_one[f1_no_one.cols - 1] = Integer::from(1);
    let f1 = concatenate_vec_row(&f1_no_one, &vec_one);

    let f2_notensor = concatenate_col(&layer2_weight, &layer2_bias);
    let f2 = f2_notensor.gen_tensored_matrix(&grp.delta);

    println!("x = {:?}", x);
    println!("f1_notensor = {}", f1_notensor);
    // println!("f1 = {:?}", f1);
    println!("f2_notensor = {}", f2_notensor);
    // println!("f2 = {:?}", f2);

    assert_eq!(x.len() * x.len(), f1.cols, "dim1 is not matched between x and f1"); // dim
    assert_eq!(f1.rows * f1.rows(), f2.cols, "dim2 is not matched between f1 and f2"); // dim2
    let dim = x.len();
    let dim2 = f1.rows;
    let dim3 = f2.rows;
    let k = 1; // fixed

    // // let b = 1;
    // let bound = 3;

    // // x = dim1
    // let x = gen_random_vector(dim, &Integer::from(bound), &mut rng);

    // let f1 = Matrix::random_quadratic_tensored(dim, dim2, &Integer::from(bound), &grp.delta, &mut rng);
    // let f2 = Matrix::random_quadratic_tensored(dim2, dim3, &Integer::from(bound), &grp.delta, &mut rng);

    // println!("x = {:?}", x);
    // println!("f1 = {:?}", f1);
    // println!("f2 = {:?}", f2);

    // Perform setup
    println!("start protocol_setup");
    let start = SystemTime::now();
    let ((dcr_sk, dcr_pk), 
        qe_sk_init, 
        qe_sk_vec,
        qe_sk_end,
    ) = protocol_setup(
        vec![dim, dim2, dim3], 
        0, 
        k, 
        &sk_bound, 
        &grp, 
        &mut rng);

    let (h0_left, h0_right) = sample_h(dim, k, &grp.delta, &mut rng);
    let (h1_left, h1_right) = sample_h(dim2, k, &grp.delta, &mut rng);
    println!("h1_left shape = {} {} ", h1_left.rows, h1_left.cols);
    let (gamma_left, gamma_right) = sample_gamma(dim, dim, &grp.delta, &mut rng);
    let end = start.elapsed();
    println!("Time elapsed in protocol_setup is: {:?}", end);
    
    println!("start protocol_keygen_switch");
    let start = SystemTime::now();
    let ((switch_key_x, switch_key_y, switch_key_h), (switch_key_dcr_x, switch_key_dcr_y, switch_key_dcr_h))
    = protocol_keygen_switch(
        &qe_sk_init,
        &dcr_sk, 
        &h0_right, 
        &gamma_left, 
        dim, 
        k, 
        &decomp, 
        &grp, 
        &mut rng);
    let end = start.elapsed();
    println!("Time elapsed in protocol_keygen_switch is: {:?}", end);

    println!("start protocol_enc_init");
    let start = SystemTime::now();
    let ctxt_x = protocol_enc_init(&dcr_pk, &gamma_right, &x, &grp);
    let end = start.elapsed();
    println!("Time elapsed in protocol_enc_init is: {:?}", end);

    println!("start protocol_keyswitch");
    let start = SystemTime::now();
    //  run keyswitch
    let (ct_out_x, ct_out_y, ct_out_h) = protocol_keyswitch(
        &ctxt_x, 
        (&switch_key_x, &switch_key_y, &switch_key_h),
        (&switch_key_dcr_x, &switch_key_dcr_y, &switch_key_dcr_h),
        &grp
    );
    let end = start.elapsed();
    println!("Time elapsed in protocol_keyswitch is: {:?}", end);
keygen_i
    println!("start protocol_");
    let start = SystemTime::now();
    let ((qe_b_x, qe_b_y, qe_b_h), 
         (sk_f_mat_x, sk_f_mat_y, sk_f_mat_h), 
         (sk_red_mat_x, sk_red_mat_y, sk_red_mat_h)
        ) = protocol_keygen_i(
            &qe_sk_vec[0], 
            &qe_sk_init, 
            &h1_right, 
            &h0_left, 
            dim2, 
            k, 
            &f1, 
            &decomp, 
            &grp, 
            &mut rng);
    let end = start.elapsed();
    println!("Time elapsed in protocol_keygen_i is: {:?}", end);

    println!("start protocol_dec_i");
    let start = SystemTime::now();
    let (ct_out_x2, ct_out_y2, ct_out_h2) = protocol_dec_i(
        (&ct_out_x, &ct_out_y, &ct_out_h),
        (&qe_b_x, &qe_b_y, &qe_b_h),
        (&sk_f_mat_x, &sk_f_mat_y, &sk_f_mat_h),
        (&sk_red_mat_x, &sk_red_mat_y, &sk_red_mat_h),
        &decomp,
        &grp
    );
    let end = start.elapsed();
    println!("Time elapsed in protocol_dec_i is: {:?}", end);

    println!("start protocol_keygen_end");
    let start = SystemTime::now();
    let (sk_f_mat, sk_f_red) = protocol_keygen_end(&qe_sk_vec[0], &h1_left, &f2, &decomp, &grp);
    let end = start.elapsed();
    println!("Time elapsed in protocol_keygen_end is: {:?}", end);

    println!("start protocol_dec_end");
    let start = SystemTime::now();
    let val_end = protocol_dec_end(
        (&ct_out_x2, &ct_out_y2, &ct_out_h2),
        (&sk_f_mat, &sk_f_red),
        &decomp,
        &grp
    );
    let end = start.elapsed();
    println!("Time elapsed in protocol_dec_end is: {:?}", end);

    println!("reprint inputs");
    println!("x: {:?}", x);
    println!("f1: {:?}", f1);
    println!("f2: {:?}", f2);
    let fx = eval_quadratic_multivariate(&x, &x, &f1);
    println!("f1(x) = {:?}", fx);

    let mut ffx = eval_quadratic_multivariate(&fx, &fx, &f2);
    println!("f2(f1(x)) = {:?}", ffx);
    vec_mod(&mut ffx, &grp.n);
    println!("real mod n = {:?}", ffx);

    println!("eval result: {:?}", val_end);

    for val in val_end {
        let val_scaled = val / Float::with_val(53, scale).pow(Integer::from(10));
    }

}
fn main() {
    // run_protocol_start_to_end();
    run_protocol_with_ml_data();
}
