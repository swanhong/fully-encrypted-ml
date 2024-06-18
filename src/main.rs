extern crate rug;

mod util;
mod parse;
mod dcr_ipe;
mod ipe;
mod qe;
mod protocol;
mod test;

use protocol::algorithm::get_sk_bound;
use rug::Integer;
use rug::ops::Pow;
use rug::rand::RandState;
use util::arguments::Argument;
use std::time::SystemTime;
use crate::protocol::scheme::{protocol_setup, protocol_keygen_switch, protocol_enc_init, protocol_keyswitch, protocol_keygen_i, protocol_dec_i, protocol_keygen_end, protocol_dec_end};
use crate::util::group::Group;
use crate::util::matrix::{Matrix, eval_quadratic_multivariate, concatenate_col, concatenate_vec_row};
use crate::util::vector::{gen_random_vector, vec_mod};
use crate::util::decomp::Decomp;
use crate::protocol::algorithm::{sample_h, sample_gamma};
use crate::parse::file_reader::read_matrix;
use clap::Parser;

fn run_protocol_start_to_end(
    bit_len: u64,
    dim_vec: &Vec<usize>, 
    bound: usize,
    n_decomp: usize,
    rng: &mut RandState
) -> Vec<Integer> {
    println!("run test_ipe_start_to_end");

    let grp = Group::new(bit_len); // Initialize the group        
    // println!("{}", grp);
    
    assert_eq!(dim_vec.len(), 3, "dim_vec must have length 3");
    let dim  = dim_vec[0];
    let dim2 = dim_vec[1];
    let dim3 = dim_vec[2];
    let k = 1;
    let sk_bound= get_sk_bound(dim, bound, 128, &grp);
    
    let decomp = Decomp::new(n_decomp, &grp);

    let x = gen_random_vector(dim, &Integer::from(bound), rng);
    let f1 = Matrix::random_quadratic_tensored(dim, dim2, &Integer::from(bound), &grp.delta, rng);
    let f2 = Matrix::random_quadratic_tensored(dim2, dim3, &Integer::from(bound), &grp.delta, rng);

    println!("x = {:?}", x);
    println!("f1 = {:?}", f1);
    println!("f2 = {:?}", f2);

    let fx = eval_quadratic_multivariate(&x, &x, &f1);
    let mut ffx = eval_quadratic_multivariate(&fx, &fx, &f2);
    for i in 0..ffx.len() {
        if ffx[i].clone().abs() > grp.n.clone() / Integer::from(2) {
            println!(" ===== abort since ffx[{}] = {} > n = {} =====", i, ffx[i], grp.n);
            return vec![Integer::from(-1); dim3]
        }
    }

    // Perform setup
    println!("start protocol_setup");
    let start = SystemTime::now();
    let ((dcr_sk, dcr_pk), 
        qe_sk_init, 
        qe_sk_vec,
        _qe_sk_end,
    ) = protocol_setup(
        vec![dim, dim2, dim3], 
        0, 
        k, 
        &sk_bound, 
        &grp, 
        rng);
    
    let (h0_left, h0_right) = sample_h(dim, k, &grp.delta, rng);
    let (h1_left, h1_right) = sample_h(dim2, k, &grp.delta, rng);

    let (gamma_left, gamma_right) = sample_gamma(dim, dim, &grp.delta, rng);
    let time_setup = start.elapsed();
    println!("Time elapsed in protocol_setup is: {:?}", time_setup);

    println!("start protocol_keygen_switch");
    let start = SystemTime::now();
    let (
        (switch_key_x, switch_key_y, switch_key_h),
        (switch_key_dcr_x, switch_key_dcr_y, switch_key_dcr_h),
    ) = protocol_keygen_switch(
        &qe_sk_init,
        &dcr_sk, 
        &h0_right, 
        &gamma_left, 
        dim, 
        k, 
        &decomp, 
        &grp, 
        rng);
    let time_keygen_switch = start.elapsed();
    println!("Time elapsed in protocol_keygen_switch is: {:?}", time_keygen_switch);

    println!("start protocol_enc_init");
    let start = SystemTime::now();
    let ctxt_x = protocol_enc_init(&dcr_pk, &gamma_right, &x, &grp, rng);
    let time_enc = start.elapsed();
    println!("Time elapsed in protocol_enc_init is: {:?}", time_enc);
    
    println!("start protocol_keyswitch");
    let start = SystemTime::now();
    //  run keyswitch
    let (ct_out_x, ct_out_y, ct_out_h) = protocol_keyswitch(
        &ctxt_x, 
        (&switch_key_x, &switch_key_y, &switch_key_h),
        (&switch_key_dcr_x, &switch_key_dcr_y, &switch_key_dcr_h),
        &decomp,
        &grp
    );
    let time_protocol_keyswitch = start.elapsed();
    println!("Time elapsed in protocol_keyswitch is: {:?}", time_protocol_keyswitch);

    println!("start protocol_keygen_i");
    let start = SystemTime::now();
    let ((qe_b_x, qe_b_y, qe_b_h), 
         (sk_f_mat_x, sk_f_mat_y, sk_f_mat_h), 
         (sk_red_mat_x, sk_red_mat_y, sk_red_mat_h),
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
            rng);
    let time_protocol_keygen_i = start.elapsed();
    println!("Time elapsed in protocol_keygen_i is: {:?}", time_protocol_keygen_i);

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
    let time_protocol_dec_i = start.elapsed();
    println!("Time elapsed in protocol_dec_i is: {:?}", time_protocol_dec_i);

    println!("start protocol_keygen_end");
    let start = SystemTime::now();
    let (sk_f_mat, sk_f_red) = protocol_keygen_end(&qe_sk_vec[0], &h1_left, &f2, &decomp, &grp);
    let time_protocol_keygen_end = start.elapsed();
    println!("Time elapsed in protocol_keygen_end is: {:?}", time_protocol_keygen_end);

    println!("start protocol_dec_end");
    let start = SystemTime::now();
    let val_end = protocol_dec_end(
        (&ct_out_x2, &ct_out_y2, &ct_out_h2),
        (&sk_f_mat, &sk_f_red),
        &decomp,
        &grp
    );
    let time_protocol_dec_end = start.elapsed();
    println!("Time elapsed in protocol_dec_end is: {:?}", time_protocol_dec_end);

    println!("reprint inputs");
    println!("x: {:?}", x);
    println!("f1: {:?}", f1);
    println!("f2: {:?}", f2);
    println!("f1(x) = {:?}", fx);
    println!("f2(f1(x)) = {:?}", ffx);
    vec_mod(&mut ffx, &grp.n);
    println!("real mod n = {:?}", ffx);
    // vec_mod(&mut val_end, &grp.n);
    println!("eval result: {:?}", val_end);
    for i in 0..val_end.len() {
        assert_eq!(val_end[i], ffx[i], "dec_end is not matched with eval result in th element");
    }

    println!(" === Time summaries ===");
    println!("protocol_setup: {:?}", time_setup);
    println!("protocol_keygen_switch: {:?}", time_keygen_switch);
    println!("protocol_enc_init: {:?}", time_enc);
    println!("protocol_keyswitch: {:?}", time_protocol_keyswitch);
    println!("protocol_keygen_i: {:?}", time_protocol_keygen_i);
    println!("protocol_dec_i: {:?}", time_protocol_dec_i);
    println!("protocol_keygen_end: {:?}", time_protocol_keygen_end);
    println!("protocol_dec_end: {:?}", time_protocol_dec_end);
    println!("");
    println!("SETUP: {:?}", time_setup.unwrap());
    println!("KEYGEN: {:?}", time_keygen_switch.unwrap() + time_protocol_keygen_i.unwrap() + time_protocol_keygen_end.unwrap());
    println!("ENC: {:?}", time_enc.unwrap());
    println!("DEC: {:?}", time_protocol_keyswitch.unwrap() + time_protocol_dec_i.unwrap() + time_protocol_dec_end.unwrap());

    val_end
}
fn run_protocol_with_ml_data(
    dataset: &str,
    bit_len: u64,
    scale: i64,
    n_decomp: usize,
) {
    println!("run test_ipe_start_to_end");

    let mut rng = RandState::new(); // Create a single RandState object
    let d = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .expect("Duration since UNIX_EPOCH failed");
    rng.seed(&Integer::from(d.as_secs()));

    let grp = Group::new(bit_len); // Initialize the group        
    let decomp = Decomp::new(n_decomp, &grp);

    println!("=== run {} dataset with bit_len {} ===", dataset, bit_len);
    let file_x_test: String = format!("plain_model/{}/X_test.csv", dataset);
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
    let scale_int_pow_1 = scale_int.clone().pow(1);
    let scale_int_pow_4 = scale_int.clone().pow(4);
    let scale_int_pow_10 = scale_int.clone().pow(10);
    layer1_bias.mul_scalar_inplace(&scale_int_pow_1); // scale * scale^1
    layer2_bias.mul_scalar_inplace(&scale_int_pow_4); // scale * scale^4

    // we only test the first row of X_test data
    let mut x = x_test.get_row(0);
    x.push(Integer::from(1));

    let f1_notensor = concatenate_col(&layer1_weight, &layer1_bias);
    let f1_no_one = f1_notensor.gen_tensored_matrix(&grp.delta);
    let mut vec_one = vec![Integer::from(0); f1_no_one.cols];
    vec_one[f1_no_one.cols - 1] = Integer::from(1);
    let f1 = concatenate_vec_row(&f1_no_one, &vec_one);

    let f2_notensor = concatenate_col(&layer2_weight, &layer2_bias);
    let f2 = f2_notensor.gen_tensored_matrix(&grp.delta);

    println!("scale = {}", scale);
    println!("x = {:?}", x);
    println!("f1_notensor = {}", f1_notensor);
    println!("f2_notensor = {}", f2_notensor);

    assert_eq!(x.len() * x.len(), f1.cols, "dim1 is not matched between x and f1"); // dim
    assert_eq!(f1.rows * f1.rows(), f2.cols, "dim2 is not matched between f1 and f2"); // dim2
    let dim = x.len();
    let dim2 = f1.rows;
    let dim3 = f2.rows;
    let k = 1; // fixed
    let sk_bound = get_sk_bound(dim, 10 * scale as usize, 128, &grp);

    // Perform setup
    println!("start protocol_setup");
    let start = SystemTime::now();
    let ((dcr_sk, dcr_pk),
        qe_sk_init,
        qe_sk_vec,
        _qe_sk_end,
    ) = protocol_setup(
        vec![dim, dim2, dim3],
        0,
        k,
        &sk_bound,
        &grp,
        &mut rng);
    println!("protocol setup done?");
    
    let (h0_left, h0_right) = sample_h(dim, k, &grp.delta, &mut rng);
    let (h1_left, h1_right) = sample_h(dim2, k, &grp.delta, &mut rng);
    let (gamma_left, gamma_right) = sample_gamma(dim, dim, &grp.delta, &mut rng);
    let time_setup = start.elapsed();
    println!("Time elapsed in protocol_setup is: {:?}", time_setup);

    println!("start protocol_keygen_switch");
    let start = SystemTime::now();
    let (
        (switch_key_x, switch_key_y, switch_key_h),
        (switch_key_dcr_x, switch_key_dcr_y, switch_key_dcr_h),
    ) = protocol_keygen_switch(
        &qe_sk_init,
        &dcr_sk,
        &h0_right,
        &gamma_left,
        dim,
        k,
        &decomp,
        &grp,
        &mut rng);
    let time_keygen_switch = start.elapsed();
    println!("Time elapsed in protocol_keygen_switch is: {:?}", time_keygen_switch);

    println!("start protocol_enc_init");
    let start = SystemTime::now();
    let ctxt_x = protocol_enc_init(&dcr_pk, &gamma_right, &x, &grp, &mut rng);
    let time_enc = start.elapsed();
    println!("Time elapsed in protocol_enc_init is: {:?}", time_enc);

    println!("start protocol_keyswitch");
    let start = SystemTime::now();
    //  run keyswitch
    let (ct_out_x, ct_out_y, ct_out_h) = protocol_keyswitch(
        &ctxt_x,
        (&switch_key_x, &switch_key_y, &switch_key_h),
        (&switch_key_dcr_x, &switch_key_dcr_y, &switch_key_dcr_h),
        &decomp,
        &grp
    );
    let time_protocol_keyswitch = start.elapsed();
    println!("Time elapsed in protocol_keyswitch is: {:?}", time_protocol_keyswitch);

    println!("start protocol_keygen_i");
    let start = SystemTime::now();
    let ((qe_b_x, qe_b_y, qe_b_h),
         (sk_f_mat_x, sk_f_mat_y, sk_f_mat_h),
         (sk_red_mat_x, sk_red_mat_y, sk_red_mat_h),
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
    let time_protocol_keygen_i = start.elapsed();
    println!("Time elapsed in protocol_keygen_i is: {:?}", time_protocol_keygen_i);

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
    let time_protocol_dec_i = start.elapsed();
    println!("Time elapsed in protocol_dec_i is: {:?}", time_protocol_dec_i);

    println!("start protocol_keygen_end");
    let start = SystemTime::now();
    let (sk_f_mat, sk_f_red) = protocol_keygen_end(&qe_sk_vec[0], &h1_left, &f2, &decomp, &grp);
    let time_protocol_keygen_end = start.elapsed();
    println!("Time elapsed in protocol_keygen_end is: {:?}", time_protocol_keygen_end);

    println!("start protocol_dec_end");
    let start = SystemTime::now();
    let val_end = protocol_dec_end(
        (&ct_out_x2, &ct_out_y2, &ct_out_h2),
        (&sk_f_mat, &sk_f_red),
        &decomp,
        &grp
    );
    let time_protocol_dec_end = start.elapsed();
    println!("Time elapsed in protocol_dec_end is: {:?}", time_protocol_dec_end);
    
    let val_scaled = val_end.iter().map(|x| x.clone() / scale_int_pow_10.clone()).collect::<Vec<Integer>>();

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
    println!("scaled result: {:?}", val_scaled);

    println!(" === Time summaries ===");
    println!("protocol_setup: {:?}", time_setup);
    println!("protocol_keygen_switch: {:?}", time_keygen_switch);
    println!("protocol_enc_init: {:?}", time_enc);
    println!("protocol_keyswitch: {:?}", time_protocol_keyswitch);
    println!("protocol_keygen_i: {:?}", time_protocol_keygen_i);
    println!("protocol_dec_i: {:?}", time_protocol_dec_i);
    println!("protocol_keygen_end: {:?}", time_protocol_keygen_end);
    println!("protocol_dec_end: {:?}", time_protocol_dec_end);
    println!("");
    println!("SETUP: {:?}", time_setup.unwrap());
    println!("KEYGEN: {:?}", time_keygen_switch.unwrap() + time_protocol_keygen_i.unwrap() + time_protocol_keygen_end.unwrap());
    println!("ENC: {:?}", time_enc.unwrap());
    println!("DEC: {:?}", time_protocol_keyswitch.unwrap() + time_protocol_dec_i.unwrap() + time_protocol_dec_end.unwrap());
}

fn main() {
    let args = Argument::parse();
    println!("{}", args);

    let mut rng = RandState::new(); // Create a single RandState object
    let d = SystemTime::now()
    .duration_since(SystemTime::UNIX_EPOCH)
    .expect("Duration since UNIX_EPOCH failed");
    rng.seed(&Integer::from(d.as_secs()));
    let bit_len = args.bit_len;
    let scale = args.scale;
    let n_decomp = args.n_decomp;

    match args.target_algorithm.as_str() {
        "protocol" => {
            let dim0 = args.dim0;
            let dim1 = args.dim1;
            let dim2 = args.dim2;
            let n_try = args.n_try;
            let bound = args.bound;
            let dim_vec = vec![dim0, dim1, dim2];
            let mut outputs = Matrix::new(n_try, dim_vec[dim_vec.len()-1]);
            for i in 0..n_try {
                println!("  ====================== i = {} ======================", i);
                let row = run_protocol_start_to_end(
                    bit_len,
                    &dim_vec, 
                    bound,
                    n_decomp,
                    &mut rng);
                outputs.set_row(i, &row);
            }
            if n_try >= 2 {
                println!("outputs = \n{}", outputs);    
            }
        },
        "ml" => {
            run_protocol_with_ml_data(
                args.data_ml.as_ref().unwrap().as_str(),
                bit_len,
                scale,
                n_decomp,
            );
        },
        _ => {
            println!("invalid target");
        }
    }
}
