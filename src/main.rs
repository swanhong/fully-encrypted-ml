extern crate rug;

mod util;
mod parse;
mod dcr;
mod ipfe;
mod qfe;
mod protocol;
mod test;

use protocol::algorithm::get_sk_bound;
use rug::Integer;
use rug::ops::Pow;
use rug::rand::RandState;

use util::arguments::Argument;
use std::time::SystemTime;
use crate::protocol::scheme::*;
use crate::util::group::Group;
use crate::util::matrix::{Matrix, eval_quadratic_multivariate, concatenate_row,  concatenate_col, concatenate_vec_row};
use crate::util::vector::{gen_random_vector, vec_mod};
use crate::util::decomp::Decomp;
use crate::protocol::algorithm::{sample_h, sample_gamma};
use crate::parse::file_reader::read_matrix;
use clap::Parser;

fn run_protocol_start_to_end(
    bit_len: u64,
    dim: &Vec<usize>, 
    bound: usize,
    n_decomp: usize,
    rng: &mut RandState
) -> Vec<Integer> {
    println!("run test_ipfe_start_to_end");

    let grp = Group::new(bit_len); // Initialize the group    
    let depth = dim.len() - 1;  
    let k = 1;
    let bound_h = Integer::from(1026); // bound_h * 2 + 1 is 11-bit prime
    let sk_bound= get_sk_bound(dim[0], bound, 128, &grp);
    let bound = Integer::from(bound);  
    
    let decomp = Decomp::new(n_decomp, &grp);

    let x = gen_random_vector(dim[0], &bound, rng);

    let f_num = dim.len() - 1;
    let mut f_origin = vec![Matrix::new(0, 0); f_num];
    let mut f = vec![Matrix::new(0, 0); f_num];
    for i in 0..f_num {
        f_origin[i] = Matrix::random_quadratic_tensored_with_one_padded(
            dim[i], 
            dim[i+1], 
            &bound,
            &grp.delta, 
            rng
        );
        if i < f_num - 1 {
            let mut unit_vector = Matrix::new(1, f_origin[i].cols);
            unit_vector.set(0, f_origin[i].cols - 1, Integer::from(1));
            f[i] = concatenate_row(&f_origin[i], &unit_vector);
        } else {
            f[i] = f_origin[i].clone();
        }
    }

    // Perform setup
    println!("start protocol_setup");
    let start = SystemTime::now();
    let ((dcr_sk, dcr_pk), 
        qfe_sk,
    ) = protocol_setup(
        &dim,
        k, 
        &sk_bound, 
        &grp, 
        rng);
    
    let mut h_left = vec![Matrix::new(0, 0); depth];
    let mut h_right = vec![Matrix::new(0, 0); depth];
    for i in 0..depth {
        (h_left[i], h_right[i]) = sample_h(dim[i], k, &bound_h, &grp.n, rng);
    }
    let (gamma_left, gamma_right) = sample_gamma(dim[0], &grp.n, rng);
    let time_setup = start.elapsed();
    println!("Time elapsed in protocol_setup is: {:?}", time_setup);

    println!("start protocol_keygen_dcr_to_qfe");
    let start = SystemTime::now();
    let (mat_ctxts, fk)
    = protocol_keygen_dcr_to_qfe(
        &dcr_sk, 
        &qfe_sk[0],
        &h_right[0], 
        &gamma_left, 
        &decomp, 
        &grp, 
        rng);
    let time_keygen_dcr_to_qfe = start.elapsed();
    println!("Time elapsed in protocol_keygen_dcr_to_qfe is: {:?}", time_keygen_dcr_to_qfe);


    println!("start protocol_keygen_qfe_to_qfe");
    let mut fk_qfe_to_qfe = vec![Matrix::new(0, 0); depth - 1];
    let start = SystemTime::now();
    for i in 0..fk_qfe_to_qfe.len() {
        fk_qfe_to_qfe[i] = protocol_keygen_qfe_to_qfe(
            &qfe_sk[i],
            &qfe_sk[i + 1],
            &h_right[i + 1],
            &h_left[i],
            &f[i],
            &decomp,
            &grp,
            rng
        );
    }
    let time_keygen_qfe_to_qfe = start.elapsed();
    println!("Time elapsed in protocol_keygen_qfe_to_qfe is: {:?}", time_keygen_qfe_to_qfe);

    
    println!("start protocol_keygen_qfe_to_plain");
    let start = SystemTime::now();
    // let mut fk_qfe_to_plain_test = vec![Matrix::new(0, 0); depth - 1];
    // for i in 0..(depth-1) {
    //     fk_qfe_to_plain_test[i] = protocol_keygen_qfe_to_plain(
    //         &qfe_sk[i], // changed
    //         &h_left[i],
    //         &f_origin[i],
    //         &grp
    //     );
    // }
    let fk_qfe_to_plain = protocol_keygen_qfe_to_plain(
        &qfe_sk[depth - 1],
        &h_left[depth - 1],
        &f_origin[depth - 1],
        &grp
    );
    let time_keygen_qfe_to_plain = start.elapsed();
    println!("Time elapsed in protocol_keygen_qfe_to_plain is: {:?}", time_keygen_qfe_to_plain);

    println!("start protocol_enc_init");
    let start = SystemTime::now();
    let ct0 = protocol_enc_init(&dcr_pk, &gamma_right, &x, &grp, rng);
    let time_enc = start.elapsed();
    println!("Time elapsed in protocol_enc_init is: {:?}", time_enc);
    
    let mut ct_vec = vec![vec![Integer::from(0);1]; depth];
    println!("start protocol_dec_dcr_to_qfe");
    let start = SystemTime::now();
    ct_vec[0] = protocol_dec_dcr_to_qfe(
        &ct0, 
        &mat_ctxts,
        &fk,
        &decomp,
        &grp
    );
    let time_dec_dcr_to_qfe = start.elapsed();
    println!("Time elapsed in protocol_dec_dcr_to_qfe is: {:?}", time_dec_dcr_to_qfe);
    
    println!("start protocol_dec_qfe_to_qfe");
    let start = SystemTime::now();
    for i in 0..depth - 1 {
        println!("ct_vec{} <- ct_vec{}", i + 1, i);
        ct_vec[i + 1] = protocol_dec_qfe(
            &ct_vec[i],
            &fk_qfe_to_qfe[i],
            dim[i] + k + 1,
            2 * (dim[i] + k + 1) + 1,
            &decomp,
            &grp,
            true,
        );
    }
    let time_dec_qfe_to_qfe = start.elapsed();
    println!("Time elapsed in protocol_dec_qfe_to_qfe is: {:?}", time_dec_qfe_to_qfe);
    
    println!("start protocol_dec_qfe_to_plain");
    let start = SystemTime::now();
    // let mut val_end = vec![vec![Integer::from(0);1]; depth];
    // for i in 0..depth - 1 {
    //     println!("i = {}", i);
    //     val_end[i] = protocol_dec_qfe(
    //         &ct_vec[i],
    //         &fk_qfe_to_plain_test[i],
    //         dim[i] + k + 1,
    //         2 * (dim[i] + k + 1) + 1,
    //         &decomp,
    //         &grp,
    //         false,
    //     );
    // }
    let val_end_final = protocol_dec_qfe(
        &ct_vec[depth - 1],
        &fk_qfe_to_plain,
        dim[depth-1] + k + 1,
        2 * (dim[depth-1] + k + 1) + 1,
        &decomp,
        &grp,
        false,
    );
    let time_dec_qfe_to_plain = start.elapsed();
    println!("Time elapsed in protocol_dec_qfe_to_plain is: {:?}", time_dec_qfe_to_plain);
    // val_end[depth - 1] = val_end_final.clone();

    println!("reprint inputs");
    println!("x: {:?}", x);
    for i in 0..f_num {
        println!("f[{}] = \n{}", i, f[i]);
    }
    // print f^i(x)
    let mut fx = x.clone();
    fx.push(Integer::from(1));
    println!("x = {:?}", fx);
    for i in 0..f_num {
        fx = eval_quadratic_multivariate(&fx, &fx, &f[i]);
        println!("f^{}(x) = {:?}", i, fx);
        // println!("eval = {:?}", val_end[i]);
    }
    println!("eval result: {:?}", val_end_final);

    println!(" === Time summaries ===");
    println!("protocol_setup: {:?}", time_setup);
    println!("protocol_keygen_dcr_to_qfe: {:?}", time_keygen_dcr_to_qfe);
    println!("protocol_enc_init: {:?}", time_enc);
    println!("protocol_keygen_qfe_to_qfe: {:?}", time_keygen_qfe_to_qfe);
    println!("protocol_keygen_qfe_to_plain: {:?}", time_keygen_qfe_to_plain);
    println!("protocol_dec_dcr_to_qfe: {:?}", time_dec_dcr_to_qfe);
    println!("protocol_dec_qfe_to_qfe: {:?}", time_dec_qfe_to_qfe);
    println!("protocol_dec_qfe_to_plain: {:?}", time_dec_qfe_to_plain);
    println!("");
    println!("SETUP: {:?}", time_setup.unwrap());
    println!("KEYGEN: {:?}", time_keygen_dcr_to_qfe.unwrap() + time_keygen_qfe_to_qfe.unwrap() + time_keygen_qfe_to_plain.unwrap());
    println!("ENC: {:?}", time_enc.unwrap());
    println!("DEC: {:?}", time_dec_dcr_to_qfe.unwrap() + time_dec_qfe_to_qfe.unwrap() + time_dec_qfe_to_plain.unwrap());

    val_end_final
}

fn run_protocol_with_ml_data(
    dataset: &str,
    bit_len: u64,
    scale: i64,
    n_decomp: usize,
) {
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
    // let scale_int_pow_10 = scale_int.clone().pow(10);
    layer1_bias.mul_scalar_inplace(&scale_int_pow_1); // scale * scale^1
    layer2_bias.mul_scalar_inplace(&scale_int_pow_4); // scale * scale^4

    // we only test the first row of X_test data
    let x = x_test.get_row(0);

    let f1_notensor = concatenate_col(&layer1_weight, &layer1_bias);
    let f1_no_one = f1_notensor.gen_tensored_matrix(&grp.delta);
    let mut vec_one = vec![Integer::from(0); f1_no_one.cols];
    vec_one[f1_no_one.cols - 1] = Integer::from(1);
    let f1 = concatenate_vec_row(&f1_no_one, &vec_one);

    let f2_notensor = concatenate_col(&layer2_weight, &layer2_bias);
    let f2 = f2_notensor.gen_tensored_matrix(&grp.delta);

    let k = 1; // fixed
    let bound_h = Integer::from(1026);
    let dim = vec![x.len(), f1.rows - 1, f2.rows];
    let depth = 2;
    let sk_bound = get_sk_bound(dim[0], 10 * scale as usize, 128, &grp);
    
    // Perform setup
    println!("start protocol_setup");
    let start = SystemTime::now();
    let ((dcr_sk, dcr_pk),
        qfe_sk,
    ) = protocol_setup(
        &dim,
        k,
        &sk_bound,
        &grp,
        &mut rng);
    println!("protocol setup done?");

    let mut h_left = vec![Matrix::new(0, 0); depth];
    let mut h_right = vec![Matrix::new(0, 0); depth];
    for i in 0..depth {
        (h_left[i], h_right[i]) = sample_h(dim[i], k, &bound_h, &grp.n, &mut rng);
    }

    let (gamma_left, gamma_right) = sample_gamma(dim[0], &grp.delta, &mut rng);
    let time_setup = start.elapsed();
    println!("Time elapsed in protocol_setup is: {:?}", time_setup);

    println!("start protocol_keygen_dcr_to_qfe");
    let (mat_ctxts, fk)
        = protocol_keygen_dcr_to_qfe(
            &dcr_sk, 
            &qfe_sk[0],
            &h_right[0], 
            &gamma_left, 
            &decomp, 
            &grp, 
            &mut rng
        );
    let time_keygen_dcr_to_qfe = start.elapsed();
    println!("Time elapsed in protocol_keygen_dcr_to_qfe is: {:?}", time_keygen_dcr_to_qfe);

    println!("start protocol_enc_init");
    let start = SystemTime::now();
    let ct0 = protocol_enc_init(&dcr_pk, &gamma_right, &x, &grp, &mut rng);
    let time_enc = start.elapsed();
    println!("Time elapsed in protocol_enc_init is: {:?}", time_enc);

    println!("start protocol_dec_dcr_to_qfe");
    let start = SystemTime::now();
    //  run keyswitch
    let ct0 = protocol_dec_dcr_to_qfe(
        &ct0, 
        &mat_ctxts,
        &fk,
        &decomp,
        &grp
    );
    let time_dec_dcr_to_qfe = start.elapsed();
    println!("Time elapsed in protocol_keyswitch is: {:?}", time_dec_dcr_to_qfe);

    println!("start protocol_keygen_qfe_to_qfe");
    let start = SystemTime::now();
    let fk_qfe_to_qfe = protocol_keygen_qfe_to_qfe(
        &qfe_sk[0],
        &qfe_sk[1],
        &h_right[1],
        &h_left[0],
        &f1,
        &decomp,
        &grp,
        &mut rng
    );
    let time_protocol_keygen_qfe_to_qfe = start.elapsed();
    println!("Time elapsed in protocol_keygen_qfe_to_qfe is: {:?}", time_protocol_keygen_qfe_to_qfe);

    println!("start protocol_dec_qfe_to_qfe");
    let start = SystemTime::now();
    let ct0 = protocol_dec_qfe(
        &ct0,
        &fk_qfe_to_qfe,
        dim[0] + k + 1,
        2 * (dim[0] + k + 1) + 1,
        &decomp,
        &grp,
        true,
    );
    let time_protocol_dec_qfe_to_qfe = start.elapsed();
    println!("Time elapsed in protocol_dec_qfe_to_qfe is: {:?}", time_protocol_dec_qfe_to_qfe);

    println!("start protocol_keygen_qfe_to_plain");
    let start = SystemTime::now();
    let fk_qfe_to_plain = protocol_keygen_qfe_to_plain(
        &qfe_sk[1], // changed
        &h_left[1],
        &f2,
        &grp
    );
    
    let time_protocol_keygen_qfe_to_plain = start.elapsed();
    println!("Time elapsed in protocol_keygen_qfe_to_plain is: {:?}", time_protocol_keygen_qfe_to_plain);

    println!("start protocol_dec_qfe_to_plain");
    let start = SystemTime::now();
    let val_end = protocol_dec_qfe(
        &ct0,
        &fk_qfe_to_plain,
        dim[1] + k + 1,
        2 * (dim[1] + k + 1) + 1,
        &decomp,
        &grp,
        false,
    );
    let time_protocol_dec_qfe_to_plain = start.elapsed();
    println!("Time elapsed in protocol_dec_qfe_to_plain is: {:?}", time_protocol_dec_qfe_to_plain);

    println!("reprint inputs");
    let mut x1 = x.clone();
    x1.push(Integer::from(1));
    println!("x: {:?}", x);
    println!("f1: {}", f1);
    println!("f2: {}", f2);
    let fx = eval_quadratic_multivariate(&x1, &x1, &f1);
    println!("f1(x) = {:?}", fx);
    let mut ffx = eval_quadratic_multivariate(&fx, &fx, &f2);
    println!("f2(f1(x)) = {:?}", ffx);
    vec_mod(&mut ffx, &grp.n);
    println!("real mod n = {:?}", ffx);
    println!("eval result: {:?}", val_end);

    println!(" === Time summaries ===");
    println!("protocol_setup: {:?}", time_setup);
    println!("protocol_keygen_dcr_to_qfe: {:?}", time_keygen_dcr_to_qfe);
    println!("protocol_enc_init: {:?}", time_enc);
    println!("protocol_keygen_qfe_to_qfe: {:?}", time_protocol_keygen_qfe_to_qfe);
    println!("protocol_keygen_qfe_to_plain: {:?}", time_protocol_keygen_qfe_to_plain);
    println!("protocol_dec_dcr_to_qfe: {:?}", time_dec_dcr_to_qfe);
    println!("protocol_dec_qfe_to_qfe: {:?}", time_protocol_dec_qfe_to_qfe);
    println!("protocol_dec_qfe_to_plain: {:?}", time_protocol_dec_qfe_to_plain);
    println!("");
    println!("SETUP: {:?}", time_setup.unwrap());
    println!("KEYGEN: {:?}", time_keygen_dcr_to_qfe.unwrap() + time_protocol_keygen_qfe_to_qfe.unwrap() + time_protocol_keygen_qfe_to_plain.unwrap());
    println!("ENC: {:?}", time_enc.unwrap());
    println!("DEC: {:?}", time_dec_dcr_to_qfe.unwrap() + time_protocol_dec_qfe_to_qfe.unwrap() + time_protocol_dec_qfe_to_plain.unwrap());
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
            let n_try = args.n_try;
            let bound = args.bound;
            let dim_vec = if args.dim.len() == 0 {
                vec![5, 4, 3]
            } else {
                args.dim
            };
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
