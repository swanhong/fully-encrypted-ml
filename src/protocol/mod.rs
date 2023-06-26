pub mod scheme;
pub mod algorithm;

#[cfg(test)]
mod test {
    use rug::Integer;
    use rug::rand::RandState;
    use std::time::SystemTime;
    use crate::protocol::scheme::{protocol_setup, protocol_keygen_switch, protocol_enc_init, protocol_keyswitch, protocol_keygen_i, protocol_dec_i, protocol_keygen_end, protocol_dec_end};
    use crate::util::group::Group;

    use crate::util::matrix::{Matrix, eval_quadratic_multivariate};
    use crate::util::vector::{gen_random_vector, vec_mod, tensor_product_vecs};
    use crate::util::decomp::Decomp;
    use super::algorithm::{sample_h, sample_gamma};

    #[test]
    fn test_protocol_start_to_end() {
        println!("run test_protocol_start_to_end");

        let mut rng = RandState::new(); // Create a single RandState object
        let d = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .expect("Duration since UNIX_EPOCH failed");
        rng.seed(&Integer::from(d.as_secs()));
        
        let grp = Group::new(100); // Initialize the group        
        // println!("{}", grp);
        
        let dim = 2;
        let k = 1;
        let f_num = dim;
        // let b = 1;
        let bound = 10;
        let sk_bound = Integer::from(10);

        let decomp = Decomp::new(4, &grp);
        // println!("decomp = {}", decomp);
        let mut x = gen_random_vector(dim, &Integer::from(bound), &mut rng);
        let mut f_mat = Matrix::new(f_num, dim * dim);
        for i in 0..f_num {
            let f = gen_random_vector(dim, &Integer::from(bound), &mut rng);
            let f_tensor = tensor_product_vecs(&f, &f, &grp.delta);
            f_mat.set_row(i, &f_tensor);
        }

        // Perform setup
        println!("start protocol_setup");
        let start = SystemTime::now();
        let ((dcr_sk, dcr_pk), 
            qe_sk, 
            qe_sk_vec,
            qe_sk_end,
        ) = protocol_setup(
            vec![dim], 
            f_num, 
            k, 
            &sk_bound, 
            &grp, 
            &mut rng);

        let (h_left, h_right) = sample_h(dim, k, &grp.delta, &mut rng);
        let (gamma_left, gamma_right) = sample_gamma(dim, dim, &grp.delta, &mut rng);
        let end = start.elapsed();
        println!("Time elapsed in protocol_setup is: {:?}", end);
        
        println!("start protocol_keygen_switch");
        let start = SystemTime::now();
        let ((switch_key_x, switch_key_y, switch_key_h), (switch_key_dcr_x, switch_key_dcr_y, switch_key_dcr_h))
        = protocol_keygen_switch(
            &qe_sk,
            &dcr_sk, 
            &h_right, 
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
        let ((qe_b_x, qe_b_y, qe_b_h), (sk_f_mat_x, sk_f_mat_y, sk_f_mat_h), (sk_red_mat_x, sk_red_mat_y, sk_red_mat_h))
            = protocol_keygen_i(&qe_sk, &qe_sk, &h_right, &h_left, dim, k, &f_mat, &decomp, &grp, &mut rng);
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
        let (sk_f_mat, sk_f_red) = protocol_keygen_end(&qe_sk, &h_left, &f_mat, &decomp, &grp);
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
        println!("f_mat: {:?}", f_mat);
        let fx = eval_quadratic_multivariate(&x, &x, &f_mat);
        println!("f(x) = {:?}", fx);

        let mut ffx = eval_quadratic_multivariate(&fx, &fx, &f_mat);
        println!("f(f(x)) = {:?}", ffx);
        println!("mod n = {:?}", grp.n);
        vec_mod(&mut ffx, &grp.n);
        println!("real mod n = {:?}", ffx);

        println!("eval result: {:?}", val_end);
    }

    #[test]
    fn test_matrix_h_and_gamma() {
        let dim = 4;
        let k = 1;
        let grp = Group::new(10);
        let n = grp.delta;

        let mut rng = rug::rand::RandState::new();

        let (h_left_1, h_right_1) = sample_h(dim, k, &n, &mut rng);

        println!("h_left_1:");
        println!("{}", h_left_1);

        println!("h_right_1:");
        println!("{}", h_right_1);

        println!("test mult");
        let mut mult = h_left_1.clone() * h_right_1.clone();
        mult.mod_inplace(&n);
        println!("{}", mult);

        let (gamma_left_1, gamma_right_1) = sample_gamma(dim, k, &n, &mut rng);

        println!("gamma_left_1:");
        println!("{}", gamma_left_1);

        println!("gamma_right_1:");
        println!("{}", gamma_right_1);

        println!("test mult");
        let mut mult = gamma_left_1.clone() * gamma_right_1.clone();
        mult.mod_inplace(&n);
        println!("{}", mult);
    }

}