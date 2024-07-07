pub mod scheme;
pub mod algorithm;

#[cfg(test)]
pub mod test {
    use rug::Integer;
    use rug::rand::RandState;
    use std::time::SystemTime;

    use crate::protocol::algorithm::{row_reduce_form, row_reduce_form_integer, sample_gamma_wo_padding, sample_h_trivial, sample_h_wo_padding};
    use crate::protocol::scheme::*;
    use crate::qfe::scheme::*;
    use crate::util::group::Group;

    use crate::util::matrix::{Matrix, eval_quadratic_multivariate};
    use crate::util::vector::{gen_random_vector, vec_mod, tensor_product_vecs};
    use crate::util::decomp::Decomp;
    use super::algorithm::{sample_h, sample_gamma, sample_gamma_trivial};

    #[test]
    fn test_protocol_start_to_end() {
        println!("run test_protocol_start_to_end");

        let mut rng = RandState::new(); // Create a single RandState object
        let d = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .expect("Duration since UNIX_EPOCH failed");
        rng.seed(&Integer::from(d.as_secs()));
        
        let bit_len = 100;
        let grp = Group::new(bit_len); // Initialize the group   
        println!("!!! WARNING !!! bit len for primes = {}: this is not secure !!!", bit_len);     
        
        let dim = 5;
        let k = 1;
        let q = 2 * dim + 1;
        let f_num = dim;
        let bound = 10;
        let sk_bound = Integer::from(10);

        let decomp = Decomp::new(3, &grp);
        // println!("decomp = {}", decomp);
        let x = gen_random_vector(dim, &Integer::from(bound), &mut rng);
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
            qfe_sk, 
            _qfe_sk_vec,
            _qfe_sk_end,
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
            &qfe_sk,
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
        let ctxt_x = protocol_enc_init(&dcr_pk, &gamma_right, &x, &grp, &mut rng);
        let end = start.elapsed();
        println!("Time elapsed in protocol_enc_init is: {:?}", end);

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
        let end = start.elapsed();
        println!("Time elapsed in protocol_keyswitch is: {:?}", end);

        // let ct_out_x = decomp.vector_inv(&ct_out_x);
        // let ct_out_y = decomp.vector_inv(&ct_out_y);
        // let ct_out_h = decomp.vector_inv(&ct_out_h);

        println!("start protocol_keygen_i");
        let start = SystemTime::now();
        let (
            (qfe_b_x, qfe_b_y, qfe_b_h),
            (sk_f_mat_x, sk_f_mat_y, sk_f_mat_h),
            (sk_red_mat_x, sk_red_mat_y, sk_red_mat_h)
        ) = protocol_keygen_i(
            &qfe_sk,
            &qfe_sk,
            &h_right,
            &h_left,
            dim,
            k,
            &f_mat,
            &decomp,
            &grp,
            &mut rng
        );
        let end = start.elapsed();
        println!("Time elapsed in protocol_keygen_i is: {:?}", end);

        println!("start protocol_dec_i");
        let start = SystemTime::now();
        let (ct_out_x2, ct_out_y2, ct_out_h2) = protocol_dec_i(
            (&ct_out_x, &ct_out_y, &ct_out_h),
            (&qfe_b_x, &qfe_b_y, &qfe_b_h),
            (&sk_f_mat_x, &sk_f_mat_y, &sk_f_mat_h),
            (&sk_red_mat_x, &sk_red_mat_y, &sk_red_mat_h),
            dim,
            q,
            &decomp,
            &grp
        );
        let end = start.elapsed();
        println!("Time elapsed in protocol_dec_i is: {:?}", end);

        println!("start protocol_keygen_end");
        let start = SystemTime::now();
        let (sk_f_mat, sk_f_red) = protocol_keygen_end(&qfe_sk, &h_left, &f_mat, &grp);
        let end = start.elapsed();
        println!("Time elapsed in protocol_keygen_end is: {:?}", end);

        println!("start protocol_dec_end");
        let start = SystemTime::now();
        let val_end = protocol_dec_end(
            (&ct_out_x2, &ct_out_y2, &ct_out_h2),
            (&sk_f_mat, &sk_f_red),
            &decomp,
            dim,
            q,
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
        vec_mod(&mut ffx, &grp.n);
        println!("real mod n = {:?}", ffx);
        println!("eval result: {:?}", val_end);
    }

    #[test]
    fn test_protocol_start_to_end_new() {
        println!("run test_protocol_start_to_end");

        for try_num in 0..20 {
            println!(" ==================== ");
            println!(" ==================== ");
            println!(" try number = {}", try_num);
            println!(" ==================== ");
            println!(" ==================== ");

        let mut rng = RandState::new(); // Create a single RandState object
        let d = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .expect("Duration since UNIX_EPOCH failed");
        rng.seed(&Integer::from(d.as_secs()));
        
        let bit_len = 100;
        let grp = Group::new(bit_len); // Initialize the group   
        println!("!!! WARNING !!! bit len for primes = {}: this is not secure !!!", bit_len);     
        
        let dim = 3;
        let k = 1;
        let q = 2 * dim + 1;
        let f_num = 2;
        let bound = 5;
        let sk_bound = Integer::from(10);
        let decomp = Decomp::new(4, &grp);
        
        let x = gen_random_vector(dim, &Integer::from(bound), &mut rng);
        // let x = vec![Integer::from(0); dim];
        let mut f_mat = Matrix::new(dim + 1, (dim + 1) * (dim + 1));
        let mut f_mat_no_padding = Matrix::new(dim, dim * dim);
        for i in 0..dim { 
            let mut f = gen_random_vector(dim, &Integer::from(bound), &mut rng);
            let f_tensor = tensor_product_vecs(&f, &f, &grp.delta);
            f_mat_no_padding.set_row(i, &f_tensor);
            f.push(Integer::from(1));
            // let f_tensor = tensor_product_vecs(&f, &f, &grp.delta);
            let f_tensor = gen_random_vector((dim + 1) * (dim + 1), &Integer::from(bound), &mut rng);
            let f_tensor = vec![Integer::from(0); (dim + 1) * (dim + 1)];
            f_mat.set_row(i, &f_tensor);
        }
        let mut unit_vector = vec![Integer::from(0); f_mat.cols];
        unit_vector[f_mat.cols - 1] = Integer::from(1);
        f_mat.set_row(dim, &unit_vector);

        // Perform setup
        println!("start protocol_setup");
        let start = SystemTime::now();
        let ((dcr_sk, dcr_pk), 
            qfe_sk_dcr_to_qe, 
            // qfe_sk_qe_to_qe,
            // qfe_sk_end,
        ) = protocol_setup_new(
            vec![dim, dim, dim], 
            f_num, 
            k,
            &sk_bound, 
            &grp, 
            &mut rng
        );

        let (h0_left, h0_right) = sample_h(dim, k, &grp.n, &mut rng);
        // let (h1_left, h1_right) = sample_h(dim, k, &grp.n, &mut rng);
        let (gamma_left, gamma_right) = sample_gamma(dim, dim, &grp.n, &mut rng);

        // remove random for test
        // println!(" !!!!! REMOVE RANDOM FOR TESTING in geenrating H, gamma for protocol_setup !!!!! ");
        // let (h0_left, h0_right) = sample_h_trivial(dim, k);
        // let (gamma_left, gamma_right) = sample_gamma_trivial(dim, dim);
        // let (h1_left, h1_right) = sample_h_trivial(dim, k);
        // let h1_left = h0_left.clone();
        // let h1_right = h0_right.clone();
        // swap gamma_left and gamma_right and transpose
        let tmp = gamma_left.clone();
        let gamma_left = gamma_right.clone();
        let gamma_right = tmp.clone();
        let gamma_left = gamma_left.transpose();
        let gamma_right = gamma_right.transpose();

        println!("gamma_left = \n{}", gamma_left);
        println!("gamma_right = \n{}", gamma_right);
        println!("gamma mult text");
        let mut mult = gamma_left.clone() * gamma_right.clone();
        mult.mod_inplace(&grp.n);
        println!("{}", mult);


        let end = start.elapsed();
        println!("Time elapsed in protocol_setup is: {:?}", end);

        println!("start protocol_keygen_dcr_to_qfe");
        let start = SystemTime::now();
        let (mat_ctxts, fk)
        = protocol_keygen_dcr_to_qfe(
            &dcr_sk, 
            &qfe_sk_dcr_to_qe,
            &h0_right, 
            &gamma_left, 
            &decomp, 
            &grp, 
            &mut rng);
        let end = start.elapsed();
        println!("Time elapsed in protocol_keygen_dcr_to_qfe is: {:?}", end);

        println!(" ======= ");
        println!("test mat_ctxts * x");
        let mut x1 = x.clone();
        x1.push(Integer::from(1));
        println!("mat_ctxts_size = {} x {}", mat_ctxts.rows, mat_ctxts.cols);
        println!("x1_size = {}", x1.len());
        let mat_ctxts_gamma = mat_ctxts.clone() * gamma_right.clone();
        let mut mat_ctxts_x = mat_ctxts_gamma * x1.clone();
        vec_mod(&mut mat_ctxts_x, &grp.n);
        let mat_ctxts_x: Vec<Integer> = decomp.vector_inv(&mat_ctxts_x);

        // println!("start protocol_keygen_qfe_to_qfe");
        // let start = SystemTime::now();
        // let fk_qfe_to_qfe = protocol_keygen_qfe_to_qfe(
        //     &qfe_sk_dcr_to_qe,
        //     &qfe_sk_dcr_to_qe,
        //     &h0_right,
        //     &h0_left,
        //     &f_mat,
        //     &decomp,
        //     &grp,
        //     &mut rng
        // );
        // let end = start.elapsed();
        // println!("Time elapsed in protocol_keygen_qfe_to_qfe is: {:?}", end);

        println!("start protocol_keygen_qfe_to_plain");
        let start = SystemTime::now();
        let fk_qfe_to_plain = protocol_keygen_qfe_to_plain(
            &qfe_sk_dcr_to_qe, // changed
            &h0_left,
            &f_mat,
            &decomp,
            &grp
        );
        println!("fk_qfe_to_plain size = {} x {}", fk_qfe_to_plain.rows, fk_qfe_to_plain.cols);
        let end = start.elapsed();
        println!("Time elapsed in protocol_keygen_qfe_to_plain is: {:?}", end);


        println!("start protocol_enc_init");
        let start = SystemTime::now();
        let ctxt_x = protocol_enc_init(&dcr_pk, &gamma_right, &x, &grp, &mut rng);
        let end = start.elapsed();
        println!("Time elapsed in protocol_enc_init is: {:?}", end);
        
        println!("start protocol_dec_dcr_to_qfe");
        let start = SystemTime::now();
        //  run keyswitch
        let ct0 = protocol_dec_dcr_to_qfe(
            &ctxt_x, 
            &mat_ctxts,
            &fk,
            &decomp,
            &grp
        );
        let end = start.elapsed();
        println!("Time elapsed in protocol_dec_dcr_to_qfe is: {:?}", end);

        println!("test ct0 vs mat_ctxts_x");
        for i in 0..mat_ctxts_x.len() {
            // println!("{}th : {} vs {}", i, ct0[i], mat_ctxts_x[i]);
            assert_eq!(ct0[i], mat_ctxts_x[i]);
        }
        println!("ct0 == mat_ctxts_x test passed");

        let f_trivial = vec![Integer::from(1);qfe_sk_dcr_to_qe.dim * qfe_sk_dcr_to_qe.dim];
        let fk = qfe_keygen(&qfe_sk_dcr_to_qe, &&f_trivial, &grp);

        let dec_test_out = qfe_dec( &fk, &ct0,qfe_sk_dcr_to_qe.dim, qfe_sk_dcr_to_qe.q, &grp);
        println!("f_trivial = {:?}", f_trivial);
        println!("dec_test_out = {:?}", dec_test_out);

        let mut f_one = vec![Integer::from(0); qfe_sk_dcr_to_qe.dim * qfe_sk_dcr_to_qe.dim];
        f_one[0] = Integer::from(1);
        let fk_one = qfe_keygen(&qfe_sk_dcr_to_qe, &&f_one, &grp);
        let dec_test_out_one = qfe_dec( &fk_one, &ct0,qfe_sk_dcr_to_qe.dim, qfe_sk_dcr_to_qe.q, &grp);
        println!("f_one = {:?}", f_one);
        println!("dec_test_out_one = {:?}", dec_test_out_one);

        // print h0_right * gamma_left * gamma_right * x
        println!("h0_right size = {} x {}", h0_right.rows, h0_right.cols);
        println!("gamma_left size = {} x {}", gamma_left.rows, gamma_left.cols);
        println!("gamma_right size = {} x {}", gamma_right.rows, gamma_right.cols);
        println!("x1 size = {}", x1.len());
        let mut mult = h0_right.clone() * gamma_left.clone();
        mult = mult * gamma_right.clone();
        let mut h0_right_x = mult * x1.clone();
        vec_mod(&mut h0_right_x, &grp.n);
        println!("h0_right_x = {:?}", h0_right_x);
        let mut h0_right_x_tensor  = tensor_product_vecs(&h0_right_x, &h0_right_x, &grp.delta);
        println!("h0_right_x_tensor = {:?}", h0_right_x_tensor);
        println!("sum = {:?}", h0_right_x_tensor.iter().sum::<Integer>());
        let ctxt_h0_right_x = qfe_enc(&qfe_sk_dcr_to_qe, &h0_right_x, &grp, &mut rng);
        println!("compare ctxt_h0_right_x vs ct0");
        for i in 0..ctxt_h0_right_x.len() {
            println!("{}th: {}", i, ctxt_h0_right_x[i] == ct0[i]);
            println!("bit lengths: {} vs {}", ctxt_h0_right_x[i].significant_bits(), ct0[i].significant_bits());
        }


        // println!("start protocol_dec_qfe_to_qfe");
        // let start = SystemTime::now();
        // let ct0 = protocol_dec_qfe(
        //     &ct0,
        //     &fk_qfe_to_qfe,
        //     dim + k,
        //     q,
        //     &decomp,
        //     &grp,
        //     true,
        // );
        // let end = start.elapsed();
        // println!("Time elapsed in protocol_dec_qfe_to_qfe is: {:?}", end);

        println!("start protocol_dec_qfe_to_plain");
        let start = SystemTime::now();
        let mut val_end = protocol_dec_qfe(
            &ct0,
            &fk_qfe_to_plain,
            dim + k,
            q,
            &decomp,
            &grp,
            false,
        );

        let mut val_end2 = protocol_dec_qfe(
            &ctxt_h0_right_x,
            &fk_qfe_to_plain,
            dim + k,
            q,
            &decomp,
            &grp,
            false,
        );
        let end = start.elapsed();
        println!("Time elapsed in protocol_dec_qfe_to_plain is: {:?}", end);

        println!("reprint inputs");
        println!("x: {:?}", x);
        println!("f_mat: {}", f_mat);
        // let fx = eval_quadratic_multivariate(&x, &x, &f_mat_no_padding);
        // println!("f(x) = {:?}", fx);
        // let mut ffx = eval_quadratic_multivariate(&fx, &fx, &f_mat_no_padding);
        // println!("f(f(x)) = {:?}", ffx);
        // vec_mod(&mut ffx, &grp.n);

        let mut x_one = x.clone();
        x_one.push(Integer::from(1));
        println!("x_one = {:?}", x_one);
        let fx_one = eval_quadratic_multivariate(&x_one, &x_one, &f_mat);
        println!("f(x) (1 padding) = {:?}", fx_one);
        let ffx_one = eval_quadratic_multivariate(&fx_one, &fx_one, &f_mat);
        // println!("f(f(x)) (1 padding) = {:?}", ffx_one);

        // println!("real mod n = {:?}", ffx);
        // println!("eval result: {:?}", val_end);
        vec_mod(&mut val_end, &grp.n);
        vec_mod(&mut val_end2, &grp.n);
        println!("eval result mod n = {:?}", val_end);
        println!("eval with h0_right_x: {:?}", val_end2);
        println!("{}th try", try_num);
        // assert ffx_one == val_end
        for i in 0..val_end.len() {
            assert_eq!(fx_one[i], val_end[i]);
        }
    }
    println!(" test end: all passed");
    }

    #[test]
    fn test_matrix_h_and_gamma() {
        let dim = 5;
        let k = 2;
        let n = Integer::from(101);

        let mut rng = RandState::new(); // Create a single RandState object
        // let d = SystemTime::now()
        // .duration_since(SystemTime::UNIX_EPOCH)
        // .expect("Duration since UNIX_EPOCH failed");
        // rng.seed(&Integer::from(d.as_secs()));

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

        println!("modulo = {}", n);
    }

    #[test]
    fn test_echelon_form() {
        let dim = 5;
        let mut rng = RandState::new(); // Create a single RandState object
        let d = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .expect("Duration since UNIX_EPOCH failed");
        rng.seed(&Integer::from(d.as_secs()));
        
        let mut rank: i32 = -1;
        let mut rref = Matrix::new(1, 1);
        let mut m = Matrix::new(1, 1);
        let mut row_op = Matrix::new(1, 1);
        while rank < (dim as i32) {
            m = Matrix::random(dim + 1, dim, &Integer::from(3), &mut rng);
            m.add_int_inplace(&Integer::from(-1));
            m.mod_inplace(&Integer::from(11));

            rref = m.clone();
            let (new_row_op, new_rank) = row_reduce_form(&mut rref, &Integer::from(11));
            rank = new_rank;
            row_op = new_row_op;
        }
        println!("m = {}", m);
        
        println!("rref = {}", rref);
        println!("row_op = {}", row_op);

        let mut m_left_inv = Matrix::new(dim, dim + 1);
        let mut counter = 0;
        let mut i = 0;
        while counter < dim {
            println!("check {} {}", i, counter);
            let val = rref.get(i, counter);
            if val != Integer::from(0) {
                let row = row_op.get_row(i);
                m_left_inv.set_row(counter, &row);
                counter += 1;
            }
            i += 1;
        } 

        println!("m_left_inv = {}", m_left_inv);

        // test mult
        let mut m2 = m_left_inv * m.clone();
        m2.mod_inplace(&Integer::from(11));
        println!("m2 = \n{}", m2);

    }

    #[test]
    fn test_row_reduce_echelon_form_integer() {
        let dim = 2;
        let mut rng = RandState::new(); // Create a single RandState object
        let d = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .expect("Duration since UNIX_EPOCH failed");
        rng.seed(&Integer::from(d.as_secs()));
        
        let mut rank: i32 = -1;
        let mut rref = Matrix::new(1, 1);
        let mut m = Matrix::new(1, 1);
        let mut row_op = Matrix::new(1, 1);
        while rank < (dim as i32) {
            m = Matrix::random(dim + 1, dim, &Integer::from(3), &mut rng);
            m.add_int_inplace(&Integer::from(-1));
            rref = m.clone();
            let (new_row_op, new_rank) = row_reduce_form_integer(&mut rref);
            rank = new_rank;
            row_op = new_row_op;
        }
        println!("m = \n{}", m);
        println!("rref = \n{}", rref);
        println!("row_op = \n{}", row_op);

        let mult = row_op.clone() * m.clone();
        println!("mult = \n{}", mult);

        let mut m_left_inv = Matrix::new(dim, dim + 1);
        let mut counter = 0;
        let mut i = 0;
        while counter < dim {
            println!("check {} {}", i, counter);
            let val = rref.get(i, counter);
            if val != Integer::from(0) {
                let row = row_op.get_row(i);
                m_left_inv.set_row(counter, &row);
                counter += 1;
            }
            i += 1;
        } 

        println!("m_left_inv = {}", m_left_inv);

        // test mult
        let m2 = m_left_inv * m.clone();
        println!("m2 = \n{}", m2);
    }
}