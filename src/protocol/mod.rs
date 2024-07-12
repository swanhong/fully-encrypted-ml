pub mod scheme;
pub mod algorithm;

#[cfg(test)]
pub mod test {
    use rug::Integer;
    use rug::rand::RandState;
    use std::time::SystemTime;

    use crate::protocol::algorithm::{row_reduce_form, row_reduce_form_integer, sample_h, sample_gamma};
    use crate::protocol::scheme::*;
    use crate::util::group::Group;

    use crate::util::matrix::{concatenate_row, eval_quadratic_multivariate, Matrix};
    use crate::util::vector::{gen_random_vector, tensor_product_vecs};
    use crate::util::decomp::Decomp;

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
        
        let dim = 3;
        let k = 1;
        let depth = dim;
        let bound = 5;
        let sk_bound = Integer::from(10);
        let decomp = Decomp::new(3, &grp);
        
        let x = gen_random_vector(dim, &Integer::from(bound), &mut rng);
        // let x = vec![Integer::from(0); dim];
        let mut f_mat = Matrix::new(dim, (dim + 1) * (dim + 1));
        let mut f_mat_no_padding = Matrix::new(dim, dim * dim);
        for i in 0..dim { 
            let mut f = gen_random_vector(dim, &Integer::from(bound), &mut rng);
            let f_tensor = tensor_product_vecs(&f, &f, &grp.delta);
            f_mat_no_padding.set_row(i, &f_tensor);
            f.push(Integer::from(1));
            let f_tensor = gen_random_vector((dim + 1) * (dim + 1), &Integer::from(bound), &mut rng);
            f_mat.set_row(i, &f_tensor);
        }
        let f_mat_origin = f_mat.clone();
        let mut unit_vector = Matrix::new(1, f_mat.cols);
        unit_vector.set(0, f_mat.cols - 1, Integer::from(1));
        let f_mat = concatenate_row(&f_mat, &unit_vector);

        // Perform setup
        println!("start protocol_setup");
        let start = SystemTime::now();
        let ((dcr_sk, dcr_pk), 
            qfe_sk, 
        ) = protocol_setup(
            &vec![dim; depth + 1],  
            k,
            &sk_bound, 
            &grp, 
            &mut rng
        );
        println!("qfe_sk.len = {}", qfe_sk.len());

        let mut h_left = vec![Matrix::new(0, 0); depth];
        let mut h_right = vec![Matrix::new(0, 0); depth];
        for i in 0..depth {
            let (h_left_i, h_right_i) = sample_h(dim, k, &grp.delta, &mut rng);
            h_left[i] = h_left_i;
            h_right[i] = h_right_i;
        }
        let (gamma_left, gamma_right) = sample_gamma(dim, &grp.n, &mut rng);

        let end = start.elapsed();
        println!("Time elapsed in protocol_setup is: {:?}", end);

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
            &mut rng);
        let end = start.elapsed();
        println!("Time elapsed in protocol_keygen_dcr_to_qfe is: {:?}", end);

        println!("start protocol_keygen_qfe_to_qfe");
        let mut fk_qfe_to_qfe = vec![Matrix::new(0, 0); depth - 1];
        let start = SystemTime::now();
        for i in 0..fk_qfe_to_qfe.len() {
            fk_qfe_to_qfe[i] = protocol_keygen_qfe_to_qfe(
                &qfe_sk[i],
                &qfe_sk[i + 1],
                &h_right[i + 1],
                &h_left[i],
                &f_mat,
                &decomp,
                &grp,
                &mut rng
            );
        }
        // let fk_qfe_to_qfe = protocol_keygen_qfe_to_qfe(
        //     &qfe_sk[0],
        //     &qfe_sk[1],
        //     &h1_right,
        //     &h0_left,
        //     &f_mat,
        //     &decomp,
        //     &grp,
        //     &mut rng
        // );
        let end = start.elapsed();
        println!("Time elapsed in protocol_keygen_qfe_to_qfe is: {:?}", end);

        println!("start protocol_keygen_qfe_to_plain");
        let start = SystemTime::now();
        let mut fk_qfe_to_plain_test = vec![Matrix::new(0, 0); depth - 1];
        for i in 0..fk_qfe_to_plain_test.len() {
            fk_qfe_to_plain_test[i] = protocol_keygen_qfe_to_plain(
                &qfe_sk[i], // changed
                &h_left[i],
                &f_mat_origin,
                &grp
            );
        }
        // let fk_qfe_to_plain_test = protocol_keygen_qfe_to_plain(
        //     &qfe_sk[0], // changed
        //     &h0_left,
        //     &f_mat,
        //     &grp
        // );
        println!("test done");
        let fk_qfe_to_plain = protocol_keygen_qfe_to_plain(
            &qfe_sk[depth - 1],
            &h_left[depth - 1],
            &f_mat_origin,
            &grp
        );
        let end = start.elapsed();
        println!("Time elapsed in protocol_keygen_qfe_to_plain is: {:?}", end);


        println!("start protocol_enc_init");
        let start = SystemTime::now();
        let ct0 = protocol_enc_init(&dcr_pk, &gamma_right, &x, &grp, &mut rng);
        let end = start.elapsed();
        println!("Time elapsed in protocol_enc_init is: {:?}", end);
        
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
        let end = start.elapsed();
        println!("Time elapsed in protocol_dec_dcr_to_qfe is: {:?}", end);

        // println!("start protocol_dec_qfe_to_plain");
        // let start = SystemTime::now();
        // let val_end_fx = protocol_dec_qfe(
        //     &ct0,
        //     &fk_qfe_to_plain_test,
        //     dim + k + 1,
        //     2 * (dim + k + 1) + 1,
        //     &decomp,
        //     &grp,
        //     false,
        // );
        // let end = start.elapsed();
        // println!("Time elapsed in protocol_dec_qfe_to_plain is: {:?}", end);


        println!("start protocol_dec_qfe_to_qfe");
        let start = SystemTime::now();
        for i in 0..depth - 1 {
            println!("ct_vec{} <- ct_vec{}", i + 1, i);
            ct_vec[i + 1] = protocol_dec_qfe(
                &ct_vec[i],
                &fk_qfe_to_qfe[i],
                dim + k + 1,
                2 * (dim + k + 1) + 1,
                &decomp,
                &grp,
                true,
            );
        }
        // let ct0 = protocol_dec_qfe(
        //     &ct0,
        //     &fk_qfe_to_qfe,
        //     dim + k + 1,
        //     2 * (dim + k + 1) + 1,
        //     &decomp,
        //     &grp,
        //     true,
        // );
        let end = start.elapsed();
        println!("Time elapsed in protocol_dec_qfe_to_qfe is: {:?}", end);

        println!("start protocol_dec_qfe_to_plain");
        let mut val_end = vec![vec![Integer::from(0);1]; depth];
        let start = SystemTime::now();
        for i in 0..depth - 1 {
            val_end[i] = protocol_dec_qfe(
                &ct_vec[i],
                &fk_qfe_to_plain_test[i],
                dim + k + 1,
                2 * (dim + k + 1) + 1,
                &decomp,
                &grp,
                false,
            );
        }
        let val_end_final = protocol_dec_qfe(
            &ct_vec[depth - 1],
            &fk_qfe_to_plain,
            dim + k + 1,
            2 * (dim + k + 1) + 1,
            &decomp,
            &grp,
            false,
        );
        let end = start.elapsed();
        println!("Time elapsed in protocol_dec_qfe_to_plain is: {:?}", end);

        val_end[depth - 1] = val_end_final.clone();
        println!("val_ewnd = {:?}", val_end_final);

        println!("reprint inputs");
        println!("x: {:?}", x);
        println!("f_mat: {}", f_mat);

        let mut x_one = x.clone();
        x_one.push(Integer::from(1));
        println!("x_one = {:?}", x_one);
        let mut fx_one = x_one.clone();
        for i in 0..depth {
            fx_one = eval_quadratic_multivariate(&fx_one, &fx_one, &f_mat);
            println!("f^{i}(x) = {:?}", fx_one);
            println!("val_end[{}]: {:?}", i, val_end[i]);
            // assert fx_one == val_end
            for j in 0..val_end[i].len() {
                assert_eq!(fx_one[j], val_end[i][j]);
            }
        }
        // let fx_one = eval_quadratic_multivariate(&x_one, &x_one, &f_mat);
        // println!("f(x) (1 padding) = {:?}", fx_one);
        // let ffx_one = eval_quadratic_multivariate(&fx_one, &fx_one, &f_mat);
        // println!("f(f(x)) (1 padding) = {:?}", ffx_one);

        // assert ffx_one == val_end
        // for i in 0..val_end.len() {
        //     assert_eq!(ffx_one[i], val_end[i]);
        // }
    }

    #[test]
    fn test_matrix_h_and_gamma() {
        let dim = 5;
        let k = 2;
        let n = Integer::from(101);

        let mut rng = RandState::new(); // Create a single RandState object
        let d = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .expect("Duration since UNIX_EPOCH failed");
        rng.seed(&Integer::from(d.as_secs()));

        let (h_left_1, h_right_1) = sample_h(dim, k, &n, &mut rng);

        println!("h_left_1:");
        println!("{}", h_left_1);

        println!("h_right_1:");
        println!("{}", h_right_1);

        println!("test mult");
        let mut mult = h_left_1.clone() * h_right_1.clone();
        mult.mod_inplace(&n);
        println!("{}", mult);

        let (gamma_left_1, gamma_right_1) = sample_gamma(dim, &n, &mut rng);

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