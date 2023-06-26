pub mod keys;
pub mod scheme;

mod tests {
    use crate::util::group::Group;

    use rug::Integer;
    use rug::rand::RandState;
    use std::time::{Duration, SystemTime};
    use crate::qe::scheme::{qe_setup, qe_keygen, qe_enc, qe_enc_matrix_same_xy, qe_dec, qe_enc_matrix_expression};
    use crate::util::vector::{gen_random_vector, vec_mod, vec_mul_scalar};
    use crate::util::decomp::Decomp;
    #[test]
    fn test_qe_start_to_end() {
        let mut rand = RandState::new(); // Create a single RandState object
        // let d = SystemTime::now()
        // .duration_since(SystemTime::UNIX_EPOCH)
        // .expect("Duration since UNIX_EPOCH failed");
        // rand.seed(&Integer::from(d.as_secs()));
        
        let grp = Group::new(10); // Initialize the group        
        // println!("{}", grp);
        
        let n_x = 6;
        let n_y = 6;
        let b = 2 * n_x + 1;
        // let b = 1;
        let bound = 3;

        let decomp = Decomp::new(4, &grp);
        // println!("decomp = {}", decomp);

        let mut x = gen_random_vector(n_x, &Integer::from(bound), &mut rand);
        let mut y = gen_random_vector(n_y, &Integer::from(bound), &mut rand);
        let mut f = gen_random_vector(n_x * n_y, &Integer::from(bound), &mut rand);
        // x[0] = Integer::from(12);
        // y[0] = Integer::from(4);
        // f[0] = Integer::from(9);

        // Perform setup
        let start = SystemTime::now();
        let qe_sk = qe_setup(&grp, n_x, n_y, b, &mut rand);
        let end = start.elapsed();
        println!("Time elapsed in qe_setup is: {:?}", end);

        // println!("{}", qe_sk);

        // Perform key generation
        let start = SystemTime::now();
        let (sk_f, sk_red) = qe_keygen(&qe_sk, &f, &grp, &decomp);
        let end = start.elapsed();
        println!("Time elapsed in qe_keygen is: {:?}", end);

        // println!("sk_f: {:?}", sk_f);
        // println!("sk_red: {:?}", sk_red);
        println!("lengths are {} and {}", sk_f.len(), sk_red.len());

        let original = false;
        
        let mut ctxt_x = Vec::new();
        let mut ctxt_y = Vec::new();
        let mut ctxt_h = Vec::new(); 
        if original {
            let start = SystemTime::now();
            let (enc_x, enc_y, enc_h) = qe_enc(&qe_sk, &x, &y, &grp, &decomp, &mut rand);
            let end = start.elapsed();
            println!("Time elapsed in qe_enc is: {:?}", end);
            ctxt_x = enc_x;
            ctxt_y = enc_y;
            ctxt_h = enc_h;
        } else {
            // Perform encryption
            let start = SystemTime::now();
            let (enc_mat_x, enc_mat_y, enc_mat_h) = qe_enc_matrix_expression(&qe_sk, n_x, n_y, &grp, &decomp, &mut rand);
            // let (enc_mat_x, enc_mat_y, enc_mat_h) = qe_enc_matrix_same_xy(&qe_sk, n_x, &grp, &decomp, &mut rand);
            let end = start.elapsed();
            println!("Time elapsed in qe_enc_matrix is: {:?}", end);

            // println!("enc_mat_x: {}", enc_mat_x);
            // println!("enc_mat_y: {}", enc_mat_y);
            // println!("enc_mat_h: {}", enc_mat_h);

            println!("sizes are {} and {}", enc_mat_x.rows, enc_mat_x.cols);
            println!("sizes are {} and {}", enc_mat_y.rows, enc_mat_y.cols);
            println!("sizes are {} and {}", enc_mat_h.rows, enc_mat_h.cols);

            // generate ctxt from matrix representation
            let x_mu: Vec<Integer> = vec_mul_scalar(&x, &grp.mu);
            let mut x_mu_1 = x_mu.clone();
            x_mu_1.push(Integer::from(1));
            let mut x1 = x.clone();
            x1.push(Integer::from(1));
            let mut y1 = y.clone();
            y1.push(Integer::from(1));
            
            // println!("x_mu_1 = {:?}", x_mu_1);
            // println!("y1 = {:?}", y1);

            // y_x_1_1 = (y || x || 1 || 1)
            let mut y_xmu_1_1 = y.clone();
            y_xmu_1_1.extend(x_mu.clone());
            y_xmu_1_1.push(Integer::from(1));
            y_xmu_1_1.push(Integer::from(1));

            ctxt_x = enc_mat_x.mul_vec(&x1);
            vec_mod(&mut ctxt_x, &grp.delta);
            ctxt_y = enc_mat_y.mul_vec(&y1);
            vec_mod(&mut ctxt_y, &grp.delta);
            ctxt_h = enc_mat_h.mul_vec(&y_xmu_1_1);
            vec_mod(&mut ctxt_h, &grp.delta);
        }
        // println!("ctxt_x (size = {}): {:?}", ctxt_x.len(), ctxt_x);
        // println!("ctxt_y (size = {}): {:?}", ctxt_y.len(), ctxt_y);
        // println!("ctxt_h (size = {}): {:?}", ctxt_h.len(), ctxt_h);
        


        let start = SystemTime::now();
        let out = qe_dec((&sk_f, &sk_red), (&ctxt_x, &ctxt_y, &ctxt_h), &grp, &decomp);
        let end = start.elapsed();
        println!("Time elapsed in qe_dec is: {:?}", end);


        // compute out2 = inner product <f, x tensor y>
        let mut out2 = Integer::from(0);
        for i in 0..n_x {
            for j in 0..n_y {
                out2 += f[i * n_y + j].clone() * x[i].clone() * y[j].clone();
            }
        }
        println!("x = {:?}", x);
        println!("y = {:?}", y);
        println!("f = {:?}", f);
        println!("out: {}", out);
        println!("out2: {}", out2);
        assert_eq!(out, out2);

    }

    #[test]
    fn test_qe_samexy_start_to_end() {
        let mut rand = RandState::new(); // Create a single RandState object
        let d = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .expect("Duration since UNIX_EPOCH failed");
        rand.seed(&Integer::from(d.as_secs()));
        
        let grp = Group::new(10); // Initialize the group        
        // println!("{}", grp);
        
        let n_x = 4;
        let b = 2 * n_x + 1;
        // let b = 1;
        let bound = 10;

        let decomp = Decomp::new(4, &grp);

        let mut x = gen_random_vector(n_x, &Integer::from(bound), &mut rand);
        let mut f = gen_random_vector(n_x * n_x, &Integer::from(bound), &mut rand);


        // compute out2 = inner product <f, x tensor y>
        let mut out2 = Integer::from(0);
        for i in 0..n_x {
            for j in 0..n_x {
                out2 += f[i * n_x + j].clone() * x[i].clone() * x[j].clone();
            }
        }

        assert!(out2 < grp.n, "out2 = {} >= n = {}", out2, grp.n);

        // Perform setup
        let start = SystemTime::now();
        let qe_sk = qe_setup(&grp, n_x, n_x, b, &mut rand);
        let end = start.elapsed();
        println!("Time elapsed in qe_setup is: {:?}", end);

        // println!("{}", qe_sk);

        // Perform key generation
        let start = SystemTime::now();
        let (sk_f, sk_red) = qe_keygen(&qe_sk, &f, &grp, &decomp);
        let end = start.elapsed();
        println!("Time elapsed in qe_keygen is: {:?}", end);

        // println!("sk_f: {:?}", sk_f);
        // println!("sk_red: {:?}", sk_red);
        // println!("lengths are {} and {}", sk_f.len(), sk_red.len());

        
        let mut ctxt_x = Vec::new();
        let mut ctxt_y = Vec::new();
        let mut ctxt_h = Vec::new(); 
        // Perform encryption
        let start = SystemTime::now();
        let (enc_mat_x, enc_mat_y, enc_mat_h) = qe_enc_matrix_same_xy(&qe_sk, n_x, &grp, &decomp, &mut rand);
        let end = start.elapsed();
        println!("Time elapsed in qe_enc_matrix is: {:?}", end);

        // println!("enc_mat_x: {}", enc_mat_x);
        // println!("enc_mat_y: {}", enc_mat_y);
        // println!("enc_mat_h: {}", enc_mat_h);

        // generate ctxt from matrix representation
        // let x_mu: Vec<Integer> = vec_mul_scalar(&x, &grp.mu);
        // let mut x_mu_1 = x_mu.clone();
        // x_mu_1.push(Integer::from(1));
        let mut x1 = x.clone();
        x1.push(Integer::from(1));
        
        // println!("x_mu_1 = {:?}", x_mu_1);
        // println!("y1 = {:?}", y1);

        // y_x_1_1 = (y || x || 1 || 1)
        // let mut y_xmu_1_1 = x.clone();
        // y_xmu_1_1.extend(x_mu.clone());
        // y_xmu_1_1.push(Integer::from(1));
        // y_xmu_1_1.push(Integer::from(1));

        ctxt_x = enc_mat_x.mul_vec(&x1);
        vec_mod(&mut ctxt_x, &grp.delta);
        ctxt_y = enc_mat_y.mul_vec(&x1);
        vec_mod(&mut ctxt_y, &grp.delta);
        ctxt_h = enc_mat_h.mul_vec(&x1);
        vec_mod(&mut ctxt_h, &grp.delta);

        // println!("ctxt_x: {:?}", ctxt_x);
        // println!("ctxt_y: {:?}", ctxt_y);
        // println!("ctxt_h: {:?}", ctxt_h);
        
        let start = SystemTime::now();
        let out = qe_dec((&sk_f, &sk_red), (&ctxt_x, &ctxt_y, &ctxt_h), &grp, &decomp);
        let end = start.elapsed();
        println!("Time elapsed in qe_dec is: {:?}", end);

        println!("x = {:?}", x);
        println!("f = {:?}", f);
        println!("out: {}", out);
        println!("out2: {}", out2);
        assert_eq!(out, out2);

    }
}