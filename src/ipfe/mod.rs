pub mod keys;
pub mod scheme;

#[cfg(test)]
mod tests {
    use crate::util::group::Group;

    use rug::Integer;
    use rug::rand::RandState;
    use crate::ipfe::scheme::{ipe_setup, ipe_keygen, ipe_enc, ipe_enc_matrix_expression, ipe_dec};
    use crate::util::vector::{gen_random_vector, vec_mod, int_mod};
    use std::time::SystemTime;
    
    #[test]
    fn test_ipe_start_to_end() {
        let mut rand = RandState::new(); // Create a single RandState object
        let d = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .expect("Duration since UNIX_EPOCH failed");
        rand.seed(&Integer::from(d.as_secs()));
        
        let grp = Group::new(100); // Initialize the group        
        
        let dim = 10;
        let q = 2 * dim + 1;
        let bound = 1024;
        // Generate random inputs for sk and y
        let mut x = gen_random_vector(dim, &Integer::from(2 * bound), &mut rand);
        let mut y = gen_random_vector(dim, &Integer::from(2 * bound), &mut rand);

        for i in 0..dim {
            x[i] -= Integer::from(bound / 2);
            y[i] -= Integer::from(bound / 2);
        }        

        // Perform setup
        println!("start ipe_setup");
        let start = SystemTime::now();
        let sk = ipe_setup(&grp, dim, q, &mut rand);
        let end = start.elapsed();
        println!("Time elapsed in ipe_setup is: {:?}", end);

        // Perform key generation
        println!("start ipe_keygen");
        let start = SystemTime::now();
        let sk_f = ipe_keygen(&sk, &y, &grp);
        let end = start.elapsed();
        println!("Time elapsed in ipe_keygen is: {:?}", end);

        // Perform encryption
        println!("start ipe_enc");
        let original = false;
        let mut ctxt: Vec<Integer>;
        let start = SystemTime::now();
        if original {
            println!("run original enc");
            ctxt = ipe_enc(&sk, &x, &grp, true, &mut rand);
        } else {
            println!("run matrix enc");
            let ipe_enc_mat = ipe_enc_matrix_expression(&sk, &grp, true, &mut rand);
            let end = start.elapsed();
            println!("Time elapsed in ipe_enc_matrix_expression is: {:?}", end);

            // println!("ipe_enc_mat: {}", ipe_enc_mat);

            let mut x1 = x.clone();
            x1.push(Integer::from(1));
            // let start = SystemTime::now();
            ctxt = ipe_enc_mat * x1.clone();
            vec_mod(&mut ctxt, &grp.delta);
            // println!("ctxt: {:?}", ctxt);
        }
        let end = start.elapsed();
        println!("Time elapsed in enc (mat_mult) is: {:?}", end);        

        // Perform decryption
        println!("start ipe_dec");
        let start = SystemTime::now();
        let out = ipe_dec(&sk_f, &ctxt, &grp, true);
        let end = start.elapsed();
        println!("Time elapsed in ipe_dec is: {:?}", end);
        
        // check inner product between x and y
        let mut out2 = Integer::from(0);
        for i in 0..x.len() {
            out2 += x[i].clone() * y[i].clone();
        }

        println!("x = {:?}", x);
        println!("y = {:?}", y);
        println!("out eval  : {}", out);
        println!("out plain : {}", out2);

        let out_mod = int_mod(&out, &grp.n);
        println!("eval mod n: {}", out_mod);
        let out2_mod = int_mod(&out2, &grp.n);
        println!("plain mod n: {}", out2_mod);

        assert_eq!(out_mod, out2_mod);

    }
    
}
