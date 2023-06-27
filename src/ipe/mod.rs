pub mod keys;
pub mod scheme;

#[cfg(test)]
mod tests {
    use crate::util::group::Group;

    use rug::Integer;
    use rug::rand::RandState;
    use crate::ipe::scheme::{ipe_setup, ipe_keygen, ipe_enc, ipe_enc_matrix_expression, ipe_dec};
    use crate::util::vector::{gen_random_vector, vec_mod};
    use std::time::SystemTime;
    
    #[test]
    fn test_ipe_start_to_end() {
        let mut rand = RandState::new(); // Create a single RandState object
        let d = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .expect("Duration since UNIX_EPOCH failed");
        rand.seed(&Integer::from(d.as_secs()));
        
        let grp = Group::new(10); // Initialize the group        
        // println!("{}", grp);
        
        let dim = 10;
        let b = dim * 2 + 1;
        let bound = 10;
        // Generate random inputs for sk and y
        let mut x = gen_random_vector(dim, &Integer::from(bound), &mut rand);
        let mut y = gen_random_vector(dim, &Integer::from(bound), &mut rand);

        for i in 0..dim {
            x[i] -= Integer::from(bound / 2);
            y[i] -= Integer::from(bound / 2);
        }        

        // Perform setup
        let start = SystemTime::now();
        let sk = ipe_setup(&grp, dim, b, &mut rand);
        let end = start.elapsed();
        println!("Time elapsed in ipe_setup is: {:?}", end);

        // Perform key generation
        let start = SystemTime::now();
        let sk_f = ipe_keygen(&sk, &y, &grp);
        let end = start.elapsed();
        println!("Time elapsed in ipe_keygen is: {:?}", end);
        
        // Print the generated secret key
        // println!("sk_f: {:?}", sk_f);

        // Perform encryption
        let original = false;
        let mut ctxt: Vec<Integer> = Vec::new();
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

        assert_eq!(out, out2);

    }
    
}
