mod keys;
mod ipe;

#[cfg(test)]
mod tests {
    use crate::util::group::Group;
    use crate::ipe::keys::IPE_sk;
    #[test]
    fn test_ipe_sk() {
        let dim = 5;
        let grp = Group::new(10);
        let b = 4;
        
        // Create an instance of IPE_sk
        let sk = IPE_sk::new(dim, b, &grp);
        
        // Print the instance
        println!("{}", sk);
        
        // Perform assertions or additional test operations
        // ...
    }

    use rug::Integer;
    use rug::rand::RandState;
    use std::time::SystemTime;
    use crate::ipe::ipe::{ipe_setup, ipe_keygen, ipe_enc_matrix_expression, ipe_dec};
    use crate::util::vector::{gen_random_vector, mod_vec};
    
    #[test]
    fn test_ipe_all() {
        let mut rand = RandState::new(); // Create a single RandState object
        // let d = SystemTime::now()
        // .duration_since(SystemTime::UNIX_EPOCH)
        // .expect("Duration since UNIX_EPOCH failed");
        // rand.seed(&Integer::from(d.as_secs()));
        
        let grp = Group::new(10); // Initialize the group        
        println!("{}", grp);
        
        let dim = 5;
        let bound = 10;
        // Generate random inputs for sk and y
        let x = gen_random_vector(dim, &Integer::from(bound), &mut rand);
        let y = gen_random_vector(dim, &Integer::from(bound), &mut rand);
        
        println!("x = {:?}", x);
        println!("y = {:?}", y);

        let sk = ipe_setup(&grp, dim, 1);

        // Print the instance
        println!("{}", sk);

        // Perform key generation
        let sk_f = ipe_keygen(&sk, &y, &grp);
        
        // Print the generated secret key
        println!("sk_f: {:?}", sk_f);

        let ipe_enc_mat = ipe_enc_matrix_expression(&sk, &grp, true, &mut rand);
        println!("ipe_enc_mat: {}", ipe_enc_mat);

        let mut x1 = x.clone();
        x1.push(Integer::from(1));
        let mut ctxt = ipe_enc_mat * x1.clone();
        mod_vec(&mut ctxt, &grp.delta);
        println!("ctxt: {:?}", ctxt);

        let out = ipe_dec(&sk_f, &ctxt, &grp, true);
        println!("out: {:?}", out);
        
        // check inner product between x and y
        let mut out2 = Integer::from(0);
        for i in 0..x.len() {
            out2 += x[i].clone() * y[i].clone();
        }
        println!("out2: {:?}", out2);

        assert_eq!(out, out2);

    }
    
}
