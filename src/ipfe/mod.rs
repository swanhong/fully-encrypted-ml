pub mod keys;
pub mod scheme;

#[cfg(test)]
mod tests {
    use crate::util::group::Group;

    use rug::Integer;
    use rug::rand::RandState;
    use crate::ipfe::scheme::{ipfe_setup, ipfe_keygen, ipfe_enc, ipfe_dec};
    use crate::util::vector::{gen_random_vector, int_mod};
    use std::time::SystemTime;
    
    #[test]
    fn test_ipfe_start_to_end() {
        let mut rand = RandState::new(); // Create a single RandState object
        let d = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .expect("Duration since UNIX_EPOCH failed");
        rand.seed(&Integer::from(d.as_secs()));
        
        let grp = Group::new(100);
        
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
        println!("start ipfe_setup");
        let start = SystemTime::now();
        let sk = ipfe_setup(&grp, dim, q, &mut rand);
        let end = start.elapsed();
        println!("Time elapsed in ipfe_setup is: {:?}", end);

        // Perform key generation
        println!("start ipfe_keygen");
        let start = SystemTime::now();
        let sk_f = ipfe_keygen(&sk, &y, &grp);
        let end = start.elapsed();
        println!("Time elapsed in ipfe_keygen is: {:?}", end);

        // Perform encryption
        println!("start ipfe_enc");
        let start = SystemTime::now();
        let ctxt = ipfe_enc(&sk, &x, &grp, true, &mut rand);
        let end = start.elapsed();
        println!("Time elapsed in enc (mat_mult) is: {:?}", end);        

        // Perform decryption
        println!("start ipfe_dec");
        let start = SystemTime::now();
        let out = ipfe_dec(&sk_f, &ctxt, &grp, true);
        let end = start.elapsed();
        println!("Time elapsed in ipfe_dec is: {:?}", end);
        
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
