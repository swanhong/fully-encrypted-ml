pub mod keys;
pub mod scheme;

mod tests {
    #![allow(unused_imports)]
    use crate::util::group::Group;

    use rug::Integer;
    use rug::rand::RandState;
    use std::time::{Duration, SystemTime};
    use crate::qfe::scheme::{qfe_setup, qfe_keygen, qfe_enc, qfe_enc_matrix_same_xy, qfe_dec, qfe_enc_matrix_expression};
    use crate::util::vector::{gen_random_vector, gen_random_vector_signed, vec_mod, vec_mul_scalar, int_mod};

    #[test]
    fn test_qfe_start_to_end() {
        let mut rand = RandState::new(); // Create a single RandState object
        let d = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .expect("Duration since UNIX_EPOCH failed");
        rand.seed(&Integer::from(d.as_secs()));
        
        let grp = Group::new(100); // Initialize the group        
        // println!("{}", grp);
        
        let dim = 6;
        let q = 2 * dim + 1;
        let bound = 1000000;

        let x = gen_random_vector_signed(dim, &Integer::from(2 * bound), &mut rand);
        let f = gen_random_vector_signed(dim * dim, &Integer::from(2 * bound), &mut rand);

        // compute out2 = inner product <f, x tensor y>
        let mut out2 = Integer::from(0);
        for i in 0..dim {
            for j in 0..dim {
                out2 += f[i * dim + j].clone() * x[i].clone() * x[j].clone();
            }
        }

        // Perform setup
        let start = SystemTime::now();
        let qfe_sk = qfe_setup(&grp, dim, q, &mut rand);
        let end = start.elapsed();
        println!("Time elapsed in qfe_setup is: {:?}", end);

        println!("qfe_sk: {}", qfe_sk);

        // Perform key generation
        let start = SystemTime::now();
        let fk = qfe_keygen(&qfe_sk, &f, &grp);
        let end = start.elapsed();
        println!("Time elapsed in qfe_keygen is: {:?}", end);
        
        let ctxt = qfe_enc(&qfe_sk, &x, &grp, &mut rand);
        println!("ctxt_size = {}", ctxt.len());
        // let mut ctxt_x: Vec<Integer>;
        // let mut ctxt_y: Vec<Integer>;
        // let mut ctxt_h: Vec<Integer>;
        // // Perform encryption
        // let start = SystemTime::now();
        // let (enc_mat_x, enc_mat_y, enc_mat_h) = qfe_enc_matrix_same_xy(&qfe_sk, dim, &grp, &mut rand);
        // let end = start.elapsed();
        // println!("Time elapsed in qfe_enc_matrix is: {:?}", end);

        // let mut x1 = x.clone();
        // x1.push(Integer::from(1));

        // ctxt_x = enc_mat_x.mul_vec(&x1);
        // vec_mod(&mut ctxt_x, &grp.delta);
        // ctxt_y = enc_mat_y.mul_vec(&x1);
        // vec_mod(&mut ctxt_y, &grp.delta);
        // ctxt_h = enc_mat_h.mul_vec(&x1);
        // vec_mod(&mut ctxt_h, &grp.delta);
        
        let start = SystemTime::now();
        let out = qfe_dec(&fk, &ctxt, dim, q, &grp);
        let end = start.elapsed();
        println!("Time elapsed in qfe_dec is: {:?}", end);

        println!("x = {:?}", x);
        println!("f = {:?}", f);
        println!("eval: {}", out);
        println!("real: {}", out2);
        let out_mod = int_mod(&out, &grp.n.clone());
        let out2_mod = int_mod(&out2, &grp.n.clone());
        println!("eval mod n: {}", out_mod);
        println!("real mod n: {}", out2_mod);
        assert_eq!(out_mod, out2_mod);

    }
}