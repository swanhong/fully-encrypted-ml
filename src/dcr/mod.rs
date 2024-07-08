pub mod scheme;

#[cfg(test)]
mod test {
    use crate::rug::Integer;
    use crate::util::group::Group;
    use crate::dcr::scheme::{dcr_setup, dcr_keygen, dcr_enc, dcr_dec};
    use std::time::{Instant, SystemTime};
    use crate::util::vector::{gen_random_vector_signed, int_mod};

    #[test]
    pub fn test_dcr_start_to_end() {
    println!(" == start test_dcr_start_to_end ==");
    
    let bit_len = 10;
    let dim: usize = 10;
    let bound = Integer::from(10000);

    println!("bit_len: {}", bit_len);
    println!("dim: {}", dim);
    println!("bound: {}", bound);
    println!(" !!! WARNING: this bound is not secure for real use case !!!");

    let grp = Group::new(bit_len);
    println!("grp: {}", grp);
    // define random Integer vector of dimension m
    let d = SystemTime::now()
    .duration_since(SystemTime::UNIX_EPOCH)
    .expect("Duration since UNIX_EPOCH failed");
    let mut rand = rug::rand::RandState::new();
    rand.seed(&Integer::from(d.as_secs()));
    let x = gen_random_vector_signed(dim, &bound, &mut rand);
    let y = gen_random_vector_signed(dim, &bound, &mut rand);
    println!("x: {:?}", x);
    println!("y: {:?}", y);
    

    let start = Instant::now();
    let (sk, pk) = dcr_setup(dim, &bound, &grp, &mut rand);
    let end = start.elapsed();
    println!("Time elapsed in dcr_setup is: {:?}", end);
    // println!("sk: {:?}", sk);
    // println!("pk: {:?}", pk);

    let start = Instant::now();
    let sk_y = dcr_keygen(&sk, &y);
    let end = start.elapsed();
    println!("Time elapsed in dcr_keygen is: {:?}", end);
    // println!("sk_y: {}", sk_y);
    
    let start = Instant::now();
    let ct_x = dcr_enc(&pk, &x, &grp, &mut rand);
    let end = start.elapsed();
    println!("Time elapsed in dcr_enc is: {:?}", end);
    // println!("ct_x: {:?}", ct_x);

    let start = Instant::now();
    let out = dcr_dec(&ct_x, &y, &sk_y, &grp);
    let end = start.elapsed();
    println!("Time elapsed in dcr_dec is: {:?}", end);
    println!("out: {}", out);   

    // check real inner product between x and y
    let mut real_ip = Integer::from(0);
    for i in 0..x.len() {
        real_ip += x[i].clone() * y[i].clone();
    }
    real_ip = int_mod(&mut real_ip, &grp.n);
    println!("real_ip: {}", real_ip);
    assert_eq!(out, real_ip);
    }
}