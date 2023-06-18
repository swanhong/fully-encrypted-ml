pub mod scheme;

#[cfg(test)]
mod test {
    use crate::rug::Integer;
    use crate::util::group::Group;
    use crate::dcr_ipe::scheme::{dcr_setup, dcr_keygen, dcr_enc, dcr_dec};
    use std::time::{Instant, SystemTime};

    #[test]
    pub fn test_dcr_ipe_start_to_end() {
    println!(" == start test_dcr_ipe_start_to_end ==");
    
    let bit_len = 3072;
    let dim: usize = 10;
    let bound = Integer::from(100);

    println!("bit_len: {}", bit_len);
    println!("dim: {}", dim);
    println!("bound: {}", bound);

    let grp = Group::new(bit_len);

    // define random Integer vector of dimension m
    let d = SystemTime::now()
    .duration_since(SystemTime::UNIX_EPOCH)
    .expect("Duration since UNIX_EPOCH failed");
    let mut rand = rug::rand::RandState::new();
    rand.seed(&Integer::from(d.as_secs()));
    let mut x = vec![Integer::from(0); dim];
    let mut y = vec![Integer::from(0); dim];
    for i in 0..y.len() {
        x[i] = bound.clone().random_below(&mut rand);
        y[i] = bound.clone().random_below(&mut rand);
    }
    // println!("x: {:?}", x);
    // println!("y: {:?}", y);
    

    let start = Instant::now();
    let (sk, pk) = dcr_setup(dim, bound.clone(), &grp);
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
    let ct_x = dcr_enc(&pk, &x, &grp);
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
    println!("real_ip: {}", real_ip);
    assert_eq!(out, real_ip);
    }
}