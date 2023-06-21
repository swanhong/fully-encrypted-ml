extern crate rug;

mod util;
mod dcr_ipe;
mod ipe;
mod qe;
mod protocol;
mod test;

#[test]
fn test_generate_matrix_with_right_inverse() {
    // let row = 3;
    // let col = 4;
    // let modulo = &Integer::from(10);

    // let (a, a_inv, n) = generate_matrix_with_right_inverse(row, col, modulo);

    // // Add assertions to check the correctness of the generated matrices
    // assert_eq!(a.len(), row);
    // assert_eq!(a[0].len(), col);
    // // ...

    // assert_eq!(a_inv.len(), col);
    // assert_eq!(a_inv[0].len(), row);
    // // ...

    // assert_eq!(n.len(), row);
    // assert_eq!(n[0].len(), col);
    // // ...

    // Add more assertions as needed
}

fn main() {
    // let grp = Group::new(3072);

    // let dim: usize = 10;
    // let bound = Integer::from(100);


    // // define random Integer vector of dimension m
    // let mut rand = rug::rand::RandState::new();
    // let mut x = vec![Integer::from(0); dim];
    // let mut y = vec![Integer::from(0); dim];
    // for i in 0..y.len() {
    //     x[i] = bound.clone().random_below(&mut rand);
    //     y[i] = bound.clone().random_below(&mut rand);
    // }
    // println!("x: {:?}", x);
    // println!("y: {:?}", y);

    // let start = Instant::now();
    // let (sk, pk) = dcr_setup(dim, bound.clone(), &grp);
    // let end = start.elapsed();
    // println!("Time elapsed in dcr_setup is: {:?}", end);
    // // println!("sk: {:?}", sk);
    // // println!("pk: {:?}", pk);

    // let start = Instant::now();
    // let sk_y = dcr_keygen(&sk, &y, &grp);
    // let end = start.elapsed();
    // println!("Time elapsed in dcr_keygen is: {:?}", end);
    // // println!("sk_y: {}", sk_y);
    
    // let start = Instant::now();
    // let ct_x = dcr_enc(&pk, &x, &grp);
    // let end = start.elapsed();
    // println!("Time elapsed in dcr_enc is: {:?}", end);
    // // println!("ct_x: {:?}", ct_x);

    // let start = Instant::now();
    // let out = dcr_dec(&ct_x, &y, &sk_y, &grp);
    // let end = start.elapsed();
    // println!("Time elapsed in dcr_dec is: {:?}", end);
    // println!("out: {}", out);   

    // // check real inner product between x and y
    // let mut real_ip = Integer::from(0);
    // for i in 0..x.len() {
    //     real_ip += x[i].clone() * y[i].clone();
    // }
    // println!("real_ip: {}", real_ip);

    // // println!("Group values:");
    // // println!("n: {}", group.n);
    // // println!("p: {}", group.p);
    // // println!("q: {}", group.q);
    // // println!("n_sq: {}", group.n_sq);
    // // println!("n_root: {}", group.n_root);
    // // println!("g: {}", group.g);
    // // println!("mu: {}", group.mu);
    // // println!("phi_n: {}", group.phi_n);
    // // println!("delta: {}", group.delta);
}
