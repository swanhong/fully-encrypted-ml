// dcr_ipe.rs
extern crate rug;
use rug::Integer;
use rug::rand::RandState;

use std::time::SystemTime;

use crate::util::group::{Group};

// takes input of dimension, secret key bound, and group
// returns a key pair (sk, pk)
pub fn dcr_setup(
    dim: usize, 
    sk_bound: &Integer, 
    grp: &Group,
    rng: &mut RandState<'_>,
) -> (Vec<Integer>, Vec<Integer>) {
    let mut sk = vec![Integer::from(0); dim];
    let mut pk = vec![Integer::from(0); dim];
        
    for i in 0..dim {
        sk[i] = sk_bound.clone().random_below(rng);
    }

    for i in 0..dim {
        pk[i] = grp.g.clone().pow_mod(&sk[i], &grp.n_sq).unwrap();
    }
    
    (sk, pk)
}

pub fn dcr_keygen(sk: &Vec<Integer>, y: &Vec<Integer>) -> Integer {
    // sk_y = inner product between sk and y
    assert_eq!(sk.len(), y.len(), "sk and y must have the same length");
    let mut sk_y = Integer::from(0);
    for i in 0..sk.len() {
        sk_y += sk[i].clone() * y[i].clone();
    }
    sk_y
}

pub fn dcr_enc(pk: &Vec<Integer>, x: &Vec<Integer>, grp: &Group) -> Vec<Integer> {
    assert!(pk.len() == x.len());
    
    let mut ct_x = vec![Integer::from(0); x.len() + 1];
    
    let mut rand = RandState::new(); // Create a single RandState object
    let d = SystemTime::now()
    .duration_since(SystemTime::UNIX_EPOCH)
    .expect("Duration since UNIX_EPOCH failed");
    rand.seed(&Integer::from(d.as_secs()));

    let r: Integer = grp.n_sq.clone().random_below(&mut rand);
    ct_x[x.len()] = grp.g.clone().pow_mod(&r, &grp.n_sq).unwrap();
    for i in 0..x.len() {
        let encode = &x[i] * &grp.n + Integer::from(1);
        ct_x[i] = pk[i].clone().pow_mod(&r, &grp.n_sq).unwrap() * encode;
        ct_x[i] = ct_x[i].clone().div_rem_euc(grp.n_sq.clone()).1;
    }
    ct_x
}

pub fn dcr_dec(ct_x: &Vec<Integer>, y: &Vec<Integer>, sk_y: &Integer, grp: &Group) -> Integer {
    assert!(ct_x.len() == y.len() + 1);

    let dim = y.len();
    let mut out = Integer::from(1);
    
    for i in 0..ct_x.len()-1 {
        out *= ct_x[i].clone().pow_mod(&y[i], &grp.n_sq).unwrap();
    }
    let num = 
            ct_x[dim].clone().pow_mod(&sk_y, &grp.n_sq).unwrap()
            .invert(&grp.n_sq).unwrap();
    out *= num;
    out = out.div_rem_euc(grp.n_sq.clone()).1;
    // out = (out - 1) / &grp.n;
    out = out - 1;
    let (quo, rem) = out.div_rem(grp.n.clone());
    assert!(rem == Integer::from(0), "Remainder is not 0");
    quo
}