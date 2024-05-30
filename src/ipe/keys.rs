// keys.rs

extern crate rug;
use rug::Integer;
use rug::rand::RandState;

use crate::util::{group::Group, matrix::*};
use crate::util::vector::vec_mod;

pub struct IpeSk {
    pub a: Vec<Integer>,
    pub u: Matrix,
    pub u_t: Matrix,
    pub d: Matrix,
    pub d_inv: Matrix,
    pub d_perp: Matrix,
    pub sk1: Vec<Integer>,
    pub sk_enc: Matrix,
}

impl IpeSk {
    pub fn new(dim: usize, b: usize, grp: &Group, rng: &mut RandState<'_>) -> IpeSk {
        // modulus = grp.delta
        // a = random vector of size 2 mod b
        // U = random matrix of size (4*dim) x 2
        // U_t = transpose of U
        // D = random matrix of size (4*dim + 2) x (4*dim + 2 + q)
        // D_inv = right inverse of D
        // D_perp = null space basis of D
        // sk1 = (D_inv_left + D_inv_right * U) * a
        // sk_enc = D_inv_right

        let mut a: Vec<Integer> = Vec::new();
        for _ in 0..2 {
            a.push(grp.delta.clone().random_below(rng));
        }

        let u = Matrix::random(dim, 2, &grp.delta, rng);
        let u_t = u.transpose();
        let (d, d_inv, d_perp) = generate_right_inverse_space(dim + 2, 4*dim + 2 + b, &grp.delta, rng);
        
        // sk1 = (D_inv_left + D_inv_right * U) * a
        let mut d_inv_left = Matrix::new(d_inv.rows, 2);
        let mut sk_enc = Matrix::new(d_inv.rows, dim);
        for i in 0..d_inv.rows {
            for j in 0..d_inv.cols {
                if j < 2 {
                    d_inv_left.set(i, j, d_inv.get(i, j));
                } else {
                    sk_enc.set(i, j - 2, d_inv.get(i, j));
                }
            }
        }
        let mut sk1 = (d_inv_left + &(sk_enc.clone() * u.clone())).mul_vec(&a);
        vec_mod(&mut sk1, &grp.delta);
        
        IpeSk {
            a,
            u,
            u_t,
            d,
            d_inv,
            d_perp,
            sk1,
            sk_enc
        }
    }
}

use std::fmt;
impl fmt::Display for IpeSk {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "a: {:?}", self.a)?;
        writeln!(f, "u: {:?}", self.u)?;
        writeln!(f, "u_t: {:?}", self.u_t)?;
        writeln!(f, "d: {:?}", self.d)?;
        writeln!(f, "d_inv: {:?}", self.d_inv)?;
        writeln!(f, "d_perp: {:?}", self.d_perp)?;
        writeln!(f, "sk1: {:?}", self.sk1)?;
        writeln!(f, "sk_enc: {:?}", self.sk_enc)?;
        Ok(())
    }
}
