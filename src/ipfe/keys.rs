// keys.rs

extern crate rug;
use rug::Integer;
use rug::rand::RandState;

use crate::util::{group::Group, matrix::*};
use crate::util::vector::vec_mod;

pub struct IpfeSk {
    pub dim: usize,
    pub a: Vec<Integer>,
    pub u: Matrix,
    pub u_t: Matrix,
    pub d: Matrix,
    pub d_inv: Matrix,
    pub d_perp: Matrix,
    pub sk1: Vec<Integer>,
    pub sk_enc: Matrix,
}

impl IpfeSk {
    pub fn new(dim: usize, q: usize, grp: &Group, rng: &mut RandState<'_>) -> IpfeSk {
        // modulus = grp.delta
        // a = random vector of size 2 mod b
        // U = random matrix of size (dim) x 2
        // U_t = transpose of U
        // D = random matrix of size (dim + 2) x (dim + 2 + q)
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
        // D = random matrix of size (dim + 2) x (dim + 2 + q)
        // D_inv = right inverse of D, size (dim + 2 + q) x (dim + 2)
        // D_perp = null space basis of D, size (dim + 2 + q) x (dim + 2)
        let (d, d_inv, d_perp) = generate_right_inverse_space(dim + 2, dim + 2 + q, &grp.delta, rng);
        println!("input dim = {}, q = {}", dim, q);
        println!("d size = {} x {}", d.rows, d.cols);
        println!("d_inv size = {} x {}", d_inv.rows, d_inv.cols);
        println!("d_perp size = {} x {}", d_perp.rows, d_perp.cols);
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
        
        IpfeSk {
            dim,
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
impl fmt::Display for IpfeSk {
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
