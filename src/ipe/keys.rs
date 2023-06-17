// keys.rs

extern crate rug;
use rug::Integer;
use rug::rand::RandState;
use std::time::{SystemTime, UNIX_EPOCH};

use crate::util::{group::Group, matrix::*};

pub struct IPE_sk {
    pub a: Vec<Integer>,
    pub u: Matrix,
    pub u_t: Matrix,
    pub d: Matrix,
    pub d_inv: Matrix,
    pub d_perp: Matrix,
    pub sk1: Vec<Integer>,
    pub sk_enc: Matrix,
}

impl IPE_sk {
    pub fn new(dim: usize, b: usize, grp: &Group) -> IPE_sk {
        // modulus = grp.delta
        // a = random vector of size 2
        // U = random matrix of size dim x 2
        // U_t = transpose of U
        // D = random matrix of size (dim + 2) x (dim + 2 + B)
        // D_inv = right inverse of D
        // D_perp = null space basis of D
        // sk1 = (D_inv_left + D_inv_right * U) * a
        // sk_enc = D_inv_right

        let mut rng = RandState::new();
        let mut rng = RandState::new();
        let d = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("Duration since UNIX_EPOCH failed")
            .as_secs();
        rng.seed(&Integer::from(d));
        
        let mut a: Vec<Integer> = Vec::new();
        for _ in 0..2 {
            a.push(grp.delta.clone().random_below(&mut rng));
        }
        
        let mut u = Matrix::random(dim, 2, &grp.delta);
        let mut u_t = transpose(&u);
        let (d, d_inv, d_perp) = generate_right_inverse_space(dim + 2, dim + 2 + b, &grp.delta);

        // D_inv = (D_inv_left (col=2) || sk_enc (col=dim))

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
        let sk1 = (d_inv_left + &(sk_enc.clone() * u.clone())).mul_vec(&a);
        
        IPE_sk {
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
impl fmt::Display for IPE_sk {
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
