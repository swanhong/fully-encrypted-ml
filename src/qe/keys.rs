extern crate rug;
use rug::rand::RandState;

use crate::util::{group::Group, matrix::*};
use crate::ipe::keys::IpeSk;

pub struct QeSk {
    pub dim: usize,
    pub q: usize,
    pub v: Matrix,
    pub w: Matrix,
    pub m: Matrix,
    pub d_x: Matrix,
    pub d_y: Matrix,
    pub d_x_null: Matrix,
    pub d_x_inv: Matrix,
    pub d_y_null: Matrix,
    pub d_y_inv: Matrix,
    pub ipe_sk: IpeSk,
}


impl QeSk {
    pub fn new(grp: &Group, dim: usize, q: usize, rng: &mut RandState<'_>) -> QeSk {
        let modulo = grp.delta.clone();
        
        let v = Matrix::random(dim, 2, &modulo, rng);
        let w = Matrix::random(dim, 2, &modulo, rng);

        let i_dim = Matrix::get_identity(dim);

        let v_t = v.transpose();
        let w_t = w.transpose();

        let m_up = Matrix::tensor_product(&v_t, &i_dim, &modulo);
        let m_down = Matrix::tensor_product(&i_dim, &w_t, &modulo);
        let m = concatenate_row(&m_up, &m_down);
        
        let (d_x, d_x_inv, d_x_null) = generate_right_inverse_space(dim, dim + q, &modulo, rng);
        let (d_y, d_y_inv, d_y_null) = generate_right_inverse_space(dim, dim + q, &modulo, rng);

        let ipe_sk = IpeSk::new(4 * dim, q, &grp, rng);

        QeSk {
            dim,
            q,
            v,
            w,
            m,
            d_x,
            d_y,
            d_x_null,
            d_x_inv,
            d_y_null,
            d_y_inv,
            ipe_sk,
        }
    }
}

use std::fmt;

impl fmt::Display for QeSk {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "QeSk:")?;
        writeln!(f, "dim: {}", self.dim)?;
        writeln!(f, "q: {}", self.q)?;
        writeln!(f, "v: \n{}", self.v)?;
        writeln!(f, "w: \n{}", self.w)?;
        writeln!(f, "m: \n{}", self.m)?;
        writeln!(f, "d_x: \n{}", self.d_x)?;
        writeln!(f, "d_y: \n{}", self.d_y)?;
        writeln!(f, "d_x_null: \n{}", self.d_x_null)?;
        writeln!(f, "d_x_inv: \n{}", self.d_x_inv)?;
        writeln!(f, "d_y_null: \n{}", self.d_y_null)?;
        writeln!(f, "d_y_inv: \n{}", self.d_y_inv)?;
        writeln!(f, "ipe_sk: \n{}", self.ipe_sk)?;
        Ok(())
    }
}
