extern crate rug;
use rug::rand::RandState;

use crate::util::{group::Group, matrix::*};
use crate::ipe::keys::IpeSk;

pub struct QeSk {
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
    pub fn new(grp: &Group, n_x: usize, n_y: usize, b: usize, rng: &mut RandState<'_>) -> QeSk {

        let modulo = grp.delta.clone();
        
        let v = Matrix::random(n_x, 2, &modulo, rng);
        let w = Matrix::random(n_y, 2, &modulo, rng);
              
        let i_n_x = Matrix::get_identity(n_x);
        let i_n_y = Matrix::get_identity(n_y);

        let v_t = v.transpose();
        let w_t = w.transpose();

        let m_up = Matrix::tensor_product(&v_t, &i_n_y, &modulo);
        let m_down = Matrix::tensor_product(&i_n_x, &w_t, &modulo);
        let m = concatenate_row(&m_up, &m_down);
        
        let (d_x, d_x_inv, d_x_null) = generate_right_inverse_space(n_x, n_x + b, &modulo, rng);
        let (d_y, d_y_inv, d_y_null) = generate_right_inverse_space(n_y, n_y + b, &modulo, rng);

        let ipe_sk = IpeSk::new(2 * (n_x + n_y), b, &grp, rng);

        QeSk {
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
