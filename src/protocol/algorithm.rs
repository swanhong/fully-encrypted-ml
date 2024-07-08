#![allow(dead_code)]

use rug::Integer;
use rug::rand::RandState;
use crate::util::matrix::{Matrix, concatenate_diag_one, matrix_inverse};
use crate::util::vector::int_mod;

pub fn row_reduce_form(m: &mut Matrix, mod_val: &Integer) -> (Matrix, i32) {
    let n_rows = m.rows();
    let n_cols = m.cols();

    let mut row_ops_matrix = Matrix::get_identity(n_rows);

    let mut pivot_row = 0;
    let mut pivot_col = 0;
    let mut rank = 0;

    while pivot_row < n_rows && pivot_col < n_cols {
        println!("pivot_row, pivot_col = {}, {}", pivot_row, pivot_col);
        println!("m = \n{}", m);
        let mut pivot = m.get(pivot_row, pivot_col);

        while pivot == 0 {
            pivot_row += 1;

            if pivot_row >= n_rows {
                return (row_ops_matrix, -1);
            }

            pivot = m.get(pivot_row, pivot_col);
        }

        let pivot_inv = match pivot.clone().invert(&mod_val) {
            Ok(x) => x,
            Err(_e) => Integer::from(-1),
        };

        // println!("pviot = {}, pivot_inv = {}", pivot, pivot_inv);

        if pivot_inv != Integer::from(-1) {
            // Scale the pivot row
            for j in 0..n_cols {
                let val = int_mod(&(m.get(pivot_row, j) * pivot_inv.clone()), mod_val);
                m.set(pivot_row, j, val.clone());
            }
            let mut m_scalar = Matrix::get_identity(n_rows);
            m_scalar.set(pivot_row, pivot_row, pivot_inv.clone());
            row_ops_matrix = m_scalar.clone() * row_ops_matrix;
            row_ops_matrix.mod_inplace(mod_val);
        }

        let mut m_scalar = Matrix::get_identity(n_rows);
        for i in 0..n_rows {
            if i != pivot_row {
                let val1 = m.get(i, pivot_col);
                let ratio = val1.clone();

                // Update the current row
                for j in 0..n_cols {
                    let val2 = m.get(pivot_row, j);
                    let entry = int_mod(&(m.get(i, j) - (val2.clone() * ratio.clone())), mod_val);
                    // println!("entry for {} {}", i, j);
                    // println!("{} - {} * {} = {}", m.get(i, j), val2.clone(), ratio.clone(), entry);
                    m.set(i, j, entry);
                }

                // Update the row operations matrix
                m_scalar.set(i, pivot_row, ratio.clone() * -1);
                // for j in 0..n_rows {
                //     let val3 = row_ops_matrix.get(pivot_row, j);
                //     let entry = (val3 * ratio.clone() - row_ops_matrix.get(i, j)) % mod_val;
                //     row_ops_matrix.set(i, j, entry);
                // }
            }
        } 
        row_ops_matrix = m_scalar.clone() * row_ops_matrix;
        row_ops_matrix.mod_inplace(mod_val);

        pivot_row += 1;
        pivot_col += 1;
        rank += 1;
    }

    (row_ops_matrix.clone(), rank)
}

#[allow(dead_code)]
pub fn row_reduce_form_integer(m: &mut Matrix) -> (Matrix, i32) {
    let n_rows = m.rows();
    let n_cols = m.cols();

    let mut row_ops_matrix = Matrix::get_identity(n_rows);

    let mut pivot_row = 0;
    let mut pivot_col = 0;
    let mut rank = 0;
    while pivot_row < n_rows && pivot_col < n_cols {
        let pivot = m.get(pivot_row, pivot_col);

        if pivot == 0 {
            let mut i = 1;
            while pivot_row + i < n_cols {
                let val = m.get(pivot_row + i, pivot_col);
                if val != 0 {
                    // swap pivot_row and pviot_row + i -th rows
                    for j in 0..n_cols {
                        let val1 = m.get(pivot_row, j);
                        let val2 = m.get(pivot_row + i, j);
                        m.set(pivot_row, j, val2);
                        m.set(pivot_row + i, j, val1);
                    }
                    let mut m_scalar = Matrix::get_identity(n_rows);
                    m_scalar.set(pivot_row, pivot_row, Integer::from(0));
                    m_scalar.set(pivot_row + i, pivot_row + i, Integer::from(0));
                    m_scalar.set(pivot_row, pivot_row + i, Integer::from(1));
                    m_scalar.set(pivot_row + i, pivot_row, Integer::from(1));
                    row_ops_matrix = m_scalar.clone() * row_ops_matrix;
                    break;
                }
                i += 1;
            }
            return (row_ops_matrix, -1);
        } else {
            for i in 0..n_rows {
                if i != pivot_row {
                    let val1 = m.get(i, pivot_col);

                    for j in 0..n_cols {
                        let val2 = m.get(pivot_row, j);
                        let val3 = m.get(i, j);

                        m.set(i, j, val3 * pivot.clone() - val1.clone() * val2);
                        row_ops_matrix.set(i, j, row_ops_matrix.get(i, j) * pivot.clone() - row_ops_matrix.get(pivot_row, j) * val1.clone());
                    }
                }
            }

            pivot_row += 1;
        }

        pivot_col += 1;
        rank += 1;
    }

    (row_ops_matrix, rank)
}

pub fn sample_h(dim: usize, k: usize, modulo: &Integer, rng: &mut RandState<'_>) -> (Matrix, Matrix) {
    // sample two matrices h_left_1 and h_right_1
    
    // h_left is dim * (dim + k ) ternary matrix
    // h_right is (dim + k) * dim matrix satisfying h_left * h_right = identity_dim

    // h_right = h^t * (h * h^t)^-1 
    // h_left_1 = (h_left, 1)
    // h_right_1 = (h_right, 1)

    let mut h_0: Matrix; // dim * (dim + k)
    let mut h_t: Matrix; // (dim + k) * dim
    let h_0_inv: Matrix; // dim * dim
    let mut h_pr_0: Matrix; // (dim + k) * dim

    loop {
        // sample h_0 from {-1, 0, 1}^{dim * (dim+k)}
        h_0 = Matrix::random(dim, dim + k, &Integer::from(3), rng);
        h_0.add_int_inplace(&Integer::from(-1));
        // h_0.mod_inplace(modulo);

        h_t = h_0.transpose();
        let mut tmp = h_0.clone() * h_t.clone();
        tmp.mod_inplace(modulo);
        match matrix_inverse(&mut tmp, modulo) {
            Ok(m_inv) => {
                h_0_inv = m_inv.clone();
                break;
            }
            Err(_rank) => {
                continue;
            }
        }

    }
    
    h_pr_0 = h_t * h_0_inv;
    h_pr_0.mod_inplace(modulo);
    

    let h = concatenate_diag_one(&h_0);
    let h_pr = concatenate_diag_one(&h_pr_0);

    (h_pr.transpose(), h.transpose())
}

pub fn sample_h_trivial(dim: usize, k: usize) -> (Matrix, Matrix) {
    // sample two matrices h_left_1 and h_right_1
    
    // h_left is dim * (dim + k ) ternary matrix
    // h_right is (dim + k) * dim matrix satisfying h_left * h_right = identity_dim

    // h_right = h^t * (h * h^t)^-1 
    // h_left_1 = (h_left, 1)
    // h_right_1 = (h_right, 1)

    let mut h_0 = Matrix::new(dim, dim + k);
    let mut h_pr_0 = Matrix::new(dim + k, dim);

    for i in 0..dim {
        h_0.set(i, i, Integer::from(1));
        h_pr_0.set(i, i, Integer::from(1));
    }


    let h = concatenate_diag_one(&h_0);
    let h_pr = concatenate_diag_one(&h_pr_0);

    (h, h_pr)
}


pub fn sample_h_wo_padding(dim: usize, k: usize, modulo: &Integer, rng: &mut RandState<'_>) -> (Matrix, Matrix) {
    // sample two matrices h_left_1 and h_right_1
    
    // h_left is dim * (dim + k ) ternary matrix
    // h_right is (dim + k) * dim matrix satisfying h_left * h_right = identity_dim

    // h_right = h^t * (h * h^t)^-1 
    // h_left_1 = (h_left, 1)
    // h_right_1 = (h_right, 1)

    let mut h_0: Matrix; // dim * (dim + k)
    let mut h_t: Matrix; // (dim + k) * dim
    let h_0_inv: Matrix; // dim * dim
    let h_pr_0: Matrix; // (dim + k) * dim

    loop {
        // sample h_0 from {-1, 0, 1}^{dim * (dim+k)}
        h_0 = Matrix::random(dim, dim + k, &Integer::from(3), rng);
        h_0.add_int_inplace(&Integer::from(-1));
        // h_0.mod_inplace(modulo);

        h_t = h_0.transpose();
        let mut tmp = h_0.clone() * h_t.clone();
        tmp.mod_inplace(modulo);
        match matrix_inverse(&mut tmp, modulo) {
            Ok(m_inv) => {
                h_0_inv = m_inv.clone();
                break;
            }
            Err(_rank) => {
                continue;
            }
        }

    }
    
    h_pr_0 = h_t * h_0_inv;

    (h_0, h_pr_0)
}

pub fn sample_gamma(
    dim: usize,
    modulo: &Integer,
    rng: &mut RandState<'_>,
) -> (Matrix, Matrix) {
    // sample two matrices gamma_left_1 and gamma_right_1
    
    // gamma_left is dim * 2dim binary matrix
    // gamma_right is 2dim * dim matrix satisfying gamma_left * gamma_right = identity_dim

    // gamma_right = gamma^t * (gamma * gamma^t)^-1 
    // gamma_left_1 = (gamma_left, 1)
    // gamma_right_1 = (gamma_right, 1)
    let mut gamma_0: Matrix;
    let mut gamma_0_t: Matrix;
    let gamma_0_inv: Matrix;

    loop {
        gamma_0 = Matrix::random( 2 * dim, dim, &Integer::from(3), rng);
        for i in 0..gamma_0.rows() {
            for j in 0..gamma_0.cols() {
                let val = gamma_0.get(i, j);
                gamma_0.set(i, j, val - 1);
            }
        }
        gamma_0_t = gamma_0.transpose();
        let mut tmp = gamma_0_t.clone() * gamma_0.clone();
        tmp.mod_inplace(modulo);
        match matrix_inverse(&mut tmp, modulo) {
            Ok(m_inv) => {
                gamma_0_inv = m_inv.clone();
                break;
            }
            Err(_rank) => {
                continue;
            }
        }
    }
    

    let mut gamma_pr_0 = gamma_0_inv * gamma_0_t;
    gamma_pr_0.mod_inplace(modulo);

    let gamma_1 = concatenate_diag_one(&gamma_0);
    let gamma_pr_1 = concatenate_diag_one(&gamma_pr_0);

    (gamma_pr_1, gamma_1)
}

pub fn sample_gamma_trivial(
    dim: usize,
    k: usize,
) -> (Matrix, Matrix) {
    // sample two matrices gamma_left_1 and gamma_right_1
    
    // gamma_left is dim * (dim + k ) binary matrix
    // gamma_right is (dim + k) * dim matrix satisfying gamma_left * gamma_right = identity_dim

    // gamma_right = gamma^t * (gamma * gamma^t)^-1 
    // gamma_left_1 = (gamma_left, 1)
    // gamma_right_1 = (gamma_right, 1)
    let mut gamma_0 = Matrix::new(dim, dim + k);
    let mut gamma_pr_0 = Matrix::new(dim + k, dim);

    for i in 0..dim {
        gamma_0.set(i, i, Integer::from(1));
        gamma_pr_0.set(i, i, Integer::from(1));
    }

    let gamma_1 = concatenate_diag_one(&gamma_0);
    let gamma_pr_1 = concatenate_diag_one(&gamma_pr_0);

    (gamma_1, gamma_pr_1)
}


pub fn sample_gamma_wo_padding(
    dim: usize,
    k: usize,
    modulo: &Integer,
    rng: &mut RandState<'_>,
) -> (Matrix, Matrix) {
    // sample two matrices gamma_left_1 and gamma_right_1
    
    // gamma_left is dim * (dim + k ) binary matrix
    // gamma_right is (dim + k) * dim matrix satisfying gamma_left * gamma_right = identity_dim

    // gamma_right = gamma^t * (gamma * gamma^t)^-1 
    // gamma_left_1 = (gamma_left, 1)
    // gamma_right_1 = (gamma_right, 1)
    let mut gamma_0: Matrix;
    let mut gamma_0_t: Matrix;
    let gamma_0_inv: Matrix;

    loop {
        gamma_0 = Matrix::random( dim, dim + k, &Integer::from(2), rng);
        gamma_0_t = gamma_0.transpose();
        let mut tmp = gamma_0.clone() * gamma_0_t.clone();
        tmp.mod_inplace(modulo);
        match matrix_inverse(&mut tmp, modulo) {
            Ok(m_inv) => {
                gamma_0_inv = m_inv.clone();
                break;
            }
            Err(_rank) => {
                continue;
            }
        }
    }
    

    let mut gamma_pr_0 = gamma_0_t * gamma_0_inv;
    gamma_pr_0.mod_inplace(modulo);

    (gamma_0, gamma_pr_0)
}

use rug::ops::Pow;
use crate::util::group::Group;
pub fn get_sk_bound(
    dim: usize,
    bound: usize,
    lambda: usize,
    grp: &Group,
) -> Integer {
    // sk_bound >= 2^{lambda + dim + 1} * (b + dim*(sqrt(dim)*b)^dim)^{dim-1} * dim * N^2
    Integer::from(2).pow(lambda as u32 + dim as u32 + 1) * Integer::from(bound + dim*(dim as f64).sqrt() as usize * bound).pow(dim as u32 - 1) * Integer::from(dim) * grp.n_sq.clone()
}