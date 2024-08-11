#![allow(dead_code)]

use rug::Integer;
use rug::rand::RandState;
use rand::seq::SliceRandom;
use crate::util::matrix::{concatenate_col, concatenate_diag_one, concatenate_row, matrix_inverse, Matrix};
use crate::util::vector::{int_mod, vec_mul_scalar};

pub fn row_reduce_form(m: &mut Matrix, mod_val: &Integer) -> (Matrix, Vec<usize>, i32) {
    let n_rows = m.rows();
    let n_cols = m.cols();

    let mut row_ops_matrix = Matrix::get_identity(n_rows);

    let mut pivot_row = 0;
    let mut pivot_col = 0;
    let mut rank = 0;

    let mut pivot_col_vec = Vec::new();

    while pivot_row < n_rows && pivot_col < n_cols {
        let mut pivot = m.get(pivot_row, pivot_col);

        while pivot == 0 {
            pivot_row += 1;

            if pivot_row >= n_rows {
                return (row_ops_matrix, pivot_col_vec, -1);
            }

            pivot = m.get(pivot_row, pivot_col);
        }

        let pivot_inv = match pivot.clone().invert(&mod_val) {
            Ok(x) => x,
            Err(_e) => Integer::from(-1),
        };

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
                    m.set(i, j, entry);
                }

                // Update the row operations matrix
                m_scalar.set(i, pivot_row, ratio.clone() * -1);
            }
        } 
        row_ops_matrix = m_scalar.clone() * row_ops_matrix;
        row_ops_matrix.mod_inplace(mod_val);

        pivot_col_vec.push(pivot_row);
        pivot_row += 1;
        pivot_col += 1;
        rank += 1;
    }

    (row_ops_matrix.clone(), pivot_col_vec, rank)
}

#[allow(dead_code)]
pub fn row_reduce_form_integer(m: &mut Matrix) -> (Matrix, Vec<usize>, i32) {
    let n_rows = m.rows();
    let n_cols = m.cols();

    let mut row_ops_matrix = Matrix::get_identity(n_rows);

    let mut pivot_row = 0;
    let mut pivot_col = 0;
    let mut rank = 0;
    let mut pivot_vec = Vec::new();
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
            return (row_ops_matrix, pivot_vec, -1);
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

        pivot_vec.push(pivot_col);
        pivot_col += 1;
        rank += 1;
    }

    (row_ops_matrix, pivot_vec, rank)
}

pub fn sample_h(dim: usize, k: usize, bound: &Integer, modulo: &Integer, rng: &mut RandState<'_>) -> (Matrix, Matrix) {
    // mat_u: k x (dim + k) matrix: id_k || 0_{dim x k}
    // let mut mat_u = Matrix::get_identity(k);
    // let mat_v = Matrix::new(k, dim);
    // mat_u = concatenate_col(&mat_u, &mat_v);

    // println!("mat_u : {}", mat_u);
    let mut mat_u: Matrix;
    let mut row_ops: Matrix;
    let mut pivot_vec: Vec<usize>;
    loop {
        // create random permutation vector of size (dim + k)
        let mut perm = (0..dim+k).collect::<Vec<usize>>();
        perm.shuffle(&mut rand::thread_rng());
        // println!("perm : {:?}", perm);

        // create permutation matrix mat_u
        // entry is +-1 in random
        mat_u = Matrix::new(k, dim + k);
        for j in 0..k {
            let val = 2 * Integer::from(2).random_below(rng) - 1;
            mat_u.set(j, perm[j], val);
        }
        
        // mat_u = Matrix::random(k, dim + k, bound, rng);
        println!("mat_u : {}", mat_u);
        (row_ops, pivot_vec, _) = row_reduce_form(&mut mat_u.transpose(), modulo);
        println!("pivot_vec.len() = {}", pivot_vec.len());
        if pivot_vec.len() == k {
            break;
        }
    }
    println!("row_ops : {}", row_ops);
    println!("pivot_vec : {:?}", pivot_vec);
    println!("mat_u = \n{}", mat_u);
    // Compute the left null space of mat_u
    let mut null_space = Matrix::new(mat_u.cols(), mat_u.cols() - mat_u.rows());

    let mut free_var_index = 0;
    for col in 0..mat_u.cols() {
        if !pivot_vec.contains(&col) {
            // Set the free variable column to the corresponding basis vector
            null_space.set(col, free_var_index, Integer::from(1));

            // Update the null space entries using the row operations
            for idx in 0..pivot_vec.len() {
                let pivot = pivot_vec[idx];
                if pivot >= col {
                    break;
                }
                println!("pivot = {}, col = {}", pivot, col);
                let val = mat_u.get(idx, col);
                let val_neg = val.clone() * Integer::from(-1);
                null_space.set(pivot, free_var_index, val_neg);
            }

            free_var_index += 1;
        }
        // println!("col = {}, free_var_index = {}", col, free_var_index);
        // println!("null_space : {}", null_space);
    }

    println!("null_space : \n{}", null_space);

    println!("test null space");
    println!("mat_u * null_space = {}", mat_u.clone() * null_space.clone());

    let mat_u = mat_u.transpose();
    let mat_v = null_space.transpose();

    println!("test v * u");
    let mut mul = mat_v.clone() * mat_u.clone();
    mul.mod_inplace(modulo);
    println!("v * u = {}", mul);

    // fix lambda = 128
    // B_h >= 2^{lambda/(dim * k) - 1}
    let bound_h = Integer::from(2).pow(128 / (dim * k) as u32 - 1);
    println!("bound_h = {}", bound_h);
    let mat_t: Matrix;
    let mut mat_h: Matrix;
    loop {
        mat_h = Matrix::random_signed(dim+k, dim, &bound_h, rng);
        let mul = mat_v.clone() * mat_h.clone();
        match matrix_inverse(&mut mul.clone(), modulo) {
            Ok(inv) => {
                mat_t = inv.clone();
                break;
            }
            Err(_rank) => {
                continue;
            }
        }
    }

    let mat_h_pr = mat_t * mat_v.clone();

    let mat_h = concatenate_diag_one(&mat_h);
    let mat_h_pr = concatenate_diag_one(&mat_h_pr);

    (mat_h_pr, mat_h)
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