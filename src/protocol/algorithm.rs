use rug::Integer;
use rug::rand::RandState;
use crate::util::matrix::{Matrix, concatenate_diag_one, concatenate_col};
use crate::util::vector::{int_mod};

pub fn echelon_form(
    m: &mut Matrix,
    mod_val: &Integer,
) -> (Vec<usize>, Vec<usize>, i32) {
    let n_rows = m.rows();
    let n_cols = m.cols();

    let mut pivot_cols: Vec<usize> = Vec::new();
    let mut free_vars: Vec<usize> = Vec::new();
    let mut rank = -1;
    let mut pivot_row = 0;
    let mut pivot_col = 0;

    while pivot_row < n_rows && pivot_col < n_cols {
        let mut pivot = m.get(pivot_row, pivot_col);
        
        while pivot == 0 {
            // if pivot_col < n_cols - 1 {
            //     free_vars.push(pivot_col);
            //     pivot_col += 1;
            // } else {
            //     free_vars.push(pivot_col);
            //     pivot_row += 1;
            //     pivot_col = 0;
            // }
            pivot_row += 1;

            if pivot_row >= n_rows {
                return (pivot_cols, free_vars, -1);
            }

            pivot = m.get(pivot_row, pivot_col);
        }

        rank += 1;
        let pivot = m.get(pivot_row, pivot_col);
        let pivot_inv = match pivot.clone().invert(&mod_val) {
            Ok(x) => x,
            Err(_e) => Integer::from(-1)
        };

        if pivot_inv != Integer::from(-1) {
            for j in pivot_col..n_cols {
                let val = m.get(pivot_row, j) * &pivot_inv % mod_val;
                m.set(pivot_row, j, val);
            }
        } else {
            return (pivot_cols, free_vars, -1);
        }

        for i in 0..n_rows {
            let pivot = m.get(pivot_row, pivot_col);
            let mut ratio = Integer::from(0);
            let pivot_inv = match pivot.clone().invert(&mod_val) {
                Ok(x) => x,
                Err(_e) => Integer::from(-1)
            };
            if pivot_inv != Integer::from(-1) {
                let val1 = m.get(i, pivot_col);
                ratio = val1 * pivot_inv % mod_val;
            } else {
                return (pivot_cols, free_vars, -1);
            }

            for j in 0..n_cols {
                if i != pivot_row {
                    let val = m.get(pivot_row, j) * ratio.clone();
                    let mut entry = m.get(i, j);
                    entry = (entry - val) % mod_val;
                    m.set(i, j, entry);
                }
            }
        }

        pivot_cols.push(pivot_col);
        pivot_row += 1;
        pivot_col += 1;
    }
    if pivot_col < n_cols {
        for j in pivot_col..n_cols {
            free_vars.push(j);
        }
    }  
    m.mod_inplace(mod_val);

    (pivot_cols, free_vars, rank)
}

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


pub fn matrix_inverse(
    m: &mut Matrix,
    mod_val: &Integer,
) -> Result<Matrix, i32> {
    assert_eq!(m.rows, m.cols);
    let n = m.rows;
    let mut m_inv = Matrix::get_identity(n);
    let mut m_aug = concatenate_col(m, &m_inv);
    let (pivot_cols, free_vars, r) = echelon_form(&mut m_aug, mod_val);
    
    if r != -1 {
        for i in 0..n {
            for j in 0..n {
                m_inv.set(i, j, m_aug.get(i, j + n));
            }
        }
    }
    
    match r {
        -1 => Err(-1),
        _ => Ok(m_inv)
    }
}

pub fn sample_h(dim: usize, k: usize, t: Integer, modulo: &Integer, rng: &mut RandState<'_>) -> (Matrix, Matrix) {
    // sample two matrices h_left_1 and h_right_1
    
    // h_left is dim * (dim + k ) ternary matrix
    // h_right is (dim + k) * dim matrix satisfying h_left * h_right = identity_dim

    // h_right = h^t * (h * h^t)^-1 
    // h_left_1 = (h_left, 1)
    // h_right_1 = (h_right, 1)

    let mut h_0: Matrix = Matrix::new(dim, dim + k);
    let mut h_t: Matrix = Matrix::new(dim + k, dim);
    let mut h_0_inv: Matrix = Matrix::new(dim, dim);
    let mut h_pr_0: Matrix = Matrix::new(dim + k, dim);

    loop {
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
                Err(rank) => {
                    continue;
                }
            }
    
        }
        
        h_pr_0 = h_t * h_0_inv;
        let mut is_h_pr_0_in_hyperball = true;
        for i in 0..dim {
            let col_i = h_pr_0.get_col(i);
            // compute norm of col_i
            let mut norm = Integer::from(0);
            for j in 0..dim + k {
                norm += col_i[j].clone().pow(2);
            }
            norm = norm.sqrt();
            if norm > t {
                is_h_pr_0_in_hyperball = false;
            }
        }

        if is_h_pr_0_in_hyperball {
            break;
        }
    }
    

    let h = concatenate_diag_one(&h_0);
    let h_pr = concatenate_diag_one(&h_pr_0);

    (h, h_pr)
}

pub fn sample_gamma(
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
    let mut gamma_0 = Matrix::new(1, 1);
    let mut gamma_0_t = Matrix::new(1, 1);
    let mut gamma_0_inv = Matrix::new(1, 1);

    while true {
        gamma_0 = Matrix::random( dim, dim + k, &Integer::from(2), rng);
        gamma_0_t = gamma_0.transpose();
        let mut tmp = gamma_0.clone() * gamma_0_t.clone();
        tmp.mod_inplace(modulo);
        match matrix_inverse(&mut tmp, modulo) {
            Ok(m_inv) => {
                gamma_0_inv = m_inv.clone();
                break;
            }
            Err(rank) => {
                continue;
            }
        }
    }
    

    let mut gamma_pr_0 = gamma_0_t * gamma_0_inv;
    gamma_pr_0.mod_inplace(modulo);

    let gamma_1 = concatenate_diag_one(&gamma_0);
    let gamma_pr_1 = concatenate_diag_one(&gamma_pr_0);

    (gamma_1, gamma_pr_1)
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

pub fn get_smallest_t(
    dim: usize,
    k: usize,
    lambda: usize,
) -> Integer {
    // t >= sqrt(dim * (dim + k)) * 2^{(lambda / dim^2)}
    Integer::from((dim * (dim + k)) as u32 + 1).sqrt() * Integer::from(2).pow((lambda / (dim * dim)) as u32 + 1)
}