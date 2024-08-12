#![allow(dead_code)]

use rug::{Integer, rand::RandState};
use std::time::{SystemTime, UNIX_EPOCH};
use std::str::FromStr;
use crate::util::vector::{gen_random_vector, tensor_product_vecs, int_mod};

#[derive(Debug)]
pub struct Matrix {
    data: Vec<Integer>,
    pub rows: usize,
    pub cols: usize,
}

impl Matrix {
    pub fn new(rows: usize, cols: usize) -> Matrix {
        Matrix {
            data: vec![Integer::from(0); rows * cols],
            rows: rows,
            cols: cols,
        }
    }

    pub fn parse_matrix(input: &str) -> Matrix {
        let lines: Vec<&str> = input.trim().split('\n').collect();
        let rows = lines.len() - 2;
        let cols = lines[1]
            .trim_matches(|c: char| !c.is_digit(10))
            .split_whitespace()
            .count();
    
        let mut matrix = Matrix::new(rows, cols);
        println!("row, col = {}, {}", rows, cols);
    
        for (i, line) in lines.iter().enumerate() {
            let values: Vec<Integer> = line
                .trim_matches(|c: char| !c.is_digit(10))
                .split_whitespace()
                .map(|value| {
                    Integer::from_str(value).unwrap_or_else(|_| {
                        panic!("Failed to parse integer value: '{}'", value);
                    })
                })
                .collect();
            if i > 0 && i < rows + 1 {
                matrix.set_row(i-1, &values);
            }
        }
    
        matrix
    }

    // random matrix with size rows x cols and elements in 0..bound
    pub fn random(rows: usize, cols: usize, bound: &Integer, rng: &mut RandState<'_>) -> Matrix {
        let mut matrix = Matrix::new(rows, cols);

        for i in 0..rows {
            for j in 0..cols {
                matrix.set(i, j, bound.clone().random_below(rng));
            }
        }
        matrix
    }

    pub fn random_signed(rows: usize, cols: usize, bound: &Integer, rng: &mut RandState<'_>) -> Matrix {
        let mut matrix = Matrix::new(rows, cols);
        let bound2: Integer = bound.clone() * 2 + 1;
        for i in 0..rows {
            for j in 0..cols {
                let val = bound2.clone().random_below(rng) - bound.clone();
                matrix.set(i, j, val);
            }
        }

        matrix
    }

    pub fn gen_tensored_matrix(&self, modulo: &Integer) -> Matrix {
        // input dim = self.rows
        // output dim = self.cols
        let mut matrix = Matrix::new(self.rows, self.cols * self.cols);
        for i in 0..self.rows {
            let f = self.get_row(i);
            let f_tensor = tensor_product_vecs(&f, &f, &modulo);
            matrix.set_row(i, &f_tensor);
        }
        matrix
    }

    pub fn random_quadratic_tensored(
        dim_input: usize,
        dim_output: usize,
        bound: &Integer,
        modulo: &Integer,
        rng: &mut RandState<'_>,
    ) -> Matrix {
        let mut matrix = Matrix::new(dim_output, dim_input * dim_input);

        for i in 0..dim_output {
            let f = gen_random_vector(dim_input, bound, rng);
            let f_tensor = tensor_product_vecs(&f, &f, &modulo);
            matrix.set_row(i, &f_tensor);
        }
        matrix
    }

    pub fn random_quadratic_tensored_with_one_padded(
        dim_input: usize,
        dim_output: usize,
        bound: &Integer,
        modulo: &Integer,
        rng: &mut RandState<'_>,
    ) -> Matrix {
        let mut matrix = Matrix::new(dim_output, (dim_input + 1) * (dim_input + 1));

        for i in 0..dim_output {
            let mut f = gen_random_vector(dim_input, bound, rng);
            f.push(Integer::from(1));
            let f_tensor = tensor_product_vecs(&f, &f, &modulo);
            matrix.set_row(i, &f_tensor);
        }
        matrix
    }

    pub fn clone(&self) -> Matrix {
        let mut matrix = Matrix::new(self.rows, self.cols);
        for i in 0..self.rows {
            for j in 0..self.cols {
                matrix.set(i, j, self.get(i, j));
            }
        }
        matrix
    }


    pub fn copy(&self) -> Matrix {
        let mut matrix = Matrix::new(self.rows, self.cols);
        for i in 0..self.rows {
            for j in 0..self.cols {
                matrix.set(i, j, self.get(i, j));
            }
        }
        matrix
    }

    pub fn get(&self, row: usize, col: usize) -> Integer {
        assert!(row < self.rows);
        assert!(col < self.cols);
        self.data[row * self.cols + col].clone()
    }

    pub fn set(&mut self, row: usize, col: usize, value: Integer) {
        assert!(row < self.rows);
        assert!(col < self.cols);
        self.data[row * self.cols + col] = value;
    }

    pub fn set_row(&mut self, row: usize, values: &Vec<Integer>) {
        assert!(row < self.rows);
        assert!(values.len() == self.cols);
        for j in 0..self.cols {
            self.set(row, j, values[j].clone());
        }
    }

    pub fn set_col(&mut self, col: usize, values: &Vec<Integer>) {
        assert!(col < self.cols);
        assert!(values.len() == self.rows);
        for i in 0..self.rows {
            self.set(i, col, values[i].clone());
        }
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }

    pub fn get_row(&self, row: usize) -> Vec<Integer> {
        assert!(row < self.rows);
        let mut result = Vec::new();
        for j in 0..self.cols {
            result.push(self.get(row, j));
        }
        result
    }

    pub fn get_col(&self, col: usize) -> Vec<Integer> {
        assert!(col < self.cols);
        let mut result = Vec::new();
        for i in 0..self.rows {
            result.push(self.get(i, col));
        }
        result
    }

    pub fn negate_inplace(&mut self) {
        for i in 0..self.rows {
            for j in 0..self.cols {
                self.data[i * self.cols + j] = -self.data[i * self.cols + j].clone();
            }
        }
    }

    pub fn negate(&self) -> Matrix {
        let mut matrix = Matrix::new(self.rows, self.cols);
        for i in 0..self.rows {
            for j in 0..self.cols {
                matrix.set(i, j, -self.get(i, j));
            }
        }
        matrix
    }

    // add
    pub fn add_inplace(&mut self, other: &Matrix) {
        assert!(self.rows == other.rows);
        assert!(self.cols == other.cols);
        for i in 0..self.rows {
            for j in 0..self.cols {
                self.data[i * self.cols + j] += other.data[i * self.cols + j].clone();
            }
        }
    }

    pub fn sub_inplace(&mut self, other: &Matrix) {
        assert!(self.rows == other.rows);
        assert!(self.cols == other.cols);
        for i in 0..self.rows {
            for j in 0..self.cols {
                self.data[i * self.cols + j] -= other.data[i * self.cols + j].clone();
            }
        }
    }

    pub fn add_int_inplace(&mut self, other: &Integer) {
        for i in 0..self.rows {
            for j in 0..self.cols {
                self.data[i * self.cols + j] += other.clone();
            }
        }
    }

    // mult
    pub fn mul_inplace(&mut self, other: &Matrix) {
        assert!(self.cols == other.rows);
        let mut data: Vec<Integer> = Vec::new();
        for i in 0..self.rows {
            for j in 0..other.cols {
                let mut val = Integer::from(0);
                for k in 0..self.cols {
                    val += self.data[i * self.cols + k].clone() * other.data[k * other.cols + j].clone();
                }
                data.push(val);
            }
        }
        self.data = data;
        self.cols = other.cols;
    }

    pub fn mul_scalar_inplace(&mut self, scalar: &Integer) {
        for i in 0..self.rows {
            for j in 0..self.cols {
                self.data[i * self.cols + j] *= scalar.clone();
            }
        }
    }

    // mod
    pub fn mod_inplace(&mut self, modulo: &Integer) {
        for i in 0..self.rows {
            for j in 0..self.cols {
                let val = self.data[i * self.cols + j].clone();
                self.data[i * self.cols + j] = int_mod(&val, modulo);
                // (_, self.data[i * self.cols + j]) = val.div_rem(modulo.clone());
                // if self.data[i * self.cols + j] < 0 {
                //     self.data[i * self.cols + j] += modulo.clone();
                // }
            }
        }
    }

    pub fn get_identity(dim: usize) -> Matrix {
        let mut matrix = Matrix::new(dim, dim);
        for i in 0..dim {
            matrix.set(i, i, Integer::from(1));
        }
        matrix
    }

    pub fn transpose(&self) -> Matrix {
        let mut matrix = Matrix::new(self.cols, self.rows);
        for i in 0..self.rows {
            for j in 0..self.cols {
                matrix.set(j, i, self.get(i, j));
            }
        }
        matrix
    }

    pub fn mul_vec(&self, vec: &Vec<Integer>) -> Vec<Integer> {
        assert!(self.cols == vec.len());
        let mut result = Vec::new();
        for i in 0..self.rows {
            let mut val = Integer::from(0);
            for j in 0..self.cols {
                val += self.data[i * self.cols + j].clone() * vec[j].clone();
            }
            result.push(val);
        }
        result
    }

    pub fn tensor_product(a: &Matrix, b: &Matrix, p: &Integer) -> Matrix {
        let mut out = Matrix::new(a.rows * b.rows, a.cols * b.cols);

        for i in 0..a.rows {
            for j in 0..a.cols {
                for m in 0..b.rows {
                    for n in 0..b.cols {
                        let val1 = a.get(i, j);
                        let val2 = b.get(m, n);
                        // let val = val1 * val2 % p;
                        let val = val1.clone() * val2.clone() % p;
                        out.set(i * b.rows + m, j * b.cols + n, val);
                    }
                }
            }
        }
        out
    }

    pub fn tensor_product_vec_left(vec: &Vec<Integer>, mat: &Matrix, p: &Integer) -> Matrix {
        let mut out = Matrix::new(vec.len() * mat.rows, mat.cols);

        for i in 0..vec.len() {
            for j in 0..mat.rows {
                for m in 0..mat.cols {
                    let val1 = vec[i].clone();
                    let val2 = mat.get(j, m);
                    let val = val1 * val2 % p;
                    out.set(i * mat.rows + j, m, val);
                }
            }
        }
        out
    }

    pub fn tensor_product_vec_right(mat: &Matrix, vec: &Vec<Integer>, p: &Integer) -> Matrix {
        let mut out = Matrix::new(mat.rows * vec.len(), mat.cols);

        for i in 0..mat.rows {
            for j in 0..mat.cols {
                let val_mat = mat.get(i, j);
                for m in 0..vec.len() {
                    let val_vec = vec[m].clone();
                    let val = val_mat.clone() * val_vec % p;
                    out.set(i * vec.len() + m, j, val);
                }
            }
        }
        out
    }
}


use std::fmt;
impl fmt::Display for Matrix {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for i in 0..self.rows {
            for j in 0..self.cols {
                write!(f, "{} ", self.get(i, j))?;
            }
            writeln!(f)?;
        }
        Ok(())
    }
}

impl Clone for Matrix {
    fn clone(&self) -> Matrix {
        self.clone()
    }
}


use std::ops::{Add, Sub, Mul, Neg, Rem, RemAssign};

use super::vector::{eval_quadratic};


impl<'a, 'b> Add<&'b Matrix> for &'a Matrix {
    type Output = Matrix;

    fn add(self, rhs: &'b Matrix) -> Matrix {
        assert_eq!(self.rows, rhs.rows);
        assert_eq!(self.cols, rhs.cols);

        let mut result = self.clone();
        result.add_inplace(&rhs);
        result
    }
}

impl Add<&Matrix> for Matrix {
    type Output = Matrix;

    fn add(self, rhs: &Matrix) -> Matrix {
        &self + rhs
    }
}

impl Add<Matrix> for Matrix {
    type Output = Matrix;

    fn add(self, rhs: Matrix) -> Matrix {
        &self + &rhs
    }
}


impl<'a, 'b> Sub<&'b Matrix> for &'a Matrix {
    type Output = Matrix;

    fn sub(self, rhs: &'b Matrix) -> Matrix {
        assert_eq!(self.rows, rhs.rows);
        assert_eq!(self.cols, rhs.cols);

        let mut result = self.clone();
        result.sub_inplace(&rhs);
        result
    }
}

impl Sub<&Matrix> for Matrix {
    type Output = Matrix;

    fn sub(self, rhs: &Matrix) -> Matrix {
        &self - rhs
    }
}

impl Sub<Matrix> for Matrix {
    type Output = Matrix;

    fn sub(self, rhs: Matrix) -> Matrix {
        &self - &rhs
    }
}

impl<'a, 'b> Mul<&'b Matrix> for &'a Matrix {
    type Output = Matrix;

    fn mul(self, rhs: &'b Matrix) -> Matrix {
        assert_eq!(self.cols, rhs.rows);

        let mut result = self.clone();
        result.mul_inplace(&rhs);
        result
    }
}

impl Mul<&Matrix> for Matrix {
    type Output = Matrix;

    fn mul(self, rhs: &Matrix) -> Matrix {
        &self * rhs
    }
}

impl Mul<Matrix> for Matrix {
    type Output = Matrix;

    fn mul(self, rhs: Matrix) -> Matrix {
        &self * &rhs
    }
}


impl Mul<Vec<Integer>> for Matrix {
    type Output = Vec<Integer>;

    fn mul(self, rhs: Vec<Integer>) -> Vec<Integer> {
        assert_eq!(self.cols, rhs.len());

        let mut result = Vec::with_capacity(self.rows);
        for i in 0..self.rows {
            let mut sum = Integer::from(0);
            for j in 0..self.cols {
                let a = self.get(i, j);
                let b = &rhs[j];
                sum += a * b;
            }
            result.push(sum);
        }
        result
    }
}

impl Mul<Integer> for Matrix {
    type Output = Matrix;

    fn mul(self, rhs: Integer) -> Matrix {
        let mut result = self.clone();
        for i in 0..self.rows {
            for j in 0..self.cols {
                result.set(i, j, self.get(i, j) * rhs.clone());
            }
        }
        result
    }
}

impl<'a> Neg for &'a Matrix {
    type Output = Matrix;

    fn neg(self) -> Matrix {
        let mut result = self.clone();
        result.negate_inplace();
        result
    }
}

impl Neg for Matrix {
    type Output = Matrix;

    fn neg(self) -> Matrix {
        -&self
    }
}

impl Rem<&Integer> for Matrix {
    type Output = Matrix;

    fn rem(self, modulo: &Integer) -> Matrix {
        let mut result = self.clone();
        result.mod_inplace(modulo);
        result
    }
}

impl RemAssign<&Integer> for Matrix {
    fn rem_assign(&mut self, modulo: &Integer) {
        self.mod_inplace(modulo);
    }
}

pub fn concatenate_row(a: &Matrix, b: &Matrix) -> Matrix {
    assert_eq!(a.cols, b.cols);
    let mut result = Matrix::new(a.rows + b.rows, a.cols);
    for i in 0..a.rows {
        for j in 0..a.cols {
            result.set(i, j, a.get(i, j));
        }
    }
    for i in 0..b.rows {
        for j in 0..b.cols {
            result.set(i + a.rows, j, b.get(i, j));
        }
    }
    result
}

pub fn concatenate_vec_row(a: &Matrix, b: &Vec<Integer>) -> Matrix {
    assert_eq!(a.cols, b.len());
    let mut result = Matrix::new(a.rows + 1, a.cols);
    for i in 0..a.rows {
        for j in 0..a.cols {
            result.set(i, j, a.get(i, j));
        }
    }
    for j in 0..b.len() {
        result.set(a.rows, j, b[j].clone());
    }
    result
}

pub fn concatenate_col(a: &Matrix, b: &Matrix) -> Matrix {
    assert_eq!(a.rows, b.rows);
    let mut result = Matrix::new(a.rows, a.cols + b.cols);
    for i in 0..a.rows {
        for j in 0..a.cols {
            result.set(i, j, a.get(i, j));
        }
    }
    for i in 0..b.rows {
        for j in 0..b.cols {
            result.set(i, j + a.cols, b.get(i, j));
        }
    }
    result
}

pub fn concatenate_vec_col(a: &Matrix, b: &Vec<Integer>) -> Matrix {
    assert_eq!(a.rows, b.len());
    let mut result = Matrix::new(a.rows, a.cols + 1);
    for i in 0..a.rows {
        for j in 0..a.cols {
            result.set(i, j, a.get(i, j));
        }
        result.set(i, a.cols, b[i].clone());
    }
    result
}

pub fn concatenate_diag_one(a: &Matrix) -> Matrix {
    let mut result = Matrix::new(a.rows + 1, a.cols + 1);
    for i in 0..a.rows {
        for j in 0..a.cols {
            result.set(i, j, a.get(i, j));
        }
    }
    for i in 0..a.rows {
        result.set(i, a.cols, Integer::from(0));
    }
    for j in 0..a.cols {
        result.set(a.rows, j, Integer::from(0));
    }
    result.set(a.rows, a.cols, Integer::from(1));
    result
}

pub fn remove_diag_one(a: &Matrix) -> Matrix {
    let mut result = Matrix::new(a.rows - 1, a.cols - 1);
    for i in 0..a.rows - 1 {
        for j in 0..a.cols - 1 {
            result.set(i, j, a.get(i, j));
        }
    }
    result
}

pub fn transpose(a: &Matrix) -> Matrix {
    let mut result = Matrix::new(a.cols, a.rows);
    for i in 0..a.rows {
        for j in 0..a.cols {
            result.set(j, i, a.get(i, j));
        }
    }
    result
}

pub fn sample_unimodular_matrix_from_perm(dim: usize, perm_loops: usize) -> (Matrix, Matrix) {
    let mut rng = RandState::new();

    let mut u = Matrix::get_identity(dim);
    let mut u_inv = Matrix::get_identity(dim);

    let d = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("Duration since UNIX_EPOCH failed")
        .as_secs();
    rng.seed(&Integer::from(d));

    for _i in 0..perm_loops {
        let x = rng.below(dim as u32) as usize;
        let y = rng.below(dim as u32) as usize;
        if x != y {
            for j in 0..dim {
                let val = u.get(x, j).clone() + u.get(y, j).clone();
                u.set(y, j, val.clone());
                let val = u_inv.get(j, x).clone() - u_inv.get(j, y).clone();
                u_inv.set(j, x, val.clone());
            }
        }
    }
    (u, u_inv)
}


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
            let ratio: Integer;
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

pub fn matrix_inverse(
    m: &mut Matrix,
    mod_val: &Integer,
) -> Result<Matrix, i32> {
    assert_eq!(m.rows, m.cols);
    let n = m.rows;
    let mut m_inv = Matrix::get_identity(n);
    let mut m_aug = concatenate_col(m, &m_inv);
    let (_pivot_cols, _free_vars, r) = echelon_form(&mut m_aug, mod_val);
    
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

pub fn generate_right_inverse_space(
    row: usize,
    col: usize,
    modulus: &Integer,
    rng: &mut RandState<'_>,
) -> (Matrix, Matrix, Matrix) {
    // u, u_inv are col x col unimodular matrices
    // r,r_inv are row x row random matrices
    // d is a row x col diagonal matrix
    // d_inv is a col x row diagonal matrix
    // n is a concatenation of a row x row zero matrix and a (col - row) x row random matrix
    // define
    //  a = r * d * u
    //  a_inv = u_inv * d_inv * r_inv
    //  null = u_inv * n

    let mod_half: Integer = modulus.clone() / 2;
    let (mut u, mut u_inv) = sample_unimodular_matrix_from_perm(col, col * 128);
    u %= modulus;
    u_inv %= modulus;

    let mut d = Matrix::new(row, col);
    let mut d_inv = Matrix::new(col, row);

    for i in 0..row {
        let mut inv = Integer::from(0);
        let mut val = Integer::from(0);
        while inv == Integer::from(0) {
            // Sample random val (odd)
            let random_val = mod_half.clone().random_below(rng);
            val = random_val.clone() * 2 + 1;
            let val_copy = val.clone();
            inv = match val_copy.invert(modulus) {
                Ok(inv) => inv,
                Err(_) => Integer::from(0),
            };
        }
        d.set(i, i, val.clone());
        d_inv.set(i, i, inv.clone());
    }

    let mut a = d * &u;
    a %= modulus;
    let mut a_inv = u_inv.clone() * d_inv;
    a_inv %= modulus;

    let n_upper = Matrix::new(row, col - row);
    let n_lower = Matrix::random(col - row, col - row, modulus, rng);
    let n = concatenate_row(&n_upper, &n_lower);
    let mut null = u_inv * n;
    null %= modulus;

    (a, a_inv, null)
}

pub fn eval_quadratic_multivariate(
    x: &Vec<Integer>,
    y: &Vec<Integer>,
    f_mat: &Matrix,
) -> Vec<Integer> {
    let mut result = vec![Integer::from(0); f_mat.rows];
    for i in 0..f_mat.rows {
        let f = f_mat.get_row(i);
        result[i] = eval_quadratic(&x, &y, &f);
    }
    result
}