use rug::{Integer, rand::RandState};
use std::time::{SystemTime, UNIX_EPOCH};
use std::str::FromStr;

#[derive(Clone, Debug)]
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

    pub fn clone(&self) -> Matrix {
        let mut matrix = Matrix::new(self.rows, self.cols);
        for i in 0..self.rows {
            for j in 0..self.cols {
                matrix.set(i, j, self.get(i, j));
            }
        }
        matrix
    }


    pub fn Copy(&self) -> Matrix {
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
                (_, self.data[i * self.cols + j]) = val.div_rem(modulo.clone());
                if self.data[i * self.cols + j] < 0 {
                    self.data[i * self.cols + j] += modulo.clone();
                }
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

use std::ops::{Add, Sub, Mul, Neg, Rem, RemAssign};


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


pub fn generate_right_inverse_space(
    row: usize,
    col: usize,
    modulus: &Integer,
    rng: &mut RandState<'_>,
) -> (Matrix, Matrix, Matrix) {
    // variables for randomness
    let mut rng = RandState::new();
    let d = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("Duration since UNIX_EPOCH failed")
        .as_secs();
   rng.seed(&Integer::from(d));

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
            let random_val = mod_half.clone().random_below(&mut rng);
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

    let n_upper = Matrix::new(row, row);
    let n_lower = Matrix::random(col - row, row, modulus, &mut rng);
    let n = concatenate_row(&n_upper, &n_lower);

    let mut null = u_inv * n;
    null %= modulus;

    (a, a_inv, null)
}
