use rug::Integer;
use super::{matrix::Matrix, group::Group};

pub struct Decomp {
    pub base: Integer,
    pub dim: usize,
    pub modulo: Integer,
    pub modulo_for_exp: Integer,
}

// decompose every Integer into base "base".
// note that base^dim = modulus

impl Decomp {
    pub fn new(dim: usize, grp: &Group) -> Decomp {
        let modulo = grp.delta.clone();
        let modulo_exp = grp.n_sq.clone();
        let (mut base, rem) = modulo.clone().root_rem(Integer::new(), dim as u32);
        if rem > 0 {
            base = base.clone() + 1;
        }
        // println!("base = {}", base);
        Decomp {
            base,
            dim,
            modulo: modulo.clone(),
            modulo_for_exp: modulo_exp.clone(),
        }
    }

    // decompose 
    pub fn int(&self, a: &Integer) -> Vec<Integer> {
        assert!(a < &self.modulo, "a = {}, self.modulo = {}", a, self.modulo);
        // println!("input : {}", a);
        let mut res: Vec<Integer> = vec![Integer::from(0); self.dim];
        let mut b = a.clone();
        let zero = Integer::from(0);
        let mut index = 0;
        while b > zero {
            let (quo, rem) = b.clone().div_rem_euc(self.base.clone());
            // println!("quo = {}, rem = {}", quo, rem);
            res[index] = rem;
            b = quo.clone();
            index += 1;
        }
        // println!("output: {:?}", res);
        assert_eq!(res.len(), self.dim, "res.len() = {}, self.dim = {}", res.len(), self.dim);
        res
    }
    
    pub fn vector(&self, a: &Vec<Integer>) -> Vec<Integer> {
        let mut res: Vec<Integer> = Vec::new();
        for i in 0..a.len() {
            res.extend(self.int(&a[i]));
        }
        res
    }

    pub fn matrix_col(&self, a: &Matrix) -> Matrix {
        let mut res = Matrix::new(a.rows * self.dim, a.cols);
        for j in 0..a.cols {
            let col_decomp = self.vector(&a.get_col(j));
            res.set_col(j, &col_decomp);
        }
        res
    }
    pub fn matrix_row(&self, a: &Matrix) -> Matrix {
        let mut res = Matrix::new(a.rows, a.cols * self.dim);
        for i in 0..a.rows {
            let row_decomp = self.vector(&a.get_row(i));
            res.set_row(i, &row_decomp);
        }
        res
    }

    pub fn int_pow(&self, a: &Integer) -> Vec<Integer> {
        let mut res: Vec<Integer> = vec![Integer::from(0); self.dim];
        let mut b = a.clone();
        for i in 0..self.dim {
            res[i] = b.clone();
            let (_, new_b) = (b.clone() * self.base.clone()).div_rem(self.modulo_for_exp.clone());
            b = new_b;
        }
        assert_eq!(res.len(), self.dim);
        res
    }

    pub fn vector_pow(&self, a: &Vec<Integer>) -> Vec<Integer> {
        let mut res: Vec<Integer> = Vec::new();
        for i in 0..a.len() {
            res.extend(self.int_pow(&a[i]));
        }
        res
    }

    pub fn int_pow_exp(&self, a: &Integer) -> Vec<Integer> {
        let mut res: Vec<Integer> = vec![Integer::from(0); self.dim];
        res[0] = a.clone();
        for i in 1..self.dim {
            res[i] = res[i-1].clone().pow_mod(&self.base, &self.modulo_for_exp).unwrap();
        }
        assert_eq!(res.len(), self.dim);
        res
    }
    pub fn vector_pow_exp(&self, a: &Vec<Integer>) -> Vec<Integer> {
        let mut res: Vec<Integer> = Vec::new();
        for i in 0..a.len() {
            res.extend(self.int_pow_exp(&a[i]));
        }
        res
    }
    

    pub fn matrix_pow_col(&self, a: &Matrix) -> Matrix {
        let mut res = Matrix::new(a.rows, a.cols * self.dim);
        for i in 0..a.rows {
            let row_pow = self.vector_pow(&a.get_row(i));
            res.set_row(i, &row_pow)
        }
        res
    }

    pub fn matrix_pow_row(&self, a: &Matrix) -> Matrix {
        let mut res = Matrix::new(a.rows * self.dim, a.cols);
        for j in 0..a.cols {
            let mut col_pow = self.vector_pow(&a.get_col(j));
            res.set_col(j, &col_pow);
        }
        res
    }
}