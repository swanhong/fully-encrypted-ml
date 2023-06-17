pub mod group;
pub mod vector;
pub mod matrix;
#[cfg(test)]
mod tests {
    use rug::Integer;
    // Matrix::mod_inplace test
    #[test]
    fn test_matrix_mod_inplace() {
        let mut mat = super::matrix::Matrix::new(2, 2);
        mat.set(0, 0, Integer::from(1));
        mat.set(0, 1, Integer::from(2));
        mat.set(1, 0, Integer::from(3));
        mat.set(1, 1, Integer::from(4));

        println!("Matrix before mod:");
        println!("{}", mat);

        mat.mod_inplace(&Integer::from(3));

        println!("Matrix after mod:");
        println!("{}", mat);
    }

    use super::matrix::Matrix;

    #[test]
    fn test_matrix_operations() {
        let bound = Integer::from(100);
        let mat1 = Matrix::random(2, 2, &bound);
        let mat2 = Matrix::random(2, 2, &bound);
        let mat3 = (mat1.clone() + &mat2) % &bound;
        let mat4 = mat1.clone() * &mat2 % &bound;

        println!("Matrix 1:");
        println!("{}", mat1);

        println!("Matrix 2:");
        println!("{}", mat2);
    
        println!("Matrix Addition Result:");
        println!("{}", mat3);
    
        println!("Matrix Multiplication Result:");
        println!("{}", mat4);
    }
    

    use super::matrix::sample_unimodular_matrix_from_perm;
    #[test]
    fn test_unimodular_sampling() {
        let dim = 3;
        let perm_loops = 100;
        let bound = Integer::from(100);
        let (mut u, mut u_inv) = sample_unimodular_matrix_from_perm(dim, perm_loops);

        u %= &bound;
        u_inv %= &bound;
        println!("Matrix U:");
        println!("{}", u);

        println!("Matrix U Inv:");
        println!("{}", u_inv);

        println!("Matrix U * U Inv:");
        let mult = (u * &u_inv) % &bound;
        println!("{}", mult);
    }
    
    use super::matrix::generate_right_inverse_space;
    #[test]
    fn test_generate_matrix_with_right_inverse() {
        let row = 3;
        let col = 4;
        let modulus = Integer::from(7);

        let (a, a_inv, null) = generate_right_inverse_space(row, col, &modulus);

        println!("Matrix A:");
        println!("{}", a);
        println!("Matrix A_inv:");
        println!("{}", a_inv);
        println!("Matrix Null Space:");
        println!("{}", null);

        println!("test A * A_inv");
        let mult = (a.clone() * &a_inv) % &modulus;
        println!("{}", mult);

        println!("test A * null");
        let mult = (a * &null) % &modulus;
        println!("{}", mult);
    }

}