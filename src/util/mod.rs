pub mod group;
pub mod vector;
pub mod matrix;
pub mod decomp;

#[cfg(test)]
mod tests {
    use rug::Integer;
    use rug::ops::Pow;
    use rug::rand::RandState;
    use std::time::SystemTime;

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

    use crate::util::vector::vec_inner_pow;

    use super::matrix::Matrix;

    #[test]
    fn test_matrix_operations() {
        let bound = Integer::from(100);
        // let mat1 = Matrix::random(2, 2, &bound);
        // let mat2 = Matrix::random(2, 2, &bound);
        // let mat3 = (mat1.clone() + &mat2) % &bound;
        // let mat4 = mat1.clone() * &mat2 % &bound;

        // println!("Matrix 1:");
        // println!("{}", mat1);

        // println!("Matrix 2:");
        // println!("{}", mat2);
    
        // println!("Matrix Addition Result:");
        // println!("{}", mat3);
    
        // println!("Matrix Multiplication Result:");
        // println!("{}", mat4);
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

        let mut rand = RandState::new(); // Create a single RandState object
        let d = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .expect("Duration since UNIX_EPOCH failed");
        rand.seed(&Integer::from(d.as_secs()));


        let (a, a_inv, null) = generate_right_inverse_space(row, col, &modulus, &mut rand);

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

    use super::decomp::Decomp;
    use super::group::Group;
    #[test]
    fn test_decompose_and_power() {
        // Create a decomposition instance
        let dim = 4; // Example dimension
        let grp = Group::new(10);
        let decomp = Decomp::new(dim, &grp);

        println!("decomp.modulo = {}", decomp.modulo);
        println!("decomp.base = {}", decomp.base);
        println!("decomp.modulo_exp = {}", decomp.modulo_for_exp);

        // Example input values
        let input_int = Integer::from(123);
        let input_vec = vec![
            Integer::from(987),
            Integer::from(654),
            Integer::from(321),
            Integer::from(123),
        ];
        let mut input_matrix = Matrix::new(2, 3);
        input_matrix.set(0, 0, Integer::from(10));
        input_matrix.set(0, 1, Integer::from(20));
        input_matrix.set(0, 2, Integer::from(30));
        input_matrix.set(1, 0, Integer::from(40));
        input_matrix.set(1, 1, Integer::from(50));
        input_matrix.set(1, 2, Integer::from(60));

        // Decompose an Integer
        let int_decomposed = decomp.int(&input_int);
        println!("input : {}", input_int);
        println!("Decomposed Integer: {:?}", int_decomposed);

        // Decompose a vector of Integers
        let vec_decomposed = decomp.vector(&input_vec);
        println!("input : {:?}", input_vec);
        println!("Decomposed Vector: {:?}", vec_decomposed);

        // Decompose a matrix by columns
        let matrix_col_decomposed = decomp.matrix_col(&input_matrix);
        println!("input : \n{}", input_matrix);
        println!("Decomposed Matrix (by columns):\n{}", matrix_col_decomposed);

        // Decompose a matrix by rows
        let matrix_row_decomposed = decomp.matrix_row(&input_matrix);
        println!("input : \n{}", input_matrix);
        println!("Decomposed Matrix (by rows):\n{}", matrix_row_decomposed);

        // Compute the power of an Integer
        let int_power = Integer::from(2);
        let int_powered = decomp.int_pow(&int_power);
        println!("input : {}", int_power);
        println!("Powered Integer: {:?}", int_powered);

        // Compute the power of a vector of Integers
        let vec_powered = decomp.vector_pow(&input_vec);
        println!("input : {:?}", input_vec);
        println!("Powered Vector: {:?}", vec_powered);

        // Compute the power of a vector of Integers using exponentiation
        let vec_powered_exp = decomp.vector_pow_exp(&input_vec); // Replace `grp` with your Group instance
        println!("input : {:?}", input_vec);
        println!("Powered Vector (using exponentiation): {:?}", vec_powered_exp);

        // Compute the power of a matrix by rows
        let matrix_row_powered = decomp.matrix_pow_row(&input_matrix);
        println!("input : \n{}", input_matrix);
        println!("Powered Matrix (by rows):\n{}", matrix_row_powered);

        // Compute the power of a matrix by columns
        let matrix_col_powered = decomp.matrix_pow_col(&input_matrix);
        println!("input : \n{}", input_matrix);
        println!("Powered Matrix (by columns):\n{}", matrix_col_powered);

        // Aggregating the results
        let aggregated_result: Integer = int_powered.iter().sum();
        println!("Aggregated Result: {}", aggregated_result);

        // generate two integer vectors
        let x = vec![
            Integer::from(123),
            Integer::from(223),
            Integer::from(312),
            Integer::from(411),
        ];
        let y = vec![
            Integer::from(511),
            Integer::from(623),
            Integer::from(723),
            Integer::from(812),
        ];

        // compute the inner product of x and y
        let mut out = Integer::from(0);
        for i in 0..x.len() {
            out += x[i].clone() * y[i].clone();
        }

        // decompose x and power y
        let x_decomposed = decomp.vector(&x);
        let y_powered = decomp.vector_pow(&y);
        let mut out2 = Integer::from(0);
        for i in 0..x_decomposed.len() {
            out2 += x_decomposed[i].clone() * y_powered[i].clone();
        }

        println!("x = {:?}", x);
        println!("y = {:?}", y);
        println!("out = {}", out);

        println!("x_decomposed = {:?}", x_decomposed);
        println!("y_powered = {:?}", y_powered);
        println!("out2 = {}", out2);

        let mut out_exp = Integer::from(1);
        for i in 0..x.len() {
            out_exp *= y[i].clone().pow_mod(&x[i], &grp.n_sq).unwrap();
            out_exp %= &grp.n_sq;
        }

        let y_powered_exp = decomp.vector_pow_exp(&y);
        let mut out2_exp = Integer::from(1);
        for i in 0..x_decomposed.len() {
            out2_exp *= y_powered_exp[i].clone().pow_mod(&x_decomposed[i], &grp.n_sq).unwrap();
            out2_exp %= &grp.n_sq;
        }

        let out3_exp = vec_inner_pow(&y_powered_exp, &x_decomposed, &grp);

        println!("y = {:?}", y);
        println!("y_powered = {:?}", y_powered);
        println!("y_powered_exp = {:?}", y_powered_exp);


        println!("out_exp = {}", out_exp);
        println!("out2_exp = {}", out2_exp);
        println!("out3_exp = {}", out3_exp);

        println!("decomp.modulo = {}", decomp.modulo);
        println!("decomp.base = {}", decomp.base);
        println!("decomp.modulo_exp = {}", decomp.modulo_for_exp);
    }

    #[test]
    fn test_parse_matrix() {
        let input = "[
            [ 92767638662  576911443709  156205951134  406459819186 ]
    [ 132160905606  406733278377  940896381053  770350121770 ]
    [ 773805851499  194799123750  974833913434  872906692024 ]
    [ 719952816917  262229557635  565821972923  410323596120 ]
    [ 817854487975  46131499252  764294592157  625909382421 ]
    [ 927585505118  81815493437  247789617931  489272859060 ]
    [ 561537790429  614934518533  642147877359  644031578406 ]
    [ 280599438082  999971579106  324364833313  570825953465 ]
    [ 195862779364  248697167504  938286923331  614741748946 ]
    ]";
        let matrix = Matrix::parse_matrix(input);

        println!("{}", matrix);
    }
}