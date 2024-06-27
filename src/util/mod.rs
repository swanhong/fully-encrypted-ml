pub mod group;
pub mod vector;
pub mod matrix;
pub mod decomp;
pub mod arguments;

#[cfg(test)]
mod tests {
    use rug::Integer;
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

    use crate::util::vector::gen_random_vector;
    use crate::util::vector::tensor_product_vecs;
    use crate::util::vector::vec_inner_pow;

    use super::matrix::Matrix;

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
        let col = 5;
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
    use super::vector::vec_exp_with_base;
    #[test]
    fn test_decompose_and_power() {
        // Create a decomposition instance
        let dim = 2; // Example dimension
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
        
        let n_x = 2;
        let b = 1;
        let rng = &mut RandState::new(); // Create a single RandState object
        let dx = Matrix::random(n_x, n_x+b, &grp.delta, rng);
        let dy = Matrix::random(n_x, n_x+b, &grp.delta, rng);
        let dxdy = Matrix::tensor_product(&dx, &dy, &grp.delta);
        let dxdyt = dxdy.transpose();
        let f = gen_random_vector(dxdyt.cols, &Integer::from(10), rng);
        let sk_red = vec_exp_with_base(&grp.g, &(dxdyt.mul_vec(&f)), &grp.n_sq);
        println!("dxdy shape = {} x {}", dxdyt.rows, dxdyt.cols);
        println!("sk_red (size = {})= {:?}", sk_red.len(), sk_red);
        println!("dx = \n{}", dx);
        println!("dy = \n{}", dy);
        let dx_power = decomp.matrix_pow_col(&dx);
        let dy_power = decomp.matrix_pow_col(&dy);
        println!("dx_power = \n{}", dx_power);
        println!("dy_power = \n{}", dy_power);
        let dxdy_power = Matrix::tensor_product(&dx_power, &dy_power, &grp.delta);
        let dxdy_power_t = dxdy_power.transpose();
        println!("dxdy_power shape = {} x {}", dxdy_power_t.rows, dxdy_power_t.cols);
        let sk_red_power = vec_exp_with_base(&grp.g, &(dxdy_power_t.mul_vec(&f)), &grp.n_sq);
        println!("sk_red_power (size = {})= {:?}", sk_red_power.len(), sk_red_power);

        let enc_x = gen_random_vector(n_x+1, &Integer::from(10), rng);
        let enc_y = gen_random_vector(n_x+1, &Integer::from(10), rng);

        let enc_x_power = decomp.vector(&enc_x);
        let enc_y_power = decomp.vector(&enc_y);

        let enc_xy = tensor_product_vecs(&enc_x, &enc_y, &grp.delta);
        let enc_xy_power = tensor_product_vecs(&enc_x_power, &enc_y_power, &grp.delta);

        let val1 = vec_inner_pow(&sk_red, &enc_xy, &grp);
        let val2 = vec_inner_pow(&sk_red_power, &enc_xy_power, &grp);

        println!("val1 = {}", val1);
        println!("val2 = {}", val2);

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

    use std::time::Instant;
    #[test]
    fn test_powm() {
        // Create a large base and exponent for testing
        let p = Integer::from_str_radix("3793070626895572467283950491156608172646280619543651052186002962401960948344595436008317725360316523061862837176776937156394891603251329682527611995573870675378784146300080405741531439158689694869658339169636599335667867300511737877416759591946775631328177983594319449174751169843162114768684946303385589552432358881413822213266477101325470478323250412128032770624674451411359365977206532598270731108378914296880563062937443011799707142015647962685071743064933626689093956215640185694480670389308276588763735300266054342210518510655503938301168859241232366499082115661704416570624348536949976639770352500574819937312763118260737810176902179509516452097179033706353327008781474150697920827389789734610788578490948322681570725850626191531687621385181189566240487471233234659263652702257011861689546257355290042178376098634268810888992322511302018026179564203024011436332852715101641086471023435541325718482408540095488497904963", 10).unwrap();
        let q = Integer::from_str_radix("3437698264967377277030634390176995339108761179197266750352568811124733920898558792968942018331196975819629005205008386113587643859808501473991043603423893930344847276501635656115934929863928487189206920868313713647633563058406832931496433057597415655426803516291613870104861183026912877156823405411321004129466972852732763205585115072430483394506661540167036456924629674894206529803341592142340709693772564216807751480082209883528322329167890554845107315685301749569796946314108563426847831535198565314370173512451464243428810143307683852607035244704650171763912281349917627223963810996092135117877342840690507979526886087476901475153679297045764348723938821743929206668894364791858405008648706877501022484609478488897352403273816492628362528291051657763913142495017627088665776991501665058363019172356379867940793605694439311020508346728498443971266879607174176083764814467005057407301185462205210470795534594752168817506103", 10).unwrap();
        let modulus = p.clone() * q.clone();

        // Measure the time taken to compute the pow_mod operation
        let start_time = Instant::now();
        let result = p.clone().pow_mod(&q, &modulus).unwrap();
        let duration = start_time.elapsed();

        // Print the result and the time taken
        println!("Result: {}", result);
        println!("Time taken: {:?}", duration);
    }

}