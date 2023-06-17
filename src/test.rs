// writh test.rs file
// here I will test generate_matrix_with_right_inverse function in src/util.rs

// here use #[cfg(test)] to only compile this module when running tests

use rug::Integer;

// use crate::util::matrix::generate_matrix_with_right_inverse;


// #[test]
// fn test_generate_matrix_with_right_inverse() {
//     let row = 3;
//     let col = 4;
//     let modulo = Integer::from(10);

//     let (a, a_inv, n) = generate_matrix_with_right_inverse(row, col, &modulo);

//     // Add assertions to check the correctness of the generated matrices
//     assert_eq!(a.len(), row);
//     assert_eq!(a[0].len(), col);
//     // ...

//     assert_eq!(a_inv.len(), col);
//     assert_eq!(a_inv[0].len(), row);
//     // ...

//     assert_eq!(n.len(), row);
//     assert_eq!(n[0].len(), col);
//     // ...

//     // Add more assertions as needed
// }