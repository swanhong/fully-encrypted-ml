pub mod file_reader;

#[cfg(test)]
mod test {
    use crate::parse::file_reader::{read_float_matrix, display_float_matrix};

    use super::file_reader::convert_matrix;
    #[test]
    fn test_read_file() {
        let filename = "plain_model/iris/input_layer.weight.csv";
        let matrix = read_float_matrix(&filename);
        display_float_matrix(&matrix);

        let matrix_new = convert_matrix(matrix, 1000);
        println!("Matrix new: {}", matrix_new);
    }
}