#![allow(dead_code)]

use std::fs::File;
use std::io::{BufRead, BufReader};
use rug::Integer;
use crate::util::matrix::Matrix;

pub fn read_file() {
    // Read the text file
    let file = File::open("plain_model/iris/input_layer.weight.csv").expect("Failed to open file");
    let reader = BufReader::new(file);

    // Initialize variables
    let mut matrix: Vec<Vec<f64>> = Vec::new();
    let mut row_count = 0;
    let mut column_count = 0;

    // Process each line of the file
    for line in reader.lines() {
        if let Ok(row_data) = line {
            let row: Vec<f64> = row_data
                .split(',')
                .map(|s| s.trim().parse().expect("Failed to parse value"))
                .collect();

            // Update column count
            column_count = row.len();

            // Add the row to the matrix
            matrix.push(row);

            // Update row count
            row_count += 1;
        }
    }

    // Display the matrix
    println!("Matrix:");
    for row in &matrix {
        println!("{:?}", row);
    }

    // Display row and column numbers
    println!("Rows: {}", row_count);
    println!("Columns: {}", column_count);
}

pub fn read_float_matrix(
    filename: &str,
) -> Vec<Vec<f64>> {
    // Read the text file
    let file = File::open(filename).expect("Failed to open file");
    let reader = BufReader::new(file);

    // Initialize variables
    let mut matrix: Vec<Vec<f64>> = Vec::new();

    // Process each line of the file
    for line in reader.lines() {
        if let Ok(row_data) = line {
            let row: Vec<f64> = row_data
                .split(',')
                .map(|s| s.trim().parse().expect("Failed to parse value"))
                .collect();

            // Add the row to the matrix
            matrix.push(row);
        }
    }    matrix
}

pub fn read_matrix(filename: &str, scale: i64) -> Matrix {
    let matrix = read_float_matrix(filename);
    convert_matrix(matrix, scale)
}

pub fn display_float_matrix(matrix: &Vec<Vec<f64>>) {
    // Display the matrix
    println!("Matrix:");
    for row in matrix {
        println!("{:?}", row);
    }

    // Display row and column numbers
    println!("Rows: {}", matrix.len());
    println!("Columns: {}", matrix[0].len());
}

pub fn convert_matrix(matrix: Vec<Vec<f64>>, scale: i64) -> Matrix {
    let mut m = Matrix::new(matrix.len(), matrix[0].len());
    for i in 0..matrix.len() {
        for j in 0..matrix[0].len() {
            let mut encode = Integer::new();
            encode.assign_f64(matrix[i][j] * scale as f64).unwrap();
            m.set(i, j, encode);
        }
    }
    m
}