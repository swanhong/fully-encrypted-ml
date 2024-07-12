<p align="center">
    <h1 align="center">ComposableFE-rs</h1>
</p>

## Overview

This repository contains the `Rust` implementation of the scheme suggested in the paper "Fully Encrypted Machine Learning Protocol using Functional Encryption," which is currently under review.

## Build
Our code is written in Rust, so it only depends on Rust. After installing Rust, run the following command:
```bash
    # cd ~
    cargo build
```

## Run

### Test
To test the code, run following command:
```bash
    ./target/debug/fe -t protocol
```
This runs a 2-level composable FE evaluation. The default dimensions are [5, 4, 3], which means that it chooses a random vector of dimension 5, evaluates a random quadratic polynomial (with constant term) to output a dimension 4 vector, and then evaluates another random quadratic polynomial to finally output a dimension 3 vector. The bit length of primes for the DCR group is set to 100, so the code shows fast (but not secure) results.

### Run Protocol

The target option `protocol` runs the protocol with a random message and evaluates the composition of two quadratic functions.
```bash
    ./target/debug/fe --t protocol --bit-len 3072 --dim 5 --dim 4 --dim 3 --dim 2 --n-decomp 4
```
The arguments are as follows:
	- `--t`: Specifies the target option. In this case, it's protocol.
	- `--bit-len`: The bit length of primes used for the DCR group. To ensure 128-bit security, this parameter should be chosen as 3072. Option: (10, 100, 3072)
	- `--dim`: The dimensions. first input set the input dimension, and the next inputs set the following dimentions
	- `--n-decomp`: Decomposition parameter

### Run ML with UCI datasets

Here, we exploit `Iris` or `breast` dataset from UCI to test evaluation of 2-layer neural network with square activation function using our composable FE. Here, the model is already trained in `plain_model/{DATA_NAME}` folder. Refer each `*.ipynb` files for the details.

In each `plain_model/{DATA_NAME}` folder, the test data is stored in `X_test.csv` file and the information of each layer are contained in `{input/output}_layer.{bias/weight}.csv` file, respectively.

Then, to test evaluation of model using "first" sample in test file, run following command:
```bash
	./target/debug/fe -t ml --bit-len 3072 -d {DATA_NAME}
```
Put `iris` or `breast` on `{DATA_NAME}`. 