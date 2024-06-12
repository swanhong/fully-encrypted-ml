#!/bin/bash
#SBATCH --job-name=bfe_ml
#SBATCH --output=out060824/ml/iris_scaled_%j.txt
#SBATCH --mem=200G
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=BEGIN,END,ERROR
#SBATCH --mail-user=shong@nygenome.org

 # ./target/debug/fe -t protocol --bit-len 3072 --dim0 6 --dim1 $1 --dim2 1
./target/debug/fe -t ml --bit-len 3072 -d iris
# ./target/debug/fe -t ml -d breast --bit-len 300
