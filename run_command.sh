#!/bin/bash

for i in 1 2 3 4 5 6; do
  echo ${i} start
  ./target/debug/fe -t protocol --dim $i --dim $i --dim 1 --bit-len 3072 > out080224/protocol/${i}.txt
  echo ${i} done
done
for d in iris breast; do 
  echo ${d} start
  ./target/debug/fe -t ml -d $d --bit-len 3072 > out080224/ml/${d}.txt
  echo ${d} done
done