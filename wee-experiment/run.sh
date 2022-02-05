#!/bin/bash

set -euo pipefail

mkdir -p results$1

time nl ../../cpaior2019-sbs-for-subgraphs-paper/experiments/instances.txt | while read n a b c d; do
    echo $n $a timeout=$1
    ../mcsp_sparse -t $1 --lad B ../../cpaior2019-sbs-for-subgraphs-paper/instances/{$b,$c} > results$1/$a.out
done
