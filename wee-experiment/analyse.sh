#!/bin/bash

set -eu

cat ../../cpaior2019-sbs-for-subgraphs-paper/experiments/instances.txt | while read a b c d; do
    head -n4 results/$a.out | awk -v f=$a '
        /CPU/ {time=$4}
        /TIMEOUT/ {time=1000000; sat="X"}
        /^SAT/ {sat=1}
        /^UNSAT/ {sat=0}
        END { print f, sat, time}
    '
done | tee results.txt

awk '
        NR==1 {print "instance", "sat", "time", $0}
        NR==FNR {r[$1] = $0}
        NR!=FNR {print $0, r[$1]}
    ' ../../cpaior2019-sbs-for-subgraphs-paper/paper/inducedruntimes.data results.txt | tee combined-results.txt
