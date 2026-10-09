#!/bin/bash
# Small cases: q=3 (n=2..12) and q=9 (n=2,3,4,6), q=27 (n=2). Output: small.log
cd "$(dirname "$0")"
gcc -O2 -march=native -fopenmp -o apnorm apnorm.c || exit 1
: > small.log
for n in 2 3 4 5 6 7 8 9 10 12; do ./apnorm 1 $n 16 2>>small.err | tee -a small.log; cp fails.txt fails_3_$n.txt; done
for n in 2 3 4 6; do ./apnorm 2 $n 16 2>>small.err | tee -a small.log; cp fails.txt fails_9_$n.txt; done
./apnorm 3 2 16 2>>small.err | tee -a small.log; cp fails.txt fails_27_2.txt
