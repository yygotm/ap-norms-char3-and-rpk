#!/usr/bin/env bash
cd "$(cd "$(dirname "$0")" && pwd)"
[ -x /tmp/brute32 ] || gcc -O2 -fopenmp -o /tmp/brute32 brute.c
export OMP_NUM_THREADS=8
for spec in "13 1 7" "5 1 12" "11 1 8" "2 4 7"; do echo "-- p e m = $spec"; /tmp/brute32 $spec; done
echo "== done"
