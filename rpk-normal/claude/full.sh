#!/usr/bin/env bash
cd "$(cd "$(dirname "$0")" && pwd)"
gcc -O2 -fopenmp -o /tmp/rsearch32 rsearch.c || exit 1
export OMP_NUM_THREADS=8
while IFS= read -r spec; do
  echo "-- $spec"; timeout 600 /tmp/rsearch32 params/$spec.txt 3e10 8 || echo "TIMEOUT_OR_FAIL $spec"
done < todo2.txt
echo "== done"
