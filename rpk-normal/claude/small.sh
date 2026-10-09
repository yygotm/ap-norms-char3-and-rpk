#!/usr/bin/env bash
# small exhaustive runs. Usage: bash rpk-normal/claude/small.sh (WSL/Linux, gcc with OpenMP)
set -u
D="$(cd "$(dirname "$0")" && pwd)"
cd "$D"
gcc -O2 -fopenmp -o /tmp/brute32 brute.c || exit 1
export OMP_NUM_THREADS=8
echo "== selftests"
/tmp/brute32 2 1 7 s
/tmp/brute32 2 2 3 s
echo "== runs"
for spec in "2 1 8" "2 1 10" "2 1 12" "2 1 14" "2 1 16" "2 1 18" "2 2 7" "2 2 8" "2 2 9" "2 2 10" "5 1 8" "7 1 7" "7 1 8" "2 2 12" "5 1 10" "2 3 8" "7 1 9"; do
  echo "-- p e m = $spec"
  /tmp/brute32 $spec
done
echo "== done"
