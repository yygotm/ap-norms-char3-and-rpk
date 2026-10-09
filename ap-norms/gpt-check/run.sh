#!/usr/bin/env bash
# Usage: DATA=<dir containing the <q>_<n> certificate folders> bash ap-norms/gpt-check/run.sh build|batch|full|sample|coverage <pair> [threads] [count]
set -euo pipefail
work="$(cd "$(dirname "$0")" && pwd)"
data=${DATA:?set DATA to the directory holding the certificate folders}
mode=${1:?mode required}
if [[ $mode == batch ]]; then
  {
    printf 'UTC batch started: '
    date -u +%FT%TZ
    sha256sum "$work/check.cpp" "$work/check" "$work/fixtures.txt"
    for pair in 3_7 3_8 3_10 3_12 9_6 9_8 27_6; do
      printf 'UTC pair started: '
      date -u +%FT%TZ
      "$work/check" full "$data/$pair" 8 2>&1 | tee "$work/${pair}-full.log"
      printf 'UTC pair finished: '
      date -u +%FT%TZ
    done
    printf 'UTC batch finished: '
    date -u +%FT%TZ
  } | tee "$work/batch.log"
  exit
fi
if [[ $mode == build ]]; then
  g++ -O3 -march=native -fopenmp -std=c++20 -Wall -Wextra "$work/check.cpp" -o "$work/check"
  "$work/check" fixtures "$work/fixtures.txt" | tee "$work/fixtures.log"
  exit
fi
pair=${2:?pair required}
case "$pair" in 3_6|3_7|3_8|3_10|3_12|9_6|9_8|27_6) ;; *) exit 2;; esac
case "$mode" in full|sample|coverage) ;; *) exit 2;; esac
threads=${3:-8}
count=${4:-30000}
printf 'UTC started: '
date -u +%FT%TZ
"$work/check" "$mode" "$data/$pair" "$threads" "$count" | tee "$work/${pair}-${mode}.log"
printf 'UTC finished: '
date -u +%FT%TZ
