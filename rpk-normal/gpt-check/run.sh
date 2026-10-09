#!/usr/bin/env bash
# Usage: [DATA=<witness dir>] bash rpk-normal/gpt-check/run.sh bench|full [pair|all]
set -euo pipefail
work="$(cd "$(dirname "$0")" && pwd)"
data=${DATA:-$work/../witnesses}
mode=${1:?mode}; pair=${2:-all}
case "$mode" in bench|full) ;; *) exit 2;; esac
g++ -O3 -march=native -fopenmp -std=c++20 -Wall -Wextra "$work/check.cpp" -o "$work/check"
{
  printf 'UTC started: '; date -u +%FT%TZ
  g++ --version | head -n 1
  sha256sum "$work/check.cpp" "$work/check" "$work/cases.txt"
  "$work/check" "$work/cases.txt" "$data" "$mode" "$pair"
  printf 'UTC completed: '; date -u +%FT%TZ
} 2>&1 | tee "$work/${pair}-${mode}.log"
