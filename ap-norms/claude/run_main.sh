#!/bin/bash
# Production: (q,n) = (9,8) and (27,6) with certificates (large: about 88 MB and 5.9 GB).
# Usage: bash ap-norms/claude/run_main.sh [test]   (WSL/Linux, gcc with OpenMP)
# Certificates go to $OUT/<q>_<n>/ (default ~/ap31); logs go to $OUT/logs/.
SRC="$(cd "$(dirname "$0")" && pwd)"
OUT=${OUT:-$HOME/ap31}
LOG=$OUT/logs; mkdir -p "$OUT" "$LOG"
cp "$SRC/apnorm.c" "$OUT/" && cd "$OUT" || exit 1
gcc -O2 -march=native -fopenmp -o apnorm apnorm.c || exit 1
run() {  # k n
  local q=$((3**$1)); local d=$OUT/${q}_$2; mkdir -p "$d"; cd "$d"
  /usr/bin/time -v "$OUT/apnorm" $1 $2 16 0 "$d" > result.txt 2> stderr.txt
  echo "exit=$? $(cat result.txt)" >> "$LOG/main.log"
  grep -E "Maximum resident|Elapsed" stderr.txt >> "$LOG/main.log"
}
if [ "$1" = "test" ]; then run 2 6; exit; fi
echo "start $(date '+%F %T')" >> "$LOG/main.log"
run 2 8
run 3 6
echo "end $(date '+%F %T')" >> "$LOG/main.log"
