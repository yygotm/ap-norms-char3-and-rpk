#!/bin/bash
# Certificates and checks for q=3 (n=6,7,8,10,12) and (9,6). Output appended to cert.log
# Usage: bash run_q3_cert.sh   (needs ./apnorm and ./checkcert built in $WORK, default ~/ap31)
WORK=${WORK:-$HOME/ap31}
cd "$WORK" || exit 1
for kn in "1 6" "1 7" "1 8" "1 10" "1 12" "2 6"; do
  set -- $kn; q=$((3**$1)); d=${q}_$2; mkdir -p $d; (cd $d && ../apnorm $1 $2 2 0 . 2>/dev/null)
  ./checkcert $d 2 | tee -a cert.log
done
