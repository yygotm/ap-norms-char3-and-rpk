#!/bin/bash
# Check the (27,6) certificate with checkcert; result appended to cert.log in $WORK (default ~/ap31)
WORK=${WORK:-$HOME/ap31}
cd "$WORK" && /usr/bin/time -f "%e s %M KB" ./checkcert 27_6 16 2>&1 | tee -a cert.log
