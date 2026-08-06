#!/bin/bash
# build.sh — build the sim-outorder simulator from ss3/.
#
# Usage:   scripts/build.sh
# Output:  simulator/ss3/sim-outorder
#
# Notes:
#  - sysprobe must be built first: the Makefile uses it to derive endianness
#    defines, and a stale/missing sysprobe silently drops them.
#  - libexo's sources are shipped pre-generated; the touch sequence below
#    keeps make from invoking flex to regenerate them.
#  - glibc 2.42+ removed <termio.h>; ss3/compat/ provides a shim.  The build
#    adds -Icompat only when the system header is absent (on older glibc the
#    shim would conflict with the system definition of struct termio).
set -e
cd "$(dirname "$0")/../simulator/ss3"

COMPAT=""
if [ ! -e /usr/include/termio.h ]; then
    COMPAT="-Icompat"
fi

gcc -O0 -o sysprobe sysprobe.c

find libexo -name '*.l' -exec touch {} + 2>/dev/null || true
sleep 1
find libexo -name '*.c' -exec touch {} +
find libexo -name '*.h' -exec touch {} +

make sim-outorder OFLAGS="-O0 -g -Wall -std=gnu99 $COMPAT \
    -Wno-implicit-function-declaration -Wno-implicit-int \
    -Wno-error=incompatible-pointer-types -Wno-error=int-conversion"

ls -la sim-outorder
echo "build complete: simulator/ss3/sim-outorder"
