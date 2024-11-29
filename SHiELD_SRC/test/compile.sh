#!/bin/bash
# Script to compile the code
# Check if three command-line arguments are provided
if [ "$#" -ne 1 ]; then
  echo "Usage: ./compile.sh COMP"
  echo "choices:  debug, repro, prod"
  exit 1
fi
COMP="$1"

cd ../../SHiELD_build/Build
./COMPILE solo sw $COMP gnu 64bit noclean

# Move the log file to test dir
mv build_solo_sw.$COMP.64bit.gnu.out ../../SHiELD_SRC/test/build_solo_sw.$COMP.64bit.gnu.out

