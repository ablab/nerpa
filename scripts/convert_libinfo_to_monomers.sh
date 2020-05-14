#!/bin/sh

set -e

# RBAN_PATH="/Data/Projects/CAB/Nerpa/rban/rBAN-1.0.jar"
SCRIPT_PATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
LIBINFO=$1
DB_PATH=$2
if [ "$#" -gt "2" ]
then
  MONOMERINFO=$3
else
  MONOMERINFO=$(basename -- "$LIBINFO").monomers
fi
RBAN_INP=$(basename -- "$LIBINFO").input.json
RBAN_OUT=$(basename -- "$LIBINFO").output.json

#echo $LIBINFO
#echo $DB_PATH
#echo $MONOMERINFO
#echo $RBAN_INP
#echo $RBAN_OUT
#echo $RBAN_PATH

# 1. make rBAN input
python3 "$SCRIPT_PATH"/convert_libinfo_to_rban_input.py "$LIBINFO" "$DB_PATH" > "$RBAN_INP"

# 2. run rBAN
java -jar "$RBAN_PATH" -inputFile "$RBAN_INP" -outputFolder ./ -outputFileName "$RBAN_OUT"

# 3. make info.monomers
python3 "$SCRIPT_PATH"/convert_rban_to_monomerinfo.py "$RBAN_OUT" "$RBAN_INP" > "$MONOMERINFO"
