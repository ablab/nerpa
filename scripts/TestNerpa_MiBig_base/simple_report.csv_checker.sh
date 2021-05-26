#!/bin/bash

report=$1
echo "specified report is $report"
if [ ! -f "$report" ]; then
  echo "Provide valid report.csv file as the first argument"
  exit 1
fi

all_matches=`cut -d',' -f1,6,7 $report | sed '1d;' | sed 's@/.*/@@g' | sed 's@,@	@g' | sort -r -g -k1,1 | sed 's@[\._][a-zA-Z0-9_]*@@g' | cut -f2,3 | awk -F"\t" '!_[$1]++' | wc -l`
correct_matches=`cut -d',' -f1,6,7 $report | sed '1d;' | sed 's@/.*/@@g' | sed 's@,@	@g' | sort -r -g -k1,1 | sed 's@[\._][a-zA-Z0-9_]*@@g' | cut -f2,3 | awk -F"\t" '!_[$1]++' | awk -F"\t" '$1==$2' | wc -l`
echo "${correct_matches} matches with correct BGC out of ${all_matches} total best matches per structure"
