#!/bin/bash

OUT=mibig_2019_new_scoring  # $1
JSONS=local_test_data/golden_79/list_jsons.txt 
SUMMARY=$OUT/summary.tsv

#rm -rf $OUT
mkdir -p $OUT
echo "BGC	CONTIG	ORF	A-ID	PRED_TOP5" > $SUMMARY

for id in `cat local_test_data/golden_79/list_ids_google_spreadsheet.txt`; do
	json=`grep $id $JSONS`;
	codes=$OUT/$id/nrpspks_predictions_txt/ctg1_nrpspredictor2_codes.txt
	if [ ! -f $codes ]; then
		./src/nerpa_pipeline/NRPSPredictor_utils/main.py -o $OUT/$id $json -m hybrid
	fi
	cat $codes | cut -f1,3 | cut -f1,2,3,4,5 -d ';' | gsed 's/_/\t/g' | gsed 's/^/'$id'\t/g' >> $SUMMARY
done

echo "Done! See results in"
echo $SUMMARY
