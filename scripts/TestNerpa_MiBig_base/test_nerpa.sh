DATE=`date +"%d-%b-%Y_%H:%M"`

mkdir -p $TEST_MIBIG_RESULT_PATH/result/res_$DATE

START_TIME=$(date +%s)
nerpa.py -a $TEST_MIBIG_DATA_PATH/prediction/ --smiles-tsv $TEST_MIBIG_DATA_PATH/mibig2019NRP.csv --col_smiles "SMILES" --col_id "Accession" --sep "," -o "$TEST_MIBIG_RESULT_PATH/result/res_$DATE"
END_TIME=$(date +%s)
DIFF_TIME=$(($END_TIME - $START_TIME))

echo "It tooks $DIFF_TIME seconds"
echo "It tooks $DIFF_TIME seconds" > $TEST_MIBIG_RESULT_PATH/result/res_$DATE/running_time

SOURCE="$( dirname ${BASH_SOURCE[0]} )"

python3 $SOURCE/test_nerpa.py res_$DATE
python3 $SOURCE/test_nerpa_details.py res_$DATE
python3 $SOURCE/merge_details_with_base.py res_$DATE
