#DATE=`date +"%d-%b-%Y_%H:%M"`
DATE='03-Jun-2021_10:54'
RES_DIR=$TEST_MIBIG_RESULT_PATH/result/res_$DATE

START_TIME=$(date +%s)
#nerpa.py -a $TEST_MIBIG_DATA_PATH/prediction/ --smiles-tsv $TEST_MIBIG_DATA_PATH/base.csv --col-smiles "SMILES" --col-id "Accession" --sep ","  -o $RES_DIR

END_TIME=$(date +%s)
DIFF_TIME=$(($END_TIME - $START_TIME))

echo "It tooks $DIFF_TIME seconds"
echo "It tooks $DIFF_TIME seconds" > $RES_DIR/running_time

SOURCE="$( dirname ${BASH_SOURCE[0]} )"

python3 $SOURCE/test_nerpa.py $RES_DIR $TEST_MIBIG_DATA_PATH
python3 $SOURCE/extract_true_align_table.py $RES_DIR $TEST_MIBIG_DATA_PATH
