DATE=`date +"%d-%b-%Y_%H:%M"`
#DATE='04-Mar-2021_16:16'
RES_DIR=$TEST_MIBIG_RESULT_PATH/result/res_$DATE

mkdir -p $RES_DIR

START_TIME=$(date +%s)
nerpa.py -a $TEST_MIBIG_DATA_PATH/prediction/ --smiles-tsv $TEST_MIBIG_DATA_PATH/base.csv --col_smiles "SMILES" --col_id "Accession" --sep "," -o $RES_DIR
END_TIME=$(date +%s)
DIFF_TIME=$(($END_TIME - $START_TIME))

echo "It tooks $DIFF_TIME seconds"
echo "It tooks $DIFF_TIME seconds" > $RES_DIR/running_time

SOURCE="$( dirname ${BASH_SOURCE[0]} )"

python3 $SOURCE/test_nerpa.py $RES_DIR $TEST_MIBIG_DATA_PATH
python3 $SOURCE/test_nerpa_details.py $RES_DIR $TEST_MIBIG_DATA_PATH
python3 $SOURCE/merge_details_with_base.py $RES_DIR $TEST_MIBIG_DATA_PATH
