DATE=`date +"%d-%b-%Y_%H:%M"`
START_TIME=$(date +%s)
nerpa.py -p prediction.info -l structure.info --predictor NRPSPREDICTOR2 --insertion --deletion --modification -o ./result/res_$DATE
END_TIME=$(date +%s)
DIFF_TIME=$(($END_TIME - $START_TIME))

echo "It tooks $DIFF_TIME seconds"
echo "It tooks $DIFF_TIME seconds" > ./result/res_$DATE/running_time

python3 ./test_nerpa.py res_$DATE
python3 ./test_nerpa_details.py res_$DATE
python3 ./merge_details_with_base.py res_$DATE
