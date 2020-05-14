OUT='res_13-May-2020_20:34'

python3 ./test_nerpa.py $OUT
python3 ./test_nerpa_details.py $OUT
python3 ./merge_details_with_base.py $OUT
