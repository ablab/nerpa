DATE=`date +"%d-%b-%Y_%H:%M"`

nerpa.py -p prediction.info -l structure.info --predictor NRPSPREDICTOR2 --insertion --deletion --modification -o ./result/res_$DATE
python3 ./test_nerpa.py res_$DATE
