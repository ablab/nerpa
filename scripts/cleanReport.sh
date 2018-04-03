#!/bin/sh

sed -i -e "s/graphs\///" $1
sed -i -e "s/.gr//" $1 
sed -i -e "s/\/media\/hosein\/My Passport\/hosein\/Desktop\/project\/sequence_data\/bacteria_complete\/antismash\///" $1
sed -i -e "s/\/nrpspks_predictions_txt\/ctg1_nrpspredictor2_codes.txt//" $1

