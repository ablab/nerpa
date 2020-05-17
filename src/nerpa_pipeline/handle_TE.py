#!/usr/bin/env python
import sys
import os
import shutil
import csv

import handle_PCP2
import handle_MT
from logger import log

def get_split_BGC(dirname):
    possible_BGC = []
    txt_folder = os.path.join(dirname, "txt")
    for filename in os.listdir(txt_folder):
        if filename.endswith("_NRPS_PKS.txt"):
            orfs_list = []
            csv_file_with_orf = os.path.join(txt_folder, filename)
            with open(csv_file_with_orf, 'r') as rf:
                csv_reader = csv.reader(rf, delimiter='\t')
                for row in csv_reader:
                    if row[1] == "NRPSPKS_ID":
                        continue

                    if row[6] == 'Thioesterase':
                        orfs_list.append([row[1], 1])
                    else:
                        orfs_list.append([row[1], 0])
            orfs_list_short = []
            for cur_orf in orfs_list:
                if len(orfs_list_short) == 0 or orfs_list_short[-1][0] != cur_orf[0]:
                    orfs_list_short.append(cur_orf)
                else:
                    orfs_list_short[-1][1] = max(orfs_list_short[-1][-1], cur_orf[-1])

            cur_parts = []
            for cur_orf in orfs_list_short:
                cur_parts.append(cur_orf[0])
                if cur_orf[-1] == 1:
                    possible_BGC.append(cur_parts)
                    cur_parts = []
            if len(cur_parts) > 0:
                if len(possible_BGC) > 0:
                    possible_BGC.append(cur_parts + possible_BGC[-1])
                possible_BGC.append(cur_parts)
    return possible_BGC

def gen_prediction_for_one_orfs_part(orf_part, prediction_dict, output_prefix, current_part, predictions_info_list):
    output_str = ""
    for current_orf in orf_part:
        if current_orf in prediction_dict:
            output_str += prediction_dict[current_orf]

    if (output_str != ""):
        output_file = output_prefix + "_part" + str(current_part)
        with open(output_file, 'w') as wf:
            wf.write(output_str)
        current_part += 1
        predictions_info_list.append(output_file)

    return current_part

def gen_predictions(bgc_orfs_parts, input_file_name, output_prefix, current_part,
                    predictions_info_list, double_orf, double_aa, mt_aa):
    prediction_dict = {}
    with open(input_file_name, 'r') as rf:
        for line in rf:
            ctgorf = '_'.join(line.split('_')[:2])

            if line.split("\t")[0] in double_aa:
                line = line.split("\t")[0] + "*\t" + '\t'.join(line.split("\t")[1:])

            if line.split("\t")[0] in mt_aa:
                line = line.split("\t")[0] + "\t" + line.split("\t")[1] + "+MT\t"  + '\t'.join(line.split("\t")[2:])

            if ctgorf in double_orf:
                line = ctgorf + "*_" + '_'.join(line.split('_')[2:])

            if (ctgorf in prediction_dict):
                prediction_dict[ctgorf] += line
            else:
                prediction_dict[ctgorf] = line

    for orf_part in bgc_orfs_parts:
        current_part = gen_prediction_for_one_orfs_part(orf_part, prediction_dict, output_prefix, current_part, predictions_info_list)
    return current_part

def create_predictions_by_antiSAMSHout(path_to_antismashouts, outdir, predictor):
    if predictor != "NRPSPREDICTOR2":
        log.err("You can provide antiSMASH output only for NRPSPREDICTOR2!")
        sys.exit()

    dir_for_predictions = os.path.join(outdir, "predictions")
    if not os.path.exists(dir_for_predictions):
        os.makedirs(dir_for_predictions)

    predictions_info_file = os.path.join(outdir, "predictions.info")
    predictions_info_list = []
    with open(path_to_antismashouts) as fr:
        for dirname in fr:
            if dirname[-1] == '\n':
                dirname = dirname[:-1]

            double_orf, double_aa = handle_PCP2.get_double_orfs_and_AA(dirname)
            mt_aa = handle_MT.get_MT_AA(dirname)
            print(mt_aa)
            bgc_orfs_parts = get_split_BGC(dirname)

            nrpspred_dir = os.path.join(dirname, "nrpspks_predictions_txt")
            if os.path.isdir(nrpspred_dir):
                for filename in os.listdir(nrpspred_dir):
                    if filename.endswith('nrpspredictor2_codes.txt'):
                        base_antismashout_name = os.path.basename(dirname)
                        base_pred_name = os.path.basename(filename)
                        #predictions_info_list.append(os.path.join(dir_for_predictions, base_antismashout_name + "_" + base_pred_name))
                        #shutil.copyfile(os.path.join(nrpspred_dir, filename), os.path.join(dir_for_predictions, base_antismashout_name + "_" + base_pred_name))
                        gen_predictions(bgc_orfs_parts, os.path.join(nrpspred_dir, filename),
                                        os.path.join(dir_for_predictions, base_antismashout_name + "_" + base_pred_name)[:-4],
                                        0, predictions_info_list, double_orf, double_aa, mt_aa)

    f = open(predictions_info_file, 'w')
    for line in predictions_info_list:
        f.write(line + "\n")
    f.close()

    return predictions_info_file

