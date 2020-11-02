#!/usr/bin/env python
import sys
import os
import shutil
import csv

import handle_PCP2
import handle_MT
import handle_E
import splitter
import handle_helper
from logger import log

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
                    predictions_info_list, double_orf, double_aa, mt_aa, d_aa, base_antismashout_name):
    prediction_dict = {}
    with open(input_file_name, 'r') as rf:
        for line in rf:
            ctgorf = '_'.join(line.split('_')[:2])

            if line.split("\t")[0] in double_aa:
                line = line.split("\t")[0] + "*\t" + '\t'.join(line.split("\t")[1:])

            if line.split("\t")[0] in mt_aa:
                line = line.split("\t")[0] + "\t" + line.split("\t")[1] + "+MT\t"  + '\t'.join(line.split("\t")[2:])

            if line.split("\t")[0] in d_aa:
                line = line.split("\t")[0] + "\td-" + line.split("\t")[1] + "\t"  + '\t'.join(line.split("\t")[2:])

            if ctgorf in double_orf:
                line = ctgorf + "*_" + '_'.join(line.split('_')[2:])

            if (ctgorf in prediction_dict):
                prediction_dict[ctgorf] += line
            else:
                prediction_dict[ctgorf] = line

    orf_ori = handle_helper.get_orf_orientation(base_antismashout_name)
    for cur_orf in orf_ori.keys():
        if orf_ori[cur_orf] == '-' and (cur_orf in prediction_dict):
            prediction_dict[cur_orf] = '\n'.join(prediction_dict[cur_orf].split('\n')[::-1])
            if prediction_dict[cur_orf][0] == '\n':
                prediction_dict[cur_orf] = prediction_dict[cur_orf][1:]
            prediction_dict[cur_orf] += '\n'

    for orf_part in bgc_orfs_parts:
        current_part = gen_prediction_for_one_orfs_part(orf_part, prediction_dict, output_prefix, current_part, predictions_info_list)
    return current_part

def create_predictions_by_antiSAMSHout(antismashouts, outdir):
    log.log("Start create predictions by antiSMASH")

    dir_for_predictions = os.path.join(outdir, "predictions")
    if not os.path.exists(dir_for_predictions):
        os.makedirs(dir_for_predictions)

    predictions_info_file = os.path.join(outdir, "predictions.info")
    predictions_info_list = []
    for dirname in antismashouts:
        if dirname[-1] == '\n':
            dirname = dirname[:-1]

        double_orf, double_aa = handle_PCP2.get_double_orfs_and_AA(dirname)
        mt_aa = handle_MT.get_MT_AA(dirname)
        d_aa = handle_E.get_D_AA(dirname)

        orf_pos = handle_helper.get_orf_position(dirname)
        orf_ori = handle_helper.get_orf_orientation(dirname)
        orf_domains = handle_helper.get_orf_domain_list(dirname)

        print("====PARTS BEFORE: ")
        parts = handle_helper.get_parts(dirname)
        handle_helper.debug_print_parts(dirname, parts, orf_domains, orf_ori, orf_pos)

        print("====SPLIT BY DIST:")
        parts = splitter.split_by_dist(parts, orf_pos)
        handle_helper.debug_print_parts(dirname, parts, orf_domains, orf_ori, orf_pos)

        print("====SPLIT BY SINGLE ORF WITH Starter-TE")
        parts = splitter.split_by_one_orf_Starter_TE(parts, orf_ori, orf_domains)
        handle_helper.debug_print_parts(dirname, parts, orf_domains, orf_ori, orf_pos)

        print("====REMOVE SINGLE DOMAINs ORFS")
        parts = splitter.split_by_single_domain_orf(parts, orf_ori, orf_domains)
        handle_helper.debug_print_parts(dirname, parts, orf_domains, orf_ori, orf_pos)

        print("====SPLIT AND REORDER")
        parts = splitter.split_and_reorder(parts, orf_ori, orf_pos, orf_domains)
        handle_helper.debug_print_parts(dirname, parts, orf_domains, orf_ori, orf_pos)

        nrpspred_dir = os.path.join(dirname, "nrpspks_predictions_txt")
        if os.path.isdir(nrpspred_dir):
            for filename in os.listdir(nrpspred_dir):
                if filename.endswith('nrpspredictor2_codes.txt'):
                    base_antismashout_name = os.path.basename(dirname)
                    base_pred_name = os.path.basename(filename)
                    #predictions_info_list.append(os.path.join(dir_for_predictions, base_antismashout_name + "_" + base_pred_name))
                    #shutil.copyfile(os.path.join(nrpspred_dir, filename), os.path.join(dir_for_predictions, base_antismashout_name + "_" + base_pred_name))
                    gen_predictions(parts, os.path.join(nrpspred_dir, filename),
                                    os.path.join(dir_for_predictions, base_antismashout_name + "_" + base_pred_name)[:-4],
                                    0, predictions_info_list, double_orf, double_aa, mt_aa, d_aa, dirname)

    f = open(predictions_info_file, 'w')
    for line in predictions_info_list:
        f.write(line + "\n")
    f.close()

    return predictions_info_file