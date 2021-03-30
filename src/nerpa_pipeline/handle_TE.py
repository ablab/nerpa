#!/usr/bin/env python

import os
import csv

import handle_PCP2
import handle_MT
import handle_E
import handle_helper


def get_split_BGC(dirname):
    possible_BGC = []
    txt_folder = os.path.join(dirname, "txt")
    if  not os.path.isdir(txt_folder):
        return possible_BGC

    orf_ori = handle_helper.get_orf_orientation(dirname)
    orf_pos = handle_helper.get_orf_position(dirname)

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
                if (len(cur_parts) > 0) and (orf_pos[cur_orf[0]][0] - orf_pos[cur_parts[-1]][1] > 20000):
                    possible_BGC.append(cur_parts)
                    cur_parts = []

                cur_parts.append(cur_orf[0])
                if cur_orf[-1] == 1:
                    if orf_ori[cur_orf[0]] == '+':
                        possible_BGC.append(cur_parts)
                        cur_parts = []
                    else:
                        cur_parts = cur_parts[:-1]
                        if len(cur_parts) > 0:
                            possible_BGC.append(cur_parts)
                        cur_parts = [cur_orf[0]]

            if len(cur_parts) > 0:
                if len(possible_BGC) > 0:
                    possible_BGC.append(cur_parts + possible_BGC[-1])
                possible_BGC.append(cur_parts)

    return possible_BGC


def reverse_prediction(pred_str):
    pred_str = '\n'.join(pred_str.split('\n')[::-1])
    pred_str = pred_str.lstrip('\n')
    pred_str += '\n'
    return pred_str

  
def gen_prediction_for_one_orfs_part(orf_part, prediction_dict, output_prefix, current_part, predictions_info_list, orf_ori):
    has_plus = False
    output_str = ""
    for current_orf in orf_part:
        if orf_ori[current_orf] == "+":
            has_plus = True

    orf_order = 1 if has_plus else -1
    for current_orf in orf_part[::orf_order]:
        if current_orf in prediction_dict:
            if orf_ori[current_orf] == "-":
                output_str += reverse_prediction(prediction_dict[current_orf])
            else:
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
        current_part = gen_prediction_for_one_orfs_part(orf_part, prediction_dict, output_prefix, current_part, predictions_info_list, orf_ori)
    return current_part


def create_predictions_by_antiSMASHout(antismashouts, outdir):
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
        # print(d_aa)  # TODO: print in debug mode
        bgc_orfs_parts = get_split_BGC(dirname)

        nrpspred_dir = os.path.join(dirname, "nrpspks_predictions_txt")
        if os.path.isdir(nrpspred_dir):
            for filename in os.listdir(nrpspred_dir):
                if filename.endswith('nrpspredictor2_codes.txt'):
                    base_antismashout_name = os.path.basename(dirname)
                    base_pred_name = os.path.basename(filename)
                    # predictions_info_list.append(os.path.join(dir_for_predictions, base_antismashout_name + "_" + base_pred_name))
                    # shutil.copyfile(os.path.join(nrpspred_dir, filename), os.path.join(dir_for_predictions, base_antismashout_name + "_" + base_pred_name))
                    gen_predictions(bgc_orfs_parts, os.path.join(nrpspred_dir, filename),
                                    os.path.join(dir_for_predictions, base_antismashout_name + "_" + base_pred_name)[:-4],
                                    0, predictions_info_list, double_orf, double_aa, mt_aa, d_aa, dirname)

    with open(predictions_info_file, 'w') as f:
        for line in predictions_info_list:
            f.write(f'{line}\n')

    return predictions_info_file
