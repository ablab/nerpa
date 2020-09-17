#!/usr/bin/env python
import sys
import os
import shutil
import csv

import handle_PCP2
import handle_MT
import handle_E
import handle_helper
from logger import log

def split_orfs_by_dist(orfs_list_short, orf_pos):
    possible_BGC = []
    cur_parts = []
    for cur_orf in orfs_list_short:
        if (len(cur_parts) > 0) and (orf_pos[cur_orf[0]][0] - orf_pos[cur_parts[-1][0]][1] > 20000):
            possible_BGC.append(cur_parts)
            cur_parts = []
        cur_parts.append(cur_orf)

    if len(cur_parts) > 0:
        possible_BGC.append(cur_parts)

    return possible_BGC


def split_orfs_by_TE(possible_BGC, orf_ori):
    res = []

    cur_parts = []
    for pos_BGC in possible_BGC:
        for cur_orf in pos_BGC:
            cur_parts.append(cur_orf[0])
            if cur_orf[-1] == 1:
                if orf_ori[cur_orf[0]] == '+':
                    res.append(cur_parts)
                    cur_parts = []
                else:
                    cur_parts = cur_parts[:-1]
                    if len(cur_parts) > 0:
                        res.append(cur_parts)
                    cur_parts = [cur_orf[0]]

            if len(cur_parts) > 0:
                res.append(cur_parts)

    return res

def preprocess_cond_start(possible_BGC, filename, orf_ori):
    orf_domain_list = {}
    with open(filename, 'r') as rf:
        csv_reader = csv.reader(rf, delimiter='\t')
        for row in csv_reader:
            if row[1] == "NRPSPKS_ID":
                continue

            if row[1] not in orf_domain_list:
                orf_domain_list[row[1]] = []

            if row[6] == 'Condensation':
                if row[7] == "Condensation_Starter":
                    orf_domain_list[row[1]].append('C_Starter')
                else:
                    orf_domain_list[row[1]].append('C_' + row[7].split('_')[-1])
            elif row[6] == 'Thioesterase':
                orf_domain_list[row[1]].append("TE")
            elif row[6] == "AMP-binding":
                orf_domain_list[row[1]].append("A")
            elif row[6] == "Epimerization":
                orf_domain_list[row[1]].append("E")
            else:
                orf_domain_list[row[1]].append(row[6])

    def is_Starter_TE(part):
        if orf_ori[part[0]] == '-':
            if orf_domain_list[part[0]][-1] == "C_Starter" and orf_domain_list[part[0]][0] == "TE":
                return True
        else:
            if orf_domain_list[part[0]][0] == "C_Starter" and orf_domain_list[part[0]][-1] == "TE":
                return True
        return False

    res = []
    for possBGC in possible_BGC:
        sep = []
        cur_i = 0
        for i in range(len(possBGC)):
            if is_Starter_TE(possBGC[i]):
                if cur_i != i:
                    sep.append(possBGC[cur_i:i])
                sep.append([possBGC[i]])
                cur_i = i + 1
        if cur_i != len(possBGC):
            sep.append(possBGC[cur_i:])
        res += sep

    return res


def check_cond_starter(possible_BGC, filename, orf_ori):
    orf_domain_list = {}
    with open(filename, 'r') as rf:
        csv_reader = csv.reader(rf, delimiter='\t')
        for row in csv_reader:
            if row[1] == "NRPSPKS_ID":
                continue

            if row[1] not in orf_domain_list:
                orf_domain_list[row[1]] = []

            if row[6] == 'Condensation':
                if row[7] == "Condensation_Starter":
                    orf_domain_list[row[1]].append('C_Starter')
                else:
                    orf_domain_list[row[1]].append('C_' + row[7].split('_')[-1])
            elif row[6] == 'Thioesterase':
                orf_domain_list[row[1]].append("TE")
            elif row[6] == "AMP-binding":
                orf_domain_list[row[1]].append("A")
            elif row[6] == "Epimerization":
                orf_domain_list[row[1]].append("E")
            else:
                orf_domain_list[row[1]].append(row[6])

    def same_ori(parts):
        et_ori = "?"
        for i in range(len(parts)):
            if et_ori == "?":
                et_ori = orf_ori[parts[i][0]]
            elif et_ori != orf_ori[parts[i][0]]:
                return False
        return True

    def print_BGC(parts):
        print("Filename: " + filename)
        for i in range(len(parts)):
            print(parts[i][0] + ":" + str(orf_domain_list[parts[i][0]]))

    for possBGC in possible_BGC:
        if not same_ori(possBGC):
            print("Different order:")
            print_BGC(possBGC)
            return


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

            cur_posBGC = split_orfs_by_dist(orfs_list_short, orf_pos)
            cur_posBGC = preprocess_cond_start(cur_posBGC, csv_file_with_orf, orf_ori)
            check_cond_starter(cur_posBGC, csv_file_with_orf, orf_ori)
            possible_BGC += split_orfs_by_TE(cur_posBGC, orf_ori)
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
        print(d_aa)
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
                                    0, predictions_info_list, double_orf, double_aa, mt_aa, d_aa, dirname)

    f = open(predictions_info_file, 'w')
    for line in predictions_info_list:
        f.write(line + "\n")
    f.close()

    return predictions_info_file