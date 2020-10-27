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
        if (len(cur_parts) > 0) and (orf_pos[cur_orf[0]][0] - orf_pos[cur_parts[-1][0]][1] > 10000):
            possible_BGC.append(cur_parts)
            cur_parts = []
        cur_parts.append(cur_orf)

    if len(cur_parts) > 0:
        possible_BGC.append(cur_parts)

    return possible_BGC

def get_domain_list(filename):
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
    return orf_domain_list

def get_dist_list(orfs_list_short, orf_pos, filename):
    orf_domain_list = get_domain_list(filename)
    dist_list = []
    cur_parts = []
    def is_big_orf(bgc_name):
        if "A" in orf_domain_list[bgc_name] and "PCP" in orf_domain_list[bgc_name]:
            for dm in orf_domain_list[bgc_name]:
                if "C_" in dm:
                    return True
        return False

    def get_str(bgc_name):
        res = ""
        for dm in orf_domain_list[bgc_name]:
            res += dm + "-"
        return res

    comment = ""
    for cur_orf in orfs_list_short:
        comment += cur_orf[0] + ": " + get_str(cur_orf[0]) + "; "
        if True:#is_big_orf(cur_orf[0]):
            if len(cur_parts) > 0:
                dist_list.append(orf_pos[cur_orf[0]][0] - orf_pos[cur_parts[-1][0]][1])
                comment += "dist=" + str(orf_pos[cur_orf[0]][0] - orf_pos[cur_parts[-1][0]][1])
            cur_parts.append(cur_orf)
        comment += "\n"

    dist_list.sort()
    return (dist_list, comment)

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

def preprocess_cond_start(possible_BGC, filename, orf_ori, orf_pos):
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
            if orf_domain_list[part[0]][-1] == "C_Starter" and (orf_domain_list[part[0]][0] == "TE" or orf_domain_list[part[0]][0] == "TD"):
                return True
        else:
            if orf_domain_list[part[0]][0] == "C_Starter" and (orf_domain_list[part[0]][-1] == "TE"  or orf_domain_list[part[0]][-1] == "TD"):
                return True
        return False

    def is_unique_A(part):
        if (len(orf_domain_list[part[0]]) == 1 and (orf_domain_list[part[0]][-1] == "A")):
            return True
        if (len(orf_domain_list[part[0]]) == 2) and ("A" in orf_domain_list[part[0]]) and ("PCP" in orf_domain_list[part[0]]):
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
            elif is_unique_A(possBGC[i]):
                if cur_i != i:
                    sep.append(possBGC[cur_i:i])
                cur_i = i + 1

        if cur_i != len(possBGC):
            sep.append(possBGC[cur_i:])
        res += sep

    def end_by_TE_TD(part):
        return orf_domain_list[part[0]][-1] in ["TE", "TD"]

    def end_by_Starter(part):
        return orf_domain_list[part[0]][-1] == "C_Starter"

    def start_by_TE_TD(part):
        return orf_domain_list[part[0]][0] in ["TE", "TD"]

    def start_by_Starter(part):
        return orf_domain_list[part[0]][0] == "C_Starter"

    def has_AA(part):
        for cur_orf in orf_domain_list[part[0]]:
            if cur_orf[0] == "A":
                return True
        return False


    def get_next_AA_id(parts, lst_id, stp):
        while lst_id < len(parts) and lst_id >= 0:
            if has_AA(parts[lst_id]):
                return lst_id
            lst_id += stp
        return lst_id


    def split_by_first_starter_te(parts):
        def separate_first_st_te(parts):
            if (orf_ori[parts[0][0]] == '+') and (orf_domain_list[parts[0][0]][0] == "C_Starter"):
                lst_id = 0
                max_dist = 0
                while (lst_id < len(parts)) and (orf_ori[parts[lst_id][0]] == '+') and (not end_by_TE_TD(parts[lst_id])):
                    lst_id += 1
                    if lst_id < len(parts):
                        max_dist = max(max_dist, orf_pos[parts[lst_id][0]][0] - orf_pos[parts[lst_id - 1][0]][1])

                if (lst_id < len(parts)) and (orf_ori[parts[lst_id][0]] == '+') and (end_by_TE_TD(parts[lst_id])):
                    #nxt_id = get_next_AA_id(parts, lst_id + 1, 1)
                    #if (nxt_id < len(parts)) and ((orf_pos[parts[nxt_id][0]][0] - orf_pos[parts[lst_id][0]][1]) > 2 * max_dist):
                    return [parts[:lst_id + 1], parts[lst_id + 1:]]
            return [parts]

        def separate_first_te_st(parts):
            if (orf_ori[parts[0][0]] == '-') and (start_by_TE_TD(parts[0])):
                lst_id = 0
                max_dist = 0
                while (lst_id < len(parts)) and (orf_ori[parts[lst_id][0]] == '-') and (not end_by_Starter(parts[lst_id])):
                    lst_id += 1
                    if lst_id < len(parts):
                        max_dist = max(max_dist, orf_pos[parts[lst_id][0]][0] - orf_pos[parts[lst_id - 1][0]][1])

                if (lst_id < len(parts)) and (orf_ori[parts[lst_id][0]] == '-') and (end_by_Starter(parts[lst_id])):
                    #nxt_id = get_next_AA_id(parts, lst_id + 1, 1)
                    #if (nxt_id < len(parts)) and ((orf_pos[parts[nxt_id][0]][0] - orf_pos[parts[lst_id][0]][1]) > 2 * max_dist):
                    return [parts[:lst_id + 1], parts[lst_id + 1:]]
            return [parts]

        def del_orfs_without_A(parts):
            without_A_cnt = 0
            for part in parts:
                has_A = False
                for cur_orf in orf_domain_list[part[0]]:
                    if cur_orf[0] == "A":
                        has_A = True
                if has_A:
                    break
                without_A_cnt += 1
            return parts[without_A_cnt:]

        sepr = []
        cur_parts = parts
        can_sep = True
        while can_sep:
            cur_parts = del_orfs_without_A(cur_parts)
            if len(cur_parts) == 0:
                break
            cur_separate = separate_first_st_te(cur_parts)
            if len(cur_separate) == 1:
                cur_parts = cur_separate[0]
                cur_separate = separate_first_te_st(cur_parts)

            if len(cur_separate) == 2:
                sepr.append(cur_separate[0])
                cur_parts = cur_separate[1]
            else:
                sepr += cur_separate
                can_sep = False
        return sepr

    def split_by_last_starter_te(parts):
        def separate_last_st_te(parts):
            if (orf_ori[parts[-1][0]] == '-') and (orf_domain_list[parts[-1][0]][-1] == "C_Starter"):
                lst_id = len(parts) - 1
                max_dist = 0
                while (lst_id >= 0) and (orf_ori[parts[lst_id][0]] == '-') and (not start_by_TE_TD(parts[lst_id])):
                    lst_id -= 1
                    if lst_id >= 0:
                        max_dist = max(max_dist, orf_pos[parts[lst_id + 1][0]][0] - orf_pos[parts[lst_id][0]][1])

                if (lst_id >= 0) and (orf_ori[parts[lst_id][0]] == '-') and (start_by_TE_TD(parts[lst_id])):
                    #nxt_id = get_next_AA_id(parts, lst_id - 1, -1)
                    #if (lst_id > 0) and ((orf_pos[parts[lst_id][0]][0] - orf_pos[parts[nxt_id][0]][1]) > 2 * max_dist):
                    return [parts[:lst_id], parts[lst_id:]]
            return [parts]

        def separate_last_te_st(parts):
            if (orf_ori[parts[-1][0]] == '+') and (end_by_TE_TD(parts[-1])):
                lst_id = len(parts) - 1
                max_dist = 0
                while (lst_id >= 0) and (orf_ori[parts[lst_id][0]] == '+') and (not start_by_Starter(parts[lst_id])):
                    lst_id -= 1
                    if lst_id >= 0:
                        max_dist = max(max_dist, orf_pos[parts[lst_id + 1][0]][0] - orf_pos[parts[lst_id][0]][1])

                if (lst_id >= 0) and (orf_ori[parts[lst_id][0]] == '+') and (start_by_Starter(parts[lst_id])):
                    #nxt_id = get_next_AA_id(parts, lst_id - 1, -1)
                    #if (nxt_id >= 0) and ((orf_pos[parts[lst_id][0]][0] - orf_pos[parts[nxt_id][0]][1]) > 2 * max_dist):
                    return [parts[:lst_id], parts[lst_id:]]
            return [parts]

        def del_orfs_without_A(parts):
            without_A_cnt = 0
            for part in parts[::-1]:
                has_A = False
                for cur_orf in orf_domain_list[part[0]]:
                    if cur_orf[0] == "A":
                        has_A = True
                if has_A:
                    break
                without_A_cnt += 1
            return parts[:len(parts) - without_A_cnt]

        sepr = []
        cur_parts = parts
        can_sep = True
        while can_sep:
            cur_parts = del_orfs_without_A(cur_parts)
            if len(cur_parts) == 0:
                break
            cur_separate = separate_last_st_te(cur_parts)
            #print("Separate: ", cur_separate)
            if len(cur_separate) == 1:
                cur_parts = cur_separate[0]
                cur_separate = separate_last_te_st(cur_parts)

            #print("Separate2: ", cur_separate)
            if len(cur_separate) == 2:
                sepr = [cur_separate[1]] + sepr
                cur_parts = cur_separate[0]
            else:
                sepr = cur_separate + sepr
                can_sep = False
        return sepr


    res2 = []
    for possBGC in res:
        sep = split_by_first_starter_te(possBGC)
        res2 += sep

    res = []
    for possBGC in res2:
        sep = split_by_last_starter_te(possBGC)
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
            if ("A" not in orf_domain_list[parts[i][0]]):
                continue
            if et_ori == "?":
                et_ori = orf_ori[parts[i][0]]
            elif et_ori != orf_ori[parts[i][0]]:
                return False
        return True

    def print_BGC(parts):
        print("Filename: " + filename)
        for i in range(len(parts)):
            print(parts[i][0] + ":" + str(orf_domain_list[parts[i][0]]))

    def has_starter(parts):
        for i in range(len(parts)):
            for orf in orf_domain_list[parts[i][0]]:
                if orf == "C_Starter":
                    return True
        return False

    def starter_count(parts):
        cnt_starter = 0
        for i in range(len(parts)):
            for orf in orf_domain_list[parts[i][0]]:
                if orf == "C_Starter":
                    cnt_starter += 1
        return cnt_starter

    def wrong_starter_pos(parts):
        for i in range(len(parts)):
            for j in range(len(orf_domain_list[parts[i][0]])):
                orf = orf_domain_list[parts[i][0]][j]
                if orf == "C_Starter":
                    if not ((i == 0 and j == 0) or (i == len(parts) - 1 and j == len(orf_domain_list[parts[i][0]]) - 1)):
                        return True
        return False


    for possBGC in possible_BGC:
        if not has_starter(possBGC):
            continue

        if wrong_starter_pos(possBGC):
            print("C_Starteron wrong poss:")
            print_BGC(possBGC)
            return

        if not same_ori(possBGC):
            print("Different order:")
            print_BGC(possBGC)
            return

        if starter_count(possBGC) > 2:
            print("Too many C_Starter:")
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

            #dist_list, comment = get_dist_list(orfs_list_short, orf_pos, csv_file_with_orf)
            #print("Comment: ", comment)
            #if (len(dist_list) == 0):
            #    dist_list.append(-1)
            #with open("dist_stat.csv", 'a') as fa:
            #    fa.write(dirname.split('/')[-1] + "\t" + str(dist_list[0]) + "\t" + str(dist_list[-1]) + "\t\"" + comment + "\"\n")

            cur_posBGC = split_orfs_by_dist(orfs_list_short, orf_pos)
            cur_posBGC = preprocess_cond_start(cur_posBGC, csv_file_with_orf, orf_ori, orf_pos)
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
        #print(d_aa)
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