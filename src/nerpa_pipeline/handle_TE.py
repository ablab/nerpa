#!/usr/bin/env python

import os
import csv

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
                if (len(cur_parts) > 0) and (orf_pos[cur_orf[0]][0] - orf_pos[cur_parts[-1]][1] > 10000):
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