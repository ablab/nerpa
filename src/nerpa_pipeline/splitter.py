#!/usr/bin/env python

import os
import csv

def split_by_dist(parts, orf_pos):
    possible_BGC = []
    for BGC in parts:
        cur_parts = []
        for cur_orf in BGC:
            if (len(cur_parts) > 0) and (orf_pos[cur_orf][0] - orf_pos[cur_parts[-1]][1] > 10000):
                possible_BGC.append(cur_parts)
                cur_parts = []
            cur_parts.append(cur_orf)

        if len(cur_parts) > 0:
            possible_BGC.append(cur_parts)

    return possible_BGC

def end_by_TE_TD(orf, orf_domains, orf_ori):
    if orf_ori[orf] == '+' and (orf_domains[orf][-1] == "TE" or orf_domains[orf][-1] == "TD"):
        return True
    if orf_ori[orf] == '-' and (orf_domains[orf][0] == "TE" or orf_domains[orf][0] == "TD"):
        return True
    return False


def start_by_Starter(orf, orf_domains, orf_ori):
    if orf_ori[orf] == '+' and (orf_domains[orf][0] == "Starter"):
        return True
    if orf_ori[orf] == '-' and (orf_domains[orf][-1] == "Starter"):
        return True
    return False


def split_by_one_orf_Starter_TE(BGCs, orf_ori, orf_domains):
    res_bgcs = []
    for BGC in BGCs:
        cur_i = 0
        for i in range(len(BGC)):
            orf = BGC[i]
            if end_by_TE_TD(orf, orf_domains, orf_ori) and start_by_Starter(orf, orf_domains, orf_ori):
                if i != cur_i:
                    res_bgcs.append(BGC[cur_i: i])
                res_bgcs.append(BGC[i])
                cur_i = i + 1
        if cur_i != len(BGC):
            res_bgcs.append(BGC[cur_i: len(BGC)])

    return res_bgcs


def split_by_single_domain_orf(BGCs, orf_ori, orf_domains):
    def is_removable(orf):
        if (len(orf_domains[orf]) == 2) and ("A" in orf_domains[orf]) and ("PCP" in orf_domains[orf]):
            return True
        if (len(orf_domains[orf]) == 1):
            return  True
        return False

    res_bgcs = []
    for BGC in BGCs:
        cur_i = 0
        for i in range(len(BGC)):
            orf = BGC[i]
            if is_removable(orf):
                if i != cur_i:
                    res_bgcs.append(BGC[cur_i:i])
                cur_i = i + 1
        if cur_i != len(BGC):
            res_bgcs.append(BGC[cur_i:len(BGC)])

    return res_bgcs