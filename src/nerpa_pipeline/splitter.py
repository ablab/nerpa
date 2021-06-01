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
    if orf_ori[orf] == '+' and (orf_domains[orf][0] == "C_Starter"):
        return True
    if orf_ori[orf] == '-' and (orf_domains[orf][-1] == "C_Starter"):
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
                res_bgcs.append([BGC[i]])
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

def C_starter_count(lft, rgh, BGC, orf_domains):
    counter = 0
    for i in range(lft, rgh):
        if "C_Starter" in orf_domains[BGC[i]]:
            counter += 1
    return counter


def TE_TD_count(lft, rgh, BGC, orf_domains):
    counter = 0
    for i in range(lft, rgh):
        if ("TE" in orf_domains[BGC[i]]) or ("TD" in orf_domains[BGC[i]]):
            counter += 1
    return counter


def has_A_domain(BGC, orf_domains):
    for orf in BGC:
        if "A" in orf_domains[orf]:
            return True
    return False


def is_correct_NCterm(BGC, orf_ori, orf_domains):
    if len(BGC) == 0:
        return True
    if orf_ori[BGC[0]] == '+' and orf_domains[BGC[0]][0] == "NRPS-COM_Nterm":
        return False
    if orf_ori[BGC[0]] == '-' and orf_domains[BGC[0]][-1] == "NRPS-COM_Nterm":
        return False
    if orf_ori[BGC[-1]] == '+' and orf_domains[BGC[-1]][-1] == "NRPS-COM_Cterm":
        return False
    if orf_ori[BGC[-1]] == '-' and orf_domains[BGC[-1]][0] == "NRPS-COM_Cterm":
        return False
    return True


def is_correct(BGC, orf_ori, orf_pos, orf_domains):
    if not is_correct_NCterm(BGC, orf_ori, orf_domains):
        return False

    if C_starter_count(0, len(BGC), BGC, orf_domains) > 1:
        return False

    if TE_TD_count(0, len(BGC), BGC, orf_domains) > 1:
        return False

    if C_starter_count(0, len(BGC), BGC, orf_domains) == 1 and (not start_by_Starter(BGC[0], orf_domains, orf_ori)):
        return False

    if TE_TD_count(0, len(BGC), BGC, orf_domains) == 1 and (not end_by_TE_TD(BGC[-1], orf_domains, orf_ori)):
        return False

    if not has_A_domain(BGC, orf_domains):
        return False

    return True


def is_correct_orfs_subseq(lft, rgh, BGC, orf_ori, orf_pos, orf_domains):
    if C_starter_count(lft, rgh + 1, BGC, orf_domains) > 1:
        return False

    if TE_TD_count(lft, rgh + 1, BGC, orf_domains) > 1:
        return False

    return True

def reverse_neg(BGC, orf_ori):
    for orf in BGC:
        if orf_ori[orf] == '+':
            return BGC

    return list(reversed(BGC))

def generate_all_perm(BGC, orf_ori, orf_pos, orf_domains):
    import itertools
    Perm_BGC = []
    row_perm_BGC = list(itertools.permutations(BGC))
    for bgc in row_perm_BGC:
        if is_correct(bgc, orf_ori, orf_pos, orf_domains):
            Perm_BGC.append(bgc)
    return Perm_BGC


def reorder(BGC, orf_ori, orf_pos, orf_domains):
    BGC = reverse_neg(BGC, orf_ori)

    C_starter_id = -1
    for i in range(len(BGC)):
        orf = BGC[i]
        if "C_Starter" in orf_domains[orf]:
            C_starter_id = i

    if C_starter_id != -1:
        BGC = [BGC[C_starter_id]] + BGC[0:C_starter_id] + BGC[C_starter_id + 1:]

    TE_TD_id = -1
    for i in range(len(BGC)):
        orf = BGC[i]
        if ("TE" in orf_domains[orf]) or ("TD" in orf_domains[orf]):
            TE_TD_id = i

    if TE_TD_id != -1:
        BGC = BGC[0:TE_TD_id] + BGC[TE_TD_id + 1:] + [BGC[TE_TD_id]]

    if not is_correct_NCterm(BGC, orf_ori, orf_domains):
        return generate_all_perm(BGC, orf_ori, orf_pos, orf_domains)

    return [BGC]


def split_and_reorder(BGCs, orf_ori, orf_pos, orf_domains):
    res_bgcs = []
    for BGC in BGCs:
        if is_correct(reverse_neg(BGC, orf_ori), orf_ori, orf_pos, orf_domains):
            res_bgcs.append(BGC)
        else:
            for lft in range(len(BGC)):
                for rgh in range(lft, len(BGC)):
                    if is_correct_orfs_subseq(lft, rgh, BGC, orf_ori, orf_pos, orf_domains):
                        res_bgcs.append(BGC[lft:rgh + 1])

    bgcs_reorder = []
    for BGC in res_bgcs:
        bgcs_reorder += reorder(BGC, orf_ori, orf_pos, orf_domains)

    res_bgcs = []
    for BGC in bgcs_reorder:
        if is_correct(BGC, orf_ori, orf_pos, orf_domains):
            res_bgcs.append(BGC)
    return res_bgcs