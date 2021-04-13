#!/usr/bin/env python

import handle_helper


def get_double_orfs_and_AA(dirname):
    double_orf_list = []
    double_AA_list = []

    # ctg_orf, ctg_orf_Aid, domain_type(AMP-binding, PCP, MT)
    domains=handle_helper.get_domains_list(dirname)
    orf_ori = handle_helper.get_orf_orientation(dirname)

    # reverse domains on - strand
    for i in range(len(domains)):
        if domains[i]:
            if orf_ori[domains[i][0][0]] == '-':
                domains[i].reverse()

    is_AA = lambda dlst : ("PKS" in dlst[2]) or ("AMP-binding" == dlst[2])
    # check dpuble PCP in the end for double orf list
    for orfds in domains:
        i = len(orfds) - 1
        while (i > 0) and (not is_AA(orfds[i])):
            if orfds[i][2] == "PCP" and orfds[i - 1][2] == "PCP":
                double_orf_list.append(orfds[i][0])
                break
            i -= 1

    # calculate PCP Condensation PCP between AMP-binding/PKS
    for orfds in domains:
        cur_i = 0
        cur_pcp_level = 0
        for i in range(len(orfds)):
            if is_AA(orfds[i]):
                if cur_pcp_level == 3:
                    double_AA_list.append(orfds[cur_i][0] + "_" + orfds[cur_i][1].split('_')[-1])
                cur_i = i
                cur_pcp_level = 0

            if ("PCP" == orfds[i][2]) and (cur_pcp_level == 0 or cur_pcp_level == 2):
                cur_pcp_level += 1

            if ("Condensation" == orfds[i][2]) and (cur_pcp_level == 1):
                cur_pcp_level += 1

        if cur_pcp_level == 3:
            double_AA_list.append(orfds[cur_i][0] + "_" + orfds[cur_i][1].split('_')[-1])

    return double_orf_list, double_AA_list
