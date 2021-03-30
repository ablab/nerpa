#!/usr/bin/env python

import handle_helper


def get_D_AA(dirname):
    D_AA_list = []
    # ctg_orf, ctg_orf_Aid, domain_type(AMP-binding, PCP, MT)
    domains=handle_helper.get_domains_list(dirname)
    orf_ori = handle_helper.get_orf_orientation(dirname)

    # reverse domains on - strand
    for i in range(len(domains)):
        if domains[i]:
            if orf_ori[domains[i][0][0]] == '-':
                domains[i].reverse()

    is_AA = lambda dlst : ("PKS" in dlst[2]) or ("AMP-binding" == dlst[2])

    for orfds in domains:
        cur_i = 0
        has_d = False
        for i in range(len(orfds)):
            if is_AA(orfds[i]):
                if has_d:
                    D_AA_list.append(orfds[cur_i][0] + "_" + orfds[cur_i][1].split('_')[-1])

                cur_i = i
                has_d = False

            if orfds[i][-1] == "Condensation_Dual" or orfds[i][-1] == "Epimerization":
                has_d = True

        if has_d:
            D_AA_list.append(orfds[cur_i][0] + "_" + orfds[cur_i][1].split('_')[-1])

    return D_AA_list
