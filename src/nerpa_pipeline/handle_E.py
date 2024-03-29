#!/usr/bin/env python

import handle_helper


def get_D_AA(dirname, orfs_order):
    D_AA_list = []
    # ctg_orf, ctg_orf_Aid, domain_type(AMP-binding, PCP, MT)
    domains=handle_helper.get_domains_list(dirname)
    orf_ori = handle_helper.get_orf_orientation(dirname)

    # reverse domains on - strand
    for i in range(len(domains)):
        if domains[i]:
            if orf_ori[domains[i][0][0]] == '-':
                domains[i].reverse()

    order_dom = []
    for corf in orfs_order:
        for j in range(len(domains)):
            if domains[j][0][0] == corf:
                order_dom.append(domains[j])
    domains = order_dom

    is_AA = lambda dlst : ("PKS" in dlst[2]) or ("AMP-binding" == dlst[2])

    has_d = False
    cur_di = 0
    cur_i = 0

    for did in range(len(domains)):
        orfds = domains[did]
        for i in range(len(orfds)):
            if is_AA(orfds[i]):
                if has_d:
                    D_AA_list.append(domains[cur_di][cur_i][0] + "_" + domains[cur_di][cur_i][1].split('_')[-1])
                cur_i = i
                cur_di = did
                has_d = False

            if orfds[i][-1] == "Condensation_Dual" or orfds[i][-1] == "Epimerization":
                has_d = True

    if has_d:
        D_AA_list.append(domains[cur_di][cur_i][0] + "_" + domains[cur_di][cur_i][1].split('_')[-1])

    return D_AA_list
