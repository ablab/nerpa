#!/usr/bin/env python3

import sys
import os
import pandas as pd

import nrp_structure_graph as nsg 

def handle_file(fname, align):
    with open(fname) as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if line.startswith('ALIGN'):
            str_id = lines[i - 3].split('_')[0].split('.')[0]
            pr_id = lines[i - 2].split('/')[-1].split('_')[0]
            if str_id == pr_id:
                j = i + 2
                
                str_id = lines[i  - 3].split('_')[0]
                stj = len(align)
                while not lines[j].startswith("GRAPH"):
                    curl = lines[j].split(' ')
                    print(curl)
                    if curl[1] == "-":
                        align.append([pr_id, "-", "-", "-", "-", str_id, "STRUCTURE", curl[-4], curl[-3], curl[-2], curl[-1].strip()])
                    else:
                        align.append([pr_id, curl[0].split("_")[0], curl[0].split("_")[1], "A" + str(int(curl[1]) + 1), curl[2] , str_id, "STRUCTURE", curl[-4], curl[-3], curl[-2], curl[-1].strip()])
                    j += 1
                
                while not lines[j].startswith("number of com"):
                    j += 1
                vcnt = int(lines[j].split(' ')[-1])
                gph = nsg.Graph(vcnt, str_id, "ctg1")

                while not lines[j].startswith("number of b"):
                    j += 1
                ecnt = int(lines[j].split(' ')[-1])
                for g in range(ecnt):
                    j += 1
                    v = int(lines[j].split(' ')[0])
                    u = int(lines[j].split(' ')[-1])
                    gph.add_e(v, u)

                strst = gph.print_structure()
                for g in range(stj, len(align)):
                    align[g][6] = strst


def main():
    res_folder = sys.argv[1]
    data_folder = sys.argv[2]
    align = []
    dpath = os.path.join(res_folder, "details_mols")
    for fname in os.listdir(dpath):
        if os.path.isfile(os.path.join(dpath, fname)):
            handle_file(os.path.join(dpath, fname), align)

    df = pd.DataFrame(align, columns=["BGC", "CONTIG", "ORF", "A-ID", "L-/D-", "STRUCTURE ID", "rBan STRUCTURE", "rBan VERTEX", "rBan AA-ID", "rBan STRUCT_CONFIGURATION", "rBan AA"])
    print(df.head())
    df.to_csv(os.path.join(res_folder, "correct_align.csv"), index=False)

    df_data = pd.read_csv(os.path.join(data_folder, "details.csv"))
    df_res = pd.merge(df_data, df,  how="outer", on=["BGC", "CONTIG", "ORF", "A-ID"])
    df_res.drop_duplicates(subset=None, keep='first', inplace=True)
    df_res.to_csv(os.path.join(res_folder, "details.csv"), index=False)


if __name__ == "__main__":
    main()
