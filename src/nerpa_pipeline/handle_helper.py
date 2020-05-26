import os
import csv

def get_domains_list(dirname):
    #ctg_orf, ctg_orf_Aid, domain_type(AMP-binding, PCP, MT)
    domains = []
    cur_orf = ""
    txt_folder = os.path.join(dirname, "txt")
    for filename in os.listdir(txt_folder):
        if filename.endswith("_NRPS_PKS.txt"):
            csv_file_with_orf = os.path.join(txt_folder, filename)
            with open(csv_file_with_orf, 'r') as rf:
                csv_reader = csv.reader(rf, delimiter='\t')
                for row in csv_reader:
                    if row[1] == "NRPSPKS_ID":
                        continue
                    if row[1] != cur_orf:
                        cur_orf = row[1]
                        domains.append([])
                    if row[6] == "Condensation":
                        row[6] = row[7]

                    domains[-1].append([row[1], row[3], row[6]])
    return domains

def get_orf_orientation(dirname):
    txt_folder = os.path.join(dirname, "txt")

    orf_ori = {}
    #get orientation of orfs
    for filename in os.listdir(txt_folder):
        if filename.endswith("_gene.txt"):
            csv_file_with_orf = os.path.join(txt_folder, filename)
            with open(csv_file_with_orf, 'r') as rf:
                csv_reader = csv.reader(rf, delimiter='\t')
                for row in csv_reader:
                    if "gene" in row[0]:
                        continue

                    orf_ori[row[0]] = row[3]
    return orf_ori