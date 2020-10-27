import os
import csv

def get_domains_list(dirname):
    #ctg_orf, ctg_orf_Aid, domain_type(AMP-binding, PCP, MT)
    domains = []
    cur_orf = ""
    txt_folder = os.path.join(dirname, "txt")
    if  not os.path.isdir(txt_folder):
        return domains

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
    if  not os.path.isdir(txt_folder):
        return orf_ori

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


def get_orf_position(dirname):
    txt_folder = os.path.join(dirname, "txt")

    orf_pos = {}
    for filename in os.listdir(txt_folder):
        if filename.endswith("_gene.txt"):
            csv_file_with_orf = os.path.join(txt_folder, filename)
            with open(csv_file_with_orf, 'r') as rf:
                csv_reader = csv.reader(rf, delimiter='\t')
                for row in csv_reader:
                    if "gene" in row[0]:
                        continue

                    orf_pos[row[0]] = (int(row[1]), int(row[2]))
    return orf_pos

def get_orf_domain_list(dirname):
    txt_folder = os.path.join(dirname, "txt")
    orf_domain_list = {}

    for filename in os.listdir(txt_folder):
        if filename.endswith("_NRPS_PKS.txt"):
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

def debug_print_parts(filename, parts, orf_domain_list, orf_ori, orf_pos):
    print("Filename: " + filename)
    for i in range(len(parts)):
        print("Part# " + str(i))
        for j in range(len(parts[i])):
            cur_orf = parts[i][j][0]
            print(cur_orf + "[" + str(orf_pos[cur_orf][0]) + "-" +
                  str(orf_pos[cur_orf][1]) + "; " + orf_ori[cur_orf] + "]: " + str(orf_domain_list[cur_orf]))