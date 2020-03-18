import os
import csv
import json

PATH = "/Bmo/kolga/runs/Nerpa/Mibig2019/GARLIC/" 
dirs = os.listdir(PATH)

def get_cnt_match(prj):
    cnt_matched = 0
    if (len(prj["garlic_results"][0]["alignment"]["subject_aliged"]) == 0):
        return 0
    lst_match = prj["garlic_results"][0]["alignment"]["subject_aliged"][0]
    for mtch in lst_match:
        if mtch is not None:
            cnt_matched += 1
    return cnt_matched

def write_to_csv(BGC, fl, csv_writer, csv_writer_per):
    with open(fl) as f:
        prj = json.load(f)

        query_aligned = prj["garlic_results"][0]["alignment"]["query_aligned"]
        subject_aligned = prj["garlic_results"][0]["alignment"]["subject_aliged"]

        #print(json.dumps(prj, indent=4))
        cluster_start = prj["garlic_results"][0]["query"]["cluster_start"]
        cluster_end = prj["garlic_results"][0]["query"]["cluster_end"]
        quer_list = []
        for cur_orf in query_aligned:
            for elem in cur_orf:
                if elem is None:
                    quer_list.append(("-", "-", "-"))
                else:
                    quer_list.append((elem["orf"], elem["module"], elem["substrate"]))
        
        subject_list = []
        for cur_sub in subject_aligned:
            for elem in cur_sub:
                if elem is None:
                    subject_list.append("-")
                else:
                    subject_list.append(elem["substrate"])

        Grape_len = 0
        CntMatch = 0
        ExactMatch = 0
        for i in range(len(subject_list)):
            csv_writer.writerow([BGC, quer_list[i][0], quer_list[i][1], quer_list[i][2], cluster_start, cluster_end, subject_list[i]])
            if (subject_list[i] not in ["", "Mal", "MeM", None]):
                Grape_len += 1
                if subject_list[i] != "-" and (quer_list[i][2] is not None):
                    CntMatch += 1
                    if (quer_list[i][2] in subject_list[i]) or (subject_list[i] in quer_list[i][2]):
                        ExactMatch += 1
        csv_writer_per.writerow([BGC, Grape_len, CntMatch, ExactMatch, CntMatch/Grape_len, ExactMatch/Grape_len])
        



def processFile(fl):
    with open(fl) as f:
        prj = json.load(f)
        cnt_matched = get_cnt_match(prj)

        if (prj["garlic_results"][0]["scores"]["self_score"] < 0):
            return (-1000, cnt_matched)
        return (prj["garlic_results"][0]["scores"]["relative_score"], cnt_matched)


def processDir(cdir, csv_writer, csv_writer_per):
    if "log" in cdir or "report.csv" in cdir:
        return
   
    genome, mol = cdir.split('_')[0], cdir.split('_')[1]
    if genome not in mol:
        return

    print(genome, mol)
    
    fls = os.listdir(PATH + cdir)
    scrs = []
    for fl in fls:
        scrs.append((processFile(PATH + cdir + "/" + fl), PATH + cdir + "/" + fl))
    scrs.sort(reverse=True)

    if (len(scrs) == 0):
        return

    print(scrs[0])
    if (scrs[0][0][0] > 0):
        write_to_csv(mol, scrs[0][1], csv_writer, csv_writer_per)
    else:
        csv_writer_per.writerow([mol, 0, 0, 0, 0, 0])
        

with open("./base_line/GARLIC/matched_percentage.csv", "w") as fwper:
    csv_writer_per = csv.writer(fwper, delimiter=',', quotechar='"')
    csv_writer_per.writerow(["Accession", "Length by GRAPE", "#MATCH", "#EXACT MATCH", "Match Percentage", "Exact Match Percentage"])
    with open("./base_line/GARLIC/details.csv", "w") as fw:
        csv_writer = csv.writer(fw, delimiter=',', quotechar='"')
        csv_writer.writerow(["BGC", "ORF", "MODULE", "PRISM-AA", "CLUSTER_START", "CLUSTER_END", "GRAPE-AA"])
        for cdir in dirs:
            processDir(cdir, csv_writer, csv_writer_per)


fw.close()
        
