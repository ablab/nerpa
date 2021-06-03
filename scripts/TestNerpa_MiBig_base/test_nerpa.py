import os
import sys
import csv
import matplotlib

import garlic_res_utils
import showFDR

matplotlib.use('Agg')

AA_info={}
list_AA=[]
res_dir = sys.argv[1]
data_path = sys.argv[2]


#all_res[0] = genome id, all_res[1] = structure id, all_res[2] = score
def getRank(all_res, BGC, score, for_structure=False):
    pos = 0
    if for_structure:
        pos = 1

    rank = 0
    for res in all_res:
        if (BGC in res[pos]) and float(res[2]) > float(score) + 0.0001:
            rank += 1
    return rank

with open(os.path.join(data_path, 'base.csv'), 'r') as f:
    csv_reader = csv.reader(f, delimiter=',')
    for row in csv_reader:
        if (row[0] == 'Accession' or row[0] == "Id"):
            continue
    
        AA_info[row[0]] = [row[2], row[3], 1, 1]
        list_AA.append(row[0])
        
all_res = []
with open( res_dir + "/report.csv", 'r') as f:
    csv_reader = csv.reader(f, delimiter=',')
    for row in csv_reader:
        if row[0] != 'score':
            all_res.append([row[6], row[5], row[0]])


#result_match[mol_id] = [score, match_cnt]
result_match = {}
cnt_all = 0
cnt_found = 0
cnt_08per = 0
cnt_05per = 0
score_iswrong = [] 
with open( res_dir + "/mibig_cmp.csv", 'w' ) as fw:
    csv_writer = csv.writer(fw, delimiter=',', quotechar='"')
    csv_writer.writerow(['Accession', 'Main Product', 'Organism', 'Structure type', 'Molecule Length', 'Has Match', 'Score', 'Matched Length', 'Matched Percentage', 'Rank for Genome', 'Rank for Structure'])
    with open( res_dir + "/report.csv", 'r') as f:
        csv_reader = csv.reader(f, delimiter=',')
        for row in csv_reader:
            row[5] = row[5].split("_")[0]
            splBGC = "_".join(row[6].split("/")[-1].split('_')[:-4])
            if (row[5].split('/')[-1].split('.')[0] in row[6]):
                score_iswrong.append((float(row[0]), False, row[5].split('/')[-1], splBGC))
                result_match[row[5].split('/')[-1]] = [row[0], row[3]]
            elif row[0] != 'score':
                score_iswrong.append((float(row[0]), True, row[5].split('/')[-1], splBGC))

    for mol_name in list_AA:
        cnt_all += 1
        if mol_name in result_match:
            cnt_found += 1
            matched_percentage=round(float(result_match[mol_name][1])/float(AA_info[mol_name][3]), 2)
            if (matched_percentage > 0.79):
                cnt_08per += 1
            if (matched_percentage > 0.49):
                cnt_05per += 1

            csv_writer.writerow([mol_name, AA_info[mol_name][0], AA_info[mol_name][1], AA_info[mol_name][2], AA_info[mol_name][3], 'True', result_match[mol_name][0], result_match[mol_name][1], matched_percentage, getRank(all_res, mol_name.split('.')[0], result_match[mol_name][0]), getRank(all_res, mol_name.split('.')[0], result_match[mol_name][0], True)])
        else:
            csv_writer.writerow([mol_name] + AA_info[mol_name] + ['False', '0', '0', '0', '-1', '-1'])


showFDR.showAllFDR(data_path, res_dir, score_iswrong)

print("Find " + str(cnt_found) + " (out of " + str(cnt_all) + ")")
print("Find " + str(cnt_05per) + " (out of " + str(cnt_all) + ") half matched")
print("Find " + str(cnt_08per) + " (out of " + str(cnt_all) + ") almost completely matched")
print("The result can be found: " + res_dir)
