import garlic_res_utils

import sys
import csv
import matplotlib
matplotlib.use('Agg')

AA_info={}
list_AA=[]
dirname = sys.argv[1]

def calcFDR(score_iswrong, top_cnt=0):
    FDR=[] # (sum wrong)/(sum all)
    scores=[]
    cnt=[]

    cnt_val = {}

    cnt_all = 0
    cnt_wrong = 0
    for x in score_iswrong:
        if x[2] not in cnt_val:
            cnt_val[x[2]] = 0

        if (top_cnt > 0 and cnt_val[x[2]] == top_cnt):
            continue

        cnt_val[x[2]] += 1
        cnt_all += 1
        cnt_wrong += x[1]
        cnt.append(cnt_all)
        scores.append(x[0])
        FDR.append(cnt_wrong/cnt_all)

    return cnt, scores, FDR

def showFDR(cnt, scores, FDR, out_prefix="FDR_"):
    import matplotlib
    import matplotlib.pyplot as plt
    
    plt.clf()
    plt.plot(cnt[:160], FDR[:160])
    plt.xlabel('number of elements')
    plt.ylabel('FDR')
    plt.savefig('result/' + dirname + "/" + out_prefix + "cnt.png") 
    plt.clf()

    with open("result/" + dirname + "/" + out_prefix + "cnt.csv", "w") as fw:
        fw.write("CNT, FDR\n")
        for i in range(len(cnt[:160])):
            fw.write(str(cnt[i]) + ", " + str(FDR[i]) + "\n")

    plt.gca().invert_xaxis()
    plt.plot(scores[:160], FDR[:160])
    plt.xlabel('scores')
    plt.ylabel('FDR')
    plt.savefig('result/' + dirname + "/" + out_prefix + "scores.png")
    plt.clf()

    with open("result/" + dirname + "/" + out_prefix + "scores.csv", "w") as fw:
        fw.write("SCORES, FDR\n")
        for i in range(len(scores[:160])):
            fw.write(str(scores[i]) + ", " + str(FDR[i]) + "\n")

def showFDRwithGARLIC(cnt_nerpa, cnt_garlic, FDR_nerpa, FDR_garlic, out_prefix="FDR_"):
    import matplotlib
    import matplotlib.pyplot as plt

    plt.plot(cnt_nerpa[:160], FDR_nerpa[:160], label="nerpa")
    plt.plot(cnt_garlic[:160], FDR_garlic[:160], label="garlic")
    plt.xlabel('number of elements')
    plt.ylabel('FDR')
    plt.legend()
    plt.savefig('result/' + dirname + "/" + out_prefix + "with_garlic.png")
    plt.clf()

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

with open('mibig_base.csv', 'r') as f:
    csv_reader = csv.reader(f, delimiter=',')
    for row in csv_reader:
        if (row[0] == 'Accession'):
            continue
        AA_info[row[0]] = [row[2], row[3], row[4], row[5]]
        list_AA.append(row[0])
        
all_res = []
with open("result/" + dirname + "/report.csv", 'r') as f:
    csv_reader = csv.reader(f, delimiter=',')
    for row in csv_reader:
        if row[0] != 'score':
            all_res.append([row[6], row[5], row[0]])


result_match = {}
cnt_all = 0
cnt_found = 0
cnt_08per = 0
cnt_05per = 0
score_iswrong = [] 
with open("result/" + dirname + "/mibig_cmp.csv", 'w' ) as fw:
    csv_writer = csv.writer(fw, delimiter=',', quotechar='"')
    csv_writer.writerow(['Accession', 'Main Product', 'Organism', 'Structure type', 'Molecule Length', 'Has Match', 'Score', 'Matched Length', 'Matched Percentage', 'Rank for Genome', 'Rank for Structure'])
    with open("result/" + dirname + "/report.csv", 'r') as f:
        csv_reader = csv.reader(f, delimiter=',')
        for row in csv_reader:
            if (row[5].split('/')[-1].split('.')[0] in row[6]):
                score_iswrong.append((float(row[0]), False, row[5].split('/')[-1][:-3], row[6]))
                result_match[row[5].split('/')[-1][:-3]] = [row[0], row[3]]
            elif row[0] != 'score':
                score_iswrong.append((float(row[0]), True, row[5].split('/')[-1][:-3], row[6]))

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

score_iswrong.sort(key=lambda x: (-x[0], x[1]))
cnt, scores, FDR = calcFDR(score_iswrong)
cntg, scoresg, FDRg = calcFDR(garlic_res_utils.get_score_iswrong("./base_line/GARLIC/report.csv"))

showFDR(cnt, scores, FDR)
showFDR(cntg, scoresg, FDR, "FDR_garlic_")
showFDRwithGARLIC(cnt, cntg, FDR, FDRg)

cnt, scores, FDR = calcFDR(score_iswrong, 1)
cntg, scoresg, FDRg = calcFDR(garlic_res_utils.get_score_iswrong("./base_line/GARLIC/report.csv"), 1)

showFDR(cnt, scores, FDR, "FDR_top1_mol_")
showFDR(cntg, scoresg, FDRg, "FDR_top1_mol_garlic_")
showFDRwithGARLIC(cnt, cntg, FDR, FDRg, "FDR_top1_mol_")

cnt, scores, FDR = calcFDR(score_iswrong, 3)
cntg, scoresg, FDRg = calcFDR(garlic_res_utils.get_score_iswrong("./base_line/GARLIC/report.csv"), 3)

showFDR(cnt, scores, FDR, "FDR_top3_mol_")
showFDR(cntg, scoresg, FDRg, "FDR_top3_mol_garlic_")
showFDRwithGARLIC(cnt, cntg, FDR, FDRg, "FDR_top3_mol_")

for i in range(len(score_iswrong)):
    score_iswrong[i]  = (score_iswrong[i][0], score_iswrong[i][1], score_iswrong[i][3], score_iswrong[i][2]) 

garlic_score = garlic_res_utils.get_score_iswrong("./base_line/GARLIC/report.csv")
for i in range(len(garlic_score)):
    garlic_score[i] = (garlic_score[i][0], garlic_score[i][1], garlic_score[i][3], garlic_score[i][2])

cnt, scores, FDR = calcFDR(score_iswrong, 1)
cntg, scoresg, FDRg = calcFDR(garlic_score, 1)

showFDR(cnt, scores, FDR, "FDR_top1_genome_")
showFDR(cntg, scoresg, FDRg, "FDR_top1_genome_garlic_")
showFDRwithGARLIC(cnt, cntg, FDR, FDRg, "FDR_top1_genome_")


cnt, scores, FDR = calcFDR(score_iswrong, 3)
cntg, scoresg, FDRg = calcFDR(garlic_score, 3)

showFDR(cnt, scores, FDR, "FDR_top3_genome_")
showFDR(cntg, scoresg, FDRg, "FDR_top3_genome_garlic_")
showFDRwithGARLIC(cnt, cntg, FDR, FDRg, "FDR_top3_genome_")


print("Find " + str(cnt_found) + " (out of " + str(cnt_all) + ")")
print("Find " + str(cnt_05per) + " (out of " + str(cnt_all) + ") half matched")
print("Find " + str(cnt_08per) + " (out of " + str(cnt_all) + ") almost completely matched")
print("The result can be found: " + "result/" + dirname)
