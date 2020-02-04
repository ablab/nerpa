import sys
import csv
import matplotlib
matplotlib.use('Agg')

AA_info={}
list_AA=[]
dirname = sys.argv[1]

def calcFDR(score_iswrong):
    FDR=[] # (sum wrong)/(sum all)
    scores=[]
    cnt=[]

    cnt_all = 0
    cnt_wrong = 0
    for x in score_iswrong:
        cnt_all += 1
        cnt_wrong += x[1]
        cnt.append(cnt_all)
        scores.append(x[0])
        FDR.append(cnt_wrong/cnt_all)

    import matplotlib
    import matplotlib.pyplot as plt
    plt.plot(cnt[:200], FDR[:200])
    plt.xlabel('number of elements')
    plt.ylabel('FDR')
    plt.savefig('result/' + dirname + "/FDR_cnt.png")
  
    plt.clf()
    plt.gca().invert_xaxis()
    plt.plot(scores[:200], FDR[:200])
    plt.xlabel('scores')
    plt.ylabel('FDR')
    plt.savefig('result/' + dirname + "/FDR_scores.png")


with open('mibig_base.csv', 'r') as f:
    csv_reader = csv.reader(f, delimiter=',')
    for row in csv_reader:
        if (row[0] == 'Accession'):
            continue
        AA_info[row[0]] = [row[2], row[3], row[4], row[5]]
        list_AA.append(row[0])
        
result_match = {}
cnt_all = 0
cnt_found = 0
score_iswrong = [] 
with open("result/" + dirname + ".csv", 'w' ) as fw:
    csv_writer = csv.writer(fw, delimiter=',', quotechar='"')
    csv_writer.writerow(['Accession', 'Main Product', 'Organism', 'Structure type', 'Molecule Length', 'Has Match', 'Score', 'Matched Length'])
    with open("result/" + dirname + "/report.csv", 'r') as f:
        csv_reader = csv.reader(f, delimiter=',')
        for row in csv_reader:
            if (row[5].split('/')[-1].split('.')[0] in row[6]):
                score_iswrong.append((float(row[0]), False))
                result_match[row[5].split('/')[-1][:-3]] = [row[0], row[3]]
            elif row[0] != 'score':
                score_iswrong.append((float(row[0]), True))

    for mol_name in list_AA:
        cnt_all += 1
        if mol_name in result_match:
            cnt_found += 1
            csv_writer.writerow([mol_name, AA_info[mol_name][0], AA_info[mol_name][1], AA_info[mol_name][2], AA_info[mol_name][3], 'True', result_match[mol_name][0], result_match[mol_name][1]])
        else:
            csv_writer.writerow([mol_name] + AA_info[mol_name] + ['False', '0', '0'])

score_iswrong.sort(key=lambda x: (-x[0], x[1]))
calcFDR(score_iswrong)


print("Find " + str(cnt_found) + " out of " + str(cnt_all))
print("The result can be found: " + "result/" + dirname)
