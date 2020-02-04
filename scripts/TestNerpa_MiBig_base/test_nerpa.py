import sys
import csv

AA_info={}
list_AA=[]
dirname = sys.argv[1]
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
with open("result/" + dirname + ".csv", 'w' ) as fw:
    csv_writer = csv.writer(fw, delimiter=',', quotechar='"')
    csv_writer.writerow(['Accession', 'Main Product', 'Organism', 'Structure type', 'Molecule Length', 'Has Match', 'Score', 'Matched Length'])
    with open("result/" + dirname + "/report.csv", 'r') as f:
        csv_reader = csv.reader(f, delimiter=',')
        for row in csv_reader:
            if (row[5].split('/')[-1].split('.')[0] in row[6]):
                result_match[row[5].split('/')[-1][:-3]] = [row[0], row[3]]

    for mol_name in list_AA:
        cnt_all += 1
        if mol_name in result_match:
            cnt_found += 1
            csv_writer.writerow([mol_name, AA_info[mol_name][0], AA_info[mol_name][1], AA_info[mol_name][2], AA_info[mol_name][3], 'True', result_match[mol_name][0], result_match[mol_name][1]])
        else:
            csv_writer.writerow([mol_name] + AA_info[mol_name] + ['False', '0', '0'])


print("Find " + str(cnt_found) + " out of " + str(cnt_all))
