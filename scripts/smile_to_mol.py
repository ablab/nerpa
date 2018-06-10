import csv
import os
smile = []
oid = []

with open('emibig.csv') as csvfile:
    csvreader = csv.reader(csvfile)
    for row in csvreader:
        if (row[1] != '-' and row[1] != 'smile'):
            smile.append(row[1])
            oid.append(row[0])

with open('emibig.info', 'w') as fw:
    for i in range(len(oid)):
        fw.write("/home/olga/bio/NRP/data/mibig2016/mol_dir/emibig_" + str(i) + ".mol " + oid[i] + "\n")

os.mkdir(os.path.dirname('./mol_dir/'))

for i in range(len(oid)):
    print("molconvert mol:V3+H --smiles \"" + smile[i] + "\" -o mol_dir/emibig_" + str(i) + ".mol")
    os.system("molconvert mol:V3+H --smiles \"" + smile[i] + "\" -o mol_dir/emibig_" + str(i) + ".mol")
