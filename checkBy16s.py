import sys
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def find_mol(cid):
    if (cid.find('/') == -1):
        return cid
    return (cid.split('/')[-1])[:-4]

mols16s_file = sys.argv[1]
report_file = sys.argv[2]
mol16s = {}
orgname = {}

fnaprefix = "/media/hosein/My Passport/hosein/Desktop/project/sequence_data/bacteria_complete/"

def run_blastn(molid, genid):
    global mol16s
    SeqIO.write(SeqRecord(Seq(mol16s[molid], generic_dna), id=molid), "tmp.fasta", "fasta")
    filebac = "\"" + fnaprefix + "fna/" + genid + ".fna\""

    os.system("blastn -query tmp.fasta -subject " + filebac + " > 16sscore")
    

def hasORG(genid, orgname):
    orgname = orgname.split(' ')[1]
    filename = fnaprefix + "fna/" + genid + ".fna"
    with open(filename) as g:
        curline = g.readline()
        print(orgname.lower())
        print(curline.lower())
        if (orgname.lower() in curline.lower()):
            return True
    return False

def getScoreId(s):
    if "inf" in s:
        return -1

    pos = s.find(".")
    print(s[:pos])
    val = int(s[:pos])
    if (val//5 >= 149):
        return -1
    return val//5

def getScor(lg):
    if "%" in lg:
        pos = lg.find("%")
        stpos = pos - 1
        while lg[stpos].isdigit():
            stpos -= 1
        val = int(lg[stpos + 1:pos])
        if (val < 85):
            return 0
        elif (val < 95):
            return 1
        else:
            return 2
        
    return -1


cntScoreValue = [[0]*5 for i in range(150)]

for record in SeqIO.parse(mols16s_file, "fasta"):
    mol16s[find_mol(record.id)] = str(record.seq)
    orgname[find_mol(record.id)] = str(record.description)

fw = open("report_16score", "w")
hf = open("Has_strep", "w")
with open(report_file) as f:
    for line in f:
        parts = line.split(' ')
        molid = find_mol(parts[0])
        parts = parts[:-1]
        print(parts)
        for j in range(2, len(parts), 2):
            print(parts[j])
            scoreId = getScoreId(parts[j + 1])
            genid = parts[j].split('/')[0] # or -1 for streptomyse
            run_blastn(molid, genid)
            fw.write(molid + " " + genid + " " + parts[j + 1] + "\n")
            bestScore = -1
            with open("16sscore") as g:
                for lg in g:
                    if ("Score =" in lg) or ("Identities" in lg):
                        fw.write(lg)
                        bestScore = max(getScor(lg), bestScore)
            if (bestScore != -1):
                cntScoreValue[scoreId][bestScore] += 1
            elif (hasORG(genid, orgname[molid])):
                hf.write("HAS ORG " + str(scoreId) + " " + str(molid) + " " + str(genid) + "\n")
                cntScoreValue[scoreId][-1] += 1
            else:
                cntScoreValue[scoreId][-2] += 1
                        
            fw.write("\n")
fw.close()
hf.close()

fw = open("report_cnt", "w")
for i in range(150):
    fw.write(str(i * 5) + "-" + str((i  + 1) * 5 - 1)  + ": " + "(0-84%) " + str(cntScoreValue[i][0]) + "; (85-94%) " +
             str(cntScoreValue[i][1]) + "; (95-100%) " + str(cntScoreValue[i][2]) + "; no data no org " +
             str(cntScoreValue[i][3]) + "; no data org " + str(cntScoreValue[i][4]) + "\n")

fw.write("inf: " + "(0-84%) " + str(cntScoreValue[-1][0]) + "; (85-94%) " +
        str(cntScoreValue[-1][1]) + "; (95-100%) " + str(cntScoreValue[-1][2]) +
        "; no data no org " + str(cntScoreValue[-1][3]) +
        "; no data org " + str(cntScoreValue[-1][4]) + "\n")


fw.write("--------------\n")
for i in range(147, -1, -1):
    for j in range(5):
        cntScoreValue[i][j] += cntScoreValue[i + 1][j]

for i in range(150):
    fw.write(">=" + str(i * 5)  + ": " + "(0-84%) " + str(cntScoreValue[i][0]) + "; (85-94%) " +
             str(cntScoreValue[i][1]) + "; (95-100%) " + str(cntScoreValue[i][2]) + "; no data no org " +
             str(cntScoreValue[i][3]) + "; no data org " + str(cntScoreValue[i][4]) + "\n")

fw.close()
