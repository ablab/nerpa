import sys
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def find_mol(cid):
    if (cid.find('/') == -1):
        return cid
    if (cid[-1] == 'l'):
        return (cid.split('/')[-1])[:-4]
    else:
        return (cid.split('/')[-1])[:-3]

def findegenid(s):
    prts = s.split('/')
    for i in range(len(prts)):
        if ("GCF" in prts[i]):
            return prts[i]
    return 

mols16s_file = sys.argv[1]
report_file = sys.argv[2]
mol16s = {}
orgname = {}

fnaprefix = "/media/hosein/My Passport/hosein/Desktop/project/sequence_data/bacteria_complete/"

def run_blastn(molid, genid):
    global mol16s
    if (molid not in mol16s):
        with open("16sscore", "w") as ffw:
            ffw.write("\n")
        return
        
    SeqIO.write(SeqRecord(Seq(mol16s[molid], generic_dna), id=molid), "tmp.fasta", "fasta")
    filebac = "\"" + fnaprefix + "fna/" + genid + ".fna\""

    print(filebac)
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

def getOrgName(genid):
    filename = fnaprefix + "fna/" + genid + ".fna"
    with open(filename) as g:
        curline = g.readline()
        curline=curline[curline.find(' ')+1:curline.find(',')]
        return curline
        

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

def getFullScore(lg):
    if "%" in lg:
        pos = lg.find("%")
        stpos = pos - 1
        while lg[stpos].isdigit():
            stpos -= 1
        val = int(lg[stpos + 1:pos])
        return val
        
    return -1



cntScoreValue = [[0]*5 for i in range(150)]

for record in SeqIO.parse(mols16s_file, "fasta"):
    mol16s[find_mol(record.id)] = str(record.seq)
    orgname[find_mol(record.id)] = str(record.description)


rw = open("newreport.csv", "w")
fw = open("report_16score", "w")
with open(report_file) as f:
    lines = f.readlines()
    header = lines[0][:-1] + ",origin organism,prediction organism,16s score"
    rw.write(header + "\n")
    for line in lines[1:]:
        parts = line.split(',')
        molid = find_mol(parts[6])
        print(molid)
        scoreId = getScoreId(parts[1])
        genid = findegenid(parts[7])
        run_blastn(molid, genid)
        fw.write(molid + " " + genid + " " + parts[1] + "\n")
        bestScore = -1
        rbstscor = -1
        with open("16sscore") as g:
            for lg in g:
                if ("Score =" in lg) or ("Identities" in lg):
                    fw.write(lg)
                    bestScore = max(getScor(lg), bestScore)
                    rbstscor = max(getFullScore(lg), rbstscor)
        if (bestScore != -1):
            cntScoreValue[scoreId][bestScore] += 1
        elif ((molid in orgname) and hasORG(genid, orgname[molid])):
            cntScoreValue[scoreId][-1] += 1
        else:
            cntScoreValue[scoreId][-2] += 1                
        fw.write("\n")
        if (molid in orgname):
            rw.write(line[:-1] + "," + orgname[molid] + "," + getOrgName(genid) + "," + str(rbstscor) + "\n")
        else:
            rw.write(line[:-1] + ",," + getOrgName(genid) + "," + str(rbstscor) + "\n")
        
        
fw.close()
rw.close()



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
