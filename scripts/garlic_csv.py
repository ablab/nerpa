import os
import json

PATH = "/home/dereplicator/kolga/data/bacteria_complete/PNP/garlic_decoy/"
#PATH = "/home/olga/" 
dirs = os.listdir(PATH)

csv = open("report.csv", "w")
csv.write("genome\tmols\ts1\ts2\ts3\ts4\ts5\ts6\ts7\ts8\ts9\ts10\n")

def isGenomeDir(cdir):
    return len(cdir) > 3 and cdir[:3] == "GCF"

def parseName(cdir):
    if ("anti" not in cdir) and ("mibig" not in cdir):
        lst = cdir.split('_')
        return '_'.join(lst[:-1]), lst[-1]
    else:
        lst = cdir.split('_')
        return '_'.join(lst[:-2]), '_'.join(lst[-2:])
    
        

def processFile(fl):
    with open(fl) as f:
        prj = json.load(f)
        if (prj["garlic_results"][0]["scores"]["self_score"] < 0):
            return -1000
        return prj["garlic_results"][0]["scores"]["relative_score"]
        
    
def processDir(cdir):
    genome, mol = parseName(cdir)
    fls = os.listdir(PATH + cdir)
    scrs = []
    for fl in fls:
        scrs.append(processFile(PATH + cdir + "/" + fl))
    scrs.sort(reverse=True)
    while len(scrs) < 10:
        scrs.append(-1000)

    if (scrs[0] > 0.3):
        csv.write(genome + "\t" + mol)
        for i in range(0, 10):
            csv.write("\t" + str(scrs[i]))
        csv.write('\n')

for cdir in dirs:
    if isGenomeDir(cdir):
        processDir(cdir)

csv.close()
        
