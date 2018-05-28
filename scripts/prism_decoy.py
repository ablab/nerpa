import os
import json
from random import *

inPATH = "/home/dereplicator/kolga/data/bacteria_complete/json/"
outPATH = "/home/dereplicator/kolga/data/bacteria_complete/decoy/"

amnNames = []
cntAmn = {}
    
def getStatisticOneFile(file):
    print(file)
    with open (inPATH + file) as f:
        djs = json.load(f)
        for si in djs["prism_results"]["clusters"]:
            for sj in si["orfs"]:
                for sg in sj["domains"]:
                    if sg["substrates"] is None:
                        continue
                    for sh in sg["substrates"]:
                        if (sh["name"] not in cntAmn):
                            cntAmn[sh["name"]] = [0]*100
                            amnNames.append(sh["name"])
                        cntAmn[sh["name"]][int(sh["score"]/100)] += 1
    return

def getStatistic():
    files = os.listdir(inPATH)
    for file in files:
        if (file[-4:] == 'json'):
            getStatisticOneFile(file)
    
    return

def genAA(scr):
    sum = 0
    for nm in amnNames:
        sum += cntAmn[nm][scr]

    val = randint(1, sum)
    cur = 0
    for nm in amnNames:
        cur += cntAmn[nm][scr]
        if (val > cur - cntAmn[nm][scr] and val <= cur):
            return nm
        
    

def rewriteOneJSON(file):
    print(file)
    with open (inPATH + file) as f:
        djs = json.load(f)
        for si in djs["prism_results"]["clusters"]:
            for sj in si["orfs"]:
                for sg in sj["domains"]:
                    if sg["substrates"] is None:
                        continue
                    for sh in sg["substrates"]:
                        sh["name"] = genAA(int(sh["score"]/100))
        with open(outPATH + file, "w") as outfile:
            json.dump(djs, outfile)

def rewriteJSON():
    files = os.listdir(inPATH)
    for file in files:
        if (file[-4:] == 'json'):
            rewriteOneJSON(file)
    

getStatistic()
rewriteJSON()

