#!/usr/bin/env python3

import argparse
import os
import math
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--details", help="path to matching csv")
    parser.add_argument("--mon", help="monomers list")
    parser.add_argument("-o", dest="output", help="path to output directory")

    return parser.parse_args()


def logPFromCnt(res):
    PT_NT, PT_NF, PF_NT, PF_NF = 0, 1, 2, 3
    sumNT = (res[PT_NT] + res[PF_NT])/sum(res)
    sumNF = (res[PT_NF] + res[PF_NF])/sum(res)

    sumPT = res[PT_NT] + res[PT_NF]
    sumPF = res[PF_NT] + res[PF_NF]

    res[PT_NT] /= (sumPT*sumNT)
    res[PF_NT] /= (sumPF*sumNT)
    res[PT_NF] /= (sumPT*sumNF)
    res[PF_NF] /= (sumPF*sumNF)
    for i in range(4):
        res[i] = math.log(res[i])


def get_MT_stats(args):
    PT_NT, PT_NF, PF_NT, PF_NF = 0, 1, 2, 3
    res = [0, 0, 0, 0]
    df = pd.read_csv(args.details, sep=",")
    for i in range(len(df)):
        mt_pred = df["MT"].iloc[i]
        mt_struct = "+MT" in str(df["rBan AA"].iloc[i])
        if mt_pred == True and mt_struct == True:
            res[PT_NT] += 1
        if mt_pred == False and mt_struct == True:
            res[PF_NT] += 1
        if mt_pred == True and mt_struct == False:
            res[PT_NF] += 1
        if mt_pred == False and mt_struct == False:
            res[PF_NF] += 1
    logPFromCnt(res)
    return res


def get_DL_stata(args):
    PT_NT, PT_NF, PF_NT, PF_NF = 0, 1, 2, 3
    res = [0, 0, 0, 0]
    df = pd.read_csv(args.details, sep=",")
    for i in range(len(df)):
        d_pred = df["L-/D-"].iloc[i]
        d_struct = str(df["rBan AA-ID"].iloc[i]).split("-")[0]
        if d_pred == "D" and d_struct == "@D":
            res[PT_NT] += 1
        if d_pred == "L" and d_struct == "@D":
            res[PF_NT] += 1
        if d_pred == "D" and d_struct == "@L":
            res[PT_NF] += 1
        if d_pred == "L" and d_struct == "@L":
            res[PF_NF] += 1
    logPFromCnt(res)
    return res


def get_modifications_csv(args):
    res = [[], []]
    res[0] = ["MT"] + get_MT_stats(args)
    res[1] = ["@D"] + get_DL_stata(args)
    df = pd.DataFrame(res, columns=["Name", "predT-nrpT", "predT-nrpF", "predF-nrpT", "predF-nrpF"])
    opath = os.path.join(args.output, "modifications.tsv")
    df.to_csv(opath, sep="\t", index=False)


def getMnIds(args):
    df = pd.read_csv(args.mon, sep="\t")
    mapMnName = {}
    for i in range(len(df)):
        mapMnName[df["Code"].iloc[i]] = df["NameID"].iloc[i]

    mnList = list({x for k, x in mapMnName.items()})
    return mapMnName, mnList


def getMonomersLogP(args, mapMnName, mnList):
    df = pd.read_csv(args.details, sep=",")
    cntMn = {mn: 0 for mn in mnList}
    for i in range(len(df)):
        aaId = df["rBan AA-ID"].iloc[i]
        if aaId != aaId:
            continue
        if aaId[0] == "@":
            aaId = aaId[3:]
        if aaId[0] == "*":
            aaId = aaId[1:]
        if aaId == "-":
            continue
        if aaId not in mapMnName:
            print(aaId, "NOT found")
            continue
        cntMn[mapMnName[aaId]] += 1

    sumCnt = sum([x for x in cntMn.values()])
    logP = [[mn, math.log(max(cntMn[mn]/sumCnt, 1/sumCnt))] for mn in mnList]
    df1 = pd.DataFrame(logP, columns=["NameID", "LogP"])
    opath = os.path.join(args.output, "monomersLogP.tsv")
    df1.to_csv(opath, sep="\t", index=False)


def getInsertion(args):
    df = pd.read_csv(args.details, sep=",")
    cntIns = 0
    cntNotIns = 0
    for i in range(len(df)):
        preds = df["PRED_TOP5"].iloc[i]
        if len(str(preds)) > 3:
            aar = str(df["rBan AA"].iloc[i]).split("+")[0]
            if aar == "-":
                cntIns += 1
                continue

            if len(aar) > 1:
                cntNotIns += 1

    print(cntIns, cntNotIns)
    logP = math.log(cntIns/(cntIns + cntNotIns))

    with open(os.path.join(args.output, "prob_gen.cfg"), "w") as fw:
        fw.write("insertion " + str(logP) + "\n")


def getDeletion(args):
    df = pd.read_csv(args.details, sep=",")
    cntDel = 0
    cntNotDel = 0
    for i in range(len(df)):
        preds = df["CONTIG"].iloc[i]
        if len(str(preds)) > 1:
            aar = str(df["rBan AA"].iloc[i]).split("+")[0]
            if len(aar) > 1:
                cntNotDel += 1
        if preds == "-":
            aar = str(df["rBan AA"].iloc[i]).split("+")[0]
            if len(aar) > 1:
                cntDel += 1


    print(cntDel, cntNotDel)
    logP = math.log(cntDel/(cntDel + cntNotDel))

    with open(os.path.join(args.output, "prob_gen.cfg"), "a") as fw:
        fw.write("deletion " + str(logP) + "\n")

def getProbGen_In_correct(args):
    df = pd.read_csv(args.details, sep=",")
    cntScore = {100: 0, 90: 0, 80: 0, 70: 0, 60: 0}
    sumScore = {100: 0, 90: 0, 80: 0, 70: 0, 60: 0}
    for i in range(len(df)):
        preds = df["PRED_TOP5"].iloc[i]
        if len(str(preds)) > 3:
            preds = preds.split(";")[0]
            aap = preds.split("(")[0]
            preds = max(60, int(float(preds.split("(")[-1].split(")")[0])))

            aar = str(df["rBan AA"].iloc[i]).split("+")[0]
            if len(aar) < 2:
                continue

            sumScore[preds] += 1
            if aap != aar:
                cntScore[preds] += 1

    keys = [100, 90, 80, 70, 60]
    logP = [str(math.log(cntScore[k]/sumScore[k])) for k in keys]
    logPCor = [str(math.log(1 - cntScore[k]/sumScore[k])) for k in keys]
    with open(os.path.join(args.output, "prob_gen.cfg"), "a") as fw:
        fw.write("ProbGenCorrect " + " ".join(logPCor) + "\n")
        fw.write("ProbGenIncorrect " + " ".join(logP) + "\n")


def getProbGetScore(args):
    df = pd.read_csv(args.details, sep=",")
    cntScore = {100: 0, 90: 0, 80: 0, 70: 0, 60: 0}
    for i in range(len(df)):
        preds = df["PRED_TOP5"].iloc[i]
        if len(str(preds)) > 3:
            preds = preds.split(";")[0]
            preds = max(60, int(float(preds.split("(")[-1].split(")")[0])))
            cntScore[preds] += 1

    sumScore = sum([x for x in cntScore.values()])
    keys = [100, 90, 80, 70, 60]
    logP = [str(math.log(cntScore[k]/sumScore)) for k in keys]
    with open(os.path.join(args.output, "prob_gen.cfg"), "a") as fw:
        fw.write("ProbGetScore " + " ".join(logP))


def main():
    args = parse_args()
    get_modifications_csv(args)
    mapMnName, mnList = getMnIds(args)
    getMonomersLogP(args, mapMnName, mnList)
    getInsertion(args)
    getDeletion(args)
    getProbGen_In_correct(args)
    getProbGetScore(args)


if __name__ == "__main__":
    main()