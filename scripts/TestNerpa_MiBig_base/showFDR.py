import os
import sys
import csv
import matplotlib

import garlic_res_utils

matplotlib.use('Agg')

MAX_ELEM=250

def FilterScore(score_iswrong, top_cnt):
    bstPred = {x[2]: [] for x in score_iswrong}
    for x in score_iswrong:
        bstPred[x[2]].append(x)

    res = []
    for bgc in bstPred.keys():
        rank = [0] * len(bstPred[bgc])
        cnt = 1
        for i in range(1, len(bstPred[bgc])):
            if bstPred[bgc][i - 1][0] != bstPred[bgc][i][0]:
                for j in range(i - cnt, i):
                    rank[j] = (2*i - cnt)/2
                cnt = 0
            cnt += 1
        
        for j in range(len(bstPred[bgc]) - cnt, len(bstPred[bgc])):
            rank[j] = (2*len(bstPred[bgc]) - cnt)/2

        rs = 0
        for i in range(0, len(bstPred[bgc])):
            if rank[i] < top_cnt:
                if bstPred[bgc][i][1] == False:
                    rs = i
        #print(rank[0])
        #print(bstPred[bgc])
        if rank[0] < top_cnt:
            res.append(bstPred[bgc][rs])
        
    res.sort(key=lambda x: (-x[0], x[1]))
    return res  
             

def calcFDR(score_iswrong, top_cnt=0, printF=False):
    if (top_cnt > 0):
        score_iswrong = FilterScore(score_iswrong, top_cnt)
    
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
        if x[1] == 1 and printF:
            print("FP: ", x)

    return cnt, scores, FDR

def showFDR(cnt, scores, FDR, res_dir, out_prefix="FDR_"):
    import matplotlib
    import matplotlib.pyplot as plt
    
    plt.clf()
    plt.plot(cnt[:MAX_ELEM], FDR[:MAX_ELEM])
    plt.xlabel('number of elements')
    plt.ylabel('FDR')
    plt.savefig(res_dir + "/" + out_prefix + "cnt.png") 
    plt.clf()

    with open(res_dir + "/" + out_prefix + "cnt.csv", "w") as fw:
        fw.write("CNT, FDR\n")
        for i in range(len(cnt[:MAX_ELEM])):
            fw.write(str(cnt[i]) + ", " + str(FDR[i]) + "\n")

    plt.gca().invert_xaxis()
    plt.plot(scores[:MAX_ELEM], FDR[:MAX_ELEM])
    plt.xlabel('scores')
    plt.ylabel('FDR')
    plt.savefig( res_dir + "/" + out_prefix + "scores.png")
    plt.clf()

    with open( res_dir + "/" + out_prefix + "scores.csv", "w") as fw:
        fw.write("SCORES, FDR\n")
        for i in range(len(scores[:MAX_ELEM])):
            fw.write(str(scores[i]) + ", " + str(FDR[i]) + "\n")

def showFDRwithGARLIC(cnt_nerpa, cnt_garlic, FDR_nerpa, FDR_garlic, res_dir, out_prefix="FDR_"):
    import matplotlib
    import matplotlib.pyplot as plt

    plt.plot(cnt_nerpa[:MAX_ELEM], FDR_nerpa[:MAX_ELEM], label="nerpa")
    plt.plot(cnt_garlic[:MAX_ELEM], FDR_garlic[:MAX_ELEM], label="garlic")
    plt.xlabel('number of elements')
    plt.ylabel('FDR')
    plt.legend()
    plt.savefig( res_dir + "/" + out_prefix + "with_garlic.png")
    plt.clf()


def showAllFDR(data_path, res_dir, score_iswrong, setAA):
    garlic_report_path = os.path.join(data_path, "base_line/GARLIC/report.csv")
    
    print(score_iswrong[:100])

    score_iswrong.sort(key=lambda x: (-x[0], x[1]))
    #cnt, scores, FDR = calcFDR(score_iswrong)
    #cntg, scoresg, FDRg = calcFDR(garlic_res_utils.get_score_iswrong(garlic_report_path))
    
    #showFDRwithGARLIC(cnt, cntg, FDR, FDRg, res_dir)
    
    #cnt, scores, FDR = calcFDR(score_iswrong, 1, True)
    #cntg, scoresg, FDRg = calcFDR(garlic_res_utils.get_score_iswrong(garlic_report_path), 1)
    
    #showFDRwithGARLIC(cnt, cntg, FDR, FDRg, res_dir, "FDR_top1_mol_")
    
    #cnt, scores, FDR = calcFDR(score_iswrong, 3)
    #cntg, scoresg, FDRg = calcFDR(garlic_res_utils.get_score_iswrong(garlic_report_path), 3)
    
    #showFDRwithGARLIC(cnt, cntg, FDR, FDRg, res_dir, "FDR_top3_mol_")
    
    for i in range(len(score_iswrong)):
        score_iswrong[i]  = (score_iswrong[i][0], score_iswrong[i][1], score_iswrong[i][3], score_iswrong[i][2]) 
        
    garlic_score = garlic_res_utils.get_score_iswrong(garlic_report_path, setAA)
    for i in range(len(garlic_score)):
        garlic_score[i] = (garlic_score[i][0], garlic_score[i][1], garlic_score[i][3], garlic_score[i][2])
        
    #cnt, scores, FDR = calcFDR(score_iswrong, 1)
    #cntg, scoresg, FDRg = calcFDR(garlic_score, 1)
    
    #showFDRwithGARLIC(cnt, cntg, FDR, FDRg, res_dir, "FDR_top1_genome_")
    
    cnt, scores, FDR = calcFDR(score_iswrong, 20, True)
    cntg, scoresg, FDRg = calcFDR(garlic_score, 20, True)
    
    showFDRwithGARLIC(cnt, cntg, FDR, FDRg, res_dir, "FDR_top20_genome_")
