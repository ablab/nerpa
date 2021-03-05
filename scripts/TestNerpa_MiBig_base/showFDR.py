import os
import sys
import csv
import matplotlib

import garlic_res_utils

matplotlib.use('Agg')

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

def showFDR(cnt, scores, FDR, res_dir, out_prefix="FDR_"):
    import matplotlib
    import matplotlib.pyplot as plt
    
    plt.clf()
    plt.plot(cnt[:160], FDR[:160])
    plt.xlabel('number of elements')
    plt.ylabel('FDR')
    plt.savefig(res_dir + "/" + out_prefix + "cnt.png") 
    plt.clf()

    with open(res_dir + "/" + out_prefix + "cnt.csv", "w") as fw:
        fw.write("CNT, FDR\n")
        for i in range(len(cnt[:160])):
            fw.write(str(cnt[i]) + ", " + str(FDR[i]) + "\n")

    plt.gca().invert_xaxis()
    plt.plot(scores[:160], FDR[:160])
    plt.xlabel('scores')
    plt.ylabel('FDR')
    plt.savefig( res_dir + "/" + out_prefix + "scores.png")
    plt.clf()

    with open( res_dir + "/" + out_prefix + "scores.csv", "w") as fw:
        fw.write("SCORES, FDR\n")
        for i in range(len(scores[:160])):
            fw.write(str(scores[i]) + ", " + str(FDR[i]) + "\n")

def showFDRwithGARLIC(cnt_nerpa, cnt_garlic, FDR_nerpa, FDR_garlic, res_dir, out_prefix="FDR_"):
    import matplotlib
    import matplotlib.pyplot as plt

    plt.plot(cnt_nerpa[:160], FDR_nerpa[:160], label="nerpa")
    plt.plot(cnt_garlic[:160], FDR_garlic[:160], label="garlic")
    plt.xlabel('number of elements')
    plt.ylabel('FDR')
    plt.legend()
    plt.savefig( res_dir + "/" + out_prefix + "with_garlic.png")
    plt.clf()


def showAllFDR(data_path, res_dir, score_iswrong):
    garlic_report_path = os.path.join(data_path, "base_line/GARLIC/report.csv")
    
    score_iswrong.sort(key=lambda x: (-x[0], x[1]))
    cnt, scores, FDR = calcFDR(score_iswrong)
    cntg, scoresg, FDRg = calcFDR(garlic_res_utils.get_score_iswrong(garlic_report_path))
    
    showFDR(cnt, scores, FDR, res_dir)
    showFDR(cntg, scoresg, FDR, res_dir, "FDR_garlic_")
    showFDRwithGARLIC(cnt, cntg, FDR, FDRg, res_dir)
    
    cnt, scores, FDR = calcFDR(score_iswrong, 1)
    cntg, scoresg, FDRg = calcFDR(garlic_res_utils.get_score_iswrong(garlic_report_path), 1)
    
    showFDR(cnt, scores, FDR, res_dir, "FDR_top1_mol_")
    showFDR(cntg, scoresg, FDRg, res_dir, "FDR_top1_mol_garlic_")
    showFDRwithGARLIC(cnt, cntg, FDR, FDRg, res_dir, "FDR_top1_mol_")
    
    cnt, scores, FDR = calcFDR(score_iswrong, 3)
    cntg, scoresg, FDRg = calcFDR(garlic_res_utils.get_score_iswrong(garlic_report_path), 3)
    
    showFDR(cnt, scores, FDR, res_dir, "FDR_top3_mol_")
    showFDR(cntg, scoresg, FDRg, res_dir, "FDR_top3_mol_garlic_")
    showFDRwithGARLIC(cnt, cntg, FDR, FDRg, res_dir, "FDR_top3_mol_")
    
    for i in range(len(score_iswrong)):
        score_iswrong[i]  = (score_iswrong[i][0], score_iswrong[i][1], score_iswrong[i][3], score_iswrong[i][2]) 
        
    garlic_score = garlic_res_utils.get_score_iswrong(garlic_report_path)
    for i in range(len(garlic_score)):
        garlic_score[i] = (garlic_score[i][0], garlic_score[i][1], garlic_score[i][3], garlic_score[i][2])
        
    cnt, scores, FDR = calcFDR(score_iswrong, 1)
    cntg, scoresg, FDRg = calcFDR(garlic_score, 1)
    
    showFDR(cnt, scores, FDR, res_dir, "FDR_top1_genome_")
    showFDR(cntg, scoresg, FDRg, res_dir, "FDR_top1_genome_garlic_")
    showFDRwithGARLIC(cnt, cntg, FDR, FDRg, res_dir, "FDR_top1_genome_")
    
    cnt, scores, FDR = calcFDR(score_iswrong, 3)
    cntg, scoresg, FDRg = calcFDR(garlic_score, 3)
    
    showFDR(cnt, scores, FDR, res_dir, "FDR_top3_genome_")
    showFDR(cntg, scoresg, FDRg, res_dir, "FDR_top3_genome_garlic_")
    showFDRwithGARLIC(cnt, cntg, FDR, FDRg, res_dir, "FDR_top3_genome_")
