import sys
import csv

def get_score_iswrong(path_to_report, setAA):
    score_iswrong=[]
    with open(path_to_report, 'r') as f:
        csv_reader = csv.reader(f, delimiter='\t')
        for row in csv_reader:
            if row[0] not in setAA:
                continue
            if (row[0] in row[1]):
                score_iswrong.append((float(row[2]), False, row[1], row[0]))
            elif row[0] != 'genome':
                score_iswrong.append((float(row[2]), True, row[1], row[0]))
        score_iswrong.sort(key=lambda x: (-x[0], x[1]))
    return score_iswrong
