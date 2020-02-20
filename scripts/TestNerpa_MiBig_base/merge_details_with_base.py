from operator import itemgetter
import os
import sys
import csv

dirname = sys.argv[1]
details_base_csv = "base_line/details_base.csv"


#(BGC, CONTIG, ORG, A-ID) -> (STRUCTURE, VERTEX, AA-ID)
mstru = dict()
BGClist = dict()


with open("result/" + dirname + "/details.csv") as f:
    csv_reader = csv.reader(f, delimiter=',')
    for row in csv_reader:
        if row[0] != "BGC":
            bgcsh = row[0].split('.')[0]
            if bgcsh not in BGClist:
                BGClist[bgcsh] = []
                
            if row[0] not in BGClist[bgcsh]:
                BGClist[bgcsh].append(row[0])

            mstru[(row[0], row[1], row[2], row[3])] = (row[4], row[5], row[6])


with open(details_base_csv) as f:
    csv_reader = csv.reader(f, delimiter=',')
    oth = next(csv_reader)[11:]

    output = []
    for row in csv_reader:
        if (row[0] not in BGClist):
            BGClist[row[0]] = [row[0]]

        for curbgc in BGClist[row[0]]:
            struct_nerpa = ['-', '-', '-']
            if (curbgc, row[1], row[2], row[3]) in mstru:
                struct_nerpa = list(mstru[(curbgc, row[1], row[2], row[3])])
            
            output.append([curbgc] + row[1:8] + struct_nerpa + row[8:])

    output.sort(key=itemgetter(0))

    with open("result/" + dirname + "/details_with_base.csv", "w") as fw:
        csv_writer = csv.writer(fw, delimiter=',')
        csv_writer.writerow(['BGC', 'CONTIG', 'ORF', 'A-ID', "PRED_TOP5", "START_POS", "END_POS", "STRAND", "STRUCTURE NERPA", "VERTEX NERPA", "AA-ID NERPA", "STRUCTURE BASE", "VERTEX BASE", "AA-ID BASE"] + oth)
        for curout in output:
            csv_writer.writerow(curout)
