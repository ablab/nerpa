GCFid = []
with open("PNP_report.csv", "r") as f:
    for line in f:
        if ('GCF' in line):
            GCFid.append((line.split('GCF_')[1]).split('.')[0])

GCFid = list(set(GCFid))

with open("bacteria_all.csv", "r") as f:
    with open("ball.csv", "w") as fw:
        for line in f:
            wasinline = False
            for gcfid in GCFid:
                if gcfid in line:
                    wasinline = True
            if ('TRUE' not in line) and ('FALSE' not in line):
                if (wasinline):
                    nl = line[:-1] + 'TRUE'
                else:
                    nl = line[:-1] + 'FALSE'
            else:
                nl = line[:-1]
            fw.write(nl + "\n")
            
        
