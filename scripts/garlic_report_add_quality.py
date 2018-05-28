import csv

qul = dict()
with open("PNP_report.csv") as f:
    pnp_report = csv.DictReader(f)
    for row in pnp_report:
        if (row["prediction id"] not in qul):
            qul[row["prediction id"]] = dict()
        qul[row["prediction id"]][row["mol id"]] = row["quality"]

with open ("qreport.csv", "w") as fw:
    fw.write("quality\tgenome\tmols\ts1\ts2\ts3\ts4\ts5\ts6\ts7\ts8\ts9\ts10\n")
    
    with open("garlic_report.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            val = -1
            if (row["genome"] in qul and row["mols"] in qul[row["genome"]]):
                val = qul[row["genome"]][row["mols"]]
            fw.write(str(val) + "\t" + row["genome"] + "\t" + row["mols"]);
            for i in range(1, 11):
                fw.write("\t"+ str(row["s" + str(i)]))
            fw.write("\n")
            
        
    
    
    
