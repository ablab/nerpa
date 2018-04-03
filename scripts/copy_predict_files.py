#!/usr/bin/python3
import os


files_to_copy = "full_list_prediction_files.txt"
write_file = "predictions_files"

fw = open(write_file, 'w')
with open(files_to_copy) as f:
    for line in f:
        if (line[-1] == '\n'):
            line=line[:-1]
        parts = line.split('/')
        shortname = ""
        for part in parts:
            if (part[:3] == "GCF"):
                shortname = part

        print("scp dereplicator@hosein.andrew.cmu.edu:\""+line + "\" predictions/" + shortname)
        os.system("scp dereplicator@hosein.andrew.cmu.edu:\""+line + "\" predictions/" + shortname)
        fw.write("predictions/" + shortname + "\n")

fw.close()