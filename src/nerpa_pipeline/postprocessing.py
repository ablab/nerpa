#!/usr/bin/env python3
import csv
import os

def originPathInReport(output_dir, origin_file):
    lines = []
    with open(os.path.join(output_dir, "report.csv"), "r") as f:
        csv_reader = csv.reader(f)
        for row in csv_reader:
            lines.append(row)

    pos_id = lines[0].index("PredictionFileName")
    lines[0] = lines[0][:pos_id + 1] + ["ContigId", "PartId"] + lines[0][pos_id + 1:]
    for i in range(1, len(lines)):
        curPath = lines[i][pos_id]
        orgPath = origin_file[curPath]
        ctgid = curPath.split('/')[-1].split('_')[-4]
        partid = curPath.split('/')[-1].split('_')[-1]

        lines[i] = lines[i][:pos_id] + [orgPath, ctgid, partid] + lines[i][pos_id + 1:]

    with open(os.path.join(output_dir, "report.csv"), "w") as fw:
        writer = csv.writer(fw, delimiter=',')
        for line in lines:
            writer.writerow(line)


def originPathInDetailReport(output_dir, origin_file):
    for filename in os.listdir(os.path.join(output_dir, "details")):
        if filename.endswith('.match'):
            lines = []
            with open(os.path.join(output_dir, "details", filename)) as f:
                for line in f:
                    if line[:-1] in origin_file:
                        orgPath = origin_file[line[:-1]]
                        ctgid = line.split('/')[-1].split('_')[-4]
                        partid = line.split('/')[-1].split('_')[-1]
                        lines.append(orgPath + "\t" + ctgid + "\t" + partid)
                    else:
                        lines.append(line)

            with open(os.path.join(output_dir, "details", filename), "w") as fw:
                for line in lines:
                    fw.write(line)


def postproccesing(output_dir, origin_file):
    originPathInReport(output_dir, origin_file)
    originPathInDetailReport(output_dir, origin_file)