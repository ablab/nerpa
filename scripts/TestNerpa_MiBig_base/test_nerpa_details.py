import os
import sys
import csv
import nrp_structure_graph
from nrp_structure_graph import Graph

dirname = sys.argv[1]
details_dir = "result/" + dirname + "/details_mols" 

graphs = []

for filename in os.listdir(details_dir):
    bgcname = filename.split('.')[0]
    with open(details_dir + "/" + filename, 'r') as f:
        lines = f.readlines()
        line_id = 0
        while (line_id < len(lines)):
            if ".gr" in lines[line_id]:
                line_id += 1
                if bgcname in lines[line_id]:
                    ctg = lines[line_id].split('/')[-1].split('_')[1]
                    line_id += 2
                    cnt_v = lines[line_id].split(' ')[-1]
                    cur_gr = Graph(cnt_v, filename, ctg)
                    line_id += 1
                    for i in range(int(cnt_v)):
                        cur_gr.add_vert_info(lines[line_id])
                        line_id += 1

                    cnt_e = int(lines[line_id].split(' ')[-1])
                    line_id += 1
                    for i in range(cnt_e):
                        cur_gr.add_e(int(lines[line_id].split()[0]), int(lines[line_id].split()[-1]))
                        line_id += 1
                    graphs.append(cur_gr)
                    break
            else:
                line_id += 1


with open("result/" + dirname + "/details.csv", "w") as fw:
    csv_writer = csv.writer(fw, delimiter=',', quotechar='"')
    csv_writer.writerow(['BGC', 'CONTIG', 'ORF', 'A-ID', 'STRUCTURE', 'VERTEX', 'AA-ID', 'rBan AA'])
    for graph in graphs:
        graph.write_to_csv(csv_writer)
