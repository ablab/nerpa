import os
import sys
import csv

class Graph:
    def __init__(self, cnt_v, bgcname, contig):
        self.gr = [[], []]
        self.str_type = "NA"
        self.first_vert = 0
        self.orientation = 0
        self.bgc = bgcname
        self.ctg = contig
        self.vert_info = []
        self.new_vert = {}

        cnt_v = int(cnt_v)
        for i in range(cnt_v):
            self.gr[0].append([])
            self.gr[1].append([])
            self.new_vert[str(i)] = str(i)

    def add_e(self, u, v):
        self.gr[0][u].append(v)
        self.gr[1][v].append(u)

    def add_vert_info(self, line):
        parts = line.split()
        if len(parts) < 7:
            return

        vert_id = parts[0]
        AA = parts[4].split('(')[0]
        Aid = 'A' + str(int(parts[-1]) + 1)
        orf = parts[-2].split('_')[-1]
        self.vert_info.append([orf, Aid, vert_id, AA])

    def is_matched(self, vid):
        for vinfo in self.vert_info:
            if str(vinfo[2]) == str(vid):
                return True
        return False

    def init_str_type(self):
        self.str_type = "cycle"
        for i in range(len(self.gr[0])):
            if (len(self.gr[0][i]) == 0) and (len(self.gr[1][i]) == 0):
                self.gr[0] = self.gr[0][:i] + self.gr[0][i + 1:]
                self.gr[1] = self.gr[1][:i] + self.gr[1][i + 1:]
                continue

            if (len(self.gr[0][i]) == 0) or (len(self.gr[1][i]) == 0):
                if (self.str_type == "cycle"):
                    self.str_type = "line"
                self.first_vert = i
                if (len(self.gr[1][i]) == 1):
                    self.orientation = 1
                if (len(self.gr[0][i]) == 1):
                    self.orientation = 0
            if (len(self.gr[0][i]) == 2) or (len(self.gr[1][i]) == 2):
                self.str_type = "branch"
        
    def get_vid_str(self, vid):
        if vid >= 0:
            return str(vid)
        return "#"

    def print_structure(self):
        if (self.str_type == "NA"):
            self.init_str_type()
        cur_v = self.first_vert
        
        id_v = 0
        if (self.str_type == "line") or (self.str_type == "cycle"):
            id_v = len(self.gr[0]) - 1

        if (not self.is_matched(self.first_vert)) and self.str_type == "branch":
            id_v -= 1

        if self.str_type == "line":
            for i in range(len(self.gr[0])):
                if ((len(self.gr[0][i]) == 1) and (len(self.gr[1][i]) == 0)) or ((len(self.gr[0][i]) == 0) and (len(self.gr[0][i]) == 1)):
                    if (i != self.first_vert):
                        if (not self.is_matched(i)):
                            id_v -= 1
                   

        structure_str = ""
       
        for i in range(len(self.gr[0])):
            if (i == 0):
                structure_str = self.get_vid_str(id_v)
            else:
                if len(self.gr[0][cur_v]) > 1 or len(self.gr[1][cur_v]) > 1:
                    structure_str = self.get_vid_str(id_v) + "*" + structure_str
                else:
                    structure_str = self.get_vid_str(id_v) + "," + structure_str
            
            self.new_vert[str(cur_v)] = self.get_vid_str(id_v)
            if (self.str_type == "branch"):
                id_v += 1
            else:
                id_v -= 1

            if len(self.gr[self.orientation][cur_v]) > 0:
                cur_v = self.gr[self.orientation][cur_v][0]

        if (self.str_type == "cycle"):
            structure_str += "*"
        return structure_str

    def write_to_csv(self, csv_writer):
        for vinfo in self.vert_info:
            csv_writer.writerow([self.bgc, self.ctg, vinfo[0], vinfo[1], self.print_structure(), self.new_vert[str(vinfo[2])], vinfo[3]])



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
                    ctg = lines[line_id].split('/')[-1].split('_')[0]
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
            else:
                line_id += 1


with open("result/" + dirname + "/details.csv", "w") as fw:
    csv_writer = csv.writer(fw, delimiter=',', quotechar='"')
    csv_writer.writerow(['BGC', 'CONTIG', 'ORF', 'A-ID', 'STRUCTURE', 'VERTEX', 'AA-ID'])
    for graph in graphs:
        graph.write_to_csv(csv_writer)
