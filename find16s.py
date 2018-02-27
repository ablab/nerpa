import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

infofile = sys.argv[1]
db = sys.argv[2]

ids = set()
paths_list = []
path_to_ids = {}
id_to_16s = {}
id_to_paths = {}

def find_best_seq(path):
    global path_to_ids
    global id_to_16s
    curids = path_to_ids[path]
    res = ""
    for cid in curids:
        if (cid in id_to_16s):
            if (len(id_to_16s[cid]) > len(res)):
                res = id_to_16s[cid]
    return res

with open(infofile) as file:
    for line in file:
        parts = line.split('\t')
        cur_ids = parts[4].split(';')
        cur_path = parts[5]
        paths_list.append(cur_path)
        path_to_ids[cur_path] = []
        for cur_id in cur_ids:
            ids.add(cur_id)
            path_to_ids[cur_path].append(cur_id)
            if (cur_id not in id_to_paths):
                id_to_paths[cur_id] = []
            id_to_paths[cur_id].append(cur_path)


for record in SeqIO.parse(db, "fasta"):
    print(record.id)
    print(record.description)
    cur_seq = str(record.seq)
    cur_seq = cur_seq.replace('U', 'T')
    cur_ids = record.description.split(';')
    cur_ids[0] = ' '.join(cur_ids[0].split(' ')[1:])
    for cur_id in cur_ids:
        cid = cur_id.replace(' ', '_')
        if (cid in ids):
            if (cid not in id_to_16s) or (len(id_to_16s[cid]) < len(cur_seq)):
                id_to_16s[cid] = cur_seq



records = []
for path in paths_list:
    seq = find_best_seq(path)
    records.append(SeqRecord(Seq(seq, generic_dna), id=path))

SeqIO.write(records, "path16s.fasta", "fasta")
