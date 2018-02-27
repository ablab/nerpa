import sys
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def find_mol(cid):
    return (cid.split('/')[-1])[:-4]

mols16s_file = sys.argv[1]
report_file = sys.argv[2]
mol16s = {}

def run_blastn(molid, genid):
    global mol16s
    SeqIO.write(SeqRecord(Seq(mol16s[molid], generic_dna), id=molid), "tmp.fasta", "fasta")
    filebac = "fna/" + genid + ".fna"

    os.system("blastn -query tmp.fasta -subject " + filebac + " > 16sscore")
    

for record in SeqIO.parse(mols16s_file, "fasta"):
    mol16s[find_mol(record.id)] = str(record.seq)


fw = open("report_16score", "w")
with open(report_file) as f:
    for line in f:
        parts = line.split(' ')
        molid = find_mol(parts[0])
        parts = parts[:-1]
        print(parts)
        for j in range(2, len(parts), 2):
            print(parts[j])
            genid = parts[j].split('/')[-1]
            run_blastn(molid, genid)
            fw.write(molid + " " + genid + " " + parts[j + 1] + "\n")
            with open("16sscore") as g:
                for lg in g:
                    if ("Score =" in lg) or ("Identities" in lg):
                        fw.write(lg)
            fw.write("\n")

fw.close()
