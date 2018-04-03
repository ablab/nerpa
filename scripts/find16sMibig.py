import sys
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

id_to_org = {}
molid_to_id = {}
org_to_16s = {}
molidlist = []
orgs = set()

with open("mibig.html") as file:
    lines = file.readlines()
    i = 0
    while (i < len(lines)):
        line = lines[i]
        if '<tr>' in line:
            i += 1
            curid = lines[i][-21:-11]
            i += 4
            org = lines[i][:-7]
            org = org[org.find('>') + 1:]
            id_to_org[curid] = org
            orgs.add(org)
        i += 1

with open("lib.info.mibig") as file:
    for line in file:
        pos = line.find('BGC')
        curid = line[pos:pos + 10]
        molid = line.split(' ')[0].split('/')[-1].split('.')[0]
        molid_to_id[molid] = curid
        molidlist.append(molid)

db = "SILVA_123_SSURef_Nr99_tax_silva.fasta"
for record in SeqIO.parse(db, "fasta"):
    print(record.id)
    print(record.description)
    cur_seq = str(record.seq)
    cur_seq = cur_seq.replace('U', 'T')
    cur_ids = record.description.split(';')
    cur_ids[0] = ' '.join(cur_ids[0].split(' ')[1:])
    for cur_id in cur_ids:
        cid = cur_id
        if (cid in orgs):
            if (cid not in org_to_16s) or (len(org_to_16s[cid]) < len(cur_seq)):
                org_to_16s[cid] = cur_seq

records = []
for molid in molidlist:
    org = "**"
    if molid_to_id[molid] in id_to_org:
        org = id_to_org[molid_to_id[molid]]
    seq = ""
    if (org in org_to_16s):
        seq = org_to_16s[org]
    records.append(SeqRecord(Seq(seq, generic_dna), id=molid, description=org))

SeqIO.write(records, "path16s.fasta", "fasta")
