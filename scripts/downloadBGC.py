from Bio import Entrez

Entrez.email = "OChe24912@yandex.ru"

handle = Entrez.esearch(db="Nucleotide", retmax=50000, term="\"biosynthetic gene cluster\"[All Fields]")
record = Entrez.read(handle)
handle.close()

ids = record['IdList']

def getFileName(s):
    return "/home/dereplicator/kolga/data/BGC/fna/" + (s.split(' ')[0])[1:] + ".fna"
    
for i in range(4521, len(ids)):
    idd = ids[i]
    handle = Entrez.efetch(db="nucleotide", id=idd, rettype="fasta", retmode="text")
    lines = handle.readlines() 
    fileName = getFileName(lines[0])
    with open(fileName, 'w') as fw:
        for i in range(len(lines)):
            fw.write(lines[i])
    
    handle.close()
