import sys

file_name = sys.argv[1]
aset = set()

with open(file_name) as f:
    print("ok")
    for line in f:
        ln = line.split()[2]
        aminoasds = ln.split(";")

        for am in aminoasds:
            s = am[:am.index('(')]
            aset.add(s)


for s in aset:
    print(s)



