import os

finfo = open('library.info')
fsmiles = open('library.smiles')

infos = finfo.readlines()
smiles = fsmiles.readlines()

for i in range(len(infos)):
    info = ((infos[i].split(' ')[0]).split('/')[-1]).split('.')[0]
    smile = smiles[i]
    os.system("java -jar grape-release.jar -s \"" + smile + "\" -img -json -txt -o json/" +  info)    

finfo.close()
fsmiles.close()
