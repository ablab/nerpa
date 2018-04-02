import os

fnaPATH = "/media/hosein/My Passport/hosein/Desktop/project/sequence_data/bacteria_complete/fna/"
prismPATH = "/home/dereplicator/kolga/soft/prism.jar"
webContentPATH = "/home/dereplicator/kolga/soft/WebContent"
outPATH = "/home/dereplicator/kolga/data/bacteria_complete/json/"

for root, dirs, files in os.walk(fnaPATH):
    print("ok")
    for file in files:
        print(file)
        os.system("java -jar " + prismPATH + " -a -p -f \"" + fnaPATH + file + "\" -tt -o "
                  + outPATH + " -w 10000 -r " + webContentPATH)
