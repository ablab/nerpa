import os

fna_folder = "/media/hosein/My Passport/Alexey_Olga/BGC/fna/"
output_folder = "/media/hosein/My Passport/Alexey_Olga/BGC/"
prediction_folder = "/media/hosein/My Passport/Alexey_Olga/BGC/prediction/"

cnt = 0
with open(output_folder + "predictions_files", "w") as fw:
    for f in os.listdir(fna_folder):
        name = f[:-4]
        fw.write("prediction/" + name + "\n")
        print("python2 /home/dereplicator/kolga/soft/antismash-3.0.5.1_CUT/run_antismash.py --outputfolder '" + output_folder + "' '" + fna_folder + f + "'")
        os.system("python2 /home/dereplicator/kolga/soft/antismash-3.0.5.1_CUT/run_antismash.py --outputfolder '" + output_folder + "' '" + fna_folder + f + "'")
        os.system("mv nrpspks_predictions_txt/ctg1_nrpspredictor2_codes.txt prediction/" + name)
        os.system("rm -r nrpspks_predictions_txt/")
        print(cnt)
        cnt += 1
    
