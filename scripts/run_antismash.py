import os

fna_folder = "/home/olga/bio/NRP/data/gene_cluster/fna/"
output_folder = "/home/olga/bio/NRP/data/gene_cluster/"
prediction_folder = "/home/olga/bio/NRP/data/prediction/"

cnt = 0
with open(output_folder + "predictions_files", "w") as fw:
    for f in os.listdir(fna_folder):
        name = f[:-4]
        fw.write("prediction/" + name + "\n")
        print("python2 /home/olga/bio/NRP/soft/antismash-3.0.5.1_CUT/run_antismash.py --outputfolder " + output_folder + " " + fna_folder + f)
        os.system("python2 /home/olga/bio/NRP/soft/antismash-3.0.5.1_CUT/run_antismash.py --outputfolder " + output_folder + " " + fna_folder + f)
        os.system("mv nrpspks_predictions_txt/ctg1_nrpspredictor2_codes.txt prediction/" + name)
        os.system("rm -r nrpspks_predictions_txt/")
        print(cnt)
        cnt += 1
    
