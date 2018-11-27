import os
from .vis_prediction import visualize_prediction

path = "/home/olga/CAB/NRP/data/serverRuns/"
genome_file = path + "genome.fna"
pathToAntismash = "/home/olga/CAB/NRP/soft/antismash-3.0.5.1_CUT/run_antismash.py"
antismashRes = path + "antismashRes/"
predictionPath = antismashRes + "nrpspks_predictions_txt/ctg1_nrpspredictor2_codes.txt"
predictionInfo = path + "predictions.info"

DBinfo = "/home/olga/CAB/NRP/data/DataBase/library.info.streptomedb"
NRPsMatcher = "/home/olga/tmp/NRP/bin/run_nrp_matcher.py"

def run_antismash():
    os.system("python2 " + pathToAntismash + " " + genome_file + " --outputfolder " + antismashRes)
    with open(predictionInfo, "w") as fw:
        fw.write(predictionPath)

def run_nrpsMatcher():
    os.system("python3 " + NRPsMatcher + " -p " + predictionInfo + " --lib_info " + DBinfo + " -o " + path)

def save_results(request_id):
    for filename in os.listdir(path + "/details_mols"):
        visualize_prediction(path + "/details_mols/" + filename, predictionPath, filename, "ctg1_nrpspredictor2_codes", request_id)

def handle_genome(f, request_id):
    with open(genome_file, "wb") as fw:
        for chunk in f.chunks():
            fw.write(chunk)

    run_antismash()
    run_nrpsMatcher()
    save_results(request_id)