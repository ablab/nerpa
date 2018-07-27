import os
from .vis_prediction import visualize_prediction
from celery import shared_task

path = "/home/olga/bio/NRP/data/serverRuns/"
genome_file = path + "genome.fna"
nrp_file = path + "nrp.mol"
smile_file = path + "nrp.sml"
predictionInfo = path + "predictions.info"
molInfo = path + "mol.info"

pathToAntismash = "/home/olga/bio/NRP/soft/antismash-3.0.5.1_CUT/run_antismash.py"
antismashRes = path + "antismashRes/"
predictionPath = antismashRes + "nrpspks_predictions_txt/ctg1_nrpspredictor2_codes.txt"

NRPsMatcher = "/home/olga/tmp/NRP/bin/run_nrp_matcher.py"

dbNRPinfo = {'streptome': "/home/olga/bio/NRP/data/DataBase/library.info.streptomedb"}
dbPredictionInfo = {'bc': "/home/olga/bio/NRP/data/DataBase/mibigNF.info"}


def run_antismash():
    os.system("python2 " + pathToAntismash + " " + genome_file + " --outputfolder " + antismashRes)
    with open(predictionInfo, "w") as fw:
        fw.write(predictionPath)


def run_nrpsMatcher(prinfo, molinfo):
    print("python3 " + NRPsMatcher + " -p " + prinfo + " --lib_info " + molinfo + " -o " + path)
    os.system("python3 " + NRPsMatcher + " -p " + prinfo + " --lib_info " + molinfo + " -o " + path)


def save_results(request_id):
    for filename in os.listdir(path + "/details_mols"):
        visualize_prediction(path + "/details_mols/" + filename, predictionPath, filename, "ctg1_nrpspredictor2_codes", request_id)


def getAllPath(filename):
    paths = []
    with open(path + "/details_mols/" + filename) as f:
        lines = f.readlines()
        cur = 0
        while (cur < len(lines)):
            paths.append(lines[cur + 1])
            while (cur < len(lines) and lines[cur] != "\n"):
                cur += 1
            while (cur < len(lines) and lines[cur] == "\n"):
                cur += 1
    return paths


def save_results_prediction(request_id):
    for filename in os.listdir(path + "/details_mols"):
        predpaths = getAllPath(filename)
        for predpath in predpaths:
            if (predpath[-1] == '\n'):
                predpath = predpath[:-1]
            visualize_prediction(path + "/details_mols/" + filename, predpath, "nrp", predpath, request_id)


def SMILE_to_MOL():
    print("molconvert mol:V3+H --smiles \"" + smile_file + "\" -o " + nrp_file)
    os.system("molconvert mol:V3+H --smiles \"" + smile_file + "\" -o " + nrp_file)


@shared_task
def handle_genome(request_id, nrpDB):
    print("START HANDLE GENOME")
    run_antismash()
    run_nrpsMatcher(predictionInfo, dbNRPinfo[nrpDB])
    save_results(request_id)


@shared_task
def handle_nrp(request_id, predictDB, is_smile=False):
    print("START HANDLE NRP")
    if (is_smile):
        SMILE_to_MOL()
    run_nrpsMatcher(dbPredictionInfo[predictDB], molInfo)
    save_results_prediction(request_id)

@shared_task
def handle_one(request_id, is_smile=False):
    print("START HANDLE ONE")
    if (is_smile):
        SMILE_to_MOL()
    run_antismash()
    run_nrpsMatcher(predictionInfo, molInfo)
    save_results(request_id)