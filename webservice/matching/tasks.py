import os
from .vis_prediction import visualize_prediction
from celery import shared_task
from .vis_prediction import DB_NONE
from .vis_prediction import DB_STREPTOME

path = "/home/olga/CAB/NRP/data/serverRuns/"
genome_file = path + "genome.fna"
nrp_file = path + "nrp.mol"
smile_file = path + "nrp.sml"
predictionInfo = path + "predictions.info"
molInfo = path + "mol.info"

pathToAntismash = "/home/olga/CAB/NRP/soft/antismash-3.0.5.1_CUT/run_antismash.py"
antismashRes = path + "antismashRes/"
predictionPath = antismashRes + "nrpspks_predictions_txt/ctg1_nrpspredictor2_codes.txt"

NRPsMatcher = "/home/olga/tmp/NRP/bin/run_nrp_matcher.py"

dbNRPinfo = {DB_STREPTOME: "/home/olga/CAB/NRP/data/DataBase/library.info.streptomedb"}
dbPredictionInfo = {'bc': "/home/olga/CAB/NRP/data/DataBase/mibigNF.info"}

def run_antismash():
    os.system("python2 " + pathToAntismash + " " + genome_file + " --outputfolder " + antismashRes)
    with open(predictionInfo, "w") as fw:
        fw.write(predictionPath)


def run_nrpsMatcher(prinfo, molinfo):
    print("python3 " + NRPsMatcher + " -p " + prinfo + " --lib_info " + molinfo + " -o " + path)
    os.system("python3 " + NRPsMatcher + " -p " + prinfo + " --lib_info " + molinfo + " -o " + path)


def save_results(request_id, nrpDB = DB_NONE):
    for filename in os.listdir(path + "/details_mols"):
        visualize_prediction(path + "/details_mols/" + filename, predictionPath, filename, "ctg1_nrpspredictor2_codes", request_id, nrpDB)


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
            visualize_prediction(path + "/details_mols/" + filename, predpath, "nrp", predpath, request_id, DB_NONE)


def SMILE_to_MOL():
    #print("molconvert mol:V3+H --smiles \"" + smile_file + "\" -o " + nrp_file)
    #os.system("molconvert mol:V3+H --smiles \"" + smile_file + "\" -o " + nrp_file)
    print("molconvert mol:V3+H " + smile_file + " -o " + nrp_file)
    os.system("molconvert mol:V3+H " + smile_file + " -o " + nrp_file)


MASS = {' C ': 12.011, ' N ': 14.007, ' O ': 16, ' Cl ': 35.453, ' P ': 30.974, ' S ': 32.065, ' H ': 1.008}

def get_mass(mol_file):
    mass = 0
    with open(mol_file) as f:
        for line in f:
            for key in MASS:
                if (key in line):
                    mass += MASS[key]
    return mass

def update_mol_info(mol_info):
    os.system("mv " + mol_info + " mol_tmp.info")
    with open(mol_info, 'w') as fw:
        with open("mol_tmp.info") as f:
            for line in f:
                fname = line.split()[0]
                ms = get_mass(fname)
                fw.write(fname + " " + str(ms) + '\n')

    os.system("rm mol_tmp.info")


def clear():
    os.system("rm -r " + path + "/details_mols/")
    os.system("rm -r " + path + "/details_predictions/")
    os.system("rm -r " + path + "/graphs/")
    os.system("rm -r " + path + "/antismashRes/nrpspks_predictions_txt/")
    os.system("rm " + path + "/genome.fna")
    os.system("rm " + path + "/nrp.mol")
    os.system("rm " + path + "/nrp.sml")
    os.system("rm " + path + "/path_to_graphs")
    os.system("rm " + path + "/report.csv")
    os.system("rm " + path + "/report_mols")

@shared_task
def handle_genome(request_id, nrpDB):
    print("START HANDLE GENOME")
    run_antismash()
    run_nrpsMatcher(predictionInfo, dbNRPinfo[nrpDB])
    save_results(request_id, nrpDB)
    clear()


@shared_task
def handle_nrp(request_id, predictDB, is_smile=False):
    print("START HANDLE NRP")
    if (is_smile):
        SMILE_to_MOL()
    update_mol_info(molInfo)
    run_nrpsMatcher(dbPredictionInfo[predictDB], molInfo)
    save_results_prediction(request_id)
    clear()

@shared_task
def handle_one(request_id, is_smile=False):
    print("START HANDLE ONE")
    if (is_smile):
        SMILE_to_MOL()
    update_mol_info(molInfo)
    run_antismash()
    run_nrpsMatcher(predictionInfo, molInfo)
    save_results(request_id, DB_NONE)
    clear()