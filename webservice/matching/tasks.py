import os
from .vis_prediction import visualize_prediction
from celery import shared_task
from .vis_prediction import DB_NONE
from .vis_prediction import DB_PNP

path = os.environ.get('DATA_PATH')
genome_file = os.path.join(path, "genome.fna")
nrp_file = os.path.join(path, "nrp.mol")
smile_file = os.path.join(path, "nrp.sml")
predictionInfo = os.path.join(path, "predictions.info")
molInfo = os.path.join(path, "mol.info")
output_folder = os.path.join(path, "out/")

pathToAntismash = os.environ.get('ANTISMASH_PATH')
antismashRes = os.path.join(path, "antismashRes/")
predictionPath = antismashRes + "nrpspks_predictions_txt/ctg1_nrpspredictor2_codes.txt"

Nerpa = "nerpa.py"

dbNRPinfo = {DB_PNP: "/home/dereplicator/kolga/data/DB/PNP/library.info"}
dbPredictionInfo = {'bc': "/home/olga/CAB/NRP/data/DataBase/mibigNF.info"}

def init_var(request_id):
    global genome_file
    global nrp_file
    global smile_file
    global predictionInfo
    global molInfo
    global antismashRes
    global predictionPath
    global output_folder

    res_folder = "res" + str(request_id)
    if not os.path.exists(os.path.join(path, res_folder)):
        os.makedirs(os.path.join(path, res_folder))

    genome_file = os.path.join(path, res_folder, "genome.fna")
    nrp_file = os.path.join(path, res_folder, "nrp.mol")
    smile_file = os.path.join(path, res_folder, "nrp.sml")
    predictionInfo = os.path.join(path, res_folder, "predictions.info")
    molInfo = os.path.join(path, res_folder, "mol.info")
    with open(molInfo, "w") as fw:
        fw.write(nrp_file)
    output_folder = os.path.join(path, res_folder, "out")
    antismashRes = os.path.join(path, res_folder, "antismashRes")
    predictionPath = os.path.join(antismashRes, "nrpspks_predictions_txt/ctg1_nrpspredictor2_codes.txt")


def run_antismash():
    print(os.environ['PATH'])
    print("python2 " + pathToAntismash + " " + genome_file + " --outputfolder "  + antismashRes)
    os.system("python2 " + pathToAntismash + " " + genome_file + " --outputfolder " + antismashRes)
    with open(predictionInfo, "w") as fw:
        fw.write(predictionPath)


def run_nrpsMatcher(prinfo, molinfo):
    print(Nerpa + " -p " + prinfo + " --lib_info " + molinfo + " --predictor NRPSPREDICTOR2 --insertion --deletion --single_match --single_match_coeff 0.2 --modification -o " + output_folder)
    os.system(Nerpa + " -p " + prinfo + " --lib_info " + molinfo + " --predictor NRPSPREDICTOR2 --insertion --deletion --single_match --single_match_coeff 0.2 --modification -o " + output_folder)


def save_results(request_id, nrpDB = DB_NONE):
    for filename in os.listdir(output_folder + "/details_mols"):
        visualize_prediction(output_folder + "/details_mols/" + filename, predictionPath, filename, "ctg1_nrpspredictor2_codes", request_id, nrpDB)


def getAllPath(filename):
    paths = []
    with open(output_folder + "/details_mols/" + filename) as f:
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
    for filename in os.listdir(output_folder + "/details_mols"):
        predpaths = getAllPath(filename)
        for predpath in predpaths:
            if (predpath[-1] == '\n'):
                predpath = predpath[:-1]
            visualize_prediction(output_folder + "/details_mols/" + filename, predpath, "nrp", predpath, request_id, DB_NONE)


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
    os.system("rm -r " + output_folder)

@shared_task
def handle_genome(request_id, nrpDB):
    print("START HANDLE GENOME")
    init_var(request_id)
    run_antismash()
    run_nrpsMatcher(predictionInfo, dbNRPinfo[nrpDB])
    save_results(request_id, nrpDB)
    #clear()


@shared_task
def handle_nrp(request_id, predictDB, is_smile=False):
    print("START HANDLE NRP")
    init_var(request_id)
    if (is_smile):
        SMILE_to_MOL()
    update_mol_info(molInfo)
    run_nrpsMatcher(dbPredictionInfo[predictDB], molInfo)
    save_results_prediction(request_id)
    #clear()

@shared_task
def handle_one(request_id, is_smile=False):
    print("START HANDLE ONE")
    init_var(request_id)
    if (is_smile):
        SMILE_to_MOL()
    update_mol_info(molInfo)
    run_antismash()
    run_nrpsMatcher(predictionInfo, molInfo)
    save_results(request_id, DB_NONE)
    #clear()
