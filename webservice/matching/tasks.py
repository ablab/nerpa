import os
from .vis_prediction import visualize_prediction
from celery import shared_task
from .vis_prediction import DB_NONE
from .vis_prediction import DB_PNP

path = os.environ.get('DATA_PATH')
pathToAntismash = os.environ.get('ANTISMASH_PATH')

Nerpa = "nerpa.py"

dbNRPinfo = {DB_PNP: "/home/dereplicator/kolga/data/DB/PNP/library.info"}
dbPredictionInfo = {'bc': "/home/olga/CAB/NRP/data/DataBase/mibigNF.info"}


class Query:
    def __init__(self, request_id):
        self.request_id = request_id

        self.res_folder = os.path.join(path, "res" + str(request_id))
        if not os.path.exists(self.res_folder):
            os.makedirs(self.res_folder)

        self.genome_file = os.path.join(self.res_folder, "genome.fna")
        self.nrp_file = os.path.join(self.res_folder, "nrp.mol")
        self.smile_file = os.path.join(self.res_folder, "nrp.sml")
        self.predictionInfo = os.path.join(self.res_folder, "predictions.info")
        self.molInfo = os.path.join(self.res_folder, "mol.info")
        with open(self.molInfo, "w") as fw:
            fw.write(self.nrp_file)
        self.output_folder = os.path.join(self.res_folder, "out")
        self.antismashRes = os.path.join(self.res_folder, "antismashRes")
        self.predictionPath = os.path.join(self.antismashRes, "nrpspks_predictions_txt/ctg1_nrpspredictor2_codes.txt")


def run_antismash(query):
    print(os.environ['PATH'])
    print("python2 " + pathToAntismash + " " + query.genome_file + " --outputfolder "  + query.antismashRes)
    os.system("python2 " + pathToAntismash + " " + query.genome_file + " --outputfolder " + query.antismashRes)
    with open(query.predictionInfo, "w") as fw:
        fw.write(query.predictionPath)


def run_nrpsMatcher(prinfo, molinfo, query):
    print(Nerpa + " -p " + prinfo + " --lib_info " + molinfo + " --predictor NRPSPREDICTOR2 --insertion --deletion --single_match --single_match_coeff 0.2 --modification -o " + query.output_folder)
    os.system(Nerpa + " -p " + prinfo + " --lib_info " + molinfo + " --predictor NRPSPREDICTOR2 --insertion --deletion --single_match --single_match_coeff 0.2 --modification -o " + query.output_folder)


def save_results(query, nrpDB = DB_NONE):
    for filename in os.listdir(query.output_folder + "/details_mols"):
        visualize_prediction(query.output_folder + "/details_mols/" + filename, query.predictionPath, filename, "ctg1_nrpspredictor2_codes", query.request_id, nrpDB)


def getAllPath(query, filename):
    paths = []
    with open(query.output_folder + "/details_mols/" + filename) as f:
        lines = f.readlines()
        cur = 0
        while (cur < len(lines)):
            paths.append(lines[cur + 1])
            while (cur < len(lines) and lines[cur] != "\n"):
                cur += 1
            while (cur < len(lines) and lines[cur] == "\n"):
                cur += 1
    return paths


def save_results_prediction(query):
    for filename in os.listdir(query.output_folder + "/details_mols"):
        predpaths = getAllPath(query, filename)
        for predpath in predpaths:
            if (predpath[-1] == '\n'):
                predpath = predpath[:-1]
            visualize_prediction(query.output_folder + "/details_mols/" + filename, predpath, "nrp", predpath, query.request_id, DB_NONE)


def SMILE_to_MOL(query):
    #print("molconvert mol:V3+H --smiles \"" + smile_file + "\" -o " + nrp_file)
    #os.system("molconvert mol:V3+H --smiles \"" + smile_file + "\" -o " + nrp_file)
    print("molconvert mol:V3+H " + query.smile_file + " -o " + query.nrp_file)
    os.system("molconvert mol:V3+H " + query.smile_file + " -o " + query.nrp_file)


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


def clear(query):
    os.system("rm -r " + query.output_folder)

@shared_task
def handle_genome(requst_id, nrpDB):
    query = Query(requst_id)
    print("START HANDLE GENOME")
    run_antismash(query)
    run_nrpsMatcher(query.predictionInfo, dbNRPinfo[nrpDB], query)
    save_results(query, nrpDB)
    #clear(query)


@shared_task
def handle_nrp(requst_id, predictDB, is_smile=False):
    query = Query(requst_id)
    print("START HANDLE NRP")
    if (is_smile):
        SMILE_to_MOL(query)
    update_mol_info(query.molInfo)
    run_nrpsMatcher(dbPredictionInfo[predictDB], query.molInfo, query)
    save_results_prediction(query)
    #clear(query)

@shared_task
def handle_one(requst_id, is_smile=False):
    query = Query(requst_id)
    print("START HANDLE ONE")
    if (is_smile):
        SMILE_to_MOL(query)
    update_mol_info(query.molInfo)
    run_antismash(query)
    run_nrpsMatcher(query.predictionInfo, query.molInfo, query)
    save_results(query, DB_NONE)
    #clear(query)
