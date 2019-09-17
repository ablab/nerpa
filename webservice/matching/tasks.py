import os
import zipfile
from .vis_prediction import visualize_prediction
from celery import shared_task
from .vis_prediction import DB_NONE
from .vis_prediction import DB_PNP

path = os.environ.get('DATA_PATH')
pathToAntismash = os.environ.get('ANTISMASH_PATH')

Nerpa = "nerpa.py"

dbNRPinfo = {DB_PNP: "/home/dereplicator/kolga/data/DB/PNP/library.info"}
dbPredictionInfo = {'bc': "/home/olga/CAB/NRP/data/DataBase/mibigNF.info"}

class GenomeQuery:
    def __init__(self, path_to_genome, genome_id, res_folder, genome_name):
        self.path_to_genome = path_to_genome
        self.genome_id = genome_id
        self.genome_name = genome_name
        self.res_folder = res_folder
        self.antismashRes = os.path.join(self.res_folder, "antismashRes", self.genome_id)
        if not os.path.exists(self.antismashRes):
            os.makedirs(self.antismashRes)
        self.predictionPath = os.path.join(self.antismashRes, "nrpspks_predictions_txt/ctg1_nrpspredictor2_codes.txt")

    def run_antismash(self, predictionInfo):
        cmdline = "python2 " + pathToAntismash + " " + self.path_to_genome + " --outputfolder " + self.antismashRes
        print(cmdline)
        os.system(cmdline)

        with open(predictionInfo, "a+") as fw:
            fw.write(self.predictionPath + "\n")


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

    def set_load_genome_filename(self, filename):
        self.load_genome_filename = filename

    def create_genomes_list(self):
        if (zipfile.is_zipfile(self.genome_file)):
            self.zip_fna_folder = os.path.join(self.res_folder, "genomes_fna")
            with zipfile.ZipFile(self.genome_file) as genomezip:
                genomezip.extractall(path=self.zip_fna_folder)

            cur_id = 0
            for filename in os.listdir(self.zip_fna_folder):
                self.genomes.append(GenomeQuery(os.path.join(self.zip_fna_folder, filename),
                                                "g" + str(cur_id), self.res_folder, os.path.splitext(filename)[0]))
                cur_id += 1
        else:
            self.genomes.append(GenomeQuery(self.genome_file, "g0", self.res_folder, self.load_genome_filename))

    def handle_mols(self):
        if (zipfile.is_zipfile(self.nrp_file)):
            self.zip_mol_folder = os.path.join(self.res_folder, "structures_mol")
            with zipfile.ZipFile(self.nrp_file) as nrpzip:
                nrpzip.extractall(path=self.zip_mol_folder)

            with open(self.molInfo, 'w') as fw:
                for filename in os.listdir(self.zip_mol_folder):
                    fw.write(os.path.join(self.zip_mol_folder, filename) + "\n")


    def handle_query(self):
        self.genomes = []
        self.create_genomes_list()
        for genome in self.genomes:
            genome.run_antismash(self.predictionInfo)

        self.handle_mols()

    def genome_name_by_path(self, path):
        for genome_query in self.genomes:
            if genome_query.predictionPath == path:
                return genome_query.genome_name


def run_nrpsMatcher(prinfo, molinfo, query):
    print(Nerpa + " -p " + prinfo + " --lib_info " + molinfo + " --predictor NRPSPREDICTOR2 --insertion --deletion --single_match --single_match_coeff 0.2 --modification -o " + query.output_folder)
    os.system(Nerpa + " -p " + prinfo + " --lib_info " + molinfo + " --predictor NRPSPREDICTOR2 --insertion --deletion --single_match --single_match_coeff 0.2 --modification -o " + query.output_folder)


def save_results(query, nrpDB = DB_NONE):
    for filename in os.listdir(query.output_folder + "/details_mols"):
        visualize_prediction(query.output_folder + "/details_mols/" + filename, query.predictionPath, filename, "ctg1_nrpspredictor2_codes", query.request_id, nrpDB, "some_genome")


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
            visualize_prediction(query.output_folder + "/details_mols/" + filename, predpath, filename, predpath, query.request_id, DB_NONE, query.genome_name_by_path(predpath))


def SMILE_to_MOL(query):
    #print("molconvert mol:V3+H --smiles \"" + smile_file + "\" -o " + nrp_file)
    #os.system("molconvert mol:V3+H --smiles \"" + smile_file + "\" -o " + nrp_file)
    res_mol = []
    with open(query.smile_file, "r") as f:
        i = 0
        query.mol_folder = os.path.join(query.res_folder, "structures_mol")
        if not os.path.exists(query.mol_folder):
            os.makedirs(query.mol_folder)

        for line in f:
            sml = line.split()[0]
            cmd = "molconvert mol:V3+H --smiles \"" + sml + "\" -o " + os.path.join(query.mol_folder, "sml" + str(i) + ".mol")
            print(cmd)
            res_mol.append(os.path.join(query.mol_folder, "sml" + str(i) + ".mol"))
            os.system(cmd)
            i += 1

    with open(query.molInfo, "w") as fw:
        for line in res_mol:
            fw.write(line + "\n")



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
    #run_antismash(query)
    query.handle_query()
    run_nrpsMatcher(query.predictionInfo, dbNRPinfo[nrpDB], query)
    #save_results(query, nrpDB)
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
def handle_one(requst_id, genome_filename, is_smile=False):
    query = Query(requst_id)
    query.set_load_genome_filename(os.path.splitext(genome_filename)[0])
    print("START HANDLE ONE")
    if (is_smile):
        SMILE_to_MOL(query)
    query.handle_query()
    update_mol_info(query.molInfo)
    #run_antismash(query)
    run_nrpsMatcher(query.predictionInfo, query.molInfo, query)
    save_results_prediction(query)
    #clear(query)
