import os
import zipfile
import tempfile
from .vis_prediction import visualize_prediction
from celery import shared_task
from .vis_prediction import DB_NONE
from .vis_prediction import DB_PNP
from nrpsmatche.settings import ANTISMASH_ROOT

path = os.environ.get('DATA_PATH')
pathToAntismash = os.environ.get('ANTISMASH_PATH')

Nerpa = "nerpa.py"

dbNRPinfo = {DB_PNP: "/home/dereplicator/kolga/data/DB/PNP/library.info"}
dbPredictionInfo = {'bc': "/home/olga/CAB/NRP/data/DataBase/mibigNF.info"}

class StructureQuery:
    def __init__(self, mol_id, path_to_mol, res_folder, smile_string="", peptide="", details=""):
        self.molId = mol_id
        self.SMILE = smile_string
        self.path_to_mol = path_to_mol
        self.peptide = peptide
        self.details = details

        if smile_string == "":
            smile_file = tempfile.NamedTemporaryFile(dir=res_folder, delete=False)
            smile_file_name = smile_file.name
            smile_file.close()
            cmd = "molconvert smiles " + self.path_to_mol + " -o " + smile_file_name
            print(cmd)
            if os.system(cmd) != 0:
                raise Exception("Convertion MOL to SMILE failed")
            self.SMILE = open(smile_file_name).read()
            os.remove(smile_file_name)
        print(self.molId, self.SMILE)


class GenomeQuery:
    def __init__(self, path_to_genome, genome_id, res_folder, genome_name, request_id):
        self.path_to_genome = path_to_genome
        self.genome_id = genome_id
        self.genome_name = genome_name
        self.res_folder = res_folder
        self.antismashRes = os.path.join(ANTISMASH_ROOT, "res" + str(request_id), self.genome_id)
        if not os.path.exists(self.antismashRes):
            os.makedirs(self.antismashRes)
        self.predictionPath = os.path.join(self.antismashRes, "nrpspks_predictions_txt/ctg1_nrpspredictor2_codes.txt")

    def run_antismash(self, predictionInfo):
        cmdline = "python2 " + pathToAntismash + " " + self.path_to_genome + " --outputfolder " + self.antismashRes
        print(cmdline)
        if os.system(cmdline) != 0:
            raise Exception("Running antismash failed")

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
        self.is_smile = False
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
                                                "g" + str(cur_id), self.res_folder, os.path.splitext(filename)[0], self.request_id))
                cur_id += 1
        else:
            self.genomes.append(GenomeQuery(self.genome_file, "g0", self.res_folder, self.load_genome_filename, self.request_id))

    def SMILE_to_MOL(self):
        # print("molconvert mol:V3+H --smiles \"" + smile_file + "\" -o " + nrp_file)
        # os.system("molconvert mol:V3+H --smiles \"" + smile_file + "\" -o " + nrp_file)
        import re
        res_mol = []
        with open(self.smile_file, "r") as f:
            i = 0
            self.mol_folder = os.path.join(self.res_folder, "structures_mol")
            if not os.path.exists(self.mol_folder):
                os.makedirs(self.mol_folder)

            for line in f:
                line_parts = re.split(r'\t+', line)
                if line_parts[0] == "SMILES":
                    continue

                name = "sml" + str(i)
                if len(line_parts) > 1:
                    name = line_parts[1]

                peptide = ""
                if len(line_parts) > 2:
                    peptide = line_parts[2]

                details = ""
                if len(line_parts) > 3:
                    details = line_parts[3]


                sml = line_parts[0]
                cmd = "molconvert mol:V3+H --smiles \"" + sml + "\" -o " + os.path.join(self.mol_folder,
                                                                                        name + ".mol")
                print(cmd)
                res_mol.append(os.path.join(self.mol_folder, name + ".mol"))


                self.structures.append(StructureQuery(name,
                                                      os.path.join(self.mol_folder, name + ".mol"),
                                                      self.res_folder, sml, peptide=peptide, details=details))
                if os.system(cmd) != 0:
                    raise Exception("Convertion SMILES to MOL failed")
                i += 1

        with open(self.molInfo, "w") as fw:
            for line in res_mol:
                fw.write(line + "\n")


    def handle_mols(self):
        if (zipfile.is_zipfile(self.nrp_file)):
            self.zip_mol_folder = os.path.join(self.res_folder, "structures_mol")
            with zipfile.ZipFile(self.nrp_file) as nrpzip:
                nrpzip.extractall(path=self.zip_mol_folder)

            with open(self.molInfo, 'w') as fw:
                for filename in os.listdir(self.zip_mol_folder):
                    self.structures.append(StructureQuery('.'.join(filename.split('.')[:-1]), os.path.join(self.zip_mol_folder, filename), self.res_folder))
                    fw.write(os.path.join(self.zip_mol_folder, filename) + "\n")
        elif (self.is_smile):
            self.SMILE_to_MOL()
        else:
            self.structures.append(StructureQuery("nrp", self.nrp_file, self.res_folder))

    def handle_query(self):
        self.genomes = []
        self.structures = []
        self.create_genomes_list()
        for genome in self.genomes:
            genome.run_antismash(self.predictionInfo)

        self.handle_mols()

    def genome_name_by_path(self, path):
        for genome_query in self.genomes:
            if genome_query.predictionPath == path:
                return genome_query.genome_name

    def genome_id_by_path(self, path):
        for genome_query in self.genomes:
            if genome_query.predictionPath == path:
                return genome_query.genome_id

    def get_SMILE_by_mol_id(self, molid):
        for structure_query in self.structures:
            if structure_query.molId == molid:
                return structure_query.SMILE

    def get_structure_name_by_path(self, path):
        return path.split('/')[-1][:-3]

    def getStructureQueryByName(self, name):
        for structure_query in self.structures:
            if structure_query.molId == name:
                return structure_query


def run_nrpsMatcher(prinfo, molinfo, query):
    cmd = Nerpa + " -p " + prinfo + " --lib_info " + molinfo + " --predictor NRPSPREDICTOR2 --insertion --deletion --single_match --single_match_coeff 0.2 --modification -o " + query.output_folder
    print(cmd)
    if os.system(cmd) != 0:
        raise Exception("Running Nerpa failed")


def save_results(query, nrpDB = DB_NONE):
    for filename in os.listdir(query.output_folder + "/details_mols"):
        visualize_prediction(query.output_folder + "/details_mols/" + filename, query.predictionPath, filename, "ctg1_nrpspredictor2_codes", query.request_id, nrpDB, "some_genome",
                             "res" + str(query.request_id), query.get_SMILE_by_mol_id(filename))


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
            structure_query = query.getStructureQueryByName(filename)
            visualize_prediction(query.output_folder + "/details_mols/" + filename, predpath, filename, predpath,
                                 query.request_id, DB_NONE, query.genome_name_by_path(predpath),
                                 os.path.join("res" + str(query.request_id), query.genome_id_by_path(predpath), "index.html"),
                                 query.get_SMILE_by_mol_id(filename), peptide=structure_query.peptide,
                                 details_structure=structure_query.details)


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
    query.is_smile = is_smile

    print("START HANDLE NRP")
    update_mol_info(query.molInfo)
    run_nrpsMatcher(dbPredictionInfo[predictDB], query.molInfo, query)
    save_results_prediction(query)
    #clear(query)

@shared_task
def handle_one(requst_id, genome_filename, is_smile=False):
    query = Query(requst_id)
    query.set_load_genome_filename(os.path.splitext(genome_filename)[0])
    query.is_smile = is_smile

    print("START HANDLE ONE")
    query.handle_query()
    update_mol_info(query.molInfo)
    #run_antismash(query)
    run_nrpsMatcher(query.predictionInfo, query.molInfo, query)
    save_results_prediction(query)
    #clear(query)
