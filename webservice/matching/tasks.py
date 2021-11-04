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
    def __init__(self, id, path_to_mol, res_folder, smile_string="", peptide="", details=""):
        self.Id = id
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
        print(self.Id, self.SMILE)


class GenomeQuery:
    def __init__(self, path_to_genome, genome_id, res_folder, genome_name, request_id, organism="", details=""):
        self.path_to_genome = path_to_genome
        self.genome_id = genome_id
        self.genome_name = genome_name
        self.res_folder = res_folder
        self.organism = organism
        self.details = details
        self.antismashRes = os.path.join(ANTISMASH_ROOT, "res" + str(request_id), self.genome_id)
        if not os.path.exists(self.antismashRes):
            os.makedirs(self.antismashRes)
        self.predictionPath = os.path.join(self.antismashRes, "nrpspks_predictions_txt/ctg1_nrpspredictor2_codes.txt")

    def generate_genome_cluster_list(self):
        import copy

        genome_cluster_list = []
        res_folder = os.path.join(self.antismashRes, "nrpspks_predictions_txt")
        for file in os.listdir(res_folder):
            if "nrpspredictor2_codes.txt" in file:
                genome_cluster_list.append(copy.copy(self))
                genome_cluster_list[-1].predictionPath = os.path.join(res_folder, file)
                genome_cluster_list[-1].cluster = int(file.split('_')[0][3:])
        return genome_cluster_list

    def run_antismash(self, predictionInfo):
        cmdline = "python2 " + pathToAntismash + " " + self.path_to_genome + " --outputfolder " + self.antismashRes
        print(cmdline)
        if os.system(cmdline) != 0:
            raise Exception("Running antismash failed")

        genome_cluster_list = self.generate_genome_cluster_list()

        with open(predictionInfo, "a+") as fw:
            for genome_cluster in genome_cluster_list:
                fw.write(genome_cluster.predictionPath + "\n")
        return genome_cluster_list


class Query:
    def __init__(self, request_id):
        self.request_id = request_id

        self.res_folder = os.path.join(path, "res" + str(request_id))
        if not os.path.exists(self.res_folder):
            os.makedirs(self.res_folder)

        self.genome_file = os.path.join(self.res_folder, "genome.fna")
        self.nrp_file = os.path.join(self.res_folder, "nrp.csv")
        self.predictionInfo = os.path.join(self.res_folder, "predictions.info")
        self.output_folder = os.path.join(self.res_folder, "out")
        self.antiSMASHout = os.path.join(ANTISMASH_ROOT, "res" + str(request_id))

    def set_load_genome_filename(self, filename):
        self.load_genome_filename = filename

    def create_genomes_list(self):
        if (zipfile.is_zipfile(self.genome_file)):
            self.zip_fna_folder = os.path.join(self.res_folder, "genomes_fna")
            with zipfile.ZipFile(self.genome_file) as genomezip:
                genomezip.extractall(path=self.zip_fna_folder)

            info_filename = ""

            cur_id = 0
            for filename in os.listdir(self.zip_fna_folder):
                if filename[-3:] == 'tsv':
                    info_filename = filename
                    continue

                self.genomes.append(GenomeQuery(os.path.join(self.zip_fna_folder, filename),
                                                "g" + str(cur_id), self.res_folder, os.path.splitext(filename)[0], self.request_id))
                cur_id += 1

            if info_filename != "":
                import re
                with open(os.path.join(self.zip_fna_folder, info_filename)) as f:
                    for line in f:
                        line_parts = re.split(r'\t+', line)
                        if line_parts[0] == 'NAME':
                            continue

                        if line_parts[0].split('.')[-1] in ["fa", "fasta", "fna"]:
                            line_parts[0] = '.'.join(line_parts.split('.')[:-1])

                        for i in range(len(self.genomes)):
                            if self.genomes[i].genome_name == line_parts[0]:
                                if len(line_parts) > 1:
                                    self.genomes[i].organism = line_parts[1]
                                if len(line_parts) > 2:
                                    self.genomes[i].details = line_parts[2]

        else:
            self.genomes.append(GenomeQuery(self.genome_file, "g0", self.res_folder, self.load_genome_filename, self.request_id))

    def handle_structurs(self):
        with open(self.nrp_file) as f:
            for line in f:
                self.structures.append(StructureQuery(line.split('\t')[1], line.split('\t')[1], self.res_folder, line.split('\t')[0]))

    def handle_query(self):
        self.genomes = []
        self.structures = []
        self.create_genomes_list()
        cluster_genomes_list = []
        for genome in self.genomes:
            cluster_genomes_list += genome.run_antismash(self.predictionInfo)

        self.genomes = cluster_genomes_list
        self.handle_structurs()

    def genome_name_by_path(self, path):
        for genome_query in self.genomes:
            if genome_query.predictionPath == path:
                return genome_query.genome_name

    def genome_id_by_path(self, path):
        for genome_query in self.genomes:
            if genome_query.predictionPath == path:
                return genome_query.genome_id

    def get_SMILE_by_id(self, id):
        for structure_query in self.structures:
            if structure_query.Id in id:
                return structure_query.SMILE

    def get_structure_name_by_path(self, path):
        return path.split('/')[-1][:-3]

    def getGenomeQueryByPath(self, path):
        for genome_query in self.genomes:
            #print("antismashREs and path", genome_query.antismashRes, path, (genome_query.antismashRes in path))
            if (str(genome_query.genome_id) + "_") in path:
                return genome_query

    def getStructureQueryByName(self, name):
        for structure_query in self.structures:
            print(structure_query.Id, name)
            if structure_query.Id in name:
                return structure_query


def run_nrpsMatcher(antiSMASHout, sturcturesCsv, query):
    cmd = Nerpa + " -a " + antiSMASHout + " --smiles-tsv " + sturcturesCsv + \
          " --col-smiles \"SMILES\" --col-id \"NAME\" -o " + query.output_folder
    print(cmd)
    if os.system(cmd) != 0:
        raise Exception("Running Nerpa failed")


def save_results(query, nrpDB = DB_NONE):
    for filename in os.listdir(query.output_folder + "/details"):
        visualize_prediction(query.output_folder + "/details/" + filename, query.predictionPath, filename, "ctg1_nrpspredictor2_codes", query.request_id, nrpDB, "some_genome",
                             "res" + str(query.request_id), query.get_SMILE_by_mol_id(filename))


def getAllPath(query, filename):
    paths = []
    with open(query.output_folder + "/details/" + filename) as f:
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
    for filename in os.listdir(query.output_folder + "/details"):
        predpaths = getAllPath(query, filename)
        for predpath in predpaths:
            if (predpath[-1] == '\n'):
                predpath = predpath[:-1]
            structure_query = query.getStructureQueryByName(filename)
            genome_query = query.getGenomeQueryByPath(predpath)
            visualize_prediction(query.output_folder + "/details/" + filename, predpath, filename, predpath,
                                 query.request_id, DB_NONE, genome_query.genome_name,
                                 os.path.join("res" + str(query.request_id), genome_query.genome_id, "index.html"),
                                 query.get_SMILE_by_id(filename), peptide=structure_query.peptide,
                                 details_structure=structure_query.details, organism=genome_query.organism,
                                 details_genome=genome_query.details, cluster=genome_query.cluster)


MASS = {' C ': 12.011, ' N ': 14.007, ' O ': 16, ' Cl ': 35.453, ' P ': 30.974, ' S ': 32.065, ' H ': 1.008}


def get_mass(mol_file):
    mass = 0
    with open(mol_file) as f:
        for line in f:
            for key in MASS:
                if (key in line):
                    mass += MASS[key]
    return mass


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
def handle_one(requst_id, genome_filename):
    query = Query(requst_id)
    query.set_load_genome_filename(os.path.splitext(genome_filename)[0])

    print("START HANDLE ONE")
    query.handle_query()
    run_nrpsMatcher(query.antiSMASHout, query.nrp_file, query)
    save_results_prediction(query)
    #clear(query)
