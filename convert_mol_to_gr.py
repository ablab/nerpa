import os

path_to_db_dir = "/home/olga/bio-project/NRP/soft/dereplicator/antimarin_mibig_streptomedb/mol_dir"
path_file = "path_to_graphs"
graphs_folder = "graphs"
path_to_soft = "/home/olga/bio-project/NRP/soft/dereplicator_build/bin/print_structure"


molfiles = [f for f in os.listdir(path_to_db_dir) if os.path.isfile(os.path.join(path_to_db_dir, f))]

f = open(path_file, 'w')

for file in molfiles:
    nfname = file[:-3] + "gr"
    os.system(path_to_soft + " " + path_to_db_dir + "/" + file + " --print_rule_fragmented_graph > graphs/" + nfname)
    f.write(("graphs/"+nfname+"\n"))


f.close()