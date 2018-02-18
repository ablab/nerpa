import os

path_to_db_dir = "/home/olga/bio-project/NRP/soft/dereplicator/"
path_file = "path_to_graphs"
graphs_folder = "graphs"
path_to_soft = "/home/olga/bio-project/NRP/soft/dereplicator_build/bin/print_structure"

lib_info_file = "library.info.combined"

#molfiles = [f for f in os.listdir(path_to_db_dir) if os.path.isfile(os.path.join(path_to_db_dir, f))]

f = open(path_file, 'w')

with open(lib_info_file) as fr:
    for line in fr:
        print(line)
        line_parts = line.split(' ')
        file = line_parts[0]
        nfname = file.split('/')[-1][:-3] + "gr"
        info = ' '.join(line_parts[1:])
        print(file)
        print(nfname)
        print(info)
        os.system(path_to_soft + " " + path_to_db_dir + "/" + file + " --print_rule_fragmented_graph > graphs/" + nfname)
        f.write(("graphs/"+nfname+ " " + info))

f.close()
