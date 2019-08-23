import os

path_to_config = "/home/olga/CAB/NRP/soft/dereplicator_build/share/npdtools/"
path_file = "/home/olga/CAB/NRP/data/DataBase/MiBig_2019/graphs.info"
graphs_folder = "/home/olga/CAB/NRP/data/DataBase/MiBig_2019/graphs/"
path_to_soft = "/home/olga/CAB/NRP/soft/dereplicator_build/bin/print_structure"

lib_info_file = "/home/olga/CAB/NRP/data/DataBase/MiBig_2019/mibig2019.info"

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
        os.system(path_to_soft + " " + file + " -C " + path_to_config + " --print_rule_fragmented_graph > " + graphs_folder + nfname)
        f.write((graphs_folder+nfname+ " " + info))

f.close()
