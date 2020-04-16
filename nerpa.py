#!/usr/bin/python3
import sys
import os
import argparse
from shutil import copyfile

#Log class, use it, not print
class Log:
    text = ""

    def log(self, s):
        self.text += s + "\n"
        print(s)

    def warn(self, s):
        msg = "WARNING: " + s
        self.text += msg + "\n"
        sys.stdout.write(msg)
        sys.stdout.flush()

    def err(self, s):
        msg = "ERROR: " + s + "\n"
        self.text += msg
        sys.stdout.write(msg)
        sys.stdout.flush()

    def print_log(self):
        print(self.text)

    def get_log(self):
        return self.text

log = Log()

path_to_exec_dir = os.path.dirname(os.path.abspath(__file__)) + "/"

def parse_args():
    global parser
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--predictions", "-p", nargs=1, dest="predictions", help="path to file with paths to prediction files", type=str)
    parser.add_argument("--lib_info", "-l", dest="lib_info", nargs=1, help="path to file with paths to mol files", type=str)
    parser.add_argument("--predictor",
                        dest="predictor",
                        default="MINOWA",
                        choices=["NRPSPREDICTOR2", "MINOWA", "PRISM", "SANDPUMA"],
                        help="AA domain predictor name [default=MINOWA]",
                        action='store')
    parser.add_argument("--antismash_output_list", dest="antismash_out", help="path to file with list of paths to antiSMASH output folders", type=str)
    parser.add_argument("--insertion", help="insertion score [default=-1]", default=-1, action="store")
    parser.add_argument("--deletion", help="deletion score [default=-1]", default=-1, action="store")
    parser.add_argument("--open_gap", default=-1, type=float, help="score for opening gap in NRP structure [default=-1]", action="store")
    parser.add_argument("--continue_gap", default=-0.1, type=float, help="score for continue gap in NRP structure [default=-0.1]", action="store")
    parser.add_argument("--mismatch", default=-1, type=float, help="mismatche score [default=-1]", action="store")
    parser.add_argument("--skip_segment", default=-1, type=float, help="score for skip orf in prediction [default=-1]", action="store")
    parser.add_argument("--modification", help="allow modification", action="store_true")
    parser.add_argument("--modification_cfg", help="path to file with modification description", action="store", type=str)
    parser.add_argument("--AAmod_cfg", help="path to file with modification for specific AA description", action="store", type=str)
    parser.add_argument("--threads", default=1, type=int, help="number of threads for running Nerpa", action="store")
    parser.add_argument("--local_output_dir", "-o", nargs=1, help="use this output dir", type=str)
    args = parser.parse_args()
    return args

def print_cfg(args, output_dir):
    cfg_file = os.path.join(output_dir, "nerpa.cfg")
    with open(cfg_file, "w") as f:
        f.write(args.predictor + "\n")
        f.write("insertion " + str(args.insertion) + "\n")
        f.write("deletion " + str(args.deletion) + "\n")
        f.write("open_gap " + str(args.open_gap) + "\n")
        f.write("continue_gap " + str(args.continue_gap) + "\n")
        f.write("mismatch " + str(args.mismatch) + "\n")
        f.write("skip_segment " + str(args.skip_segment) + "\n")

        if args.modification:
            f.write("modification on\n")
        else:
            f.write("modification off\n")

        f.write(os.path.abspath(os.path.join(output_dir, "modifications.tsv")) + "\n")
        f.write(os.path.abspath(os.path.join(output_dir, "AAmod.tsv")) + "\n")
        f.write("threads " + str(args.threads) + "\n")

    return cfg_file

def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)

    cur_path = os.path.dirname(os.path.abspath(__file__))
    if is_exe(os.path.join(cur_path, fname)):
        return os.path.join(cur_path, fname)

    if fpath:
        if is_exe(program):
            return program
    elif "PATH" in os.environ:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

def gen_graphs_by_mol(args, main_out_dir):
    if not os.path.exists(os.path.dirname(main_out_dir + 'graphs/')):
        os.makedirs(os.path.dirname(main_out_dir + 'graphs/'))

    path_file = main_out_dir + "/path_to_graphs"
    files_list = []

    f = open(path_file, 'w')
    with open(args.lib_info[0]) as fr:
        for line in fr:
            line_parts = line.split()
            file = line_parts[0]
            if (file[0] != '/'):
                file = '/'.join(os.path.abspath(args.lib_info[0]).split('/')[:-1]) + "/" + file

            nfname = file.split('/')[-1][:-3] + "gr"
            info = ' '.join(line_parts[1:])

            prfix_to_search = ["./", "../", "../../", "../share/", "../share/npdtools/", "../share/nerpa/" ]
            config_folder = "configs"
            fragmintation_rule_folder = "Fragmentation_rule"
            path_to_program = which("print_structure")[:-15]


            for prefix in prfix_to_search:
                if (os.path.exists(path_to_program + prefix + fragmintation_rule_folder)):
                    config_folder = path_to_program + prefix

            print(os.path.join(path_to_program, "print_structure") + " " + file + " --print_rule_fragmented_graph -C "+ config_folder + " > " + main_out_dir + "graphs/" + nfname)
            os.system(os.path.join(path_to_program, "print_structure") + " " + file + " --print_rule_fragmented_graph -C "+ config_folder + " > " + main_out_dir + "graphs/" + nfname)
            f.write((main_out_dir + "graphs/" + nfname + " " + info + "\n"))
            files_list.append(main_out_dir + "graphs/" + nfname)

    f.close()
    return path_file, files_list

def copy_prediction_list(args, main_out_dir):
    new_prediction_path = os.path.join(main_out_dir, "prediction.info")
    f = open(new_prediction_path, 'w')
    with open(args.predictions[0]) as fr:
        for line in fr:
            line_parts = line.split()
            file = line_parts[0]
            if (file[0] != '/'):
                file = os.path.join(os.path.dirname(os.path.abspath(args.predictions[0])), file)
            f.write(file + "\n")
    f.close()
    return new_prediction_path

def run(args):
    path_to_cur = os.path.dirname(os.path.abspath(__file__))

    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    if (args.predictions is None) and (args.antismash_out is None):
        log.err("None prediction info file provide")
        parser.print_help()
        sys.exit()

    if (args.predictions is not None) and (args.antismash_out is not None):
        log.err("You cann't use --predictions and --antismash_output_list simultaneously")
        parser.print_help()
        sys.exit()

    if (args.lib_info == None):
        log.err("None NRP structure info file provide")
        parser.print_help()
        sys.exit()
    if (which("print_structure") == None):
        log.err("dereplicator not found. Please install dereplicator and add it to PATH.")
        sys.exit()

    main_out_dir = os.path.abspath(".") + "/"
    if args.local_output_dir is not None:
        main_out_dir = os.path.abspath(args.local_output_dir[0]) + "/"

    path_to_graphs, files_list = gen_graphs_by_mol(args, main_out_dir)
    path_predictions = os.path.abspath(copy_prediction_list(args, main_out_dir))

    directory = os.path.dirname(main_out_dir)
    if not os.path.exists(directory):
        os.makedirs(directory)
    os.chdir(directory)

    path_to_cfg = print_cfg(args, main_out_dir)

    if not os.path.exists(os.path.dirname('details_mols/')):
        os.makedirs(os.path.dirname('details_mols/'))

    path_to_AA = "./resources/aminoacids.tsv"
    if (os.path.exists(os.path.join(path_to_cur, 'NRPsMatcher'))):
        path_to_AA = "../share/nerpa/aminoacids.tsv"
    path_to_AA = os.path.join(path_to_cur, path_to_AA)

    path_to_modification_cfg = "./resources/modifications.tsv"
    if (os.path.exists(os.path.join(path_to_cur, 'NRPsMatcher'))):
        path_to_modification_cfg = "../share/nerpa/modifications.tsv"
    path_to_modification_cfg = os.path.join(path_to_cur, path_to_modification_cfg)
    if args.modification_cfg is not None:
        path_to_modification_cfg = os.path.abspath(args.modification_cfg)

    local_modifications_cfg = os.path.join(main_out_dir, "modifications.tsv")
    copyfile(path_to_modification_cfg, local_modifications_cfg)

    path_to_AAmod_cfg = "./resources/AAmod.tsv"
    if (os.path.exists(os.path.join(path_to_cur, 'NRPsMatcher'))):
        path_to_AAmod_cfg = "../share/nerpa/AAmod.tsv"
    path_to_AAmod_cfg = os.path.join(path_to_cur, path_to_AAmod_cfg)
    if args.AAmod_cfg is not None:
        path_to_AAmod_cfg = os.path.abspath(args.AAmod_cfg)

    local_AAmod_cfg = os.path.join(main_out_dir, "AAmod.tsv")
    copyfile(path_to_AAmod_cfg, local_AAmod_cfg)

    comand = path_to_exec_dir + "/NRPsMatcher \"" +  path_predictions + "\" \"" + path_to_graphs + "\" \"" + path_to_AA + "\" \"" + path_to_cfg + "\"\n"
    print(comand)
    os.system(comand)
    return

args = parse_args()
run(args)