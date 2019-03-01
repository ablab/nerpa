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
    parser.add_argument("--lib_info", dest="lib_info", nargs=1, help="path to file with paths to mol files", type=str)
    parser.add_argument("--predictor",
                        dest="predictor",
                        default="NRPSPREDICTOR2",
                        choices=["NRPSPREDICTOR2", "MINOWA"],
                        help="AA domain predictor name (NRPSPREDICTOR2, MINOWA) [default=NRPSPREDICTOR2]",
                        action='store')
    parser.add_argument("--local_output_dir", "-o", nargs=1, help="use this output dir", type=str)
    args = parser.parse_args()
    return args

def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
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

            prfix_to_search = ["./", "../", "../../", "../share/", "../share/npdtools/" ]
            config_folder = "configs"
            fragmintation_rule_folder = "Fragmentation_rule"
            path_to_program = which("print_structure")[:-15]


            for prefix in prfix_to_search:
                if (os.path.exists(path_to_program + prefix + fragmintation_rule_folder)):
                    config_folder = path_to_program + prefix

            print("print_structure " + file + " --print_rule_fragmented_graph -C "+ config_folder + " > " + main_out_dir + "graphs/" + nfname)
            os.system("print_structure " + file + " --print_rule_fragmented_graph -C "+ config_folder + " > " + main_out_dir + "graphs/" + nfname)
            f.write((main_out_dir + "graphs/" + nfname + " " + info + "\n"))
            files_list.append(main_out_dir + "graphs/" + nfname)

    f.close()
    return path_file, files_list


def gen_abs_paths_to_prediction(args):
    predictions = []
    with open(args.predictions[0]) as f:
        for line in f:
            if (line[-1] == '\n'):
                line = line[:-1]
            if line[0] == '/':
                predictions.append(line)
            else:
                line = '/'.join(os.path.abspath(args.predictions[0]).split('/')[:-1]) + "/" + line
                predictions.append(line)

    return predictions

def run(args):
    if (len(sys.argv) == 1):
        parser.print_help()
        sys.exit()

    if (args.predictions == None):
        log.err("None prediction info file provide")
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
    if args.local_output_dir != None:
        main_out_dir = os.path.abspath(args.local_output_dir[0]) + "/"

    path_to_graphs, files_list = gen_graphs_by_mol(args, main_out_dir)
    predictions = gen_abs_paths_to_prediction(args)

    path_to_pred = os.path.abspath(args.predictions[0])
    directory = os.path.dirname(main_out_dir)
    if not os.path.exists(directory):
        os.makedirs(directory)
    os.chdir(directory)

    if not os.path.exists(os.path.dirname('details_predictions/')):
        os.makedirs(os.path.dirname('details_predictions/'))

    if not os.path.exists(os.path.dirname('details_mols/')):
        os.makedirs(os.path.dirname('details_mols/'))

    print(path_to_exec_dir + "/NRPsMatcher \"" +  path_to_pred + "\" \"" + path_to_graphs + "\" " + args.predictor + "\n")
    os.system(path_to_exec_dir + "/NRPsMatcher \"" +  path_to_pred + "\" \"" + path_to_graphs + "\" " + args.predictor + "\n")
    return

args = parse_args()
run(args)