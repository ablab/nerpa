#!/usr/bin/python3
import sys
import os
import argparse
import csv
from shutil import copyfile

import nerpa_init

nerpa_init.init()
import handle_TE
from logger import log

path_to_exec_dir = os.path.dirname(os.path.abspath(__file__)) + "/"

def parse_args():
    global parser
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--predictions", "-p", nargs=1, dest="predictions",
                        help="path to file with paths to prediction files", type=str)
    parser.add_argument("--smiles", dest="smiles", nargs=1,
                        help="string with smiles of a single compound", type=str)
    parser.add_argument("--smiles-tsv", dest="smiles_tsv", nargs=1,
                        help="csv with named columns", type=str)
    parser.add_argument("--col_smiles", dest="col_smiles", nargs=1,
                        help="name of column with smiles", type=str, default='SMILES')
    parser.add_argument("--col_id", dest="col_id", nargs=1,
                        help="name of col with ids (if not provided, id=[number of row])", type=str)
    parser.add_argument("--sep", dest="sep", nargs=1,
                        help="separator in smiles-tsv", type=str, default='\t')

    parser.add_argument("--antismash_output_list", dest="antismash_out", help="path to file with list of paths to antiSMASH output folders", type=str)
    parser.add_argument("--insertion", help="insertion score [default=-1]", default=-1, action="store")
    parser.add_argument("--deletion", help="deletion score [default=-5]", default=-5, action="store")
    parser.add_argument("--modification_cfg", help="path to file with modification description", action="store", type=str)
    parser.add_argument("--monomer_cfg", help="path to file with monomer description", action="store", type=str)
    parser.add_argument("--threads", default=1, type=int, help="number of threads for running Nerpa", action="store")
    parser.add_argument("--local_output_dir", "-o", nargs=1, help="use this output dir", type=str)
    args = parser.parse_args()
    return args


def check_tsv_ids_duplicates(reader, col_id):
    ids = [row[col_id] for row in reader]
    return len(ids) == len(set(ids))

def check_if_only_one_is_present(arg1, arg2):
    return not (arg1 is not None and arg2 is not None)

def check_if_at_least_one_is_present(arg1, arg2):
    return arg1 is not None or arg2 is not None

class ValidationError(Exception):
    pass

def validate(expr, msg=''):
    if not expr:
        raise ValidationError(msg)

def validate_arguments(args):
    try:
        validate(len(sys.argv) > 1)

        validate(check_if_at_least_one_is_present(args.predictions, args.antismash_out),
            "None prediction info file provide")

        validate(check_if_only_one_is_present(args.predictions, args.antismash_out),
            "You cann't use --predictions and --antismash_output_list simultaneously")

        validate(check_if_at_least_one_is_present(args.smiles, args.smiles_tsv),
            "No structures provided")

        validate(check_if_only_one_is_present(args.smiles, args.smiles_tsv),
            "Simultaneous use of --smiles and --smiles-tsv is not allowed.")

        try:
            with open(args.smiles_tsv, newline='') as f_in:
                reader = csv.DictReader(f_in, delimiter=args.sep, quoting=csv.QUOTE_NONE)
                validate(args.col_smiles in reader.fieldnames,
                    f'Column "{args.col_smiles}" was specified but does not exist in {args.smiles_tsv}.')
                if args.col_id:
                    validate(args.col_id in reader.fieldnames,
                        f'Column "{args.col_id}" was specified but does not exist in {args.smiles_tsv}.')
                    validate(check_tsv_ids_duplicates(reader, args.col_id),
                        f'Duplicate IDs are found.')
        except FileNotFoundError:
            raise ValidationError(f'No such file: "{args.smiles_tsv}".')
        except csv.Error as e:
            raise ValidationError(f'Cannot parse "{args.smiles_tsv}": {e}.')
        except Exception as e:
            raise ValidationError(f'Invalid input file "{args.smiles_tsv}": {e}.')

    except ValidationError as e:
        if str(e):
            log.err(e)
        parser.print_help()
        sys.exit()


def print_cfg(args, output_dir):
    cfg_file = os.path.join(output_dir, "nerpa.cfg")
    with open(cfg_file, "w") as f:
        f.write("insertion " + str(args.insertion) + "\n")
        f.write("deletion " + str(args.deletion) + "\n")
        f.write(os.path.abspath(os.path.join(output_dir, "modifications.tsv")) + "\n")
        f.write(os.path.abspath(os.path.join(output_dir, "monomers.tsv")) + "\n")
        f.write(os.path.abspath(os.path.join(output_dir, "monomersLogP.tsv")) + "\n")
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

    validate_arguments(args)

    main_out_dir = os.path.abspath(".") + "/"
    if args.local_output_dir is not None:
        main_out_dir = os.path.abspath(args.local_output_dir[0]) + "/"

    if not os.path.exists(main_out_dir):
        os.makedirs(main_out_dir)

    path_to_graphs = os.path.join(main_out_dir, 'path_to_graphs')
    copyfile(args.smiles_tsv[0], path_to_graphs)


    if (args.predictions is not None):
        path_predictions = os.path.abspath(copy_prediction_list(args, main_out_dir))
    else:
        path_predictions = handle_TE.create_predictions_by_antiSAMSHout(args.antismash_out, main_out_dir)

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


    path_to_monomer_logP = "./resources/monomersLogP.tsv"
    path_to_monomer_cfg = "./resources/monomers.tsv"
    if (os.path.exists(os.path.join(path_to_cur, 'NRPsMatcher'))):
        path_to_monomer_logP =  "../share/nerpa/monomersLogP.tsv"
        path_to_monomer_cfg = "../share/nerpa/monomers.tsv"
    path_to_monomer_logP =  os.path.join(path_to_cur, path_to_monomer_logP)
    path_to_monomer_cfg = os.path.join(path_to_cur, path_to_monomer_cfg)
    if args.monomer_cfg is not None:
        path_to_monomer_cfg = os.path.abspath(args.monomer_cfg)

    local_monomers_logP = os.path.join(main_out_dir, "monomersLogP.tsv")
    local_monomers_cfg = os.path.join(main_out_dir, "monomers.tsv")
    copyfile(path_to_monomer_cfg, local_monomers_cfg)
    copyfile(path_to_monomer_logP, local_monomers_logP)

    comand = path_to_exec_dir + "/NRPsMatcher \"" +  path_predictions + "\" \"" + path_to_graphs + "\" \"" + path_to_AA + "\" \"" + path_to_cfg + "\"\n"
    print(comand)
    os.system(comand)
    return

args = parse_args()
run(args)