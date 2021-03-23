#!/usr/bin/env python3

import sys
import os
import argparse
import csv
import site
from shutil import copyfile
import subprocess

import nerpa_init
nerpa_init.init()

site.addsitedir(nerpa_init.python_modules_dir)

import nerpa_utils
import nerpa_config
import handle_TE
import handle_rban
import logger
# from logger import log  # TODO FIXME

path_to_exec_dir = os.path.dirname(os.path.abspath(__file__)) + "/"


def parse_args():
    global parser
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    genomic_group = parser.add_argument_group('Genomic input', 'description')
    genomic_group.add_argument("--predictions", "-p", nargs=1, dest="predictions",
                                     help="path to file with paths to prediction files", type=str)
    genomic_group.add_argument("--antismash_output_list", dest="antismash_out",
                               help="path to file with list of paths to antiSMASH output folders", type=str)
    genomic_group.add_argument("--antismash", "-a", dest="antismash", help="path to antiSMASH output or to directory with many antiSMASH outputs", type=str)

    struct_group = parser.add_argument_group('Chemical input', 'description')
    struct_input_group = struct_group.add_mutually_exclusive_group(required=True)
    struct_input_group.add_argument("--smiles", dest="smiles", nargs='*',
                        help="string (or several strings) with smiles", type=str)
    struct_input_group.add_argument("--smiles-tsv", dest="smiles_tsv",
                        help="csv with named columns", type=str)
    struct_input_group.add_argument("--graphs", dest="graphs",
                                    help="", type=str)
    struct_group.add_argument("--col_smiles", dest="col_smiles",
                        help="name of column with smiles [default='SMILES']", type=str, default='SMILES')
    struct_group.add_argument("--col_id", dest="col_id",
                        help="name of col with ids (if not provided, id=[number of row])", type=str)
    struct_group.add_argument("--sep", dest="sep",
                        help="separator in smiles-tsv", type=str, default='\t')

    # parser.add_argument("--insertion", help="insertion score [default=-2.8]", default=-2.8, action="store")
    # parser.add_argument("--deletion", help="deletion score [default=-5]", default=-5, action="store")
    parser.add_argument("--scoring_cfg", help="path to file with scoring weights", action="store", type=str)
    parser.add_argument("--modification_cfg", help="path to file with modification description", action="store", type=str)
    parser.add_argument("--monomer_cfg", help="path to file with monomer description", action="store", type=str)
    parser.add_argument("--threads", default=1, type=int, help="number of threads for running Nerpa", action="store")
    parser.add_argument("--output_dir", "-o", help="output dir [default: nerpa_results/results_<datetime>]",
                        type=str, default=None)
    args = parser.parse_args()
    return args


def check_tsv_ids_duplicates(reader, col_id):
    ids = [row[col_id] for row in reader]
    return len(ids) == len(set(ids))

class ValidationError(Exception):
    pass

def validate(expr, msg=''):
    if not expr:
        raise ValidationError(msg)

def validate_arguments(args):
    try:
        if not (args.predictions or args.antismash or args.antismash_out):
            raise ValidationError(f'one of the arguments --predictions --antismash/-a --antismash_output_list is required')
        if args.predictions and (args.antismash or args.antismash_out):
            raise ValidationError(f'argument --predictions: not allowed with argument --antismash/-a or --antismash_output_list')

        if args.smiles_tsv:
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
            logger.error(f'{e}\n')
        parser.print_help()
        sys.exit()


def print_cfg(args, output_dir):
    cfg_file = os.path.join(output_dir, "nerpa.cfg")
    with open(cfg_file, "w") as f:
        f.write(os.path.abspath(os.path.join(output_dir, "modifications.tsv")) + "\n")
        f.write(os.path.abspath(os.path.join(output_dir, "monomers.tsv")) + "\n")
        f.write(os.path.abspath(os.path.join(output_dir, "monomersLogP.tsv")) + "\n")
        f.write(os.path.abspath(os.path.join(output_dir, "prob_gen.cfg")) + "\n")
        f.write("threads " + str(args.threads) + "\n")

        # TODO: add to cmd parameters
        f.write(f'min_score 0.05\n')
        f.write(f'min_explain_part 0\n')
        f.write(f'default_monomer_logp -6.2\n')
        f.write(f'default_aminoacid_logp -6.64\n')
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


def gen_graphs_from_smiles_tsv(args, main_out_dir, path_to_monomers_tsv, path_to_graphs, path_to_rban_jar, path_to_monomers):
    path_to_rban_input = os.path.join(main_out_dir, 'rban.input.json')
    if args.smiles_tsv:
        handle_rban.generate_rban_input_from_smiles_tsv(
            args.smiles_tsv, path_to_rban_input, sep=args.sep, id_col_name=args.col_id, smi_col_name=args.col_smiles)
    else:
        handle_rban.generate_rban_input_from_smiles_string(args.smiles, path_to_rban_input)

    path_to_rban_output = os.path.join(main_out_dir, 'rban.output.json')
    logger.info('Running rBAN...')
    command = ['java', '-jar', path_to_rban_jar,
               '-monomersDB', path_to_monomers,
               '-inputFile', path_to_rban_input,
               '-outputFolder', main_out_dir + '/',  # rBAN specifics
               '-outputFileName', os.path.basename(path_to_rban_output)]
    try:
        p = subprocess.run(command,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE,  # instead of Python 3.7+ only: capture_output=True
                           universal_newlines=True,  # instead of Python 3.7+ only: text=True,
                           # check=True
                           )
        if p.stderr:
            logger.error(p.stderr)
        for line in p.stdout.split('\n'):
            if line.startswith('WARNING'):
                logger.warning('rBAN ' + line)
    except subprocess.CalledProcessError as e:
        logger.error(str(e))
        sys.exit()
    logger.info('Done.')

    handle_rban.generate_graphs_from_rban_output(path_to_rban_output, path_to_monomers_tsv, path_to_graphs, main_out_dir, path_to_rban_jar)


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


def get_antismash_list(args):
    aouts = []
    if (args.antismash_out is not None):
        with open(args.antismash_out) as f:
            for dir in f:
                aouts.append(dir)

    is_one = False
    if (args.antismash is not None):
        for filename in os.listdir(args.antismash):
            if filename.endswith(".json"):
                is_one = True
            if filename == "nrpspks_predictions_txt":
                is_one = True
        if is_one:
            aouts.append(args.antismash)
        else:
            for filename in os.listdir(args.antismash):
                aouts.append(os.path.join(args.antismash, filename))

    return aouts

def run(args):
    logger.init(nerpa_config.LOGGER_NAME)

    validate_arguments(args)

    output_dir = nerpa_utils.set_up_output_dir(output_dirpath=args.output_dir)
    log_fpath = os.path.join(output_dir, nerpa_config.LOGGER_NAME + '.log')
    logger.add_file_handler(log_fpath)

    if args.predictions is not None:
        path_predictions = copy_prediction_list(args, output_dir)
    else:
        path_predictions = handle_TE.create_predictions_by_antiSAMSHout(get_antismash_list(args), output_dir)

    path_to_cfg = print_cfg(args, output_dir)

    path_to_AA = os.path.join(nerpa_init.configs_dir, "aminoacids.tsv")

    path_to_modification_cfg = os.path.join(nerpa_init.configs_dir, "modifications.tsv")
    if args.modification_cfg is not None:
        path_to_modification_cfg = os.path.abspath(args.modification_cfg)

    local_modifications_cfg = os.path.join(output_dir, "modifications.tsv")
    copyfile(path_to_modification_cfg, local_modifications_cfg)

    if args.scoring_cfg is not None:
        path_to_scoring_cfg = os.path.abspath(args.scoring_cfg)
    else:
        path_to_scoring_cfg = os.path.join(nerpa_init.configs_dir, "prob_gen.cfg")
    local_scoring_cfg = os.path.join(output_dir, "prob_gen.cfg")
    copyfile(path_to_scoring_cfg, local_scoring_cfg)

    path_to_monomer_logP = os.path.join(nerpa_init.configs_dir, "monomersLogP.tsv")
    path_to_monomer_cfg = os.path.join(nerpa_init.configs_dir, "monomers.tsv")
    if args.monomer_cfg is not None:
        path_to_monomer_cfg = os.path.abspath(args.monomer_cfg)

    local_monomers_logP = os.path.join(output_dir, "monomersLogP.tsv")
    local_monomers_cfg = os.path.join(output_dir, "monomers.tsv")
    copyfile(path_to_monomer_cfg, local_monomers_cfg)
    copyfile(path_to_monomer_logP, local_monomers_logP)

    path_to_graphs = os.path.join(output_dir, 'path_to_graphs')
    if args.graphs:
        copyfile(args.graphs, path_to_graphs)
    else:
        path_to_rban = os.path.join(nerpa_init.external_tools_dir, 'rBAN', 'rBAN-1.0.jar')
        path_to_monomers = os.path.join(nerpa_init.external_tools_dir, 'rBAN', 'nrproMonomers_nerpa.json')
        gen_graphs_from_smiles_tsv(args, output_dir, local_monomers_cfg, path_to_graphs, path_to_rban, path_to_monomers)

    # FIXME: let NRPsMatcher work normally without changing dir!
    cwd = os.getcwd()
    os.chdir(output_dir)

    if not os.path.exists(os.path.dirname('details_mols/')):
        os.makedirs(os.path.dirname('details_mols/'))

    comand = os.path.join(nerpa_init.bin_dir, "NRPsMatcher") + " \"" + path_predictions + "\" \"" + path_to_graphs + "\" \"" + path_to_AA + "\" \"" + path_to_cfg + "\"\n"
    logger.info(comand)
    os.system(comand)
    logger.info("Nerpa results and log are saved to " + output_dir)

    os.chdir(cwd)


if __name__ == "__main__":
    try:
        args = parse_args()
        run(args)
    except Exception:
        _, exc_value, _ = sys.exc_info()
        logger.exception(exc_value)
        logger.error('Exception caught!')
    finally:
        # TODO: clean up: remove all intermediate files
        logger.cleanup()
