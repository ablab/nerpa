#!/usr/bin/env python3

import sys
import os
import argparse
import csv
import site
import shutil
import json
import yaml

import nerpa_init
nerpa_init.init()


site.addsitedir(nerpa_init.python_modules_dir)

import predictions_preprocessor
import nerpa_utils
import handle_rban
import logger

from src.nerpa_pipeline.rban_names_helper import rBAN_Names_Helper
from pathlib import Path

# for detecting and processing antiSMASH v.5 output
site.addsitedir(os.path.join(nerpa_init.python_modules_dir, 'NRPSPredictor_utils'))
from NRPSPredictor_utils.json_handler import get_main_json_fpath
from NRPSPredictor_utils.main import main as convert_antiSMASH_v5
from src.NewMatcher.scoring_helper import ScoringHelper
from src.NewMatcher.scoring_config import load_config as load_scoring_config
from src.NewMatcher.matcher import get_matches

from src.data_types import (
    BGC_Variant,
    NRP_Variant,
    UNKNOWN_RESIDUE
)

from src.write_results import write_results, write_nrp_variants, write_bgc_variants


def parse_args(log):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    genomic_group = parser.add_argument_group('Genomic input', 'Genomes of NRP-producing organisms (i.e. BGC predictions)')
    genomic_group.add_argument("--antismash_output_list", dest="antismash_out",
                               help="file with list of paths to antiSMASH output directories", type=str)
    genomic_group.add_argument("--antismash", "-a", dest="antismash", action='append',
                               help="single antiSMASH output directory or directory with many antiSMASH outputs")
    genomic_group.add_argument("--sequences", dest="seqs", help="GenBank/EMBL/FASTA file containing DNA sequences")

    struct_group = parser.add_argument_group('Chemical input', 'Structures of NRP molecules')
    struct_input_group = struct_group.add_mutually_exclusive_group()
    struct_input_group.add_argument("--rban-json", dest="rban_output",
                                    help="json file with rBAN-preprocessed NRP structures", type=str)
    struct_input_group.add_argument("--smiles", dest="smiles", nargs='*',
                        help="string (or several strings) with structures in the SMILES format", type=str)
    struct_input_group.add_argument("--smiles-tsv", dest="smiles_tsv",
                        help="multi-column file containing structures in the SMILES format", type=str)
    struct_group.add_argument("--col-smiles", dest="col_smiles",
                        help="column name in smiles-tsv for structures in the SMILES format [default: 'SMILES']",
                        type=str, default='SMILES')
    struct_group.add_argument("--col-id", dest="col_id",
                        help="column name in smiles-tsv for structure identifier (if not provided, row index will be used)",
                        type=str)
    struct_group.add_argument("--sep", dest="sep",
                        help="column separator in smiles-tsv", type=str, default='\t')

    advanced_input_group = parser.add_argument_group('Advanced input', 'Preprocessed BGC predictions and NRP structures in custom Nerpa-compliant formats')
    advanced_input_group.add_argument("--predictions", "-p", nargs=1, dest="predictions",
                                          help="file with paths to preprocessed BGC prediction files", type=str)
    advanced_input_group.add_argument("--structures", "-s", dest="structures",
                                          help="file with Nerpa-preprocessed NRP structures", type=str)
    advanced_input_group.add_argument("--configs_dir", help="custom directory with adjusted Nerpa configs", action="store", type=str)
    advanced_input_group.add_argument("--force-existing-outdir", dest="output_dir_reuse", action="store_true", default=False,
                                      help="don't crash if the output dir already exists")

    # parser.add_argument("--insertion", help="insertion score [default=-2.8]", default=-2.8, action="store")
    # parser.add_argument("--deletion", help="deletion score [default=-5]", default=-5, action="store")
    parser.add_argument('--rban-monomers-db', dest='rban_monomers', type=str, default=None,
                        help='file with custom monomers in rBAN compatible format')
    parser.add_argument("--process-hybrids", dest="process_hybrids", action="store_true", default=False,
                        help="process NRP-PK hybrid monomers (requires use of rBAN)")
    parser.add_argument('--antismash-path', dest='antismash_path', type=str, default=None,
                        help='path to antismash source directory')
    parser.add_argument("--threads", default=1, type=int, help="number of threads for running Nerpa", action="store")
    parser.add_argument("--min-score", default=0, type=float, help="minimum score to report a match", action="store")
    parser.add_argument("--num-matches", default=None, type=int, help="maximum number of matches to report", action="store")
    parser.add_argument("--debug", action="store_true", default=False,
                        help="run in the debug mode (keep intermediate files)")
    parser.add_argument("--output_dir", "-o", help="output dir [default: nerpa_results/results_<datetime>]",
                        type=str, default=None)

    parsed_args = parser.parse_args()

    validate_arguments(parsed_args, parser, log)
    return parsed_args


def check_tsv_ids_duplicates(reader, col_id):
    from collections import defaultdict
    counts = defaultdict(int)
    for row in reader:
        counts[row[col_id]] += 1
    duplicates = [(k, v) for k,v in counts.items() if v > 1]
    return duplicates


class ValidationError(Exception):
    pass


def validate(expr, msg=''):
    if not expr:
        raise ValidationError(msg)


def validate_arguments(args, parser, log):
    try:
        if not (args.predictions or args.antismash or args.antismash_out or args.seqs):
            raise ValidationError(f'one of the arguments --predictions --antismash/-a --antismash_output_list '
                                  f'--sequences is required')
        if args.predictions and (args.antismash or args.antismash_out or args.seqs):
            raise ValidationError(f'argument --predictions: not allowed with argument --antismash/-a '
                                  f'or --antismash_output_list or --sequences')
        if not (args.structures or args.smiles or args.smiles_tsv or args.rban_output):
            raise ValidationError(f'one of the arguments --rban-json --smiles-tsv --smiles --structures/-s'
                                  f'is required')
        if args.structures and (args.smiles or args.smiles_tsv or args.rban_output):
            raise ValidationError('argument --structures/-s: not allowed with argument --rban-json or --smiles '
                                  'or --smiles-tsv')
        if args.smiles_tsv:
            try:
                with open(args.smiles_tsv, newline='') as f_in:
                    reader = csv.DictReader(f_in, delimiter=args.sep, quoting=csv.QUOTE_NONE)
                    validate(args.col_smiles in reader.fieldnames,
                        f'Column "{args.col_smiles}" was specified but does not exist in {args.smiles_tsv}.')
                    if args.col_id:
                        validate(args.col_id in reader.fieldnames,
                            f'Column "{args.col_id}" was specified but does not exist in {args.smiles_tsv}.')
                        duplicates = check_tsv_ids_duplicates(reader, args.col_id)
                        validate(len(duplicates) == 0, f'Duplicate IDs are found: {duplicates}')
            except FileNotFoundError:
                raise ValidationError(f'No such file: "{args.smiles_tsv}".')
            except csv.Error as e:
                raise ValidationError(f'Cannot parse "{args.smiles_tsv}": {e}.')
            except Exception as e:
                raise ValidationError(f'Invalid input file "{args.smiles_tsv}": {e}.')

    except ValidationError as e:
        parser.print_help()
        error_msg = f'{e}\n' if str(e) else 'Options validation failed!'
        log.error(error_msg, to_stderr=True)


def validate_rban_output(path_to_rban_input, path_to_rban_output):
    """
    Checks whether the output is produced. Needed since rBAN silently crashes on certain inputs.

    :param path_to_rban_input: rban input json
    :param path_to_rban_output: rban output json
    :return:
    """
    with open(path_to_rban_input) as f:
        ids_in = set(x['id'] for x in json.load(f))
    with open(path_to_rban_output) as f:
        ids_out = set(x['id'] for x in json.load(f))
    for idx in ids_in - ids_out:
        log.warning(f'No rBAN output for structure "{idx}"')


def run_rban_on_smiles(args, main_out_dir, path_to_rban_jar, path_to_monomers_db, log):
    path_to_rban_input = os.path.join(main_out_dir, 'rban.input.json')
    if args.smiles_tsv:
        handle_rban.generate_rban_input_from_smiles_tsv(
            args.smiles_tsv, path_to_rban_input, sep=args.sep, id_col_name=args.col_id, smi_col_name=args.col_smiles)
    else:
        handle_rban.generate_rban_input_from_smiles_strings(args.smiles, path_to_rban_input)

    path_to_rban_output = os.path.join(main_out_dir, 'rban.output.json')
    log.info('\n======= Structures preprocessing with rBAN')
    handle_rban.run_rban(path_to_rban_jar, path_to_rban_input, path_to_rban_output, path_to_monomers_db, main_out_dir, log)
    validate_rban_output(path_to_rban_input, path_to_rban_output)
    log.info("\n======= Done with Structures preprocessing with rBAN")
    return path_to_rban_output


def create_merged_monomers_db(path_to_rban, path_to_nerpa_monomers, path_to_user_monomers, output_dir):
    from zipfile import ZipFile
    import json
    with ZipFile(path_to_rban) as zf:
        default_db = json.loads(zf.read('molecules/monomer/nrproMonomers.json'))

    def _append_db(path):
        if not path:
            return []
        start_id = 1 + max(m['id'] for m in default_db)
        with open(path_to_nerpa_monomers) as f:
            custom_db = json.loads(f.read())
        return [{**m, 'id':i} for i, m in enumerate(custom_db, start=start_id)]

    default_db += _append_db(path_to_nerpa_monomers)
    default_db += _append_db(path_to_user_monomers)

    path_to_merged_db = os.path.join(output_dir, 'rban_monomers_db.json')
    with open(path_to_merged_db, 'w') as f:
        f.write(json.dumps(default_db))
    return path_to_merged_db


def copy_prediction_list(args, main_out_dir):
    new_prediction_path = os.path.join(main_out_dir, "prediction.info")
    with open(new_prediction_path, 'w') as f:
        with open(args.predictions[0]) as fr:
            for line in fr:
                line_parts = line.split()
                file = line_parts[0]
                if (file[0] != '/'):
                    file = os.path.join(os.path.dirname(os.path.abspath(args.predictions[0])), file)
                f.write(file + "\n")
    return new_prediction_path


def get_antismash_v3_compatible_input_paths(listing_fpath, list_of_paths, output_dir, log):
    '''
    Parses all antiSMASH-related options,
    detects all relevant output dirs (either with .json [aS v.5] or with ./txt/ & ./nrpspks_predictions_txt [aS v.3],
    converts aS v.5 to aS v.3-compliants if needed,
    returns list of paths to each v.3-compliant directory
    :param args:
    :return:
    '''

    def _get_input_antiSMASH_paths(lookup_paths):
        def _is_antiSMASHv3_path(path):
            if os.path.isdir(path) and \
               os.path.isdir(os.path.join(path, 'txt')) and \
               os.path.isdir(os.path.join(path, 'nrpspks_predictions_txt')):
                return True
            return False

        def _is_antiSMASHv5_path(path):
            if os.path.isfile(path) and path.endswith('.json'):
                return True
            if os.path.isdir(path) and get_main_json_fpath(dirpath=path) is not None:
                return True
            return False

        antiSMASHv3_paths = []
        antiSMASHv5_paths = []
        for entry in lookup_paths:
            if _is_antiSMASHv3_path(entry):
                antiSMASHv3_paths.append(entry)
            elif _is_antiSMASHv5_path(entry):
                antiSMASHv5_paths.append(entry)
            elif os.path.isdir(entry):
                # excluding dirs in runtime in os.walk: https://stackoverflow.com/questions/19859840/excluding-directories-in-os-walk
                for root, dirs, files in os.walk(entry, topdown=True):
                    # always ignore files since a single json in a dir should be caught one step before when path was the dir
                    # (see _is_antiSMASHv5_path() ) while multiple jsons in a dir probably means a false positive
                    dirs_to_keep_walking = []
                    for dir in dirs:
                        full_dir_path = os.path.join(root, dir)
                        if _is_antiSMASHv3_path(full_dir_path):
                            antiSMASHv3_paths.append(full_dir_path)
                        elif _is_antiSMASHv5_path(full_dir_path):
                            antiSMASHv5_paths.append(full_dir_path)
                        else:
                            dirs_to_keep_walking.append(dir)
                    dirs[:] = dirs_to_keep_walking
        return antiSMASHv3_paths, antiSMASHv5_paths

    lookup_locations = []
    if listing_fpath is not None:
        with open(listing_fpath) as f:
            for path in f:
                lookup_locations.append(path.strip())

    if list_of_paths:
        lookup_locations += list_of_paths

    antiSMASHv3_paths, antiSMASHv5_paths = _get_input_antiSMASH_paths(lookup_locations)
    log.info("\n=== Genome predictions found: %d antiSMASH v3 inputs; %d antiSMASH v5 inputs" %
             (len(antiSMASHv3_paths), len(antiSMASHv5_paths)))
    if antiSMASHv5_paths:
        log.info("\n======= Preprocessing antiSMASH v5 inputs")
        converted_antiSMASH_v5_outputs_dir = os.path.join(output_dir, "converted_antiSMASH_v5_outputs")
        log.info(f'results will be in {converted_antiSMASH_v5_outputs_dir}', indent=1)
        converted_antiSMASH_v5_paths = convert_antiSMASH_v5(antiSMASHv5_paths +
                                                            ['-o', converted_antiSMASH_v5_outputs_dir, '-m', 'hybrid', "-n", "v3" ])
        antiSMASHv3_paths += converted_antiSMASH_v5_paths
        log.info("\n======= Done with Preprocessing antiSMASH v5 inputs")

    return antiSMASHv3_paths


def run(args, log):
    output_dir = nerpa_utils.set_up_output_dir(output_dirpath=args.output_dir,
                                               crash_if_exists=not args.output_dir_reuse, log=log)
    log.set_up_file_handler(output_dir)
    log.start()

    if args.predictions is not None:
            bgc_variants = []
            for path_to_predictions in args.predictions:
                for file_with_bgc_variants in filter(lambda f: f.suffix in ('.yml', '.yaml'),
                                                     Path(path_to_predictions).iterdir()):
                    bgc_variants.extend(BGC_Variant.from_yaml_dict(yaml_record)
                                        for yaml_record in yaml.safe_load(file_with_bgc_variants.read_text()))
    else:
        antismash_out_dirs = args.antismash if args.antismash is not None else []
        if args.seqs:
            cur_antismash_out = os.path.join(output_dir, 'antismash_output')
            if args.antismash_path:
                antismash_exe = nerpa_utils.get_path_to_program('run_antismash.py', dirpath=args.antismash_path, min_version='5.0')
            else:
                antismash_exe = nerpa_utils.get_path_to_program('antismash', min_version='5.0')
            if antismash_exe is None:
                log.error("Can't find antismash 5.x executable. Please make sure that you have antismash 5.x installed "
                          "in your system or provide path to antismash source directory via --antismash-path option.")
            command = [antismash_exe,
                       '--genefinding-tool', 'prodigal',
                       '--output-dir', cur_antismash_out,
                       '--minimal', '--skip-zip-file', '--enable-nrps-pks',
                       '--cpus', str(args.threads), args.seqs]
            nerpa_utils.sys_call(command, log, cwd=output_dir)
            antismash_out_dirs.append(cur_antismash_out)

        bgc_variants = predictions_preprocessor.parse_antismash_output(get_antismash_v3_compatible_input_paths(
                listing_fpath=args.antismash_out, list_of_paths=antismash_out_dirs,
                output_dir=output_dir, log=log), output_dir, args.debug, log)
        if not args.debug: # in the debug mode all the variants are written during generation
            write_bgc_variants(bgc_variants, output_dir)

    input_configs_dir = args.configs_dir if args.configs_dir else nerpa_init.configs_dir
    current_configs_dir = os.path.join(output_dir, "configs")
    # remember shutil.copytree caveat (compared to dir_util.copy_tree):
    # directory metadata will be copied that may cause potential problems
    # if src dir is too old and there is a cluster cronjob which
    # automatically remove old files from the temporary workspace
    shutil.copytree(input_configs_dir, current_configs_dir, copy_function=shutil.copy)

    path_to_graphs = os.path.join(output_dir, 'structures.info')
    local_monomers_cfg = os.path.join(current_configs_dir, "monomers.tsv")
    path_to_rban = os.path.join(nerpa_init.external_tools_dir, 'rBAN', 'rBAN-1.0.jar')
    path_to_monomers_db = create_merged_monomers_db(
        path_to_rban, os.path.join(current_configs_dir, "monomersdb_nerpa.json"), args.rban_monomers, output_dir)

    if args.structures is not None:
        nrp_variants = []
        for file_with_nrp_variants in filter(lambda f: f.suffix in ('.yml', '.yaml'),
                                             Path(args.structures).iterdir()):
            nrp_variants.extend(NRP_Variant.from_yaml_dict(yaml_record)
                                for yaml_record in yaml.safe_load(file_with_nrp_variants.read_text()))
        graph_records = []  # TODO: should we load them as well?
    else:
        if args.rban_output:
            path_rban_output = args.rban_output
        else:
            path_rban_output = run_rban_on_smiles(args, output_dir, path_to_rban, path_to_monomers_db, log)

        rban_names_helper = rBAN_Names_Helper(Path(local_monomers_cfg))
        graph_records, nrp_variants = handle_rban.parse_rban_output(
            path_rban_output, local_monomers_cfg, path_to_graphs, output_dir, path_to_rban, path_to_monomers_db, log,
            names_helper=rban_names_helper,
            process_hybrids=args.process_hybrids
        )
        write_nrp_variants(nrp_variants, output_dir, graph_records)

    scoring_config = load_scoring_config(Path(__file__).parent / Path('src/NewMatcher/scoring_config.yaml'))  # TODO: this is ugly
    scoring_helper = ScoringHelper(scoring_config)

    log.info("\n======= Nerpa matching")
    matches = get_matches(bgc_variants, nrp_variants, scoring_helper,
                          min_score=args.min_score,
                          max_num_matches=args.num_matches,
                          num_threads=args.threads,
                          log=log)

    write_results(matches, output_dir)
    log.info("RESULTS:")
    log.info("Main report is saved to " + str(output_dir / Path('report.tsv')), indent=1)
    log.info("Detailed reports are saved to " + str(output_dir / Path('matches_details')), indent=1)
    log.finish()


if __name__ == "__main__":
    log = logger.NerpaLogger()
    try:
        args = parse_args(log)
        run(args, log)
    except Exception:
        _, exc_value, _ = sys.exc_info()
        log.exception(exc_value)
    finally:
        # TODO: clean up: remove all intermediate files
        pass
