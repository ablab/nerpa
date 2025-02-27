#!/usr/bin/env python

import argparse
import os
import sys
from pathlib import Path
import json

import config
import codes_handler
import json_handler
from log_utils import error, info


def main(args):
    parser = argparse.ArgumentParser(description='Parser of antiSMASH v5 output and interpreter of Stachelhaus codes')
    # TODO: add argument for classical output (as much NRPSpredictor2-compatibe as possible) and advanced (our improved style)
    parser.add_argument(
        'inputs',
        metavar='FILE/DIR',
        type=str,
        nargs='+',
        help='paths to antiSMASH v.5 JSONs or its output directories containing JSONs'
    )
    parser.add_argument(
        '-o', '--output-dir',
        metavar='DIR',
        type=str,
        default=None,
        help='Output dir ("nrpspks_predictions_txt" and "txt" subdirs will be created inside). '
             'If not specified: use parent dir of each input JSON file. '
             'If specified and multiple inputs are provided, the dir is used as a root dir where each input'
             ' processing results are stored separately (subdir name is a basename of input JSON parent dir)'
    )
    parser.add_argument(
        '-c', '--codes',
        metavar='FILE',
        type=str,
        default=None,
        help='Path to either a file with preprocessed codes (custom JSON format) '
             'or to a file with raw known codes (TXT or FASTA format). '
             'Default: %s (if exists) or %s (otherwise)' % (config.DEFAULT_OUTPUT_CODES, config.DEFAULT_INPUT_CODES)
    )
    parser.add_argument(
        '--preprocess',
        metavar='FILE/"default"',
        type=str,
        default=None,
        help='Preprocess raw known codes to a custom format and save in the specified file '
             '(or in the default location if "default" is specified: %s)' % config.DEFAULT_OUTPUT_CODES
    )
    parser.add_argument(
        '-n', '--naming-style',
        choices=['v3', 'v5', 'mix'],
        default='v3',
        help='Naming style for A domains and TXT files.'
             '"v3" results in names like "ctg34_orf00040_A1" and "ctg34_nrpspredictor2_codes.txt";'
             '"v5" results in names like "nrpspksdomains_ctg34_39_AMP-binding.1" and "JOGG01000034.1_nrpspredictor2_codes.txt"'
             '"mix" results in names like "nrpspksdomains_ctg34_39_AMP-binding.1" and "ctg34_nrpspredictor2_codes.txt". '
             'Default: "%(default)s"'
    )
    parser.add_argument(
        '-m', '--mode',
        choices=['stachelhaus', 'hybrid'],
        default='stachelhaus',
        help='Scoring mode (hybrid takes into account both SVM and Stachelhaus code), default: "%(default)s"'
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Produce verbose output. Otherwise only the most important messages are printed.'
    )

    args = parser.parse_args(args)

    # special case: determining --codes default value if nothing is specified
    if args.codes is None:
        if os.path.exists(config.DEFAULT_OUTPUT_CODES):
            args.codes = config.DEFAULT_OUTPUT_CODES
        else:  # elif os.path.exists(config.DEFAULT_INPUT_CODES):
            args.codes = config.DEFAULT_INPUT_CODES

    codes_ftype = codes_handler.get_codes_file_type(args.codes)
    if codes_ftype == 'known':
        info('Preprocessing known codes into a custom format', verbose=args.verbose)
        known_codes = codes_handler.parse_known_codes(args.codes)
        if args.preprocess is not None:
            preprocessed_fpath = args.preprocess if args.preprocess != 'default' else config.DEFAULT_OUTPUT_CODES
            info('Saving the preprocessed codes into ' + preprocessed_fpath, verbose=args.verbose)
            with open(preprocessed_fpath, 'w') as f:
                json.dump(known_codes, f)
    elif codes_ftype == 'preprocessed':
        info('Reusing already preprocessing known codes', verbose=args.verbose)
        with open(args.codes) as f:
            known_codes = json.load(f)
    else:
        error('File with codes (--code) does not exist or its format cannot be recognized. Aborting..', exit=True)

    is_root_outdir = True if (args.output_dir is not None and len(args.inputs) > 1) else False
    processed_output_dirs = []
    for input_path in args.inputs:
        try:
            processed_output_dirs.append(json_handler.handle_single_input(
                Path(input_path), args.output_dir, is_root_outdir, args.naming_style,
                known_codes, scoring_mode=args.mode, verbose=args.verbose))
        except KeyboardInterrupt as e:
            raise e
        except RuntimeError as e:
            info(f'ERROR: Unable to parse the input at "{input_path}": {e}')
        except Exception as e:
            info(f'ERROR: Unmanaged Exception while parsing the input at "{input_path}": {e}')    
    return processed_output_dirs


if __name__ == "__main__":
    main(sys.argv[1:])
