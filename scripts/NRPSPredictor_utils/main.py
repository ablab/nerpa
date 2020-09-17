#!/usr/bin/env python

import argparse
import os
import json

import config
import codes_handler
import json_handler
from log_utils import error, info


def main():
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
             'If not specified: use parent dir of each input JSON file'
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
        '-m', '--mode',
        choices=['classic', 'weighted', 'hybrid'],
        default='classic',
        help='Scoring mode, default: "%(default)s"'
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Produce verbose output. Otherwise only the most important messages are printed.'
    )

    args = parser.parse_args()

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

    for input_path in args.inputs:
        json_handler.handle_single_input(input_path, args.output_dir, known_codes, args.verbose)


if __name__ == "__main__":
    main()
