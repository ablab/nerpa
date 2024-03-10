#!/usr/bin/env python
import sys
import os
import shutil
import csv

import handle_PCP2
import handle_MT
import handle_E
import splitter
import handle_helper
from collections import defaultdict

def gen_prediction_dict(orf_modules_names, input_file_name, dirname):
    '''
    dirname: converted_antiSMASH_v5_outputs
    input_file_name: CONTIGNAME_nrpspredictor2_codes.txt
    orf_part: list of the form ['ctg1_trsI_A1', 'ctg1_trsI_A2']
    '''
    double_orf, double_aa = handle_PCP2.get_double_orfs_and_AA(dirname, orf_modules_names)  # the lists of orfs and modules which can be duplicated
    mt_aa = handle_MT.get_MT_AA(dirname, orf_modules_names)  # the list of modules with methylation
    d_aa = handle_E.get_D_AA(dirname, orf_modules_names)  # the list of modules with epimerization

    prediction_dict = defaultdict(lambda: '')
    with open(input_file_name, 'r') as rf:
        for line in rf:
            module_name, aa_prediction, scores = line.split('\t')
            ctgorf, domain_id = module_name.rsplit('_', 1)

            new_module_name = '_'.join([ctgorf + ('*' if ctgorf in double_orf else ''),  # indicate that orf can be duplicated
                                        domain_id + ('*' if module_name in double_aa else '')])  # indicate that modules can be duplicated

            new_aa_prediction = ''.join(['d-' if module_name in d_aa else '',  # label indicating chirality
                                         aa_prediction,
                                         '+MT' if module_name in mt_aa else ''])  # label indicating methylation

            prediction_dict[ctgorf] += '\t'.join([new_module_name, new_aa_prediction, scores])

    return prediction_dict


def gen_prediction_for_one_orfs_part(orf_part, input_file_name, output_prefix,
                                     current_part, predictions_info_list, base_antiSMASH_out_name):
    prediction_dict = gen_prediction_dict(orf_part, input_file_name, base_antiSMASH_out_name)

    output_str = ""
    for current_orf in orf_part:
        if current_orf in prediction_dict:
            output_str += prediction_dict[current_orf]

    if (output_str != ""):
        output_file = output_prefix + "_part" + str(current_part)
        with open(output_file, 'w') as wf:
            wf.write(output_str)
        current_part += 1
        predictions_info_list.append(output_file)

    return current_part

def gen_predictions(bgc_orfs_parts, input_file_name, output_prefix, current_part,
                    predictions_info_list, base_antiSMASH_out_name):
    for orf_part in bgc_orfs_parts:
        current_part = gen_prediction_for_one_orfs_part(orf_part, input_file_name, output_prefix, current_part, predictions_info_list, base_antiSMASH_out_name)
    return current_part

def create_predictions_by_antiSMASH_out(antiSMASH_outs, outdir, log):
    log.info("Start create predictions by antiSMASH")

    if not antiSMASH_outs:
        log.info("Error: no antiSMASH results found")
        raise ValueError("Could not find antiSMASH output")

    dir_for_predictions = os.path.join(outdir, "predictions")
    if not os.path.exists(dir_for_predictions):
        os.makedirs(dir_for_predictions)

    predictions_info_file = os.path.join(outdir, "predictions.info")
    predictions_info_list = []
    for dirname in antiSMASH_outs:
        try:
            if dirname[-1] == '\n':
                dirname = dirname[:-1]

            orf_pos = handle_helper.get_orf_position(dirname)
            orf_ori = handle_helper.get_orf_orientation(dirname)
            orf_domains = handle_helper.get_orf_domain_list(dirname)

            print("====PARTS BEFORE: ")
            parts = handle_helper.get_parts(dirname)
            handle_helper.debug_print_parts(dirname, parts, orf_domains, orf_ori, orf_pos)

            #print("====SPLIT BY DIST:")
            parts = splitter.split_by_dist(parts, orf_pos)
            #handle_helper.debug_print_parts(dirname, parts, orf_domains, orf_ori, orf_pos)

            #print("====SPLIT BY SINGLE ORF WITH Starter-TE")
            parts = splitter.split_by_one_orf_Starter_TE(parts, orf_ori, orf_domains)
            #handle_helper.debug_print_parts(dirname, parts, orf_domains, orf_ori, orf_pos)

            #print("====REMOVE SINGLE DOMAINs ORFS")
            parts = splitter.split_by_single_domain_orf(parts, orf_ori, orf_domains)
            #handle_helper.debug_print_parts(dirname, parts, orf_domains, orf_ori, orf_pos)

            print("====SPLIT AND REORDER")
            parts = splitter.split_and_reorder(parts, orf_ori, orf_pos, orf_domains)
            handle_helper.debug_print_parts(dirname, parts, orf_domains, orf_ori, orf_pos)

            if len(parts) > 100:
                raise RuntimeError(f'Too many parts: {len(parts)}')


            nrpspred_dir = os.path.join(dirname, "nrpspks_predictions_txt")
            if os.path.isdir(nrpspred_dir):
                for filename in os.listdir(nrpspred_dir):
                    if filename.endswith('nrpspredictor2_codes.txt'):
                        base_antiSMASHout_name = os.path.basename(dirname)
                        base_pred_name = os.path.basename(filename)
                        #predictions_info_list.append(os.path.join(dir_for_predictions, base_antiSMASHout_name + "_" + base_pred_name))
                        #shutil.copyfile(os.path.join(nrpspred_dir, filename), os.path.join(dir_for_predictions, base_antiSMASHout_name + "_" + base_pred_name))
                        gen_predictions(parts, os.path.join(nrpspred_dir, filename),
                                        os.path.join(dir_for_predictions, base_antiSMASHout_name + "_" + base_pred_name)[:-4],
                                        0, predictions_info_list, dirname)
        except KeyboardInterrupt as e:
            raise e
        except Exception as e:
            print(f'Error: {type(e).__name__}: {e}')
            print(f'Skipping {dirname}')

    f = open(predictions_info_file, 'w')
    for line in predictions_info_list:
        f.write(line + "\n")
    f.close()

    return predictions_info_file