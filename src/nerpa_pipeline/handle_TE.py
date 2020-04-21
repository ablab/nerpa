#!/usr/bin/env python
import sys
import os
import shutil
from logger import log

def create_predictions_by_antiSAMSHout(path_to_antismashouts, outdir, predictor):
    if predictor != "NRPSPREDICTOR2":
        log.err("You can provide antiSMASH output only for NRPSPREDICTOR2!")
        sys.exit()

    dir_for_predictions = os.path.join(outdir, "predictions")
    if not os.path.exists(dir_for_predictions):
        os.makedirs(dir_for_predictions)

    predictions_info_file = os.path.join(outdir, "predictions.info")
    predictions_info_list = []
    with open(path_to_antismashouts) as fr:
        for dirname in fr:
            if dirname[-1] == '\n':
                dirname = dirname[:-1]
            nrpspred_dir = os.path.join(dirname, "nrpspks_predictions_txt")
            if os.path.isdir(nrpspred_dir):
                for filename in os.listdir(nrpspred_dir):
                    if filename.endswith('nrpspredictor2_codes.txt'):
                        base_antismashout_name = os.path.basename(dirname)
                        base_pred_name = os.path.basename(filename)
                        predictions_info_list.append(os.path.join(dir_for_predictions, base_antismashout_name + "_" + base_pred_name))
                        shutil.copyfile(os.path.join(nrpspred_dir, filename), os.path.join(dir_for_predictions, base_antismashout_name + "_" + base_pred_name))

    f = open(predictions_info_file, 'w')
    for line in predictions_info_list:
        f.write(line + "\n")
    f.close()

    return predictions_info_file

