#!/usr/bin/env python
import sys
from logger import log

def create_predictions_by_antiSAMSHout(path_to_antismashouts, outdir, predictor):
    if predictor != "NRPSPREDICTOR2":
        log.err("You can provide antiSMASH output only for NRPSPREDICTOR2!")
        sys.exit()