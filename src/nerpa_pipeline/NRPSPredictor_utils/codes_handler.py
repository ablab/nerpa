import os
import math
import pandas as pd
from collections import Counter, OrderedDict, defaultdict
from log_utils import error

from typing import Dict, List, Tuple
from src.data_types import (
    MonomerResidue,
    ResidueScores,
)
import config

ResidueSignatures = List[str]
ResidueSignaturesDict = Dict[MonomerResidue, Tuple[ResidueSignatures, ResidueSignatures]]


def get_codes_file_type(fpath):
    '''
    :param fpath:
    :return: 'known_fasta', 'known_txt', 'preprocessed_tsv' or None
    '''
    if not os.path.isfile(fpath):
        return None

    with open(fpath) as f:
        first_line = f.readline()
    if len(first_line.split()) == 1:
        return 'known'
    if first_line.startswith('[{'):
        return 'preprocessed'
    return None


def parse_known_codes(fpath):
    raw_codes = __parse_raw_codes(fpath)
    all_parsed_codes = []
    for (group_label, code), count in Counter(raw_codes).items():
        group_label = group_label.replace('/', '|')  # alternatives separator could be both '/' and '|'
        if group_label == 'ser-thr':  # super special case in nrpspredictor2_labeled_sigs.txt
            group_label = 'ser|thr'
        for label in group_label.split('|'):
            aa, mods, confs = __parse_aa_and_mod(label)
            if aa:
                all_parsed_codes.append({'code': code, 'label': label,
                                         'aa': aa, 'mods': mods, 'confs': confs,
                                         'count': count})
    return all_parsed_codes


def __parse_raw_codes(fpath):
    raw_codes = []
    with open(fpath) as f:
        prev_line = None
        for line in f:
            if prev_line is not None:
                label = prev_line.strip()
                if label.startswith('>'):
                    label = label[1:].split('_')[-1]
                code = line.strip()
                raw_codes.append((label, code))
                prev_line = None
            else:
                prev_line = line
    return raw_codes


def __parse_aa_and_mod(label):
    '''
    :param label: raw signature
    :return: parsed: AA, list of MODs, list of CONFs
    '''
    aa = None
    mods = []
    confs = []
    orig_label = label
    # label = label.lower() ## think about it

    ### parsing modifications ###

    # we suppose that D/L-configurations are determined by E or C/E-domains, not by Stachelhaus code
    # still, lets store them and allow the Nerpa main module to decide whether to use this info or not
    if 'd-' in label:
        label = label.replace('d-', '')
        confs.append('D')
    if label.endswith('-d'):
        label = label[:-2]
        confs.append('D')

    # beta-ala/lys things
    if label.endswith('-b'):
        label = 'b-' + label[:-2]

    if label in config.SPECIAL_CASES.keys():
        mods += config.SPECIAL_CASES[label][1]
        label = config.SPECIAL_CASES[label][0]
    if label in config.SYNONYMS_AA_SIGNATURES.keys():
        label = config.SYNONYMS_AA_SIGNATURES[label]

    if 'allo-' in label:
        label = label.replace('allo-', '')
        confs.append('allo')

    if label in config.KNOWN_AA_SIGNATURES:
        aa = label
    else:
        for known_aa in config.KNOWN_AA_SIGNATURES:
            if label.endswith(known_aa):
                aa = known_aa
                label = label[:-len(known_aa)]
                break
        if aa is not None:
            if label:
                for mod in config.MODS.keys():
                    if mod in label:
                        mods.append(mod)
                        label = label.replace(mod, '')
                if label:
                    error('failed to parse modifications: %s (full label: %s)' % (label, orig_label), exit=False)
        else:
            error('failed to parse core AA (%s) (full label: %s)!' %(label, orig_label), exit=False)

    return aa, mods, confs


def __get_aa_fullname(code_metadata, mode='classic'):
    '''
    code_metadata: {'code': code, 'label': label,
                    'aa': aa, 'mods': mods, 'confs': confs,
                    'count': count})
    '''
    if mode == 'classic':
        return code_metadata['label']
    elif mode == 'detailed':
        return '@'.join(code_metadata['confs']) + code_metadata["aa"] + '+'.join(code_metadata["mods"])
    elif mode == 'raw':
        return code_metadata["aa"]
    else:
        error('Internal bug: not defined mode (%s) for aa fullname pretty print!' % mode, exit=True)


def __get_aa_score(code: str, known_codes: ResidueSignatures) -> float:
    if len(known_codes):
        assert len(code) == len(next(iter(known_codes)))
    max_matched_letters = 0
    for c in known_codes:
        max_matched_letters = max(max_matched_letters, sum(code[i] == c[i] for i in range(len(c))))
    return float(max_matched_letters) / len(code)


def __get_svm_score(residue: str, svm_prediction: Tuple[float, List[MonomerResidue]]) -> float:
    if residue in config.SVM_SUBSTRATES:
        if residue in svm_prediction[1]:
            return svm_prediction[0]
        return 0.0
    return -1.0


def dummy_model(scoring_table: pd.DataFrame, model_params=None) -> ResidueScores:
    result: ResidueScores = OrderedDict()

    # Extracting and processing data
    substrates = scoring_table.index
    aa10_scores = scoring_table["aa10_score"].apply(lambda x: float(x))

    # Computing logarithm and populating results
    for substrate, aa10_score in zip(substrates, aa10_scores):
        if aa10_score is not None and aa10_score > 0:
            log_aa10_score = round(math.log(aa10_score), 4)
            result[substrate] = log_aa10_score

    return result


def get_prediction_from_signature(nrpys_prediction: dict, known_codes_dict: ResidueSignaturesDict,
                                  model=dummy_model, model_params=None) -> Tuple[str, str]:
    '''

    :param nrpys_prediction: includes aa10/34 and four level SVM predictions
    :param known_codes_dict: known aa10/34 signatures for all Nerpa-supported monomers (sorted!)
    :param model: stub for Azat's classifier
    :return: in the future: ResidueScores (with the default dummy model it is just log from aa10 score)
             right now (for legacy reasons): the most probably residue name and
                                             the sorted list of all residues with their log probs in '()', separated by ';'
    '''

    aa10_code = nrpys_prediction["aa10"]
    aa34_code = nrpys_prediction["aa34"]

    # parsing SVM prediction
    Level = str
    SvmScore = float
    svm_prediction: Dict[Level, Tuple[SvmScore, List[MonomerResidue]]] = dict()
    for level in config.SVM_LEVELS:
        score = nrpys_prediction[level]['score']
        substrate_short_names = [substrate['short'] for substrate in nrpys_prediction[level]['substrates']]
        substrate_nerpa_names = [config.antismash_substrate_to_nerpa_substrate(substrate) for substrate in substrate_short_names]
        svm_prediction[level] = (score, substrate_nerpa_names)

    scoring_table = pd.DataFrame([], columns=config.SCORING_TABLE_COLUMNS).set_index(config.SCORING_TABLE_INDEX)
    for residue, (known_aa10_codes, known_aa34_codes) in known_codes_dict.items():
        scoring_table_row = [__get_aa_score(aa10_code, known_aa10_codes),
                             __get_aa_score(aa34_code, known_aa34_codes)]

        for level in config.SVM_LEVELS:
            scoring_table_row.append(__get_svm_score(residue, svm_prediction[level]))

        scoring_table.loc[residue] = scoring_table_row

    residue_log_probs: ResidueScores = model(scoring_table, model_params)

    residue_with_highest_prob = max(residue_log_probs, key=residue_log_probs.get)  # for legacy reasons and just for debug, never used in practice
    residue_log_probs_string = ";".join(f"{key}({value})" for key, value in residue_log_probs.items())

    return residue_with_highest_prob, residue_log_probs_string