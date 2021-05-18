import os
import itertools
from collections import Counter, OrderedDict
from log_utils import error

import config


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


def __get_stachelhaus_score(signature, code, mode='classic'):
    if mode == 'classic':
        score = 0  # always making +1 for the last 'K' in 10aa code
        num_chars = min(len(signature), len(code))
        for i in range(num_chars):
            score += int(signature[i] == code[i])
        # special case: sometimes code (or signature) has 9 chars instead of the default length 10.
        # in this case, we assume the last char is 'K' and it always matches.
        if num_chars == config.STANDARD_STACHELHAUS_CODE_LENGTH - 1:
            score += 1
        score *= 10
    else:
        error('Internal bug: score mode %s is not supported yet!' % mode, exit=True)
        score = None
    return float(score)


def __get_svm_score(aa, svm):
    # TODO -- check the 'aa' is always in proper format, e.g. just 'cys' not '@Lcys+me+h'
    if aa not in config.SVM_DETECTABLE_AA:
        return None
    if aa == svm.single_amino_pred:
        return config.SVM_LEVEL_TO_SCORE['single_amino_pred']
    elif aa in svm.small_cluster_pred:
        return config.SVM_LEVEL_TO_SCORE['small_cluster_pred']
    elif aa in svm.large_cluster_pred:
        return config.SVM_LEVEL_TO_SCORE['large_cluster_pred']
    elif svm.physicochemical_class in config.PHYSICOCHEMICAL_CLASSES and \
            aa in config.PHYSICOCHEMICAL_CLASSES[svm.physicochemical_class]:
        return config.SVM_LEVEL_TO_SCORE['physicochemical_class']
    return config.SVM_LEVEL_TO_SCORE['not_matched']


def get_prediction_from_signature(signature, known_codes, svm_prediction, scoring_mode):
    # TODO: actually signature == svm_prediction['stachelhaus_seq'].upper() -- refactor out the excessive parameter
    scores = dict()
    aa_name_to_sorting_index = dict()  # to define some sorting order for amino acids with exactly the same score
    aa_name_to_raw_aa = dict()

    for code_metadata in known_codes:
        aa_name = __get_aa_fullname(code_metadata, mode='classic')
        if aa_name in config.FORBIDDEN_AA_SIGNATURES:
            continue
        if aa_name not in aa_name_to_sorting_index:
            aa_name_to_sorting_index[aa_name] = config.KNOWN_AA_SIGNATURES.index(code_metadata['aa'])
        if aa_name not in aa_name_to_raw_aa:
            aa_name_to_raw_aa[aa_name] = __get_aa_fullname(code_metadata, mode='raw')
        code = code_metadata['code']
        score = __get_stachelhaus_score(signature, code, mode='classic')
        if aa_name not in scores or score > scores[aa_name]:  # FIXME: should we count how many times the top score per AA was reached?
            scores[aa_name] = score

    if scoring_mode == 'hybrid':
        for aa_name in scores.keys():
            stachelhaus_score = scores[aa_name]
            svm_score = __get_svm_score(aa_name_to_raw_aa[aa_name], svm_prediction)
            if svm_score is not None:
                hybrid_score = (stachelhaus_score + svm_score) / 2.0
            else:
                hybrid_score = stachelhaus_score  # TODO: think of it: shouldn't we penalize the score in this case?
            scores[aa_name] = hybrid_score

    # TODO sort scores appropriately
    sorted_scores = OrderedDict(sorted(scores.items(),
                                       key=lambda x: (-x[1], aa_name_to_sorting_index[x[0]])))  # x = (aa_name, score)
    first_prediction = next(iter(sorted_scores.items()))
    if first_prediction[1] > config.MIN_SCORE_TO_REPORT:
        main_aa_pred = first_prediction[0]
    else:
        main_aa_pred = 'N/A'
    aa_pred_list = ';'.join(map(lambda x: '%s(%.1f)' % (x[0], x[1]), sorted_scores.items()))  # x = (aa_name, score)
    return main_aa_pred, aa_pred_list

    # how it works in antiSMASH v.5
    #     parts = line.split("\t")
    #     # EXAMPLE: ctg1_orf00006_A1        L--SFDASTLEGWLLTGGDTNGYGPTENTTFTTT      DALWLGGTFK      hydrophobic-aliphatic   gly,ala,val,leu,ile,abu,iva     val,leu,ile,abu,iva     val     N/A     orn,lys,arg     orn,horn        0       0:0     0.000000e+00
    #     # 0: sequence-id
    #     # 1: 8A-signature
    #     # 2: stachelhaus-code:
    #     # 3: 3class-pred
    #     # 4: large-class-pred
    #     # 5: small-class-pred
    #     # 6: single-class-pred
    #     # 7: nearest stachelhaus code
    #     # 8: NRPS1pred-large-class-pred
    #     # 9: NRPS2pred-large-class-pred
    #     # 10: outside applicability domain (1 or 0)
    #     # 11: coords
    #     # 12: pfam-score
    #     if not len(parts) == 13:
    #         raise ValueError("Invalid SVM result line: %s" % line)
    #     query_stach = parts[2]
    #     pred_from_stach = parts[7]
    #     best_stach_match = query_stach.lower()
    #     stach_count = 0
    #     for possible_hit in KNOWN_STACH_CODES[pred_from_stach]:
    #         # the datafile sometimes has - for the trailing char, but not all the time
    #         matches = [int(a == b) for a, b in list(zip(query_stach, possible_hit))[:9]] + [1]
    #         count = sum(matches)
    #         if count > stach_count:
    #             stach_count = count
    #             best_stach_match = "".join(c if match else c.lower() for (c, match) in zip(query_stach, matches))
    #
    # how it works in antiSMASH v.3
    # #Compare NRPS signature with database of signatures and write output to txt file
    # for k in querysignames:
    #   querysigseq = querysigseqs[querysignames.index(k)]
    #   scoredict = {}
    #   for i in signaturenames:
    #     sigseq = signatureseqs[signaturenames.index(i)]
    #     positions  = range(len(querysigseq))
    #     score = 0
    #     for j in positions:
    #       if querysigseq[j] == sigseq[j]:
    #         score += 1
    #     score = ((float(score) / 10) * 100)
    #     scoredict[i] = score
    #   sortedhits = sortdictkeysbyvalues(scoredict)
    #   sortedhits = sortedhits
    #   sortedscores = []
    #   sortedhits2 = []
    #   for i in sortedhits:
    #     score = scoredict[i]
    #     if score > 40:
    #       score = "%.0f"%(score)
    #       sortedscores.append(score)
    #       sortedhits2.append(i)
    #   allsortedhits = sortedhits
    #   sortedhits = sortedhits2
    #   #Find all other scores after best score
    #   nextbesthitsdict = {}
    #   nextbesthits = []
    #   for hit in allsortedhits:
    #     aa = hit.split("__")[-1]
    #     score = scoredict[hit]
    #     if not nextbesthitsdict.has_key(aa) and aa != "xxx":
    #       nextbesthitsdict[aa] = score
    #       nextbesthits.append(aa)
    #   #Write output to txt file
    #   outfile1.write(querysig34codes[querysignames.index(k)] + "\t" + k + "\n")
    #   if len(sortedhits) > 0:
    #     outfile2.write(k + "\t" + sortedhits[0].split("__")[-1] + "\t")
    #   else:
    #     outfile2.write(k + "\t" + "N/A" + "\t")
    #   outfile2.write(";".join(["%s(%s)" % (aa, nextbesthitsdict[aa]) for aa in nextbesthits]) + "\n")
