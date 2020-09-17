import os
import json
from log_utils import error, info
from codes_handler import get_prediction_from_signature


SVM_HEADER = '#sequence-id<tab>8A-signature<tab>stachelhaus-code<tab>3class-pred<tab>large-class-pred<tab>small-class-pred<tab>single-class-pred<tab>nearest stachelhaus code<tab>NRPS1pred-large-class-pred<tab>NRPS2pred-large-class-pred<tab>outside applicability domain?<tab>coords<tab>pfam-score\n'
GENE_HEADER = '\t'.join(["gene ID", "gene start", "gene end", "gene strand", "smCOG", "locus_tag", "annotation"]) + '\n'
NRPS_PKS_HEADER = ''
PREDICTIONS_TXT_DIR = 'nrpspks_predictions_txt'
FEATURES_TXT_DIR = 'txt'


def __get_svm_results(domain_prediction):
    '''
    # example1 (domain_prediction):
   {'method': 'NRPSPredictor2',
    'angstrom_code': 'FEQSFDASIFEQFIVFGGECNAYGPTEATIMANI',
    'physicochemical_class': 'hydrophobic-aliphatic',
    'large_cluster_pred': ['gly', 'ala', 'val', 'leu', 'ile', 'abu', 'iva'],
    'small_cluster_pred': ['N/A'],
    'single_amino_pred': 'N/A',
    'stachelhaus_predictions': ['phe'],
    'uncertain': False,
    'stachelhaus_seq': 'DAFfVgAimK',
    'stachelhaus_match_count': 6}

    # example2 (domain_prediction):
   {'method': 'NRPSPredictor2',
    'angstrom_code': 'RWMTFDVSVWEWHFFVSGEHNLYGPTEASVDVTS',
    'physicochemical_class': 'hydrophobic-aliphatic',
    'large_cluster_pred': ['ser', 'thr', 'dhpg', 'hpg'],
    'small_cluster_pred': ['ser'],
    'single_amino_pred': 'ser',
    'stachelhaus_predictions': ['ser'],
    'uncertain': False,
    'stachelhaus_seq': 'DVWHFSLVDK',
    'stachelhaus_match_count': 10}

    # example (output):
#sequence-id<tab>8A-signature<tab>stachelhaus-code<tab>3class-pred<tab>large-class-pred<tab>small-class-pred<tab>single-class-pred<tab>nearest stachelhaus code<tab>NRPS1pred-large-class-pred<tab>NRPS2pred-large-class-pred<tab>outside applicability domain?<tab>coords<tab>pfam-score
ctg1_orf00225_A1	FEQSFDASIFEQFIVFGGECNAYGPTEATIMANI	DAFFVGAIMK	hydrophobic-aliphatic	gly,ala,val,leu,ile,abu,iva	N/A	N/A	N/A	gly,ala,val,leu,ile,abu,iva	N/A	0	0:0	0.000000e+00
ctg1_orf00225_A2	RWMTFDVSVWEWHFFVSGEHNLYGPTEASVDVTS	DVWHFSLVDK	hydrophobic-aliphatic	ser,thr,dhpg,hpg	ser	ser	N/A	ser,thr,dhpg,hpg	ser	0	0:0	0.000000e+00

    '''
    return '\t'.join([domain_prediction['angstrom_code'],
                      domain_prediction['stachelhaus_seq'].upper(),
                      domain_prediction['physicochemical_class'],
                      ','.join(domain_prediction['large_cluster_pred']),
                      ','.join(domain_prediction['small_cluster_pred']),
                      domain_prediction['single_amino_pred'],
                      'N/A', 'N/A', 'N/A',  # can't determine from antiSMASH v.5 output
                      str(int(bool(domain_prediction['uncertain']))),
                      '0:0', '0.000000e+00'  # seems that these columns are always the same
                      ])


def handle_single_input(path, output_dir, known_codes, verbose=False):
    info('Processing ' + path, verbose=verbose)
    path = os.path.abspath(path)
    main_json_path = path
    if os.path.isdir(path):
        # TODO: if there is a single JSON inside the dir, just take it independently of its name
        main_json_path = os.path.join(path, os.path.basename(os.path.normpath(path)) + '.json')
    if not os.path.isfile(main_json_path):
        error('Main antiSMASH v.5 JSON file not found: %s. Skipping this input..' % main_json_path)
        return
    if output_dir is None:
        output_dir = os.path.dirname(main_json_path)
    output_dir = os.path.abspath(output_dir)
    __create_output_dirs(output_dir)
    info('Processing JSON %s, saving results to %s' % (main_json_path, output_dir), verbose=verbose)
    with open(main_json_path, 'r') as f:
        data = json.load(f)

    for contig_data in data["records"]:
        # TODO: process features as well and output them in antiSMASH v.3-compatible style in ./txt/ subdirectory
        # example:
        # 'features' (list of len 658)
        #     000: {'location': '[0:40388]', 'type': 'protocluster', 'id': '<unknown id>', 'qualifiers': {'aStool': ['rule-based-clusters'], 'contig_edge': ['True'], 'core_location': ['[3484:20388]'], 'cutoff': ['20000'], 'detection_rule': ['(cds(Condensation and (AMP-binding or A-OX)) or (Condensation and AMP-binding))'], 'neighbourhood': ['20000'], 'product': ['NRPS'], 'protocluster_number': ['1'], 'tool': ['antismash']}}
        #     001: {'location': '[3484:20388]', 'type': 'proto_core', 'id': '<unknown id>', 'qualifiers': {'aStool': ['rule-based-clusters'], 'tool': ['antismash'], 'cutoff': ['20000'], 'detection_rule': ['(cds(Condensation and (AMP-binding or A-OX)) or (Condensation and AMP-binding))'], 'neighbourhood': ['20000'], 'product': ['NRPS'], 'protocluster_number': ['1']}}

        if "antismash.modules.nrps_pks" in contig_data["modules"]:
            # part 1: parsing domain predictions and writing _codes.txt and _svm.txt files
            parsed_predictions = []
            if "domain_predictions" in contig_data["modules"]["antismash.modules.nrps_pks"]:
                domain_predictions = contig_data["modules"]["antismash.modules.nrps_pks"]["domain_predictions"]
                for prediction in domain_predictions.keys():
                    if "AMP-binding" in prediction:
                        prefix, ctg_id, orf_idx, amp_binding = prediction.split('_')
                        a_idx = amp_binding.split('.')[1]
                        stachelhaus_seq = domain_predictions[prediction]["NRPSPredictor2"]["stachelhaus_seq"].upper()
                        parsed_predictions.append({"orf": int(orf_idx), "A": int(a_idx),
                                                   "signature": stachelhaus_seq,
                                                   "svm": __get_svm_results(domain_predictions[prediction]["NRPSPredictor2"])})

                        # example:
                        # {'nrpspksdomains_ctg1_5_AMP-binding.1':
                        #      {'NRPSPredictor2':
                        #           {'method': 'NRPSPredictor2',
                        #             'angstrom_code': 'L--SFDASTLEGWLLTGGDTNGYGPTENTTFTTT',
                        #             'physicochemical_class': 'hydrophobic-aliphatic',
                        #             'large_cluster_pred': ['gly', 'ala',
                        #                                    'val', 'leu',
                        #                                    'ile', 'abu',
                        #                                    'iva'],
                        #             'small_cluster_pred': ['val', 'leu',
                        #                                    'ile', 'abu',
                        #                                    'iva'],
                        #             'single_amino_pred': 'val',
                        #             'stachelhaus_predictions': ['val'],
                        #             'uncertain': False,
                        #             'stachelhaus_seq': 'DALWLGGTFK',
                        #             'stachelhaus_match_count': 10}},

            if parsed_predictions:
                # FIXME: remove: for debugging purposes only
                # if ctg_id == 'ctg12':
                #     t = 1
                info('\tprocessing contig: %s' % ctg_id, verbose=verbose)
                parsed_predictions.sort(key=lambda x: (x["orf"], x["A"]))
                cur_contig_codes_output_fpath = __get_contig_output_fpath(output_dir, ctg_id, type='codes')
                cur_contig_svm_output_fpath = __get_contig_output_fpath(output_dir, ctg_id, type='svm')
                with open(cur_contig_codes_output_fpath, 'w') as codes_f:
                    with open(cur_contig_svm_output_fpath, 'w') as svm_f:
                        svm_f.write(SVM_HEADER)
                        for prediction in parsed_predictions:
                            entry_id = __get_entry_id(ctg_id, prediction["orf"], prediction["A"])
                            main_aa_pred, aa_pred_list = get_prediction_from_signature(prediction["signature"], known_codes)
                            codes_f.write('\t'.join([entry_id, main_aa_pred, aa_pred_list]) + '\n')
                            svm_f.write('\t'.join([entry_id, prediction["svm"]]) + '\n')
                            info('\t\tprocessed ORF: %s, A-domain: %s, Stachelhaus code: %s' %
                                 (prediction["orf"], prediction["A"], prediction["signature"]), verbose=verbose)

            # part 2: parsing features and writing _genes.txt and _NRPS_PKS.txt files
            seq_record_id = contig_data['id'].split('.')[0]  # for consistency with antiSMASH v.3 naming logic, e.g. 'JNWS01000001.1' --> 'JNWS01000001'
            cur_contig_gene_output_fpath = __get_contig_output_fpath(output_dir, seq_record_id, type='gene')
            cur_contig_NRPS_PKS_output_fpath = __get_contig_output_fpath(output_dir, seq_record_id, type='NRPS_PKS')

            regions_of_interest = []
            for feature in contig_data['features']:
                if feature['type'] == 'region':
                    location = feature['location']  # e.g. 'location' = '[351:486](+)'
                    coords = location.split(']')[0][1:]
                    start, end = map(int, coords.split(':'))
                    regions_of_interest.append((start, end))

            with open(cur_contig_gene_output_fpath, 'w') as gene_f:
                gene_f.write(GENE_HEADER)
                # features1 = contig_data['features'][:300]
                # features2 = contig_data['features'][300:600]
                # features3 = contig_data['features'][600:900]
                # content example: ctg1_orf00189	127377	128739	-		ctg1_orf00189	unannotated orf
                cur_region = (0, 0)
                for feature in contig_data['features']:
                    if feature['type'] == 'CDS':
                        location = feature['location']  # e.g. 'location' = '[351:486](+)'
                        coords = location.split(']')[0][1:]
                        start, end = map(int, coords.split(':'))
                        while start > cur_region[1]:
                            if regions_of_interest:
                                cur_region = regions_of_interest.pop(0)
                            else:
                                cur_region = None
                                break
                        if cur_region is None:
                            break
                        if end < cur_region[0]:
                            continue

                        strand = location.split('](')[1][0]
                        ctg_id, orf_idx = feature['qualifiers']['locus_tag'][0].split('_')  # e.g. 'locus_tag' = ['ctg1_1']
                        orf_id = __get_entry_id(ctg_id, orf_idx)
                        gene_f.write('\t'.join(map(str,
                                                   [orf_id, start, end, strand, '', orf_id, 'unannotated orf']))
                                     + '\n')

    info('Done with %s, see results in %s' % (main_json_path, output_dir), verbose=verbose)


def __get_entry_id(ctg_id, orf_idx, a_idx=None):
    orf_id ="%s_orf%05d" % (str(ctg_id), int(orf_idx))
    if a_idx is not None:
        return orf_id + "_A%d" % a_idx
    return orf_id


def __get_contig_output_fpath(output_dir, ctg_id, type='codes'):
    if type in ['codes', 'svm']:
        return os.path.join(output_dir, PREDICTIONS_TXT_DIR, ctg_id + "_nrpspredictor2_%s.txt" % type)
    elif type in ['gene', 'NRPS_PKS']:
        return os.path.join(output_dir, FEATURES_TXT_DIR, ctg_id + "_%s.txt" % type)
    return None


def __create_output_dirs(base_output_dir):
    for subdir in [PREDICTIONS_TXT_DIR, FEATURES_TXT_DIR]:
        output_dir = os.path.join(base_output_dir, subdir)
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

