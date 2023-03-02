import os
import json
import glob
import re
from log_utils import error, info
from codes_handler import get_prediction_from_signature


SVM_HEADER = '#sequence-id<tab>8A-signature<tab>stachelhaus-code<tab>3class-pred<tab>large-class-pred<tab>small-class-pred<tab>single-class-pred<tab>nearest stachelhaus code<tab>NRPS1pred-large-class-pred<tab>NRPS2pred-large-class-pred<tab>outside applicability domain?<tab>coords<tab>pfam-score\n'  #FIXME: convert into '<tab>'.join([..])
GENE_HEADER = '\t'.join(["gene ID", "gene start", "gene end", "gene strand", "smCOG", "locus_tag", "annotation"]) + '\n'
NRPS_PKS_HEADER = '\t'.join(["Cluster_ID", "NRPSPKS_ID", "annotation", "aSDomain", "score", "evalue", "domain_type",
                             "subtype", "domain_start", "domain_end", "KR activity", "KR stereochemistry",
                             "NRPSPredictor2", "Stachelhaus", "Minowa", "pkssignature", "consensus"]) + '\n'
PREDICTIONS_TXT_DIR = 'nrpspks_predictions_txt'
FEATURES_TXT_DIR = 'txt'


def __get_svm_results(domain_prediction):

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


class SVM_entry:
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

    def __init__(self, domain_prediction):
        self.angstrom_code = domain_prediction['angstrom_code']
        self.stachelhaus_seq = domain_prediction['stachelhaus_seq']
        self.physicochemical_class = domain_prediction['physicochemical_class']
        self.large_cluster_pred = domain_prediction['large_cluster_pred']
        self.small_cluster_pred = domain_prediction['small_cluster_pred']
        self.single_amino_pred = domain_prediction['single_amino_pred']
        self.uncertain = domain_prediction['uncertain']

    def __str__(self):
        return '\t'.join([self.angstrom_code,
                          self.stachelhaus_seq.upper(),
                          self.physicochemical_class,
                          ','.join(self.large_cluster_pred),
                          ','.join(self.small_cluster_pred),
                          self.single_amino_pred,
                          'N/A', 'N/A', 'N/A',  # can't determine from antiSMASH v.5 output
                          str(int(bool(self.uncertain))),
                          '0:0', '0.000000e+00'  # seems that these columns are always the same
                          ])

class NRPS_PKS_entry:
    def __init__(self):
        self.Cluster_ID = ''
        self.NRPSPKS_ID = ''
        self.annotation = ''
        self.aSDomain = ''
        self.score = ''
        self.evalue = ''
        self.domain_type = ''
        self.subtype = ''
        self.domain_start = ''
        self.domain_end = ''
        self.KR_activity = ''
        self.KR_stereochemistry = ''
        self.NRPSPredictor2 = ''
        self.Stachelhaus = ''
        self.Minowa = ''
        self.pkssignature = ''
        self.consensus = ''

    def __str__(self):
        return '\t'.join(map(str,
                             [self.Cluster_ID, self.NRPSPKS_ID, self.annotation,
                              self.aSDomain, self.score, self.evalue, self.domain_type,
                              self.subtype, self.domain_start, self.domain_end,
                              self.KR_activity, self.KR_stereochemistry,
                              self.NRPSPredictor2, self.Stachelhaus, self.Minowa,
                              self.pkssignature, self.consensus]))


def __parse_location(location):
    # e.g. 'location' = '[351:486](+)'
    match = re.match("\\[<?(?P<start>[0-9]+):>?(?P<end>[0-9]+)\\](\\((?P<strand>[+-])\\))?", location)
    start = int(match.group("start"))
    end = int(match.group("end"))
    strand = match.group("strand") or ''
    return start, end, strand


def __parse_amp_binding_domain(prediction):
    match = re.match("nrpspksdomains_(?P<orf_id>.*)_AMP-binding\\.(?P<a_idx>[0-9]+)$", prediction)
    return match.group("orf_id"), match.group("a_idx")


def __parse_locus_tag(locus_tag):
    match = re.match("^(?P<ctg_id>.*)_(?P<orf_idx>[^_]*)$", locus_tag)
    return match.group("ctg_id"), match.group("orf_idx")


def get_main_json_fpath(dirpath):
    if not os.path.isdir(dirpath):
        return None
    tentative_json_fpath = os.path.join(dirpath, os.path.basename(os.path.normpath(dirpath)) + '.json')
    if os.path.isfile(tentative_json_fpath):
        return tentative_json_fpath
    all_jsons_in_dir = list(glob.glob(os.path.join(dirpath, "*.json")))
    if len(all_jsons_in_dir) != 1:
        return None
    return all_jsons_in_dir[0]


def handle_single_input(path, output_dir, is_root_outdir, naming_style, known_codes, scoring_mode, verbose=False):
    info('Processing ' + path, verbose=verbose)
    path = os.path.abspath(path)
    main_json_path = path
    if os.path.isdir(path):
        main_json_path = get_main_json_fpath(path)
    if main_json_path is None or not os.path.isfile(main_json_path):
        error('Main antiSMASH v.5 JSON file not found in %s. Skipping this input..' % path)
        return
    if output_dir is None:
        output_dir = os.path.dirname(main_json_path)
    elif is_root_outdir:
        output_dir = os.path.join(output_dir, os.path.splitext(os.path.basename(main_json_path))[0])
    output_dir = os.path.abspath(output_dir)
    output_dir = __create_output_dirs(output_dir)
    info('Processing JSON %s, saving results to %s' % (main_json_path, output_dir), verbose=verbose)
    with open(main_json_path, 'r') as f:
        data = json.load(f)
    for ctg_idx, contig_data in enumerate(data["records"]):
        ctg_id = 'ctg%d' % (ctg_idx + 1)
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
                        # ctg_id, orf_idx, a_idx = __parse_amp_binding_domain(prediction)
                        orf_id, a_idx = __parse_amp_binding_domain(prediction)
                        stachelhaus_seq = domain_predictions[prediction]["NRPSPredictor2"]["stachelhaus_seq"].upper()
                        parsed_predictions.append({"v5_name": "%s_AMP-binding.%s" % (orf_id, a_idx),
                                                   "locus_tag": orf_id, "A": int(a_idx),
                                                   "signature": stachelhaus_seq,
                                                   "svm": SVM_entry(domain_predictions[prediction]["NRPSPredictor2"])})
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
                seq_entry_id = contig_data['id']
                info('\tprocessing contig (%s): %s' % (seq_entry_id, ctg_id), verbose=verbose)
                parsed_predictions.sort(key=lambda x: (x["locus_tag"], x["A"]))  # FIXME: "locus_tag" could be unsortable, e.g., ctg1_19 vs ctg1_2; TODO: check whether we can use coordinates (what happens if multi-entry fasta is on input)
                cur_contig_codes_output_fpath = __get_contig_output_fpath(output_dir,
                                                seq_entry_id if naming_style == 'v5' else ctg_id, type='codes')
                cur_contig_svm_output_fpath = __get_contig_output_fpath(output_dir,
                                              seq_entry_id if naming_style == 'v5' else ctg_id, type='svm')
                with open(cur_contig_codes_output_fpath, 'w') as codes_f:
                    with open(cur_contig_svm_output_fpath, 'w') as svm_f:
                        svm_f.write(SVM_HEADER)
                        for prediction in parsed_predictions:
                            entry_id = __get_entry_id(ctg_id, prediction["locus_tag"], prediction["A"]) \
                                if naming_style == 'v3' else prediction["v5_name"]
                            main_aa_pred, aa_pred_list = get_prediction_from_signature(prediction["signature"],
                                                                                       known_codes, prediction["svm"],
                                                                                       scoring_mode)
                            codes_f.write('\t'.join([entry_id, main_aa_pred, aa_pred_list]) + '\n')
                            svm_f.write('\t'.join([entry_id, str(prediction["svm"])]) + '\n')
                            info('\t\tprocessed (%s) ORF: %s, A-domain: %s, Stachelhaus code: %s' %
                                 (prediction["v5_name"], prediction["locus_tag"], prediction["A"], prediction["signature"]),
                                 verbose=verbose)

            # part 2: parsing features and writing _genes.txt and _NRPS_PKS.txt files
            if not parsed_predictions:
                continue  # TODO: check whether we need "empty" files for entries without NRPS/PKS just for consistency
            seq_record_id = contig_data['id'].split('.')[0]  # for consistency with antiSMASH v.3 naming logic, e.g. 'JNWS01000001.1' --> 'JNWS01000001'
            # TODO: check whether it is imporant to keep seq_record_id style in Nerpa;
            # otherwise it is better to have only seq_entry_id or ctg_id (depending on the naming style)
            cur_contig_gene_output_fpath = __get_contig_output_fpath(output_dir,
                                           seq_entry_id if naming_style == 'v5' else ctg_id, type='gene')
            cur_contig_NRPS_PKS_output_fpath = __get_contig_output_fpath(output_dir,
                                               seq_entry_id if naming_style == 'v5' else ctg_id, type='NRPS_PKS')

            regions_of_interest = []
            for feature in contig_data['features']:
                if feature['type'] == 'region':
                    start, end, _ = __parse_location(feature['location'])
                    regions_of_interest.append((start, end))

            num_genes_without_locus_tag = 0
            with open(cur_contig_gene_output_fpath, 'w') as gene_f:
                gene_f.write(GENE_HEADER)
                # content example: ctg1_orf00189	127377	128739	-		ctg1_orf00189	unannotated orf
                cur_reg_idx = 0
                for feature in contig_data['features']:
                    if feature['type'] == 'CDS':
                        start, end, strand = __parse_location(feature['location'])
                        while cur_reg_idx < len(regions_of_interest) and start > regions_of_interest[cur_reg_idx][1]:
                            cur_reg_idx += 1
                        if cur_reg_idx == len(regions_of_interest):
                            break
                        if end < regions_of_interest[cur_reg_idx][0]:
                            continue

                        if 'locus_tag' in feature['qualifiers']:  # antiSMASH v.5 and earlier
                            locus_tag = feature['qualifiers']['locus_tag'][0]
                        # antiSMASH v.6 and (hopefully) newer
                        elif 'gene' in feature['qualifiers']:
                            locus_tag = feature['qualifiers']['gene'][0]
                        elif 'protein_id' in feature['qualifiers']:
                            locus_tag = feature['qualifiers']['protein_id'][0]
                        else:
                            num_genes_without_locus_tag += 1
                            locus_tag = 'None_%d' % num_genes_without_locus_tag
                        orf_id = __get_entry_id(ctg_id, locus_tag) if naming_style == 'v3' else locus_tag
                        gene_f.write('\t'.join(map(str,
                                                   [orf_id, start, end, strand, '', orf_id, 'unannotated orf']))
                                     + '\n')
            with open(cur_contig_NRPS_PKS_output_fpath, 'w') as nrps_pks_f:
                nrps_pks_f.write(NRPS_PKS_HEADER)
                # features1 = contig_data['features'][:300]
                # features2 = contig_data['features'][300:600]
                # features3 = contig_data['features'][600:900]
                # # content example 1: JNWS01000001.1_c1	ctg1_orf00198	PKS/NRPS-like protein	nrpspksdomains_ctg1_orf00198_KR1	21.1	6.10E-07	PKS_KR		134108	134471	inactive	C1
                # # content example 1: JNWS01000001.1_c1	ctg1_orf00225	NRPS	nrpspksdomains_ctg1_orf00225_Xdom17	104.6	2.70E-32	Thioesterase		147794	148490
                # t = 1

                cur_reg_idx = 0
                last_CDS_before_aSDomain = None
                for feature in contig_data['features']:
                    if feature['type'] == 'CDS':
                        last_CDS_before_aSDomain = feature
                    elif feature['type'] == 'aSDomain':
                        start, end, strand = __parse_location(feature['location'])
                        while cur_reg_idx < len(regions_of_interest) and start > regions_of_interest[cur_reg_idx][1]:
                            cur_reg_idx += 1
                        if cur_reg_idx == len(regions_of_interest):
                            break
                        if end < regions_of_interest[cur_reg_idx][0]:
                            continue

                        locus_tag = feature['qualifiers']['locus_tag'][0]
                        entry = NRPS_PKS_entry()
                        entry.Cluster_ID = contig_data['id'] + '_c%d' % (cur_reg_idx + 1)
                        entry.NRPSPKS_ID = __get_entry_id(ctg_id, locus_tag) if naming_style == 'v3' else locus_tag
                        if last_CDS_before_aSDomain is not None:
                            if 'qualifiers' in last_CDS_before_aSDomain and 'NRPS_PKS' in last_CDS_before_aSDomain['qualifiers']:
                                NRPS_PKS_type = last_CDS_before_aSDomain['qualifiers']['NRPS_PKS'][-1]
                                entry.annotation = NRPS_PKS_type.split(' ')[1]
                        entry.aSDomain = feature['qualifiers']['domain_id'][0].replace(locus_tag, entry.NRPSPKS_ID)
                        if 'AMP-binding.' in entry.aSDomain and naming_style == 'v3':
                            entry.aSDomain = entry.aSDomain.replace('AMP-binding.', 'A')
                        entry.score = feature['qualifiers']['score'][0]
                        entry.evalue = feature['qualifiers']['evalue'][0]
                        entry.domain_type = feature['qualifiers']['aSDomain'][0]
                        # processing special cases with subtypes
                        if entry.domain_type.startswith('Condensation'):
                            entry.subtype = entry.domain_type
                            entry.domain_type = entry.domain_type.split('_')[0]
                        elif entry.domain_type.endswith('MT'):
                            entry.subtype = entry.domain_type
                            entry.domain_type = 'MT'
                        if 'domain_subtype' in feature['qualifiers']:  # for antiSMASH v.6 and (hopefully) newer
                            entry.subtype = feature['qualifiers']['domain_subtype'][0]
                        entry.domain_start = start
                        entry.domain_end = end
                        # TODO: process the rest fields (KR_activity, etc for PK and NRPSPredictor2, etc for AMP-binding)
                        nrps_pks_f.write(str(entry) + '\n')

    info('Done with %s, see results in %s' % (main_json_path, output_dir), verbose=verbose)
    return output_dir


def __get_entry_id(ctg_id, locus_tag, a_idx=None):
    if locus_tag.startswith(ctg_id):
        orf_id = locus_tag
    else:
        orf_idx = locus_tag
        # ctg_id, orf_idx = __parse_locus_tag(locus_tag)
        try:
            orf_id ="%s_orf%05d" % (str(ctg_id), int(orf_idx))
        except ValueError:
            orf_id = "%s_%s" % (str(ctg_id), str(orf_idx))
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
    # rather rare and trange case, but we should be ready to change the dir name if it is already occupied
    if os.path.isdir(base_output_dir):
        i = 2
        base_dirpath = base_output_dir
        while os.path.isdir(base_output_dir):
            base_output_dir = str(base_dirpath) + '__' + str(i)
            i += 1

    for subdir in [PREDICTIONS_TXT_DIR, FEATURES_TXT_DIR]:
        output_dir = os.path.join(base_output_dir, subdir)
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

    return base_output_dir
