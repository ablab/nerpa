import os
import json
import glob
import re
from log_utils import error, info
from codes_handler import get_prediction_from_signature
from pathlib import Path
from typing import Union, Literal, TypedDict, List

SVM_HEADER = '#sequence-id<tab>8A-signature<tab>stachelhaus-code<tab>3class-pred<tab>large-class-pred<tab>small-class-pred<tab>single-class-pred<tab>nearest stachelhaus code<tab>NRPS1pred-large-class-pred<tab>NRPS2pred-large-class-pred<tab>outside applicability domain?<tab>coords<tab>pfam-score\n'  #FIXME: convert into '<tab>'.join([..])
GENE_HEADER = '\t'.join(["gene ID", "gene start", "gene end", "gene strand", "smCOG", "locus_tag", "annotation"]) + '\n'
NRPS_PKS_HEADER = '\t'.join(["Cluster_ID", "NRPSPKS_ID", "annotation", "aSDomain", "score", "evalue", "domain_type",
                             "subtype", "domain_start", "domain_end", "KR activity", "KR stereochemistry",
                             "NRPSPredictor2", "Stachelhaus", "Minowa", "pkssignature", "consensus"]) + '\n'
PREDICTIONS_TXT_DIR = 'nrpspks_predictions_txt'
FEATURES_TXT_DIR = 'txt'


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
        if 'nrpys' in domain_prediction:  # antismash >=7
            prediction_data = domain_prediction['nrpys']

            self.angstrom_code = prediction_data['aa34']
            self.stachelhaus_seq = prediction_data['aa10']
            self.physicochemical_class = prediction_data['physiochemical_class']['name']
            self.large_cluster_pred = [substrate['short'].lower()
                                       for substrate in prediction_data['large_cluster']['substrates']]
            self.small_cluster_pred = [substrate['short'].lower()
                                       for substrate in prediction_data['small_cluster']['substrates']]

            substrates = prediction_data['single_amino']['substrates']
            self.single_amino_pred = substrates[0]['short'].lower() if substrates else 'N/A'

            stachelhaus_match_count = max([round(stachelhaus_match['aa10_score'] * 10)
                                           for stachelhaus_match in prediction_data['stachelhaus_matches']],
                                           default=0)
            self.uncertain = stachelhaus_match_count < 7  # not so sure about this
        elif 'NRPSPredictor2' in domain_prediction:  # older version of antismash
            prediction_data = domain_prediction['NRPSPredictor2']

            self.angstrom_code = prediction_data['angstrom_code']
            self.stachelhaus_seq = prediction_data['stachelhaus_seq']
            self.physicochemical_class = prediction_data['physicochemical_class']
            self.large_cluster_pred = prediction_data['large_cluster_pred']
            self.small_cluster_pred = prediction_data['small_cluster_pred']
            self.single_amino_pred = prediction_data['single_amino_pred']
            self.uncertain = prediction_data['uncertain']
        else:
            raise RuntimeError('Neither "nrpys" nor "NRPSPredictor2" in domain prediction.') 

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
    def parsed_block(block):
        match = re.match("\\[<?(?P<start>[0-9]+):>?(?P<end>[0-9]+)\\](\\((?P<strand>[+-])\\))?", block)
        return {'start': int(match.group('start')),
                'end': int(match.group('end')),
                'strand': match.group("strand") or ''}

    location_trimmed = location[len('join{') : -1] if location.startswith('join{') else location[:]

    blocks = sorted([parsed_block(block.strip())
                     for block in location_trimmed.split(',')],
                    key=lambda block: block['start'])

    if not len(set(block['strand'] for block in blocks)) == 1:
        raise  # for testing: all blocks are oriented in the same way

    # merge all blocks
    start = blocks[0]['start']
    end = blocks[-1]['end']
    strand = blocks[0]['strand']
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


# TODO: I don't know what the data below is about, need to consult with someone to come up with better names
class KnownCode(TypedDict):
    code: str  # stachelhaus code
    label: str  # label and aa are both about amino acid but sometimes differ: ('dpg', 'dhpg') or ('allothr', 'thr')
    aa: str
    mods: list  # don't know what it is about: in the data I see it is always empty
    confs: list  # same as mods
    count: int  # count of what?

class ParsedPrediction(TypedDict):
    v5_name: str  # domain id in the antiSMASH 5 format
    locus_tag: str  # module id
    A: int  # index of the A-domain in the module
    signature: str  # stachelhaus code
    svm: SVM_entry  # parsed prediction data


def contig_domain_predictions(ctg_id: str, contig_data: dict) -> List[ParsedPrediction]:
    if "antismash.modules.nrps_pks" not in contig_data["modules"]:
        return []

    def parsed_prediction(domain_id, prediction) -> ParsedPrediction:
        orf_id, a_idx = __parse_amp_binding_domain(domain_id)
        svm = SVM_entry(prediction)
        return ParsedPrediction({"v5_name": "%s_AMP-binding.%s" % (orf_id, a_idx),
                                 "locus_tag": orf_id,
                                 "A": int(a_idx),
                                 "signature": svm.stachelhaus_seq,
                                 "svm": svm})

    return [parsed_prediction(domain_id, prediction)
            for domain_id, prediction in contig_data["modules"]["antismash.modules.nrps_pks"]["domain_predictions"].items()
            if 'AMP-binding' in domain_id]


def handle_single_input(antismash_results: Path,  # path to either the folder with antismash results or to the json file TODO: refactor this: a variable should not be 'either this or that'
                        maybe_output_dir,  # folder to save parsed antismash results, however, it can be None or something else TODO: this is some mess, definitely needs refactoring
                        is_root_outdir: bool,
                        naming_style: Literal['v3', 'v5', 'mix'],
                        known_codes: List[KnownCode],
                        scoring_mode: Literal['hybrid', 'stachelhaus'],
                        verbose=False):
    # TODO: this should not be a responsibility of this function
    def get_antismash_json(antismash_results: Path) -> Path:
        return Path(get_main_json_fpath(antismash_results.resolve())) \
               if antismash_results.is_dir() else \
               antismash_results

    # TODO: that should not be a responsibility of this function: after all, the output folder should be the same for all antismash records
    def get_output_dir(maybe_output_dir,
                       is_root_outdir: bool,
                       antismash_results_json: Path) -> Path:
        if maybe_output_dir is None:
            output_dir = antismash_results_json
        elif is_root_outdir:
            output_dir = Path(maybe_output_dir) / antismash_results_json.stem
        else:
            output_dir = Path(maybe_output_dir)
        return Path(__create_output_dirs(str(output_dir.resolve())))


    info('Processing ' + str(antismash_results), verbose=verbose)
    try:
        antismash_results_json: Path = get_antismash_json(antismash_results)
    except TypeError:
        error('Main antiSMASH v.5 JSON file not found in %s. Skipping this input..' % antismash_results)
        return
    output_dir: Path = get_output_dir(maybe_output_dir, is_root_outdir, antismash_results_json)

    info('Processing JSON %s, saving results to %s' % (antismash_results_json, output_dir), verbose=verbose)
    with antismash_results_json.open('r') as f:
        data = json.load(f)

    for ctg_idx, contig_data in enumerate(data["records"]):
        ctg_id = 'ctg%d' % (ctg_idx + 1)

        # part 1: parsing domain predictions and writing _codes.txt and _svm.txt files
        parsed_predictions = contig_domain_predictions(ctg_id, contig_data)
        if not parsed_predictions:
            continue  # TODO: check whether we need "empty" files for entries without NRPS/PKS just for consistency

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
        seq_record_id = contig_data['id'].split('.')[0]  # for consistency with antiSMASH v.3 naming logic, e.g. 'JNWS01000001.1' --> 'JNWS01000001'
        # TODO: check whether it is imporant to keep seq_record_id style in Nerpa;
        # otherwise it is better to have only seq_entry_id or ctg_id (depending on the naming style)
        cur_contig_gene_output_fpath = __get_contig_output_fpath(output_dir,
                                       seq_entry_id if naming_style == 'v5' else ctg_id, type='gene')
        cur_contig_NRPS_PKS_output_fpath = __get_contig_output_fpath(output_dir,
                                           seq_entry_id if naming_style == 'v5' else ctg_id, type='NRPS_PKS')

        # TODO: comment to explain what is this variable about
        regions_of_interest = [__parse_location(feature['location'])[:-1]  # omit orientation
                               for feature in contig_data['features']
                               if feature['type'] == 'region']

        features_iter = iter((feature, __parse_location(feature['location']))
                             for feature in contig_data['features']
                             if feature['type'] == 'CDS')

        num_genes_without_locus_tag = 0
        with open(cur_contig_gene_output_fpath, 'w') as gene_f:
            gene_f.write(GENE_HEADER)
            # content example: ctg1_orf00189	127377	128739	-		ctg1_orf00189	unannotated orf
            regions_of_interest_iter = iter(regions_of_interest)
            features_iter = iter((feature, __parse_location(feature['location']))
                                 for feature in contig_data['features']
                                 if feature['type'] == 'CDS')

            region_start, region_end = 0, -1  # placeholder to enter the while loop
            while True:
                try:
                    feature, (cds_start, cds_end, cds_strand) = next(features_iter)
                    # find the next intersecting (cds, region) pair
                    while region_end < cds_start or cds_end < region_start:
                        if region_end < cds_start:
                            region_start, region_end = next(regions_of_interest_iter)
                        if cds_end < region_start:
                            feature, (cds_start, cds_end, cds_strand) = next(features_iter)
                except StopIteration:
                    break  # no more intersecting (cds, region) pairs

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
                                           [orf_id, cds_start, cds_end, cds_strand, '', orf_id, 'unannotated orf']))
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
                    cds_start, cds_end, cds_strand = __parse_location(feature['location'])
                    while cur_reg_idx < len(regions_of_interest) and cds_start > regions_of_interest[cur_reg_idx][1]:
                        cur_reg_idx += 1
                    if cur_reg_idx == len(regions_of_interest):
                        break
                    if cds_end < regions_of_interest[cur_reg_idx][0]:
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
                    if 'domain_subtypes' in feature['qualifiers']:  # for antiSMASH v.7 (and hopefully newer)
                        entry.subtype = feature['qualifiers']['domain_subtypes'][0]
                    elif 'domain_subtype' in feature['qualifiers']:  # for antiSMASH v.6
                        entry.subtype = feature['qualifiers']['domain_subtype'][0]
                    else:  # manually for antiSMASH v.5 and older
                        if entry.domain_type.startswith('Condensation'):
                            entry.subtype = entry.domain_type
                            entry.domain_type = entry.domain_type.split('_')[0]
                        elif entry.domain_type.endswith('MT'):
                            entry.subtype = entry.domain_type
                            entry.domain_type = 'MT'
                    entry.domain_start = cds_start
                    entry.domain_end = cds_end
                    # TODO: process the rest fields (KR_activity, etc for PK and NRPSPredictor2, etc for AMP-binding)
                    nrps_pks_f.write(str(entry) + '\n')

    info('Done with %s, see results in %s' % (antismash_results_json, output_dir), verbose=verbose)
    return str(output_dir)


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
    # rather rare and strange case, but we should be ready to change the dir name if it is already occupied
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
