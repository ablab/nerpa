import itertools
from collections import OrderedDict
from os.path import abspath, dirname, realpath, join

root_dir = abspath(dirname(realpath(__file__)))

DEFAULT_INPUT_CODES = join(root_dir, 'codes', 'nrpspredictor2_knowncodes.fasta')  # or 'nrpspredictor2_labeled_sigs.txt') or 'sandpuma_mibig_based.fna'
DEFAULT_OUTPUT_CODES = join(root_dir, 'codes', 'parsed_codes.json')
STANDARD_STACHELHAUS_CODE_LENGTH = 10

KNOWN_AA_SIGNATURES = [
    # 20 standard amino acids
    'ala', 'arg', 'asn', 'asp', 'cys', 'gln', 'glu', 'gly', 'his', 'ile',
    'leu', 'lys', 'met', 'phe', 'pro', 'ser', 'thr', 'trp', 'tyr', 'val',
    # beta-versions
    'b-ala', 'b-lys',
    # common non-proteinogenic amino acids
    'orn', 'aad', 'abu', 'hiv', 'pip', 'dab', 'dhb', 'dhpg', 'hpg', 'dht', 'bht', 'sal',
    # acceptable but very rare (appears only once in NRPSpredictor codes)
    'aeo', 'lyserg', 'cap', 'LDAP', '4ppro',
    '2-oxo-isovaleric-acid', 'alaninol',
    'iva', 'phg', 'tcl', 'vol', 'HICA',
    # additional from Sandpuma -- need to carefully re-process manually their codes!
    'bmt', 'piperazic', 'hyv', 'hse', 'kyn', 'hty', 'glycolic-acid',
    'indole-3-carboxylic', 'phe-ac', 'cysa', 'dpr', '33p-l-ala',
    # ERRONEOUS ENTRIES IN SANDPUMA LIST, added just for the sake of completeness
    'glycy',
    # additional from nrpspredictor2 KNOWNCODES
    'apa', 'apc', 'gua', 'cit', 'end', 'uda', 'cha', 'ahp',
    'xxx', # special case -- non-functional A-domain prediction
    ]

FORBIDDEN_AA_SIGNATURES = ['xxx']  # will be ignored

SPECIAL_CASES = {'N-(1,1-dimethyl-1-allyl)Trp': ('trp', ['N-(1,1-dimethyl-1-allyl)']),
                 '3-me-glu': ('glu', ['me']),
                 'hpg2cl': ('hpg', ['cl']),
                 'pro-4p': ('4ppro', []),
                 'sar': ('gly', ['me'])}

SYNONYMS_AA_SIGNATURES = {'dhab': 'dht', 'dpg': 'dhpg', 'dhp': 'dhpg', 'beta-ala': 'b-ala', 'blys': 'b-lys',
                          'athr': 'allo-thr', 'allothr': 'allo-thr', 'aile': 'allo-ile', 'alloile': 'allo-ile',
                          'alpha-hydroxy-isocaproic-acid': 'HICA', 'gly-ph': 'phg',
                          'd-lyserg': 'lyserg', 'lys-erg': 'lyserg', 'val-vol': 'vol'}

MODS = OrderedDict(
    {'dh': 'dehydro', 'cl': 'chloro', 'oh': 'hydroxy', # note! two keys for 'hydroxy'
     'do': 'deoxy',
     'h': 'hydroxy', 'a': 'acetyl', 'f': 'formyl', 'me': 'methyl',
     'N-(1,1-dimethyl-1-allyl)': 'N-(1,1-dimethyl-1-allyl)'})
# technically, Oxidation ('o': 'oxy') is also possible

CONF_MODS = ['allo', 'D', 'L']

## GENERAL AA DESCRIPTION FORMAT:
# [@CONF]BASE[+MOD], base may include numbers and "-", both CONF and MOD may be a list of more than 1 item
# Examples:
# @allo@Dthr+me+h  # allo-threonine in D-configuration with methylation and hydroxylation
# thr              # just threonine
# @Lb-ala          # beta-alanine in L-configuration
# b-lys+me+me      # beta-lysine with dimethylation

MIN_SCORE_TO_REPORT = 40.5  # for classic NRPSPredictor2 representation of a single AA in front of the ranked list of AAs

# NRPSPredictor2 SVM classes and clusters (from Roettig et al, 2011)
PHYSICOCHEMICAL_CLASSES = {'hydrophobic-aliphatic': ['ala', 'gly', 'val', 'leu', 'ile', 'abu', 'iva', 'ser',
                                                     'thr', 'hpg', 'dhpg', 'cys', 'pro', 'pip'],
                           'hydrophilic': ['arg', 'asp', 'glu', 'his', 'asn', 'lys', 'gln', 'orn', 'aad'],
                           'hydrophobic-aromatic': ['phe', 'tyr', 'trp', 'dhb', 'phg', 'bht']}

LARGE_CLUSTERS = {'Hydroxy-benzoic acid derivates': ['dhb', 'sal'],
                  'Polar, uncharged (aliphatic with -SH)': ['cys'],
                  'Aliphatic chain or phenyl group with -OH': ['ser', 'thr', 'dhpg', 'hpg'],
                  'Aliphatic chain with H-bond donor': ['asp', 'asn', 'glu', 'gln', 'aad'],
                  'Apolar, aliphatic': ['gly', 'ala', 'val', 'leu', 'ile', 'abu', 'iva'],
                  'Aromatic side chain': ['phe', 'trp', 'phg', 'tyr', 'bht'],
                  'Cyclic aliphatic chain (polar NH2 group)': ['pro', 'pip'],
                  'Long positively charged side chain': ['orn', 'lys', 'arg']}

SMALL_CLUSTERS = {'2-amino-adipic acid': ['aad'],
                  'Dhb, Sal': ['dhb', 'sal'],
                  'Polar, uncharged (hydroxy-phenyl)': ['dhpg', 'hpg'],
                  'Cys': ['cys'],
                  'Serine-specific': ['ser'],
                  'Threonine-specific': ['thr'],
                  'Asp-Asn': ['asp', 'asn'],
                  'Orn and hydroxy-Orn specific': ['orn'],
                  'Aliphatic, branched hydrophobic': ['val', 'leu', 'ile', 'abu', 'iva'],
                  'Tiny, hydrophilic, transition to aliphatic': ['gly', 'ala'],
                  'Pro-specific': ['pro'],
                  'Polar aromatic ring': ['tyr', 'bht'],
                  'Glu-Gln': ['glu', 'gln'],
                  'Arg-specific': ['arg'],
                  'Unpolar aromatic ring': ['phe', 'trp']}

SVM_LEVEL_TO_SCORE = {'single_amino_pred': 100.0,
                      'small_cluster_pred': 90.0,
                      'large_cluster_pred': 80.0,
                      'physicochemical_class': 50.0,
                      'not_matched': 0.0}

SVM_DETECTABLE_AA = set(itertools.chain.from_iterable(PHYSICOCHEMICAL_CLASSES.values()))
