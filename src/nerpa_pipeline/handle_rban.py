import sys
import os
import json
import csv
from collections import defaultdict
import subprocess
from logger import log

import handle_DL


# TODO: move this section to some sort of config file
# see rban.src.main.java.molecules.bond.BondType for the full list
PNP_BONDS = ['AMINO', 'AMINO_TRANS', 'ESTER', 'THIOETHER', 'OXAZOLE_CYCLE', 'THIAZOLE_CYCLE']
O_BONDS = ['ESTER', 'THIOETHER']
UNDEFINED_NAME = 'NA'


class NumNodesError(Exception):
    def __init__(self, identifier, coverage):
        self.identifier = identifier
        self.coverage = coverage

    def __str__(self):
        return f'Poor monomer coverage for id="{self.identifier}": only {self.coverage} monomers identified. ' \
               f'Skipping.'


def process_single_record(rban_record, nerpa_monomers, allowed_bond_types, hybrid_monomers_dict, na=UNDEFINED_NAME, min_recognized_nodes=1):
    g = rban_record['monomericGraph']['monomericGraph']
    nodes = []
    nodes = [m['monomer']['monomer']['monomer'] for m in g['monomers']]
    is_pk = [False] * len(nodes)

    try:
        is_d = handle_DL.get_monomers_chirality(rban_record)
    except Exception as e:
        log.warn(f'Unable to get isomeric information for {rban_record["id"]} : {e}')
        is_d = [None] * len(nodes)


    for k, v in hybrid_monomers_dict.items():
        nodes[k] = v
        is_d[k] = None
        is_pk[k] = True

    pnp_edges = []
    oedges = []
    other_edges = defaultdict(list)
    for bond in g['bonds']:
        s, e = bond['bond']['monomers']
        s -= 1
        e -= 1
        if bond['bond']['bondTypes'][0] in O_BONDS:
            oedges.append((s, e))
        elif bond['bond']['bondTypes'][0] in allowed_bond_types:
            pnp_edges.append((s, e))
        else:
            other_edges[s].append(e)
            other_edges[e].append(s)

    def dfs(v, e, mask):
        mask[v] = True
        for v_ in e[v]:
            if not mask[v_]:
                mask = dfs(v_,e,mask)
        return mask

    mask = [False] * len(nodes)
    node_to_component = [0] * len(nodes)
    components = []
    for i in range(len(nodes)):
        if not mask[i]:
            mask_ = mask[:]
            mask = dfs(i, other_edges, mask)
            comp = []
            for i in range(len(nodes)):
                if mask[i] and not mask_[i]:
                    comp.append(i)
                    node_to_component[i] = len(components)
            components.append(comp)

    registered_set = set(nerpa_monomers)
    def get_comp_name(comp):
        namelist = sorted([(nodes[n], n) for n in comp])
        for _name, idx in namelist:
            if _name in registered_set:
                name = _name
                if is_d[idx] is not None:
                    name = f'@D-{name}' if is_d[idx] else f'@L-{name}'
                if is_pk[idx]:
                    name = f'@PK-{name}'
                if len(comp) > 1:
                    name = f'*{name}'
                break
        else:
            cat_name = "|".join([n for n, _ in namelist if not n.startswith("X")])
            name = f'<{cat_name}>' if cat_name else na
        return name

    new_nodes = [get_comp_name(comp) for comp in components]
    num_recognized_nodes = sum(n != na and not n.startswith('<') for n in new_nodes)
    new_edges = [(node_to_component[s], node_to_component[e]) for s,e in pnp_edges]
    new_oedges = [(node_to_component[s], node_to_component[e]) for s,e in oedges]

    if num_recognized_nodes < min_recognized_nodes:
        raise NumNodesError(rban_record['id'], num_recognized_nodes)

    result = ','.join(new_nodes) + ';' \
             + ';'.join([f'{s},{e}' for s,e in new_edges]) + ";" \
             + ';'.join([f'{s},o,{e}' for s,e in new_oedges])
    return result


def rban_postprocessing(path_to_rban_output, main_out_dir, path_to_rban):

    # check all unrecognized monomers for PK involvement
    new_monomers = []
    with open(path_to_rban_output) as f_in:
        for struct_id, rban_record in enumerate(json.load(f_in)):
            for monomer_id, monomer in enumerate(rban_record["monomericGraph"]["monomericGraph"]['monomers']):
                if monomer['monomer']['monomer']['monomer'].startswith('X'):
                    smi = monomer['monomer']['monomer']['smiles']
                    try:
                        aa_smi, pk_smi, _ = handle_DL.split_aa_pk_hybrid(smi)
                        new_id = f'{struct_id}_{monomer_id}'
                        new_monomers.append((aa_smi, new_id))
                    except handle_DL.PKError:
                        pass
    new_rban_input = os.path.join(main_out_dir, 'rban-putative-hybrids.input.json')
    generate_rban_input_from_list(new_monomers, new_rban_input)
    new_rban_output = os.path.join(main_out_dir, 'rban-putative-hybrids.output.json')
    try:
        run_rban(path_to_rban, new_rban_input, new_rban_output, main_out_dir)
    except subprocess.CalledProcessError as e:
        raise e

    hybrid_monomers_dict = defaultdict(dict)
    new_monomers_processed = json.loads(open(new_rban_output).read())
    for rban_record in new_monomers_processed:
        struct_id, monomer_id = map(int, rban_record['id'].split('_'))
        aa_smi = rban_record["isomericSmiles"]
        aa_code = rban_record["monomericGraph"]["monomericGraph"]['monomers'][0]['monomer']['monomer']['monomer']
        if not aa_code.startswith('X'):
            hybrid_monomers_dict[struct_id][monomer_id] = aa_code

    return hybrid_monomers_dict


def generate_graphs_from_rban_output(path_to_rban_output, path_to_monomers_tsv, path_to_graphs, main_out_dir, path_to_rban):
    recognized_monomers = [x.split()[0] for x in open(path_to_monomers_tsv)]
    recognized_monomers = recognized_monomers[1:]
    hybrid_monomers_dict = rban_postprocessing(path_to_rban_output, main_out_dir, path_to_rban)

    with open(path_to_graphs, 'w') as f_out:
        with open(path_to_rban_output) as f_in:
            for i, rban_record in enumerate(json.load(f_in)):
                try:
                    graph = process_single_record(rban_record, recognized_monomers, PNP_BONDS, hybrid_monomers_dict[i],
                                                  UNDEFINED_NAME, min_recognized_nodes=1)
                    f_out.write(f'{rban_record["id"]} {graph}\n')
                except NumNodesError as e:
                    log.warn(e)


def generate_rban_input_from_smiles_string(smiles, path_to_rban_input):
    with open(path_to_rban_input, 'w') as f:
        json.dump([{'id': i, 'smiles': smi.strip()} for i, smi in enumerate(smiles)], f, indent=2)


def generate_rban_input_from_smiles_tsv(path_to_csv, path_to_rban_input, sep='\t',
                                 id_col_name=None, smi_col_name='SMILES'):
    result = []
    with open(path_to_csv, newline='') as f_in:
        reader = csv.DictReader(f_in, delimiter=sep, quoting=csv.QUOTE_NONE)
        for i, row in enumerate(reader, 1):
            idx = str(i) if id_col_name is None else row[id_col_name]
            smi = row[smi_col_name]
            result.append({'id': idx, 'smiles': smi.strip()})

    with open(path_to_rban_input, 'w') as f:
        json.dump(result, f, indent=2)


def generate_rban_input_from_list(lst, path_to_rban_input):
    """
    :param lst: list of tuples (smiles, id)
    :param path_to_rban_input: filename
    :return:
    """
    dicts = [{'id': idx, 'smiles': smi.strip()} for smi, idx in lst]
    with open(path_to_rban_input, 'w') as f:
        json.dump(dicts, f, indent=2)


def run_rban(path_to_rban, path_to_rban_input, path_to_rban_output, main_out_dir, log=None):
    """
    raises: subprocess.CalledProcessError

    :param path_to_rban:
    :param path_to_rban_input:
    :param path_to_rban_output:
    :param main_out_dir:
    :param log:
    :return:
    """
    command = ['java', '-jar', path_to_rban,
               # '-monomersDB', '',
               '-inputFile', path_to_rban_input,
               '-outputFolder', main_out_dir,
               '-outputFileName', os.path.basename(path_to_rban_output)]
    p = subprocess.run(command, text=True, capture_output=True,
                       # check=True
                       )
    if p.stderr and log:
        log.err(p.stderr)
    for line in p.stdout.split('\n'):
        if line.startswith('WARNING') and log:
            log.warn('rBAN ' + line)