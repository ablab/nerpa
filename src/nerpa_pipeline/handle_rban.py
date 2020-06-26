import sys
import os
import json
import csv
from collections import defaultdict
from logger import log
try:
    import handle_DL
    USE_DL=True
except ImportError:
    log.warn(f'Unable to import RDKit. Stereoisometry information will not be used.')
    USE_DL=False


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


def process_single_record(rban_record, nerpa_monomers, allowed_bond_types, na=UNDEFINED_NAME, min_recognized_nodes=3):
    g = rban_record['monomericGraph']['monomericGraph']
    nodes = [m['monomer']['monomer']['monomer'] for m in g['monomers']]
    if USE_DL:
        try:
            is_d = handle_DL.get_monomers_chirality(rban_record)
        except Exception as e:
            log.warn(f'Unable to get isomeric information for {rban_record["id"]} : {e}')
            is_d = [None] * len(nodes)
    else:
        is_d = [None] * len(nodes)

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
        name = na
        for _name, idx in sorted([(nodes[n], n) for n in comp]):
            if _name in registered_set:
                name = _name
                if is_d[idx] is not None:
                    name = f'@D-{name}' if is_d[idx] else f'@L-{name}'
                if len(comp) > 1:
                    name = f'*{name}'
                break
        return name

    new_nodes = [get_comp_name(comp) for comp in components]
    num_recognized_nodes = sum(n != na for n in new_nodes)
    new_edges = [(node_to_component[s], node_to_component[e]) for s,e in pnp_edges]
    new_oedges = [(node_to_component[s], node_to_component[e]) for s,e in oedges]

    if num_recognized_nodes < min_recognized_nodes:
        raise NumNodesError(rban_record['id'], num_recognized_nodes)

    result = ','.join(new_nodes) + ';' \
             + ';'.join([f'{s},{e}' for s,e in new_edges]) + ";" \
             + ';'.join([f'{s},o,{e}' for s,e in new_oedges])
    return result


def generate_graphs_from_rban_output(path_to_rban_output, path_to_monomers_tsv, path_to_graphs):
    recognized_monomers = [x.split()[0] for x in open(path_to_monomers_tsv)]
    recognized_monomers = recognized_monomers[1:]
    with open(path_to_graphs, 'w') as f_out:
        with open(path_to_rban_output) as f_in:
            for i, rban_record in enumerate(json.load(f_in)):
                try:
                    graph = process_single_record(rban_record, recognized_monomers, PNP_BONDS, UNDEFINED_NAME)
                    f_out.write(f'{rban_record["id"]} {graph}\n')
                except NumNodesError as e:
                    log.warn(e)


def generate_rban_input_from_smiles_string(smiles, path_to_rban_input):
    with open(path_to_rban_input, 'w') as f:
        json.dump([{'id': i, 'smiles': smi} for i, smi in enumerate(smiles)], f, indent=2)


def generate_rban_input_from_smiles_tsv(path_to_csv, path_to_rban_input, sep='\t',
                                 id_col_name=None, smi_col_name='SMILES'):
    result = []
    with open(path_to_csv, newline='') as f_in:
        reader = csv.DictReader(f_in, delimiter=sep, quoting=csv.QUOTE_NONE)
        for i, row in enumerate(reader, 1):
            idx = str(i) if id_col_name is None else row[id_col_name]
            smi = row[smi_col_name]
            result.append({'id': idx, 'smiles': smi})

    with open(path_to_rban_input, 'w') as f:
        json.dump(result, f, indent=2)