import os
import json
import csv
from collections import defaultdict
from itertools import permutations
import networkx as nx

import handle_monomers
import nerpa_utils


# TODO: move this section to some sort of config file
# see rban.src.main.java.molecules.bond.BondType for the full list
PNP_BONDS = ['AMINO', 'AMINO_TRANS', 'OXAZOLE_CYCLE', 'THIAZOLE_CYCLE', 'PYRIMIDINE_CYCLE']
UNDEFINED_NAME = 'NA'


class NumNodesError(Exception):
    def __init__(self, identifier, coverage):
        self.identifier = identifier
        self.coverage = coverage

    def __str__(self):
        return f'Poor monomer coverage for id="{self.identifier}": only {self.coverage} monomers identified. ' \
               f'Skipping.'


def build_nx_graph(rban_record, backbone_bonds, recognized_monomers, cut_lipids=True):
    monomeric_graph = rban_record['monomericGraph']['monomericGraph']

    nodes = []
    for monomer in monomeric_graph['monomers']:
        idx = monomer['monomer']['index']
        name =  monomer['monomer']['monomer']['monomer']
        attr = {
            'name': name,
            'isIdentified': name in recognized_monomers
            # 'isIdentified': monomer['monomer']['monomer']['isIdentified']
        }
        nodes.append((idx, attr))

    lipid_nodes = set()
    if cut_lipids:
        lipid_nodes = set(i for i, attr in nodes if ':' in attr['name'])

    G = nx.DiGraph()
    G.add_nodes_from(nodes)
    for bond in monomeric_graph['bonds']:
        # rBAN: N -> C
        s, e = bond['bond']['monomers']

        if s in lipid_nodes or e in lipid_nodes:
            continue

        bond_type = bond['bond']['bondTypes'][0]
        if bond_type in backbone_bonds:
            # BGC: upstream_C -> N_downstream
            G.add_edge(e, s, type=bond_type)

    return G


def hamiltonian_path(G, s, skip=None):
    node_to_id = {n:i for i, n in enumerate(G.nodes())}
    if skip is None:
        skip = []
    mask = [idx in skip for idx in G.nodes]
    path = []

    def dfs(G, v):
        mask[node_to_id[v]] = True
        path.append(v)
        if all(mask):
            return True
        for u in G.successors(v):
            if not mask[node_to_id[u]]:
                if dfs(G, u):
                    return True
        del path[-1]
        mask[node_to_id[v]] = False

    if dfs(G, s):
        return path


def putative_backbones(G, min_nodes=2):
    singular_nodes = []
    cycles = []
    paths = []
    for component in nx.connected_components(nx.Graph(G)):
        if len(component) < min_nodes:
            singular_nodes.append(component)
            continue

        Gs = G.subgraph(component)
        try:
            ced = nx.find_cycle(Gs)
            if set(ced) == set(Gs.edges()):
                cycles.append(list(nx.simple_cycles(Gs))[0])
                continue
        except:
            pass

        sources = []
        sinks = []
        for node in Gs.nodes():
            if Gs.in_degree(node) == 0 and Gs.out_degree(node) > 0:
                sources.append(node)
            elif Gs.in_degree(node) > 0 and Gs.out_degree(node) == 0:
                sinks.append(node)

        if len(sources) > 1 or len(sinks) > 1:
            raise NotImplementedError

        if len(sources) == 1:
            path = hamiltonian_path(Gs, sources[0])
            if path:
                paths.append(path)
                continue

        if len(sinks) == 1:
            path = hamiltonian_path(Gs.reverse(), sinks[0])
            if path:
                paths.append(path[::-1])
                continue

        raise NotImplementedError

    return cycles, paths, singular_nodes


def process_single_record(log, rban_record, recognized_monomers, backbone_bond_types,
                          hybrid_monomers_dict, na=UNDEFINED_NAME, min_recognized_nodes=2):
    G = build_nx_graph(rban_record, backbone_bond_types, recognized_monomers)
    structure_id = rban_record['id']

    try:
        chirality = handle_monomers.get_monomers_chirality(rban_record)
        for i, d in chirality.items():
            G.nodes[i]['isD'] = d
    except Exception as e:
        log.warning(f'Unable to get isomeric information for {rban_record["id"]}')

    for i, name in hybrid_monomers_dict.items():
        G.nodes[i]['name'] = name
        G.nodes[i]['isPKHybrid'] = True
        G.nodes[i]['isIdentified'] = True
    try:
        cycles, paths, _ = putative_backbones(G, min_nodes=2)
    except Exception as e:
        log.warning(f'Unable to linearize {rban_record["id"]}')
        return []

    def gen_nerpa_input(path, cyclic=False):
        inp_nodes = []
        for node in path:
            attr = G.nodes[node]
            name = attr['name']
            # if not attr.get('isIdentified', False):
            #     name = f'<{name}>'
            is_d = attr.get('isD', None)
            if is_d is not None:
                name = f'@D-{name}' if is_d else f'@L-{name}'
            if attr.get('isPKHybrid', False):
                name = f'*{name}'
            inp_nodes.append(name)
        inp_edges = [f'{i},{i+1}' for i in range(len(inp_nodes)-1)]
        if cyclic:
            inp_edges.append(f'{len(inp_nodes)-1},{0}')
        res = ','.join(inp_nodes)
        res += ';'
        res += ';'.join(inp_edges)
        return res

    def sufficiently_covered(path):
        n = sum(G.nodes[node]['isIdentified'] for node in path)
        return n >= min_recognized_nodes

    structures = {}
    for path in cycles:
        if sufficiently_covered(path):
            gr = gen_nerpa_input(path, cyclic=True)
            structures[gr] = path
    for path in paths:
        if sufficiently_covered(path):
            gr = gen_nerpa_input(path)
            structures[gr] = path
    if 1 < len(paths) < 4:
        for x in permutations(paths):
            path = [node for comp in x for node in comp]
            if sufficiently_covered(path):
                gr = gen_nerpa_input(path)
                structures[gr] = path
    structures = [(f'{structure_id}_variant{i}', gr, structures[gr]) for i, gr in enumerate(sorted(structures.keys()))]

    return structures


def rban_postprocessing(path_to_rban_output, main_out_dir, path_to_rban, log):

    # check all unrecognized monomers for PK involvement
    new_monomers = []
    with open(path_to_rban_output) as f_in:
        for struct_id, rban_record in enumerate(json.load(f_in)):
            for monomer in rban_record["monomericGraph"]["monomericGraph"]['monomers']:
                if monomer['monomer']['monomer']['monomer'].startswith('X'):
                    smi = monomer['monomer']['monomer']['smiles']
                    monomer_id = monomer['monomer']['index']
                    try:
                        aa_smi, pk_smi, _ = handle_monomers.split_aa_pk_hybrid(smi)
                        new_id = f'{struct_id}_{monomer_id}'
                        new_monomers.append((aa_smi, new_id))
                    except handle_monomers.PKError:
                        pass
    new_rban_input = os.path.join(main_out_dir, 'rban-putative-hybrids.input.json')
    generate_rban_input_from_list(new_monomers, new_rban_input)
    new_rban_output = os.path.join(main_out_dir, 'rban-putative-hybrids.output.json')
    run_rban(path_to_rban, new_rban_input, new_rban_output, main_out_dir, log)

    hybrid_monomers_dict = defaultdict(dict)
    new_monomers_processed = json.loads(open(new_rban_output).read())
    for rban_record in new_monomers_processed:
        struct_id, monomer_id = map(int, rban_record['id'].split('_'))
        aa_smi = rban_record["isomericSmiles"]
        aa_code = rban_record["monomericGraph"]["monomericGraph"]['monomers'][0]['monomer']['monomer']['monomer']
        if not aa_code.startswith('X'):
            hybrid_monomers_dict[struct_id][monomer_id] = aa_code

    return hybrid_monomers_dict


def generate_graphs_from_rban_output(path_to_rban_output, path_to_monomers_tsv, path_to_graphs,
                                     main_out_dir, path_to_rban, log):
    recognized_monomers = [x.split()[0] for x in open(path_to_monomers_tsv)]
    recognized_monomers = set(recognized_monomers[1:])
    hybrid_monomers_dict = rban_postprocessing(path_to_rban_output, main_out_dir, path_to_rban, log)

    with open(path_to_graphs, 'w') as f_out:
        with open(path_to_rban_output) as f_in:
            for i, rban_record in enumerate(json.load(f_in)):
                try:
                    graphs = process_single_record(log, rban_record, recognized_monomers, PNP_BONDS, hybrid_monomers_dict[i],
                                                   UNDEFINED_NAME, min_recognized_nodes=2)
                    for i, gr, pt in graphs:
                        f_out.write(f'{i} {gr} {pt}\n')
                except NumNodesError as e:
                    log.warning(e)


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


def run_rban(path_to_rban, path_to_rban_input, path_to_rban_output, main_out_dir, log):
    """
    :param path_to_rban:
    :param path_to_rban_input:
    :param path_to_rban_output:
    :param main_out_dir:
    :return:
    """
    command = ['java', '-jar', path_to_rban,
               # '-monomersDB', '',
               '-inputFile', path_to_rban_input,
               '-outputFolder', main_out_dir + '/',  # rBAN specifics
               '-outputFileName', os.path.basename(path_to_rban_output)]
    nerpa_utils.sys_call(command, log)
