import os
import json
import csv
from collections import defaultdict
from itertools import permutations
from rdkit.Chem import Descriptors
from rdkit.Chem import MolFromSmiles
import networkx as nx

import handle_monomers
import nerpa_utils


# TODO: move this section to some sort of config file
# see rban.src.main.java.molecules.bond.BondType for the full list
PNP_BONDS = ['AMINO', 'AMINO_TRANS', 'OXAZOLE_CYCLE', 'THIAZOLE_CYCLE', 'PYRIMIDINE_CYCLE', 'HETEROCYCLE']
UNDEFINED_NAME = 'NA'


class NumNodesError(Exception):
    def __init__(self, identifier, coverage):
        self.identifier = identifier
        self.coverage = coverage

    def __str__(self):
        return f'Poor monomer coverage for id="{self.identifier}": only {self.coverage} monomers identified. ' \
               f'Skipping.'


def build_nx_graph(rban_record, backbone_bonds, recognized_monomers, cut_lipids=True):
    '''
    This function creates a networkx graph from rban_record (monomer graph in json format)
    All non-backbone edges are removed.
     With cut_lipids flag all edges incident to the lipid node are also removed
    '''
    monomeric_graph = rban_record['monomericGraph']['monomericGraph']

    nodes = []
    for monomer in monomeric_graph['monomers']:
        idx = monomer['monomer']['index']
        name =  monomer['monomer']['monomer']['monomer']
        attr = {
            'name': name.replace('C10:0-NH2(2)-Ep(9)-oxo(8)', 'Aeo'),
            'isIdentified': name in recognized_monomers
            # 'isIdentified': monomer['monomer']['monomer']['isIdentified']
        }
        nodes.append((idx, attr))
    '''the list of nodes is created'''

    lipid_nodes = set()
    if cut_lipids:
        lipid_nodes = set(i for i, attr in nodes if ':' in attr['name'])

    G = nx.DiGraph()
    G.add_nodes_from(nodes)
    for bond in monomeric_graph['bonds']:
        # rBAN: N -> C
        s, e = bond['bond']['monomers']

        if s in lipid_nodes or e in lipid_nodes:  # cut edges incident to lipid nodes
            continue

        bond_type = bond['bond']['bondTypes'][0]
        if bond_type in backbone_bonds:  # keeps only backbone edges
            # BGC: upstream_C -> N_downstream
            '''ok, so in rBAN the edges are oriented in the opposite way. 
            Does that hold for all types of bonds? I hope so'''
            G.add_edge(e, s, type=bond_type)

    return G
def my_build_nx_graph(rban_record, recognized_monomers):
    '''
    The same a build_nx_graph, but does not cut any edges
    and keeps in nodes 'smiles' and 'mass'
    '''
    monomeric_graph = rban_record['monomericGraph']['monomericGraph']

    nodes = []
    for monomer in monomeric_graph['monomers']:
        idx = monomer['monomer']['index']
        name =  monomer['monomer']['monomer']['monomer']
        # TODO this is mass of the WHOLE aminoacid, the mass of monomer can differ
        mass = Descriptors.ExactMolWt(MolFromSmiles(monomer['monomer']['monomer']['smiles']))
        attr = {
            'name': name.replace('C10:0-NH2(2)-Ep(9)-oxo(8)', 'Aeo'),
            'isIdentified': name in recognized_monomers,
            'mass': mass,
            'smiles': monomer['monomer']['monomer']['smiles']
            # 'isIdentified': monomer['monomer']['monomer']['isIdentified']
        }
        nodes.append((idx, attr))

    G = nx.DiGraph()
    G.add_nodes_from(nodes)
    for bond in monomeric_graph['bonds']:
        # rBAN: N -> C
        s, e = bond['bond']['monomers']
        bond_type = bond['bond']['bondTypes'][0]
        G.add_edge(e, s, type=bond_type)

    return G

def hamiltonian_path(G, s, skip=None):
    '''
    This function just tries a single traversal of the graph, beginning at s.
    If it does not succeed, it reports a failure
    I feel that the same could be done with just a few lines of code...
    '''
    node_to_id = {n:i for i, n in enumerate(G.nodes())}
    if skip is None:  # why not just set skip=[] be default? and why do we need skip in the first place?
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
    '''
    Split the given graph G into linear and cyclic fragments
    '''
    singular_nodes = []
    cycles = []
    paths = []
    '''
    The graph is, in general, not connected because of the removal non-backbone edges in the build_nx_graph function
    '''
    for component in nx.connected_components(nx.Graph(G)):
        if len(component) < min_nodes:
            singular_nodes.append(component)
            continue

        Gs = G.subgraph(component)
        '''
         try to parse the component as a single simple cycle
        '''
        try:
            ced = nx.find_cycle(Gs)
            if set(ced) == set(Gs.edges()):
                cycles.append(list(nx.simple_cycles(Gs))[0])
                continue
        except:
            pass

        '''
        try to parse the component as a chupa-choops graph (a path attached to a simple cycle)
        '''
        sources = []
        sinks = []
        for node in Gs.nodes():
            if Gs.in_degree(node) == 0 and Gs.out_degree(node) > 0:
                sources.append(node)
            elif Gs.in_degree(node) > 0 and Gs.out_degree(node) == 0:
                sinks.append(node)

        if len(sources) > 1 or len(sinks) > 1:  # graph is definitely not a chupa-choops
            raise NotImplementedError

        '''try to traverse graph in forward direction'''
        if len(sources) == 1:
            '''
            hamiltonan_path is a simple dfs. It works with chupa-choops graphs 
            and may accidentally work with some other graphs
            '''
            path = hamiltonian_path(Gs, sources[0])
            if path:
                paths.append(path)
                continue

        '''try to traverse graph in reverse direction'''
        if len(sinks) == 1:
            path = hamiltonian_path(Gs.reverse(), sinks[0])
            if path:
                paths.append(path[::-1])
                continue

        raise NotImplementedError  # hamiltonian path not found

    return cycles, paths, singular_nodes  # singular nodes are never used!


def write_initial_graph(G, rban_record, struct_id, main_out_dir):
    '''
   Write a monomer graph in json format to the 'main_out_dir/initial_monomer_graphs' folder
    '''
    initial_graphs_dir = os.path.join(main_out_dir, 'initial_monomer_graphs')
    if not os.path.exists(initial_graphs_dir):
        os.mkdir(initial_graphs_dir)
    out_file_rban_postprocessed = os.path.join(initial_graphs_dir, struct_id + '.json')
    with open(out_file_rban_postprocessed, 'w') as out:
        contents = {'nodes': G.nodes._nodes,
                    'edges': G.adj._atlas}
        json.dump(contents, out)


def write_rban_record(rban_record, main_out_dir):
    initial_graphs_dir = os.path.join(main_out_dir, 'initial_monomer_graphs')
    if not os.path.exists(initial_graphs_dir):
        os.mkdir(initial_graphs_dir)
    struct_id = rban_record['id']
    out_file_rban_initial = os.path.join(initial_graphs_dir, struct_id + '_rban.json')
    with open(out_file_rban_initial, 'w') as out:
        json.dump(rban_record, out)


def process_single_record(log, rban_record, recognized_monomers, backbone_bond_types,
                          hybrid_monomers_dict, main_out_dir, na=UNDEFINED_NAME, min_recognized_nodes=2):
    '''
    build_nx_graphs creates a networkx DiGraph from a monomer graph in rBAN format.
    In doing so all edges which are not of backbone_bond_types are removed.
    Also, by default all edges incident to lipid monomers are removed
    '''

    def add_chirality(graph):
        try:
            chirality = handle_monomers.get_monomers_chirality(rban_record)
            for i, d in chirality.items():
                graph.nodes[i]['isD'] = d
        except Exception as e:
            log.warning(f'Structure "{rban_record["id"]}": unexpected error while determining stereoisomeric configuration '
                        f'of monomers. Stereoisomeric configuration will be omitted.')

    def add_hybrid_monomers(graph):
        for i, name in hybrid_monomers_dict.items():
            G.nodes[i]['name'] = name
            G.nodes[i]['isPKHybrid'] = True
            G.nodes[i]['isIdentified'] = True

    '''
    Try to recognize monomers chirality and
    add newly identified monomers from hybrid monomers
    '''

    structure_id = rban_record['id']
    print(f'Processing {structure_id}')
    G = build_nx_graph(rban_record, backbone_bond_types, recognized_monomers)
    add_chirality(G)
    add_hybrid_monomers(G)

    write_rban_record(rban_record, main_out_dir)  # output rban_record for Louis feature
    '''
    # I used to need these nerpa processed monomer graphs at some point but now I work directly from rBAN output
    # Nevertheless, sometime in the future I might need this newly identified monomers or chirality 
    try:
        # My version of build_nx_graph, which does not remove any edges and keeps some additional information in the nodes
        my_G = my_build_nx_graph(rban_record, recognized_monomers)
        add_chirality(my_G)
        add_hybrid_monomers(my_G)
        # I output the initial monomer graph in networkx format to use it later in Louis feature
        write_initial_graph(my_G, rban_record, structure_id, main_out_dir)
    except:
        # sometimes Chem fails to parse individual monomers and calculate their weight
        pass
    '''
    '''
    Split the graph into paths and simple cycles
    '''
    try:
        cycles, paths, _ = putative_backbones(G, min_nodes=2)
        if not cycles and not paths:
            raise
    except Exception as e:
        log.warning(f'Structure "{rban_record["id"]}": unable to determine backbone sequence. '
                    f'Skipping "{rban_record["id"]}".')
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

    '''create a dictionary with records of all fragments
    structures.keys() are the str records of fragments in Nerpa format
    structures.values() are these fragments in the original graph (as List[int])'''
    structures = {}
    for path in cycles:
        if sufficiently_covered(path):  # path/cycle contains enough recognized monomers
            gr = gen_nerpa_input(path, cyclic=True)
            structures[gr] = path
    for path in paths:
        if sufficiently_covered(path):
            gr = gen_nerpa_input(path)
            structures[gr] = path
    if 1 < len(paths) < 4:
        # try to join linear fragments in all possible ways
        for x in permutations(paths):
            path = [node for comp in x for node in comp]
            if sufficiently_covered(path):  # kind of strange to check it inside the cycle, as 'coverage' does not depend on the permutation, but here we are
                gr = gen_nerpa_input(path)
                structures[gr] = path
    structures = [(f'{structure_id}_variant{i}', gr, structures[gr]) for i, gr in enumerate(sorted(structures.keys()))]

    return structures


def rban_postprocessing(path_to_rban_output, main_out_dir, path_to_rban, path_to_monomers_db, log):
    # check all unrecognized monomers for PK involvement
    '''
    This function tries to remove a lipid tail from each of the unrecognized monomers and
    then runs rBAN on these individual trimmed monomers. If rBAN recognizes a trimmed monomer,
    this information is added to the hybrid_monomers_dict
    '''
    new_monomers = []
    with open(path_to_rban_output) as f_in:
        for struct_id, rban_record in enumerate(json.load(f_in)):
            for monomer in rban_record["monomericGraph"]["monomericGraph"]['monomers']:
                if monomer['monomer']['monomer']['monomer'].startswith('X'):  # monomer was not recognized
                    smi = monomer['monomer']['monomer']['smiles']
                    monomer_id = monomer['monomer']['index']
                    try:
                        aa_smi, pk_smi, _ = handle_monomers.split_aa_pk_hybrid(smi)
                        new_id = f'{struct_id}_{monomer_id}'
                        new_monomers.append((aa_smi, new_id))
                    except handle_monomers.PKError:
                        pass # it's okay
                    except:
                        log.warn(f'\nStructure "{struct_id}": unexpected error while resolving NRP-PK hybrid monomer '
                                 f'candidate. Skipping.')

    '''
    new_monomers is the list of trimmed monomers
    '''
    if not new_monomers:
        return defaultdict(dict)

    '''run rBAN on the trimmed monomers'''
    log.info('\n=== Resolving NRP-PK hybrid monomers candidates')
    new_rban_input = os.path.join(main_out_dir, 'rban-putative-hybrids.input.json')
    generate_rban_input_from_list(new_monomers, new_rban_input)
    new_rban_output = os.path.join(main_out_dir, 'rban-putative-hybrids.output.json')
    run_rban(path_to_rban, new_rban_input, new_rban_output, path_to_monomers_db, main_out_dir, log)

    '''
    Add newly recognized monomers to the dictionaty
    '''
    hybrid_monomers_dict = defaultdict(dict)
    new_monomers_processed = json.loads(open(new_rban_output).read())
    for rban_record in new_monomers_processed:
        struct_id, monomer_id = map(int, rban_record['id'].rsplit('_', 1))
        aa_smi = rban_record["isomericSmiles"]
        aa_code = rban_record["monomericGraph"]["monomericGraph"]['monomers'][0]['monomer']['monomer']['monomer']
        if not aa_code.startswith('X'):  # rBAN managed to recognize the trimmed monomer
            hybrid_monomers_dict[struct_id][monomer_id] = aa_code

    return hybrid_monomers_dict


'''
This was my function to retrieve monomer masses, but it is no longer needed

def write_masses(rban_record, main_out_dir):
    struct_id = rban_record['id']
    out_dir = os.path.join(main_out_dir, 'structures_masses')
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    filename = os.path.join(out_dir, struct_id)
    file_corrupted = False
    with open(filename, 'w') as out:
        records = ['index', 'monomer', 'smiles', 'mass_heavy_atoms', 'mass_full']
        out.write(','.join(records) + '\n')
        for monomer in rban_record['monomericGraph']['monomericGraph']['monomers']:
            m = monomer['monomer']
            try:
                to_write = ','.join(map(str, [m['index'],
                                              m['monomer']['monomer'],
                                              m['monomer']['smiles'],
                                              m['monomer']['mwHeavyAtoms'],
                                              Descriptors.ExactMolWt(MolFromSmiles(m['monomer']['smiles']))]))
            except:
                file_corrupted = True
                print(f'Exception while handling record: {rban_record}')
            else:
                out.write(to_write + '\n')
    if file_corrupted:
        os.remove(filename)
'''

def generate_info_from_rban_output(path_to_rban_output, path_to_monomers_tsv, path_to_graphs,
                                   main_out_dir, path_to_rban, path_to_monomers_db, log, process_hybrids=False):
    """
    path_to_rban_output --- output of the rBAN ('rban.output.json'). Used as the main argument for this function
    path_to_monomers_tsv --- the file with the list of monomers supported by Nerpa
    path_to_graphs --- output of this function (processed rBAN output) will be stored in this file ('structures.info')
    path_to_monomers_db --- path to the database of monomers (didn't yet quite understand this mess with monomer databases)
"""
    # takes rBAN output and 1. tries to annotate some unrecognized monomers, if process_hybrids=True 2. extracts monomer graphs from it, producing the file 'structures.info' (path_to_graphs)
    log.info('\n======= Processing rBAN output')
    log.info(f'results will be in {path_to_graphs}', indent=1)
    recognized_monomers = [x.split()[0] for x in open(path_to_monomers_tsv)]
    recognized_monomers = set(recognized_monomers[1:])  # reading the set of monomers supported by nerpa from the config file ([1:] because first row is annotations)
    hybrid_monomers_dict = defaultdict(dict)
    if process_hybrids:
        '''
        try to recognize some unidentified monomers by removing a lipid tail from them and then running rban again on the trimmed monomers
        the initial rban results are unaffected (and used further), only the dictionary of newly recognized monomers is created
        '''
        hybrid_monomers_dict = rban_postprocessing(path_to_rban_output, main_out_dir, path_to_rban, path_to_monomers_db, log)

    with open(path_to_graphs, 'w') as f_out:
        with open(path_to_rban_output) as f_in:
            for i, rban_record in enumerate(json.load(f_in)):
                try:
                    '''
                    process_single_record performs:
                    1. Naming newly identified monomers from hybrid_monomer_dict
                    2. Identifying chirality of each monomer (if possible)
                    3. Splitting graph into linear (or simple-cyclic) fragments
                    4. Translating each fragment into Nerpa format (in doing so original node indexes are lost)
                    '''
                    graphs = process_single_record(log, rban_record, recognized_monomers, PNP_BONDS, hybrid_monomers_dict[i],
                                                   main_out_dir, na=UNDEFINED_NAME, min_recognized_nodes=2)  # these are not graphs! As I understood, these are different linear representations of the same structure
                    '''
                    Each element of graphs is of form i, gr, pt, where:
                    'i' is the name of the fragment in the form structId_variantX
                    'gr' is the representation of the fragment in Nerpa format
                    'pt' is the fragment (path or cycle) in the original monomer graph (with original node indexes)
                    '''

                    '''
                    All edges of the original monomer graph. I don't know why those are needed
                    '''
                    rban_graph_edges = [
                        ','.join(map(str, b['bond']['monomers']))
                        for b in rban_record['monomericGraph']['monomericGraph']['bonds']
                    ]
                    rban_graph_edges_str = ';'.join(rban_graph_edges)
                    for i, gr, pt in graphs:
                        pt_str = ','.join(map(str, pt))
                        f_out.write(f'{i} {gr} {pt_str};{rban_graph_edges_str}\n')
                except NumNodesError as e:
                    log.warning(e)
    log.info('\n======= Done with Processing rBAN output')


def generate_rban_input_from_smiles_strings(smiles, path_to_rban_input):
    with open(path_to_rban_input, 'w') as f:
        json.dump([{'id': f'compound_{i:06d}', 'smiles': smi.strip()} for i, smi in enumerate(smiles)], f, indent=2)


def generate_rban_input_from_smiles_tsv(path_to_csv, path_to_rban_input, sep='\t',
                                        id_col_name=None, smi_col_name='SMILES'):
    result = []
    with open(path_to_csv, newline='') as f_in:
        reader = csv.DictReader(f_in, delimiter=sep, quoting=csv.QUOTE_NONE)
        for i, row in enumerate(reader, 1):
            idx = f'compound_{i:06d}' if id_col_name is None else row[id_col_name]
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


def run_rban(path_to_rban_jar, path_to_rban_input, path_to_rban_output, path_to_monomers, main_out_dir, log):
    """
    executes rBAN on the structures from 'rban.input.json' (path_to_rban_input). The results are written to 'rban.output.json' (path_to_rban_output)
    path_to_monomers -- the path to the monomer database
    :param path_to_rban_jar:
    :param path_to_rban_input:
    :param path_to_rban_output:
    :param main_out_dir:
    :return:
    """
    command = ['java', '-jar', path_to_rban_jar,
               '-inputFile', path_to_rban_input,
               '-outputFolder', main_out_dir + '/',  # rBAN specifics
               '-outputFileName', os.path.basename(path_to_rban_output)]
    if path_to_monomers:
        command += ['-monomersDB', path_to_monomers]
    nerpa_utils.sys_call(command, log)
