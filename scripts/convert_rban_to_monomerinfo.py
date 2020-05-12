import json
import sys
import os
from collections import defaultdict


# see rban.src.main.java.molecules.bond.BondType for the full list
PNP_BONDS = ['AMINO', 'AMINO_TRANS', 'ESTER', 'THIOETHER', 'OXAZOLE_CYCLE', 'THIAZOLE_CYCLE']
UNDEFINED_NAME = 'NA'
MONOMERS_TSV = '//Data/Projects/CAB/Nerpa/nerpa/resources/monomers.tsv'


def monomericgraph_to_string(g):
    monomers = [(m['monomer']['index'], m['monomer']['monomer']['monomer']) for m in g['monomers']]
    monomers = [m[1] for m in sorted(monomers)]
    bonds = [b['bond']['monomers'] for b in g['bonds']]
    result = ','.join(monomers) + ';' + ';'.join([f'{b[0]-1},{b[1]-1}' for b in bonds])
    return result


def monomericgraph_components_to_string(g, nerpa_monomers, allowed_bond_types, na='NA'):

    nodes = [m['monomer']['monomer']['monomer'] for m in g['monomers']]
    pnp_edges = []
    other_edges = defaultdict(list)
    for bond in g['bonds']:
        s, e = bond['bond']['monomers']
        s -= 1
        e -= 1
        if bond['bond']['bondTypes'][0] in allowed_bond_types:
            pnp_edges.append((s, e))
        else:
            other_edges[s].append(e)
            other_edges[e].append(s)

    def dfs(v,e,mask):
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
            mask = dfs(i,other_edges,mask)
            comp = []
            for i in range(len(nodes)):
                if mask[i] and not mask_[i]:
                    comp.append(i)
                    node_to_component[i] = len(components)
            components.append(comp)

    #     print(components)
    #     print(node_to_component)
    #     for s, e in pnp_edges:
    #         print(node_to_component[s], node_to_component[e])

    registered_set = set(nerpa_monomers)
    def get_comp_name(comp):
        names = sorted(list(set(nodes[i] for i in comp) & registered_set))
        name = names[0] if names else na
        if len(comp) > 1 and names:
            name = f'*{name}'
        return name

    new_nodes = [get_comp_name(comp) for comp in components]
    new_edges = [(node_to_component[s], node_to_component[e]) for s,e in pnp_edges]

    result = ','.join(new_nodes) + ';' + ';'.join([f'{s},{e}' for s,e in new_edges])
    return result


def main():
    graphs = json.load(open(sys.argv[1]))
    rban_input_json = json.load(open(sys.argv[2])) if len(sys.argv) > 2 else None

    monomers_tsv = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../resources/monomers.tsv')
    if not os.path.exists(monomers_tsv):
        monomers_tsv = MONOMERS_TSV

    nerpa_monomers = [x.split()[0] for x in open(monomers_tsv)]
    nerpa_monomers = nerpa_monomers[1:]

    out = []
    for i, g in enumerate(graphs):
        # gr = monomericgraph_to_string(g['monomericGraph']['monomericGraph'])
        gr = monomericgraph_components_to_string(g['monomericGraph']['monomericGraph'],
                                                 nerpa_monomers, PNP_BONDS, UNDEFINED_NAME)
        out.append(f'{g["id"]} {gr}')
        if rban_input_json and 'extra' in rban_input_json[i].keys():
            out[-1] += f' {rban_input_json[i]["extra"]}'

    print('\n'.join(out))


if __name__ == '__main__':
    main()