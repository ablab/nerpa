import json
import sys


def monomericgraph_to_string(g):
    monomers = [(m['monomer']['index'], m['monomer']['monomer']['monomer']) for m in g['monomers']]
    monomers = [m[1] for m in sorted(monomers)]
    bonds = [b['bond']['monomers'] for b in g['bonds']]
    result = ','.join(monomers) + ';' + ';'.join([f'{b[0]-1},{b[1]-1}' for b in bonds])
    return result


def main():
    graphs = json.load(open(sys.argv[1]))
    extras = None
    if len(sys.argv) > 2:
        extras = open(sys.argv[2]).read().strip().split()
    out = []
    for i, g in enumerate(graphs):
        gr = monomericgraph_to_string(g['monomericGraph']['monomericGraph'])
        extra = '' if extras is None else extras[i]
        idx = g['id']
        out.append(f'{idx} {gr} {extra}')

    print('\n'.join(out))


if __name__ == '__main__':
    main()