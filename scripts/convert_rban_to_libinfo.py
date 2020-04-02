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
    rban_input_json = json.load(open(sys.argv[2])) if len(sys.argv) > 2 else None

    out = []
    for i, g in enumerate(graphs):
        gr = monomericgraph_to_string(g['monomericGraph']['monomericGraph'])
        out.append(f'{g["id"]} {gr}')
        if rban_input_json and 'extra' in rban_input_json[i].keys():
            out[-1] += f' {rban_input_json[i]["extra"]}'

    print('\n'.join(out))


if __name__ == '__main__':
    main()