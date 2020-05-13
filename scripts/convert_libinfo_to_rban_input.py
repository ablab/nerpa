import os
import sys
import json
from rdkit import Chem as rdc


def main():
    filename = sys.argv[1]
    dbpath = sys.argv[2] if len(sys.argv) > 2 else None
    j = []
    for line in open(filename):
        if line.strip():
            spl = line.split()
            molfile = spl[0]
            if dbpath:
                molfile = os.path.join(dbpath, molfile)
            if not os.path.exists(molfile):
                raise FileNotFoundError(molfile)
            idx, _ = os.path.splitext(os.path.basename(molfile))
            idx = f'{idx}.gr'
            extra = ' '.join(spl[1:]) if len(spl) > 1 else None
            mol = rdc.MolFromMolBlock(open(molfile).read())
            if mol is None:
                mol = rdc.MolFromMolBlock('\n' + open(molfile).read())
            smi = rdc.MolToSmiles(mol)
            if not smi:
                raise ValueError(f'Failed to process {molfile}')
                # print(f'Failed to process {molfile}')
                # continue
            j.append({'id': idx, 'smiles': smi})
            if extra:
                j[-1]['extra'] = extra

    print(json.dumps(j, indent=2))


if __name__ == '__main__':
    main()