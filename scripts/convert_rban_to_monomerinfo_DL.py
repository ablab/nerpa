import json
import sys
import os
from collections import defaultdict
from rdkit import Chem as rdc
from collections import defaultdict
import json

# see rban.src.main.java.molecules.bond.BondType for the full list
PNP_BONDS = ['AMINO', 'AMINO_TRANS', 'ESTER', 'THIOETHER', 'OXAZOLE_CYCLE', 'THIAZOLE_CYCLE']
UNDEFINED_NAME = 'NA'


class NoChiralCenters(Exception):
    pass
class NoCarboxyl(Exception):
    pass
class WTF(Exception):
    pass


def get_aa_chirality(smi):
    mol = rdc.MolFromSmiles(smi)
    chiral_centers = rdc.FindMolChiralCenters(mol)
    if not chiral_centers:
        raise NoChiralCenters
    RCN_pattern = 'RCN' * 2
    NCR_pattern = 'NCR' * 2
    def _get_chirality(match, q):
        for idx, rs in chiral_centers:
            # Sadly, we cannot make use of rs here, see
            # https://chemistry.stackexchange.com/questions/90826/difference-between-ld-and-sr-in-naming
            if idx in match:
                ch_center = mol.GetAtomWithIdx(idx)
                neighbors = [a.GetIdx() for a in ch_center.GetNeighbors()]
                for match_cooh in mol.GetSubstructMatches(q):
                    neighbor_tags = []
                    for n in neighbors:
                        if n in match_cooh:
                            neighbor_tags.append('C')
                        elif mol.GetAtomWithIdx(n).GetSymbol() == 'N':
                            neighbor_tags.append('N')
                        else:
                            neighbor_tags.append('R')
                    neighbor_tags = ''.join(neighbor_tags)

                    # since glutamic acid has several COOH matches
                    if neighbor_tags in RCN_pattern or neighbor_tags in NCR_pattern:
                        if ch_center.GetChiralTag() == rdc.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
                            return neighbor_tags in RCN_pattern
                        else:
                            return neighbor_tags in NCR_pattern
        return None

    # try beta
    query_beta = rdc.MolFromSmarts('NCC([C](=O)O)')
    match_beta = mol.GetSubstructMatch(query_beta)
    res_beta = None
    if match_beta:
        q = rdc.MolFromSmarts('CC([C](=O)O)')
        res_beta = _get_chirality(match_beta, q)

    # try alpha
    query_alpha = rdc.MolFromSmarts('NC([C](=O)O)')
    match_alpha = mol.GetSubstructMatch(query_alpha)
    res_alpha = None
    if match_alpha:
        q = rdc.MolFromSmarts('C([C](=O)O)')
        res_alpha = _get_chirality(match_alpha, q)

    if not match_alpha and not match_beta:
        raise NoCarboxyl

    if res_alpha is not None:
        return res_alpha
    if res_beta is not None:
        return res_beta
    raise WTF


def get_monomer_bonds(mol, rban_record):
    atomic_bonds = rban_record["atomicGraph"]["atomicGraph"]['bonds']
    monomer_bonds_list = rban_record["monomericGraph"]["monomericGraph"]['bonds']
    monomer_bonds = []
    monomer_bond_types = []
    for bond_rec in monomer_bonds_list:
        for x in bond_rec['bond']['atomicIndexes']:
            b_ = mol.GetBondBetweenAtoms(*atomic_bonds[x]['atoms'])
            if b_:
                monomer_bonds.append(b_.GetIdx())
                monomer_bond_types.append(bond_rec['bond']["bondTypes"][0])
    return monomer_bonds, monomer_bond_types


def split_amino(rwmol, bonds, trans=False):
    for b_idx in sorted(bonds, reverse=True):
        # reverse because when a bond or atom is deleted, the bond or
        # atom indices are remapped. If you remove bonds with a higher index first,
        # bonds with lower indices will not be remapped.
        b = rwmol.GetBondWithIdx(b_idx)
        st_ = b.GetBeginAtomIdx()
        end_ = b.GetEndAtomIdx()
        rwmol.RemoveBond(st_, end_)
        new_bond_type = rdc.BondType.DOUBLE if trans else rdc.BondType.SINGLE
        if rwmol.GetAtomWithIdx(end_).GetSymbol() == 'C':
            new_id = rwmol.AddAtom(rdc.Atom(8))
            rwmol.AddBond(end_, new_id, new_bond_type)
        if rwmol.GetAtomWithIdx(st_).GetSymbol() == 'C':
            new_id = rwmol.AddAtom(rdc.Atom(8))
            rwmol.AddBond(st_, new_id, new_bond_type)
    return rwmol


def split_hetero(rwmol, bonds):
    carbon_ends = set()
    for b_idx in sorted(bonds, reverse=True):
        b = rwmol.GetBondWithIdx(b_idx)
        st_ = b.GetBeginAtomIdx()
        end_ = b.GetEndAtomIdx()
        rwmol.RemoveBond(st_, end_)
        for a in (st_, end_):
            if rwmol.GetAtomWithIdx(a).GetSymbol() == 'C':
                carbon_ends.add(a)

    for a in carbon_ends:
        atom = rwmol.GetAtomWithIdx(a)
        atom.SetNumExplicitHs(0)
        atom.UpdatePropertyCache()
        valence = atom.GetExplicitValence()

        if valence == 1:
            new_id1 = rwmol.AddAtom(rdc.Atom(8))
            rwmol.AddBond(a, new_id1, rdc.BondType.SINGLE)
            new_id2 = rwmol.AddAtom(rdc.Atom(8))
            rwmol.AddBond(a, new_id2, rdc.BondType.DOUBLE)
        elif valence == 2:
            new_id2 = rwmol.AddAtom(rdc.Atom(8))
            rwmol.AddBond(a, new_id2, rdc.BondType.DOUBLE)
        elif valence == 3:
            new_id2 = rwmol.AddAtom(rdc.Atom(8))
            rwmol.AddBond(a, new_id2, rdc.BondType.SINGLE)
    return rwmol


def split_by_monomer_bonds(rban_record):
    HYDROX = ['AMINO', 'ESTER', 'THIOETHER', 'CARBON_CARBON', 'NITROGEN_CARBON2',
              'NITROGEN_CARBON',
              'CARBON_CARBOXIL', 'CARBOXIL_CARBOXIL', 'GLYCOSIDIC']
    HYDROX_TRANS = ['AMINO_TRANS']
    HETERO = ['HETEROCYCLE', 'THIAZOLE_CYCLE', 'PYRIMIDINE_CYCLE']

    mol = rdc.MolFromSmiles(rban_record['isomericSmiles'])
    # kekulize so that we would not get this annoying Aromatic flag exception
    rdc.Kekulize(mol, clearAromaticFlags=True)

    # we will need editable mol instance
    rwmol = rdc.RWMol(mol)
    monomers = rban_record['monomericGraph']['monomericGraph']['monomers']

    bonds, types = get_monomer_bonds(rwmol, rban_record)
    amino_bonds = [b for b, t in zip(bonds, types) if t in HYDROX]
    rwmol = split_amino(rwmol, amino_bonds, trans=False)

    bonds, types = get_monomer_bonds(rwmol, rban_record)
    amino_trans_bonds = [b for b, t in zip(bonds, types) if t in HYDROX_TRANS]
    if amino_trans_bonds:
        rwmol = split_amino(rwmol, amino_trans_bonds, trans=True)
    # TODO: consider other(?) bond types

    bonds, types = get_monomer_bonds(rwmol, rban_record)
    hetero_bonds = [b for b, t in zip(bonds, types) if t in HETERO]
    if hetero_bonds:
        rwmol = split_hetero(rwmol, hetero_bonds)

    frag_ids = rdc.GetMolFrags(rwmol)
    frag_smi = [rdc.MolToSmiles(x) for x in rdc.GetMolFrags(rwmol, asMols=True)]
    monomer_smiles_dict = defaultdict(list)
    for k, monomer in enumerate(monomers):
        monomer_atoms = set(monomer['monomer']['atoms'])
        for i, frag_atoms in enumerate(map(set, frag_ids)):
            if len(monomer_atoms & frag_atoms) > 0:
                monomer_smiles_dict[k].append(frag_smi[i])
    return monomer_smiles_dict


def get_monomers_chirality(rban_record, debug=False):
    monomers = rban_record['monomericGraph']['monomericGraph']['monomers']
    monomer_smiles_dict = split_by_monomer_bonds(rban_record)
    res = []
    for i in range(len(monomers)):
        try:
            res.append(get_aa_chirality(monomer_smiles_dict[i][0]))
        except NoChiralCenters:
            if debug: print(f'{i}: NoChiralCenters')
            res.append(False)
        except NoCarboxyl:
            if debug: print(f'{i}: NoCarboxyl')
            res.append(False)
        except WTF:
            if debug: print(f'{i}: WTF:', monomer_smiles_dict[i][0])
            res.append(False)
    return res


def monomericgraph_to_string(g):
    monomers = [(m['monomer']['index'], m['monomer']['monomer']['monomer']) for m in g['monomers']]
    monomers = [m[1] for m in sorted(monomers)]
    bonds = [b['bond']['monomers'] for b in g['bonds']]
    result = ','.join(monomers) + ';' + ';'.join([f'{b[0]-1},{b[1]-1}' for b in bonds])
    return result


def monomericgraph_components_to_string(rban_record, nerpa_monomers, allowed_bond_types, na='NA'):
    g = rban_record['monomericGraph']['monomericGraph']
    nodes = [m['monomer']['monomer']['monomer'] for m in g['monomers']]
    is_d = get_monomers_chirality(rban_record)
    assert len(is_d) == len(nodes)
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

    registered_set = set(nerpa_monomers)
    def get_comp_name(comp):
        name = na
        for _name, idx in sorted([(nodes[n], n) for n in comp]):
            if _name in registered_set:
                name = _name
                if is_d[idx]:
                    name = f'@{name}'
                if len(comp) > 1:
                    name = f'*{name}'
                break
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
    for i, rban_record in enumerate(graphs):
        # gr = monomericgraph_to_string(g['monomericGraph']['monomericGraph'])
        gr = monomericgraph_components_to_string(rban_record,
                                                 nerpa_monomers, PNP_BONDS, UNDEFINED_NAME)
        out.append(f'{rban_record["id"]} {gr}')
        if rban_input_json and 'extra' in rban_input_json[i].keys():
            out[-1] += f' {rban_input_json[i]["extra"]}'

    print('\n'.join(out))


if __name__ == '__main__':
    main()