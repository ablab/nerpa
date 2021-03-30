from collections import defaultdict
from io import StringIO
from contextlib import redirect_stderr
from functools import wraps

from rdkit import Chem as rdc
import networkx as nx


class NoChiralCenters(Exception):
    pass
class NoCarboxyl(Exception):
    pass
class UnexpectedError(Exception):
    pass
class RDKitError(Exception):
    pass
class PKError(Exception):
    pass


def handle_RDKit_errors(suppress=True):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            sio = StringIO()
            with redirect_stderr(sio):
                res = func(*args, **kwargs)
            errors = sio.getvalue()
            if not suppress and errors:
                raise RDKitError(errors)
            return res
        return wrapper
    return decorator


def split_amino_bonds(rwmol, bonds):
    """
        TODO: Incorrect results for CARBOXYL_CARBOXYL. It looks pretty symmetrical, but
        rBAN adds -OH only to one of the ends of the bond (and it's not like to start/end
        every time)
    :param rwmol:
    :param bonds:
    :return:
    """

    for b_idx in sorted(bonds, reverse=True):
        # reverse because when a bond or atom is deleted, the bond or
        # atom indices are remapped. If you remove bonds with a higher index first,
        # bonds with lower indices will not be remapped.
        b = rwmol.GetBondWithIdx(b_idx)
        st_ = b.GetBeginAtomIdx()
        end_ = b.GetEndAtomIdx()

        rwmol.RemoveBond(st_, end_)
        for idx in (st_, end_):
            atom_ = rwmol.GetAtomWithIdx(idx)
            atom_.SetNumExplicitHs(0)
            atom_.UpdatePropertyCache()
            if atom_.GetSymbol() == 'C':
                neigh_oxygen = [x.GetIdx() for x in atom_.GetNeighbors() if x.GetSymbol() == 'O']
                bond_types = [rwmol.GetBondBetweenAtoms(idx, x).GetBondType()
                              for x in neigh_oxygen ]
                for bt in bond_types:
                    if (bt == rdc.BondType.SINGLE and atom_.GetExplicitValence() < 3):
                        new_idx = rwmol.AddAtom(rdc.Atom(8))
                        rwmol.AddBond(idx, new_idx, rdc.BondType.DOUBLE)
                    elif bt == rdc.BondType.DOUBLE:
                        new_idx = rwmol.AddAtom(rdc.Atom(8))
                        rwmol.AddBond(idx, new_idx, rdc.BondType.SINGLE)
    return rwmol


def split_heterocycle_bonds(rwmol, bonds):
    # TODO: sometimes adds OH when rBAN doesn't (when valence == 1 and type == HETEROCYCLE),
    # but in other cases (e.g. val == 1, type == Oxazole), the result matches w/ rBAN
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


def split_other_bonds(rwmol, bonds):
    for b_idx in sorted(bonds, reverse=True):
        b = rwmol.GetBondWithIdx(b_idx)
        st_ = b.GetBeginAtomIdx()
        end_ = b.GetEndAtomIdx()
        rwmol.RemoveBond(st_, end_)
    return rwmol


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


def split_by_monomer_bonds(rban_record):

    # see rban.src.main.java.molecules.bond.BondType for the full list
    HYDROX = [
        'AMINO', 'AMINO_TRANS', 'ESTER', 'THIOETHER',  'CARBON_CARBOXIL', 'GLYCOSIDIC', 'CARBOXIL_CARBOXIL'
    ]
    HETERO = [
        'HETEROCYCLE', 'THIAZOLE_CYCLE', 'PYRIMIDINE_CYCLE', 'OXAZOLE_CYCLE'
    ]
    OTHER = [
        'ARYL_ETHER', 'ARYL_ETHER2', 'DISULFIDE', 'NITROGEN_CARBON', 'CARBON_CARBON',
        'NITROGEN_CARBON2', 'SULFUR_CARBON',
    ]

    mol = rdc.MolFromSmiles(rban_record['isomericSmiles'])
    # kekulize so that we would not get this annoying Aromatic flag error
    rdc.Kekulize(mol, clearAromaticFlags=True)

    # we need editable mol instance
    rwmol = rdc.RWMol(mol)
    monomers = rban_record['monomericGraph']['monomericGraph']['monomers']

    bonds, types = get_monomer_bonds(rwmol, rban_record)
    amino_bonds = [b for b, t in zip(bonds, types) if t in HYDROX]
    rwmol = split_amino_bonds(rwmol, amino_bonds)

    bonds, types = get_monomer_bonds(rwmol, rban_record)
    hetero_bonds = [b for b, t in zip(bonds, types) if t in HETERO]
    rwmol = split_heterocycle_bonds(rwmol, hetero_bonds)

    bonds, types = get_monomer_bonds(rwmol, rban_record)
    rwmol = split_other_bonds(rwmol, bonds)

    frags = rdc.GetMolFrags(rwmol)
    frag_smi = [rdc.MolToSmiles(x) for x in rdc.GetMolFrags(rwmol, asMols=True)]

    monomer_smiles_dict = {}
    for monomer in monomers:
        monomer_atoms = set(monomer['monomer']['atoms'])
        monomer_idx = monomer['monomer']['index']
        for frag_idx, frag_atoms in enumerate(map(set, frags)):
            if len(monomer_atoms & frag_atoms) > 0:
                monomer_smiles_dict[monomer_idx] = frag_smi[frag_idx]
    return monomer_smiles_dict


def get_aa_chirality(smi):
    mol = rdc.MolFromSmiles(smi)
    chiral_centers = rdc.FindMolChiralCenters(mol)
    if not chiral_centers:
        raise NoChiralCenters

    def _get_chirality(match, q):
        RCN_pattern = 'RCN' * 2
        NCR_pattern = 'NCR' * 2
        for idx, rs in chiral_centers:
            # Sadly, we cannot make use of rs here, see
            # https://rdcistry.stackexchange.com/questions/90826/difference-between-ld-and-sr-in-naming
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
    raise UnexpectedError


def has_chiral_centers(smi):
    mol = rdc.MolFromSmiles(smi)
    chiral_centers = rdc.FindMolChiralCenters(mol)
    if not chiral_centers:
        return False
    return True


@handle_RDKit_errors(suppress=True)
def get_monomers_chirality(rban_record, debug=False):
    """
    Returns list of True/False/None for each monomer in rban's monomeric graph.
    True ~ D-, False ~ L-, None ~ unrecognized

    :param rban_record:
    :param debug:
    :return:
    """
    if not has_chiral_centers(rban_record['isomericSmiles']):
        return {}

    monomer_smiles_dict = split_by_monomer_bonds(rban_record)
    res = {}
    for monomer_idx, smi in monomer_smiles_dict.items():
        try:
            res[monomer_idx] = get_aa_chirality(smi)
        except NoChiralCenters:
            if debug: print(f'{monomer_idx}: NoChiralCenters')
        except NoCarboxyl:
            if debug: print(f'{monomer_idx}: NoCarboxyl')
        except UnexpectedError:
            if debug: print(f'{monomer_idx}: UnexpectedError:', smi)

    return res


def mol_to_nx(mol):
    G = nx.Graph()
    for atom in mol.GetAtoms():
        G.add_node(atom.GetIdx(),
                   atomic_num=atom.GetAtomicNum(),
                   formal_charge=atom.GetFormalCharge(),
                   chiral_tag=atom.GetChiralTag(),
                   hybridization=atom.GetHybridization(),
                   num_explicit_hs=atom.GetNumExplicitHs(),
                   is_aromatic=atom.GetIsAromatic())
    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(),
                   bond.GetEndAtomIdx(),
                   bond_type=bond.GetBondType())
    return G


def find_aa_pk_bond(mol):
    """
        returns pair of atom indices (a1, a2)
        a1: nrp end
        a2: pk end
    """
    G = mol_to_nx(mol)

    cooh = rdc.MolFromSmarts('[CX3](=[OX1])O')
    cooh_match = mol.GetSubstructMatches(cooh)
    if not cooh_match:
        raise PKError

    # backbone starting carbon
    starts = [x[0] for x in cooh_match]

    # identify all nitrogens
    nitrogens = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if not nitrogens:
        raise PKError

    # identify all putative backbone ends (carbons connected to nitrogens)
    ends = sorted(list(
        set(c for n in nitrogens for c in G.neighbors(n)
            if G.nodes[c]['atomic_num'] == 6)
    ))

    cyclic_nodes = set(node for cycle in nx.cycle_basis(G) for node in cycle)

    # list of all backbone candidates (no cycles allowed!)
    paths = []
    for s in starts:
        for e in ends:
            for pth in nx.shortest_simple_paths(G, s, e):
                all_atomic_nums_are_six = True
                for x in pth:
                    if G.nodes[x]['atomic_num'] != 6:
                        all_atomic_nums_are_six = False
                        break
                if all_atomic_nums_are_six:
                    # exclude last carbon to support Proline & Homoproline
                    for x in pth[:-1]:
                        if x in cyclic_nodes:
                            break
                if len(pth) > 3:
                    paths.append((len(pth), e, pth))
    if not paths:
        raise PKError
    paths = sorted(paths)
    shortest_path = paths[0][2]
    if len(shortest_path) % 2:
        bond_break = (shortest_path[-3], shortest_path[-4])
    else:
        bond_break = (shortest_path[-2], shortest_path[-3])
    return bond_break


def split_aa_pk_hybrid(smi):
    mol = rdc.MolFromSmiles(smi)
    if not mol:
        raise PKError
    bond_break = find_aa_pk_bond(mol)

    mw = rdc.RWMol(mol)
    oxs = [a.GetIdx() for a in mw.GetAtomWithIdx(bond_break[0]).GetNeighbors()
           if a.GetAtomicNum() == 8]
    if len(oxs) > 1:
        raise PKError

    if oxs:
        bond_type = mw.GetBondBetweenAtoms(bond_break[0], oxs[0]).GetBondType()
        if bond_type == rdc.BondType.DOUBLE:
            cond_type = 'KS'
            new_ox = mw.AddAtom(rdc.Atom(8))
            mw.AddBond(bond_break[0], new_ox, rdc.BondType.SINGLE)
        elif bond_type == rdc.BondType.SINGLE:
            cond_type = 'KS+KR'
            new_ox = mw.AddAtom(rdc.Atom(8))
            mw.AddBond(bond_break[0], new_ox, rdc.BondType.DOUBLE)
        else:
            raise PKError
    else:
        bond_type = mw.GetBondBetweenAtoms(*bond_break).GetBondType()
        if bond_type == rdc.BondType.DOUBLE:
            cond_type = 'KS+KR+DH'
        elif bond_type == rdc.BondType.SINGLE:
            cond_type = 'KS+KR+DH+ER'
        else:
            raise PKError
        new_ox1 = mw.AddAtom(rdc.Atom(8))
        new_ox2 = mw.AddAtom(rdc.Atom(8))
        mw.AddBond(bond_break[0], new_ox1, rdc.BondType.SINGLE)
        mw.AddBond(bond_break[0], new_ox2, rdc.BondType.DOUBLE)
    mw.RemoveBond(*bond_break)

    aa_part_smi = None
    pk_part_smi = None
    for frag_atoms in rdc.GetMolFrags(mw):
        if bond_break[0] in frag_atoms:
            aa_part_smi = rdc.MolFragmentToSmiles(mw, frag_atoms)
        else:
            pk_part_smi = rdc.MolFragmentToSmiles(mw, frag_atoms)
    if (not aa_part_smi) or (not pk_part_smi):
        raise PKError
    return aa_part_smi, pk_part_smi, cond_type


