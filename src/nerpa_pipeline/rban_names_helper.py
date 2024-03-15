from __future__ import annotations
from typing import (
    Dict,
    NamedTuple,
    Tuple
)
from src.data_types import (
    MonomerResidue,
    NRP_Monomer_Modification,
    UNKNOWN_RESIDUE
)
import csv
from collections import defaultdict
from pathlib import Path


class Parsed_rBAN_Name(NamedTuple):
    residue: MonomerResidue
    modifications: Tuple[NRP_Monomer_Modification, ...]


def parse_modifications(mods_str: str) -> Tuple[NRP_Monomer_Modification, ...]:
    def parse_mod(mod: str) -> NRP_Monomer_Modification:
        match mod:
            case 'MT': return NRP_Monomer_Modification.METHYLATION
            case 'unk': return NRP_Monomer_Modification.UNKNOWN
            case other: raise ValueError(f"Unknown modification")

    return () if mods_str == '-' else tuple(map(parse_mod, mods_str.split('+')))


class rBAN_Names_Helper:
    parsed_rban_name: Dict[str, Parsed_rBAN_Name]


    def __init__(self, names_data_file: Path):  # receives contents of file monomers.tsv
        self.parsed_rban_name = defaultdict(lambda: Parsed_rBAN_Name(residue=UNKNOWN_RESIDUE,
                                                                     modifications=()))
        for row in csv.DictReader(names_data_file.open('r'), delimiter='\t'):
            self.parsed_rban_name[row['Code']] = Parsed_rBAN_Name(residue=row['NameID'],
                                                                  modifications=parse_modifications(row['Modification']))
