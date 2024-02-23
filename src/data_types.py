from __future__ import annotations
from typing import (
    Dict,
    List,
    Tuple)
from enum import auto, Enum
import yaml
from dataclasses import asdict, dataclass

LogProb = float
MonomerResidue = str


def enum_representer(dumper, e: Enum):
    return dumper.represent_scalar(f'!{e.__class__.__name__}', e.name)


class Chirality(Enum):
    L = auto()
    D = auto()
    UNKNOWN = auto()

yaml.add_representer(Chirality, enum_representer)


class NRP_Monomer_Modification(Enum):  # post-translational modification
    METHYLATION = auto()
    UNKNOWN = auto()

yaml.add_representer(NRP_Monomer_Modification, enum_representer)


@dataclass
class NRP_Monomer:
    residue: MonomerResidue
    modifications: Tuple[NRP_Monomer_Modification, ...]
    chirality: Chirality
    rban_name: str
    rban_idx: int


class BGC_Module_Modification(Enum):
    EPIMERIZATION = auto()
    METHYLATION = auto()

yaml.add_representer(BGC_Module_Modification, enum_representer)


@dataclass
class BGC_Module_Prediction:
    gene_id: str
    module_idx: int
    residue_score: Dict[MonomerResidue, LogProb]
    modifications: Tuple[BGC_Module_Modification, ...]
    iterative_module: bool
    iterative_gene: bool


@dataclass
class BGC_Variant:
    genome_id: str
    bgc_id: str
    tentative_assembly_line: List[BGC_Module_Prediction]


@dataclass
class NRP_Fragment:
    monomers: List[NRP_Monomer]
    is_cyclic: bool


@dataclass
class NRP_Variant:
    nrp_id: str
    fragments: List[NRP_Fragment]
