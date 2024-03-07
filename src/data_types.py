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
UNKNOWN_RESIDUE = 'unk'
rBAN_Residue_Name = str


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

    @classmethod
    def from_yaml_dict(cls, data: dict) -> NRP_Monomer:
        return cls(residue=MonomerResidue(data['residue']),
                   modifications=tuple(NRP_Monomer_Modification[mod]
                                       for mod in data['modifications']),
                   chirality=Chirality(data['chirality']),
                   rban_name=data['rban_name'],
                   rban_idx=data['rban_idx'])


class BGC_Module_Modification(Enum):
    EPIMERIZATION = auto()
    METHYLATION = auto()

yaml.add_representer(BGC_Module_Modification, enum_representer)


@dataclass
class BGC_Module:
    gene_id: str
    module_idx: int
    residue_score: Dict[MonomerResidue, LogProb]
    modifications: Tuple[BGC_Module_Modification, ...]
    iterative_module: bool
    iterative_gene: bool

    @classmethod
    def from_yaml_dict(cls, data: dict) -> BGC_Module:
        return cls(gene_id=data['gene_id'],
                   module_idx=data['module_idx'],
                   residue_score=data['residue_score'],
                   modifications=tuple(BGC_Module_Modification[mod]
                                       for mod in data['modifications']),
                   iterative_module=data['iterative_module'],
                   iterative_gene=data['iterative_gene'])

@dataclass
class BGC_Variant:
    genome_id: str
    bgc_id: str
    tentative_assembly_line: List[BGC_Module]

    @classmethod
    def from_yaml_dict(cls, data: dict) -> BGC_Variant:
        return cls(genome_id=data['genome_id'],
                   bgc_id=data['bgc_id'],
                   tentative_assembly_line=list(map(BGC_Module.from_yaml_dict,
                                                    data['tentative_assembly_line'])))


@dataclass
class NRP_Fragment:
    monomers: List[NRP_Monomer]
    is_cyclic: bool

    @classmethod
    def from_yaml_dict(cls, data: dict) -> NRP_Fragment:
        return cls(is_cyclic=data['is_cyclic'],
                   monomers=list(map(NRP_Monomer.from_yaml_dict, data['monomers'])))


@dataclass
class NRP_Variant:
    nrp_id: str
    fragments: List[NRP_Fragment]

    @classmethod
    def from_yaml_dict(cls, data: dict) -> NRP_Variant:
        return cls(nrp_id=data['nrp_id'],
                   fragments=list(map(NRP_Fragment.from_yaml_dict, data['fragments'])))
