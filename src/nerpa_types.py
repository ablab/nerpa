from __future__ import annotations
from typing import (
    Dict,
    List,
    Tuple)
from enum import auto, Enum
import yaml
from dataclasses import asdict, dataclass

PredictionScore = float
MonomerResidue = str


def enum_representer(dumper, e: Enum):
    return dumper.represent_scalar(f'!{e.__class__.__name__}', e.name)


class Chirality(Enum):
    L = auto()
    D = auto()
    Unknown = auto()

yaml.add_representer(Chirality, enum_representer)


class PTM(Enum):  # post-translational modification
    Methylation = auto()
    Unknown = auto()

yaml.add_representer(PTM, enum_representer)


@dataclass
class NRP_Monomer:
    residue: MonomerResidue
    ptms: Tuple[PTM, ...]
    chirality: Chirality
    rBAN_name: str


class BGC_Module_Modification(Enum):
    Epimerization = auto()
    Methylation = auto()

yaml.add_representer(BGC_Module_Modification, enum_representer)


@dataclass
class BGC_Module_Prediction:
    gene_id: str
    module_idx: int
    residue_score: Dict[MonomerResidue, PredictionScore]
    modifications: Tuple[BGC_Module_Modification, ...]
    iterative_module: bool
    iterative_gene: bool


@dataclass
class BGC_variant:
    genome_id: str
    bgc_id: str
    predictions: List[BGC_Module_Prediction]


@dataclass
class NRP_fragment:
    monomers: List[NRP_Monomer]
    is_cyclic: bool
    to_rban_indexes: List[int]


@dataclass
class NRP_variant:
    nrp_name: str
    fragments: List[NRP_fragment]


def main():
    nrp_fragments = [NRP_fragment(monomers=[NRP_Monomer(residue='arg',
                                                        ptms=(PTM.Methylation,),
                                                        chirality=Chirality.Unknown,
                                                        rBAN_name='Arg'),
                                            NRP_Monomer(residue='leu',
                                                        ptms=(PTM.Unknown,),
                                                        chirality=Chirality.L,
                                                        rBAN_name='Leu')],
                                  is_cyclic=False,
                                  to_rban_indexes=[0,1]),
                     NRP_fragment(monomers=[NRP_Monomer(residue='ile',
                                                        ptms=(),
                                                        chirality=Chirality.Unknown,
                                                        rBAN_name='Ile'),
                                            NRP_Monomer(residue='pro',
                                                        ptms=(PTM.Unknown, PTM.Methylation),
                                                        chirality=Chirality.L,
                                                        rBAN_name='Pro'),
                                            NRP_Monomer(residue='val',
                                                        ptms=(PTM.Methylation,),
                                                        chirality=Chirality.D,
                                                        rBAN_name='Val')],
                                  is_cyclic=True,
                                  to_rban_indexes=[3, 7, 2])]
    nrp_variants = [NRP_variant(fragments=[nrp_fragments[0]], nrp_name='terrific_nrp#42'),
                    NRP_variant(fragments=[nrp_fragments[1]], nrp_name='terrific_nrp#48')]
    bgc_preds = [BGC_Module_Prediction(residue_score={'arg': 0.8, 'leu': 0.7},
                                       modifications=(),
                                       iterative_module=False,
                                       iterative_gene=False,
                                       gene_id='gene1',
                                       module_idx=0),
                 BGC_Module_Prediction(residue_score={'trp': 0.8, 'phe': 0.7},
                                       modifications=(BGC_Module_Modification.Epimerization, BGC_Module_Modification.Methylation,),
                                       iterative_module=True,
                                       iterative_gene=False,
                                       gene_id='gene1',
                                       module_idx=1)]
    bgc_variant = BGC_variant(predictions=bgc_preds,
                              genome_id='genome#189',
                              bgc_id='bgc#39')

    # dirty hack to erase information about types and make output less verbose
    # https://github.com/yaml/pyyaml/issues/408
    yaml.emitter.Emitter.prepare_tag = lambda self, tag: ''

    # another hack (albeit less dirty) to forbid yaml creating references
    # https://stackoverflow.com/questions/13518819/avoid-references-in-pyyaml
    yaml.Dumper.ignore_aliases = lambda *args: True
    with open('test_bgc_predictions_output.yaml', 'w') as out:
        yaml.dump(asdict(bgc_variant), out,
                  default_flow_style=None, sort_keys=False)
    with open('test_nrp_variants_output.yaml', 'w') as out:
        yaml.dump(list(map(asdict, nrp_variants)), out,
                  default_flow_style=None, sort_keys=False)


if __name__ == "__main__":
    main()

