from typing import Callable, Dict, NamedTuple, Tuple
from src.data_types import (
    Chirality,
    LogProb,
    MonomerResidue,
    NRP_Monomer_Modification,
    UNKNOWN_RESIDUE
)
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
import yaml


class ChiralityMatch(NamedTuple):
    bgc_epim: bool
    nrp_chr: Chirality


class ModMatch(NamedTuple):
    mod: NRP_Monomer_Modification
    bgc_mod: bool
    nrp_mod: bool


@dataclass
class ScoringConfig:
    bgc_module_skip_score: LogProb
    nrp_mon_skip_score: Dict[MonomerResidue, LogProb]
    num_unknown_residues: int
    mod_score: Dict[ModMatch, LogProb]
    chirality_score: Dict[ChiralityMatch, LogProb]

    null_hypothesis_residue_score: Dict[MonomerResidue, LogProb]
    null_hypothesis_mod_score: Dict[Tuple[NRP_Monomer_Modification, bool], LogProb]
    null_hypothesis_chirality_score: Dict[Chirality, LogProb]

    max_gene_reps: int
    max_module_reps: int


def load_config(path_to_config: Path) -> ScoringConfig:
    cfg = yaml.safe_load(path_to_config.open('r'))

    mod_score = {ModMatch(mod=mod, bgc_mod=bgc_mod, nrp_mod=nrp_mod):
                     cfg['modification_score'][f'Mod_{mod.name}_BGC_{bgc_mod}_NRP_{nrp_mod}']
                 for mod in NRP_Monomer_Modification
                 for bgc_mod in (False, True)
                 for nrp_mod in (False, True)
                 if mod != NRP_Monomer_Modification.UNKNOWN}
    chirality_score = {ChiralityMatch(bgc_epim=bgc_epim, nrp_chr=nrp_chr):
                           cfg['chirality_score'][f'BGC_{bgc_epim}_NRP_{nrp_chr.name}']
                       for bgc_epim in (False, True)
                       for nrp_chr in Chirality}

    nrp_mon_default_skip_score = cfg['nrp_monomer_default_skip_score']
    nrp_mon_skip_score = defaultdict(lambda: nrp_mon_default_skip_score + cfg['nrp_monomer_frequencies'][UNKNOWN_RESIDUE],
                                     {residue: nrp_mon_default_skip_score + residue_frequency
                                      for residue, residue_frequency in cfg['nrp_monomer_frequencies'].items()})

    null_hypothesis_residue_score = defaultdict(lambda: cfg['null_hypothesis_residue_score'][UNKNOWN_RESIDUE],
                                                cfg['null_hypothesis_residue_score'])
    def get_null_hyp_mod_score(mod: NRP_Monomer_Modification,
                               mod_in_nrp: bool) -> LogProb:
        mod_freq = cfg['module_modification_frequency'][mod.name]
        return mod_freq * mod_score[ModMatch(mod=mod, bgc_mod=True, nrp_mod=mod_in_nrp)] \
            + (1-mod_freq) * mod_score[ModMatch(mod=mod, bgc_mod=False, nrp_mod=mod_in_nrp)]

    null_hypothesis_mod_score = {(mod, mod_in_nrp): get_null_hyp_mod_score(mod, mod_in_nrp)
                                 for mod in NRP_Monomer_Modification
                                 for mod_in_nrp in (False, True)
                                 if mod != NRP_Monomer_Modification.UNKNOWN}

    def get_null_hyp_chr_score(nrp_chr: Chirality) -> LogProb:
        epim_freq = cfg['module_epimerization_frequency']
        return epim_freq * chirality_score[ChiralityMatch(bgc_epim=True, nrp_chr=nrp_chr)] \
            + (1-epim_freq) * chirality_score[ChiralityMatch(bgc_epim=False, nrp_chr=nrp_chr)]

    null_hypothesis_chirality_score = {chr: get_null_hyp_chr_score(chr)
                                       for chr in Chirality}
    return ScoringConfig(bgc_module_skip_score=cfg['bgc_module_skip_score'],
                         nrp_mon_skip_score=nrp_mon_skip_score,
                         num_unknown_residues=cfg['num_unknown_residues'],
                         mod_score=mod_score,
                         chirality_score=chirality_score,
                         null_hypothesis_residue_score=null_hypothesis_residue_score,
                         null_hypothesis_mod_score=null_hypothesis_mod_score,
                         null_hypothesis_chirality_score=null_hypothesis_chirality_score,
                         max_gene_reps=cfg['max_gene_repetitions'],
                         max_module_reps=cfg['max_module_repetitions'])

