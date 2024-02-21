from typing import Callable, Dict, NamedTuple
from src.data_types import (
    Chirality,
    LogProb,
    MonomerResidue,
    PTM)
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
import yaml


class ChiralityMatch(NamedTuple):
    bgc_epim: bool
    nrp_chr: Chirality

class PTM_Match(NamedTuple):
    ptm: PTM
    bgc_mod: bool
    nrp_ptm: bool

@dataclass
class DP_Config:
    bgc_module_remove_score: LogProb
    nrp_mon_remove_score: Dict[MonomerResidue, LogProb]
    ptm_score: Dict[PTM_Match, LogProb]
    chirality_score: Dict[ChiralityMatch, LogProb]

    null_hypothesis_residue_score: Dict[MonomerResidue, LogProb]
    null_hypothesis_ptm_score: Dict[PTM, LogProb]
    null_hypothesis_chirality_score: Dict[Chirality, LogProb]


def load_config(path_to_config: Path) -> DP_Config:
    cfg = yaml.safe_load(path_to_config.open('r'))

    ptm_score = {PTM_Match(ptm=ptm, bgc_mod=bgc_mod, nrp_ptm=nrp_ptm):
                     cfg['ptm_score'][f'PTM_{ptm.name}_BGC_{bgc_mod}_NRP_{nrp_ptm}']
                 for ptm in PTM
                 for bgc_mod in (False, True)
                 for nrp_ptm in (False, True)}
    chirality_score = {ChiralityMatch(bgc_epim=bgc_epim, nrp_chr=nrp_chr):
                           cfg.chirality_score[f'BGC_{bgc_epim}_NRP_{nrp_chr.name}']
                       for bgc_epim in (False, True)
                       for nrp_chr in Chirality}

    nrp_mon_remove_score = defaultdict(lambda: cfg['nrp_mon_remove_score']['unk'],
                                       cfg['nrp_mon_remove_score'])
    null_hypothesis_residue_score = defaultdict(lambda: cfg['null_hypothesis_residue_score']['unk'],
                                                cfg['null_hypothesis_residue_score'])
    null_hypothesis_ptm_score = {PTM[ptm]: score
                                 for ptm, score in cfg.null_hypothesis_ptm_score.items()}
    null_hypothesis_chirality_score = {Chirality[chr]: score
                                       for chr, score in cfg.null_hypothesis_chirality_score.items()}
    return DP_Config(bgc_module_remove_score=cfg['bgc_module_remove_score'],
                     nrp_mon_remove_score=nrp_mon_remove_score,
                     ptm_score=ptm_score,
                     chirality_score=chirality_score,
                     null_hypothesis_residue_score=null_hypothesis_residue_score,
                     null_hypothesis_ptm_score=null_hypothesis_ptm_score,
                     null_hypothesis_chirality_score=null_hypothesis_chirality_score)

