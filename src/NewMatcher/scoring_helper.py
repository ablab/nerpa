from src.data_types import (
    BGC_Module,
    BGC_Module_Modification,
    NRP_Monomer,
    LogProb,
    NRP_Monomer_Modification,
    UNKNOWN_RESIDUE
)
from src.NewMatcher.scoring_config import ScoringConfig, ModMatch, ChiralityMatch
from dataclasses import dataclass
from enum import Enum, auto


@dataclass
class ScoringHelper:
    scoring_config: ScoringConfig

    def bgc_module_skip(self, bgc_pred: BGC_Module) -> LogProb:
        return self.scoring_config.bgc_module_skip_score

    def nrp_mon_skip(self, mon: NRP_Monomer) -> LogProb:
        return self.scoring_config.nrp_mon_skip_score[mon.residue]  # should we take into account methylation and chirality as well?

    def match(self, bgc_pred: BGC_Module,
              nrp_mon: NRP_Monomer) -> LogProb:
        residue_score = bgc_pred.residue_score[nrp_mon.residue]
        if nrp_mon.residue == UNKNOWN_RESIDUE:
            residue_score /= self.scoring_config.num_unknown_residues
        return residue_score \
            + self.scoring_config.mod_score[ModMatch(mod=NRP_Monomer_Modification.METHYLATION,
                                                     bgc_mod=BGC_Module_Modification.METHYLATION in bgc_pred.modifications,
                                                     nrp_mod=NRP_Monomer_Modification.METHYLATION in nrp_mon.modifications)] \
            + self.scoring_config.chirality_score[ChiralityMatch(bgc_epim=BGC_Module_Modification.EPIMERIZATION in bgc_pred.modifications,
                                                                 nrp_chr=nrp_mon.chirality)]

    def iterate_module(self) -> LogProb:
        return 0

    def iterate_gene(self) -> LogProb:
        return 0

    def null_hypothesis_score(self, nrp_monomer: NRP_Monomer) -> LogProb:
        return self.scoring_config.null_hypothesis_residue_score[nrp_monomer.residue] \
            + self.scoring_config.null_hypothesis_chirality_score[nrp_monomer.chirality] \
            + self.scoring_config.null_hypothesis_mod_score[(NRP_Monomer_Modification.METHYLATION,
                                                             NRP_Monomer_Modification.METHYLATION in nrp_monomer.modifications)]

