from src.data_types import (
    BGC_Module_Prediction,
    BGC_Module_Modification,
    NRP_Monomer,
    LogProb,
    NRP_Monomer_Modification)
from src.NewMatcher.dp_config import DP_Config, ModMatch, ChiralityMatch
from dataclasses import dataclass
from enum import Enum, auto


@dataclass
class DP_Helper:
    dp_config: DP_Config

    def bgc_module_remove(self, bgc_pred: BGC_Module_Prediction) -> LogProb:
        return self.dp_config.bgc_module_remove_score

    def nrp_mon_remove(self, mon: NRP_Monomer) -> LogProb:
        return self.dp_config.nrp_mon_remove_score[mon.residue]  # should we take into account methylation and chirality as well?

    def match(self, bgc_pred: BGC_Module_Prediction,
              nrp_mon: NRP_Monomer) -> LogProb:

        return bgc_pred.residue_score[nrp_mon.residue] \
            + self.dp_config.mod_score[ModMatch(mod=NRP_Monomer_Modification.METHYLATION,
                                                bgc_mod=BGC_Module_Modification.METHYLATION in bgc_pred.modifications,
                                                nrp_mod=NRP_Monomer_Modification.METHYLATION in nrp_mon.modifications)] \
            + self.dp_config.chirality_score[ChiralityMatch(bgc_epim=BGC_Module_Modification.EPIMERIZATION in bgc_pred.modifications,
                                                            nrp_chr=nrp_mon.chirality)]

    def iterate_module(self) -> LogProb:
        return 0

    def iterate_gene(self) -> LogProb:
        return 0

    def null_hypothesis_score(self, nrp_monomer: NRP_Monomer) -> LogProb:
        return self.dp_config.null_hypothesis_residue_score[nrp_monomer.residue] \
            + self.dp_config.null_hypothesis_chirality_score[nrp_monomer.chirality] \
            + self.dp_config.null_hypothesis_mod_score[(NRP_Monomer_Modification.METHYLATION,
                                                        NRP_Monomer_Modification.METHYLATION in nrp_monomer)]

