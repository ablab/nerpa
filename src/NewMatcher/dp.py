from typing import (
    Any,
    Dict,
    List,
    NamedTuple,
    Tuple)
from src.data_types import (
    BGC_Module_Prediction,
    Chirality,
    LogProb,
    NRP_Monomer,
)
from src.NewMatcher.dp_helper import DP_Helper
from src.NewMatcher.alignment_types import AlignmentStep, Alignment
import numpy as np
from itertools import groupby
from functools import partial


DPTable = Any  # TODO: use proper numpy.typing types

class DPState(NamedTuple):
    module_pos: int
    monomer_pos: int
    num_gene_reps: int
    num_module_reps: int

def valid_state(state: DPState) -> bool:
    return all(idx >= 0 for idx in state)

class DP_Value(NamedTuple):
    score: LogProb
    parent: DPState
    action: str  # TODO: make it a proper type


def dp_recalc(dp_table: DPTable,
              transitions: List[Tuple[DPState, LogProb, str]]) -> Tuple[LogProb, DPState, str]:
    return max((dp_table[state] + score, state, transition_name)
               for state, score, transition_name in transitions)


def compute_gene_lengths(bgc_predictions: List[BGC_Module_Prediction]) -> Dict[str, int]:
    return {gene_id: len(modules)
            for gene_id, modules in groupby(bgc_predictions, lambda module: module.gene_id)}


def calculate_dp(assembly_line: List[BGC_Module_Prediction],
                 nrp_monomers: List[NRP_Monomer],
                 dp_helper: DP_Helper) -> DPTable:  # functions computing scores and other parameters
    gene_lengths = compute_gene_lengths(assembly_line)
    max_gene_reps = dp_helper.dp_config.max_gene_reps if any(module.iterative_gene for module in assembly_line) \
        else 1
    max_module_reps = dp_helper.dp_config.max_module_reps if any(module.iterative_module for module in assembly_line) \
        else 1

    dp_table = np.empty((len(assembly_line), len(nrp_monomers), max_gene_reps, max_module_reps),
                        dtype=DP_Value)

    pred_dummy = BGC_Module_Prediction('', -1, {}, (), False, False)  # placeholders for more concise code
    mon_dummy = NRP_Monomer('', (), Chirality.UNKNOWN, '')  # dummies are not actually used in the alignment

    recalc = partial(dp_recalc, dp_table)
    for gene_reps in range(max_gene_reps + 1):
        for module_reps in range(max_module_reps + 1):
            for i, bgc_pred in enumerate([pred_dummy] + assembly_line):
                for j, nrp_mon in enumerate([mon_dummy] + nrp_monomers):
                    if j == 0 and i == 0:
                        continue

                    transitions = [(DPState(i-1, j, gene_reps, module_reps), dp_helper.bgc_module_remove, (bgc_pred,)),
                                   (DPState(i, j-1, gene_reps, module_reps), dp_helper.nrp_mon_remove, (nrp_mon,)),
                                   (DPState(i-1, j-1, gene_reps, module_reps), dp_helper.match, (bgc_pred, nrp_mon))]

                    if i < len(assembly_line) and assembly_line[i + 1].iterative_module:
                        transitions.append((DPState(i+1, j, gene_reps, module_reps-1),
                                            dp_helper.iterate_module, ()))

                    if i < len(assembly_line) and assembly_line[i + 1].iterative_gene:
                        gene_len = gene_lengths[assembly_line[i + 1].gene_id]
                        transitions.append((DPState(i+gene_len, j, gene_reps-1, module_reps),
                                            dp_helper.iterate_gene, ()))
                    dp_table[DPState(i, j, gene_reps, module_reps)] = recalc((dp_state, transition(*args), transition.__name__)
                                                                             for dp_state, transition, args in transitions
                                                                             if valid_state(dp_state))

    return dp_table


def retrieve_alignment(dp_table: DPTable, state: DPState) -> Alignment:
    if state == DPState(0, 0, 0, 0):
        return []

    parent = dp_table[state].parent
    bgc_module_pos = parent.module_pos if parent.module_pos < state.module_pos else None
    nrp_monomer_pos = parent.monomer_pos if parent.monomer_pos < state.monomer_pos else None
    score = dp_table[state].score - dp_table[parent].score
    return (retrieve_alignment(dp_table, parent)
            + [AlignmentStep(bgc_module_pos=bgc_module_pos,
                             nrp_monomer_pos=nrp_monomer_pos,
                             score=score,
                             action=dp_table[state].action)])


def get_alignment(assembly_line: List[BGC_Module_Prediction],
                  nrp_monomers: List[NRP_Monomer],
                  dp_helper: DP_Helper) -> Alignment:
    dp_table = calculate_dp(assembly_line, nrp_monomers, dp_helper)
    final_state = max([DPState(len(assembly_line), len(nrp_monomers), gene_reps, module_reps)
                       for gene_reps in range(dp_table.shape[2])
                       for module_reps in range(dp_table.shape[3])],
                      key=lambda dp_state: dp_table[dp_state])
    return retrieve_alignment(dp_table, final_state)
