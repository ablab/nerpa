from typing import (
    Any,
    Dict,
    List,
    NamedTuple,
    Tuple,
    Union)

from src.NewMatcher.dp_types import DPValue
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


DP_Table = Any  # TODO: use proper numpy.typing types

class DP_State(NamedTuple):
    module_pos: int
    monomer_pos: int
    num_gene_reps: int
    num_module_reps: int

def valid_state(dp_table: DP_Table, state: DP_State) -> bool:
    return all(0 <= idx < dp_table.shape[i]
               for i, idx in enumerate(state)) and \
        dp_table[state] is not None

START_STATE = DP_State(0,0,0,0)

class DP_Value(NamedTuple):
    score: LogProb
    parent: DP_State
    action: str  # TODO: make it a proper type


def dp_recalc(dp_table: DP_Table,
              transitions: List[Tuple[DP_State, LogProb, str]]) -> Union[DP_Value, None]:
    return max((DP_Value(dp_table[state].score + score, state, transition_name)
                for state, score, transition_name in transitions),
               default=None)


def compute_gene_lengths(bgc_predictions: List[BGC_Module_Prediction]) -> Dict[str, int]:
    return {gene_id: len(list(modules))
            for gene_id, modules in groupby(bgc_predictions, lambda module: module.gene_id)}


def calculate_dp(assembly_line: List[BGC_Module_Prediction],
                 nrp_monomers: List[NRP_Monomer],
                 dp_helper: DP_Helper) -> DP_Table:  # functions computing scores and other parameters
    gene_lengths = compute_gene_lengths(assembly_line)
    max_gene_reps = dp_helper.dp_config.max_gene_reps if any(module.iterative_gene for module in assembly_line) \
        else 0
    max_module_reps = dp_helper.dp_config.max_module_reps if any(module.iterative_module for module in assembly_line) \
        else 0

    dp_table = np.empty((len(assembly_line)+1, len(nrp_monomers)+1, max_gene_reps+1, max_module_reps+1),
                        dtype=DP_Value)
    dp_table[START_STATE] = DP_Value(0, START_STATE, '')

    pred_dummy = BGC_Module_Prediction('', -1, {}, (), False, False)  # placeholders for more concise code
    mon_dummy = NRP_Monomer('', (), Chirality.UNKNOWN, '')  # dummies are not actually used in the alignment

    recalc = partial(dp_recalc, dp_table)
    for gene_reps in range(max_gene_reps + 1):
        for module_reps in range(max_module_reps + 1):
            for i, bgc_pred in enumerate([pred_dummy] + assembly_line):
                for j, nrp_mon in enumerate([mon_dummy] + nrp_monomers):
                    if i == 0 and j == 0:
                        continue

                    transitions = [(DP_State(i - 1, j, gene_reps, module_reps), dp_helper.bgc_module_remove, (bgc_pred,)),
                                   (DP_State(i, j - 1, gene_reps, module_reps), dp_helper.nrp_mon_remove, (nrp_mon,)),
                                   (DP_State(i - 1, j - 1, gene_reps, module_reps), dp_helper.match, (bgc_pred, nrp_mon))]

                    if i < len(assembly_line) - 1 and assembly_line[i + 1].iterative_module:
                        transitions.append((DP_State(i + 1, j, gene_reps, module_reps - 1),
                                            dp_helper.iterate_module, ()))

                    if i < len(assembly_line) - 1 and assembly_line[i + 1].iterative_gene:
                        gene_len = gene_lengths[assembly_line[i + 1].gene_id]
                        transitions.append((DP_State(i + gene_len, j, gene_reps - 1, module_reps),
                                            dp_helper.iterate_gene, ()))
                    dp_table[DP_State(i, j, gene_reps, module_reps)] = recalc((dp_state, transition(*args), transition.__name__)
                                                                              for dp_state, transition, args in transitions
                                                                              if valid_state(dp_table, dp_state))

    return dp_table


def retrieve_alignment(dp_table: DP_Table, state: DP_State) -> Alignment:
    if state == START_STATE:
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
    final_state = max([DP_State(len(assembly_line), len(nrp_monomers), gene_reps, module_reps)
                       for gene_reps in range(dp_table.shape[2])
                       for module_reps in range(dp_table.shape[3])],
                      key=lambda dp_state: dp_table[dp_state])
    return retrieve_alignment(dp_table, final_state)
