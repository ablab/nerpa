from typing import (
    Any,
    Dict,
    List,
    NamedTuple,
    Tuple,
    Union)

from src.data_types import (
    BGC_Module,
    Chirality,
    LogProb,
    NRP_Monomer,
)
from src.NewMatcher.dp_helper import DP_Helper
from src.NewMatcher.alignment_types import AlignmentStep, Alignment
import numpy as np
import numpy.typing as npt
from itertools import groupby
from functools import partial


class DP_State(NamedTuple):
    module_pos: int
    monomer_pos: int
    num_gene_reps: int
    num_module_reps: int


class DP_Value(NamedTuple):
    score: LogProb
    parent: DP_State
    action: str  # TODO: make it a proper type

DP_Table = npt.NDArray[Union[DP_Value, None]]

def valid_state(dp_table: DP_Table, state: DP_State) -> bool:
    return all(0 <= idx < dp_table.shape[i]
               for i, idx in enumerate(state)) and \
        dp_table[state] is not None

START_STATE = DP_State(0,0,0,0)


def dp_recalc(dp_table: DP_Table,
              transitions: List[Tuple[DP_State, LogProb, str]]) -> Union[DP_Value, None]:
    return max((DP_Value(dp_table[state].score + score, state, transition_name)
                for state, score, transition_name in transitions),
               default=None)


def calculate_dp(assembly_line: List[BGC_Module],
                 nrp_monomers: List[NRP_Monomer],
                 dp_helper: DP_Helper) -> DP_Table:  # functions computing scores and other parameters
    gene_lengths = {gene_id: len(list(modules))
                    for gene_id, modules in groupby(assembly_line, lambda module: module.gene_id)}
    max_gene_reps = dp_helper.dp_config.max_gene_reps if any(module.iterative_gene for module in assembly_line) \
        else 0
    max_module_reps = dp_helper.dp_config.max_module_reps if any(module.iterative_module for module in assembly_line) \
        else 0

    dp_table = np.empty((len(assembly_line)+1, len(nrp_monomers)+1, max_gene_reps+1, max_module_reps+1),
                        dtype=DP_Value)
    dp_table[START_STATE] = DP_Value(0, START_STATE, '')

    pred_dummy = BGC_Module('', -1, {}, (), False, False)  # placeholders for more concise code
    mon_dummy = NRP_Monomer('', (), Chirality.UNKNOWN, '', 0)  # dummies are not actually used in the alignment

    recalc = partial(dp_recalc, dp_table)
    for gene_reps in range(max_gene_reps + 1):
        for module_reps in range(max_module_reps + 1):
            for i, bgc_pred in enumerate([pred_dummy] + assembly_line):
                for j, nrp_mon in enumerate([mon_dummy] + nrp_monomers):
                    if i == 0 and j == 0:
                        continue

                    transitions = [(DP_State(i - 1, j, gene_reps, module_reps), dp_helper.bgc_module_skip, (bgc_pred,)),
                                   (DP_State(i, j - 1, gene_reps, module_reps), dp_helper.nrp_mon_skip, (nrp_mon,)),
                                   (DP_State(i - 1, j - 1, gene_reps, module_reps), dp_helper.match, (bgc_pred, nrp_mon))]

                    if i < len(assembly_line) and assembly_line[i].iterative_module:  # note that assembly_line[i] is right AFTER bgc_pred
                        transitions.append((DP_State(i + 1, j, gene_reps, module_reps - 1),
                                            dp_helper.iterate_module, ()))

                    if i < len(assembly_line) and assembly_line[i].iterative_gene:
                        gene_len = gene_lengths[assembly_line[i].gene_id]
                        transitions.append((DP_State(i + gene_len, j, gene_reps - 1, module_reps),
                                            dp_helper.iterate_gene, ()))
                    dp_table[DP_State(i, j, gene_reps, module_reps)] = recalc((dp_state, transition(*args), transition.__name__)
                                                                              for dp_state, transition, args in transitions
                                                                              if valid_state(dp_table, dp_state))

    return dp_table


def retrieve_alignment(dp_table: DP_Table, state: DP_State,
                       assembly_line: List[BGC_Module],
                       nrp_monomers: List[NRP_Monomer]) -> Alignment:
    if state == START_STATE:
        return []

    parent = dp_table[state].parent
    bgc_module = assembly_line[parent.module_pos] if parent.module_pos < state.module_pos else None
    nrp_monomer = nrp_monomers[parent.monomer_pos] if parent.monomer_pos < state.monomer_pos else None
    score = dp_table[state].score - dp_table[parent].score
    return (retrieve_alignment(dp_table, parent, assembly_line, nrp_monomers)
            + [AlignmentStep(bgc_module=bgc_module,
                             nrp_monomer=nrp_monomer,
                             score=score,
                             action=dp_table[state].action)])


def get_alignment(assembly_line: List[BGC_Module],
                  nrp_monomers: List[NRP_Monomer],
                  dp_helper: DP_Helper) -> Alignment:
    dp_table = calculate_dp(assembly_line, nrp_monomers, dp_helper)

    def last_state(gene_reps: int, module_reps: int) -> DP_State:
        return DP_State(len(assembly_line), len(nrp_monomers), gene_reps, module_reps)

    final_state = max([last_state(gene_reps, module_reps)
                       for gene_reps in range(dp_table.shape[2])
                       for module_reps in range(dp_table.shape[3])
                       if dp_table[last_state(gene_reps, module_reps)] is not None],
                      key=lambda dp_state: dp_table[dp_state])
    return retrieve_alignment(dp_table, final_state,
                              assembly_line, nrp_monomers)
