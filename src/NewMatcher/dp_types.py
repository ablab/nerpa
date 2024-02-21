from typing import Any, NamedTuple, List, Tuple
from src.data_types import LogProb


DPTable = Any  # TODO: use proper numpy.typing types


class DPState(NamedTuple):
    num_gene_reps: int
    num_module_reps: int
    module_pos: int
    monomer_pos: int


def valid_state(state: Tuple[int, ...]) -> bool:
    return all(idx >=0 for idx in state)


class DPValue(NamedTuple):
    score: LogProb
    parent: DPState
    action: str  # TODO: make it a proper type


def dp_recalc(dp_table: DPTable,
              transitions: List[Tuple[DPState, LogProb, str]]) -> Tuple[LogProb, DPState, str]:
    return max((dp_table[state] + score, state, transition_name)
               for state, score, transition_name in transitions)
