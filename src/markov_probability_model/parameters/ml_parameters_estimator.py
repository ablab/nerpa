import numpy as np
import os

from src.markov_probability_model.parameters.parameters_calculator import ParametersCalculator, \
    PairwiseAlignmentHMMParameters
from src.markov_probability_model.data_loader.alignments_loader import PairwiseAlignmentOutputWithLogs
from src.markov_probability_model.data_loader.data_loader import TwoSequenceListsData
from src.markov_probability_model.base.alphabet import Gap, Aminoacid, ScoredAminoacid, AlignedAminoacid, \
    AlignedScoredAminoacid
from src.markov_probability_model.parameters.utils import get_alphabets_from_data, parse_nerpa_config, estimate_p, \
    estimate_qa_qb, array_to_dict
from typing import List, Dict, Optional, Tuple, Callable
from collections import defaultdict


class MaxLikelihoodParametersEstimator(ParametersCalculator):
    def __init__(self, alignments: List[PairwiseAlignmentOutputWithLogs],
                 data: TwoSequenceListsData, nerpa_cfg_path: str,
                 estimate_transition_probs: bool = False,
                 log_dir: Optional[str] = None):
        self._alignments = alignments
        self._data = data
        self._nerpa_cfg = parse_nerpa_config(nerpa_cfg_path)
        self._estimate_transition_probs = estimate_transition_probs
        self._log_dir = log_dir

    def calculate_parameters(self) -> PairwiseAlignmentHMMParameters:
        if self._estimate_transition_probs:
            tau, mu = estimate_tau_mu(self._alignments)
        else:
            tau, mu = get_tau_mu_from_nerpa_cfg(self._nerpa_cfg)
        omega_a, omega_b = get_alphabets_from_data(self._data, self._alignments)
        p = self._estimate_p(omega_a, omega_b)
        q_a, q_b = estimate_qa_qb(p)
        p, q_a, q_b = array_to_dict(omega_a, omega_b, p, q_a, q_b)
        res = PairwiseAlignmentHMMParameters(omega_a, omega_b, mu, tau, p, q_a, q_b)
        if self._log_dir is not None:
            res.log_to(os.path.join(self._log_dir, 'max_likelihood_parameters.txt'))
        return res

    def _estimate_p(self, omega_a: List[Aminoacid], omega_b: List[ScoredAminoacid]) -> np.ndarray:
        observations = [(alignment.aligned_sequence1.symbols[i],
                         alignment.aligned_sequence2.symbols[i])
                        for alignment in self._alignments
                        for i in range(len(alignment.aligned_sequence1))]
        # 1. Estimate g(a) = P(a | mismatch)
        g: Dict[str, float] = estimate_p_on_condition(omega_a, omega_b, observations,
                                                      condition=lambda x, y: x.name != y.name)
        # 2. Estimate f(a) = P(a | match)
        f: Dict[str, float] = estimate_p_on_condition(omega_a, omega_b, observations,
                                                      condition=lambda x, y: x.name == y.name)
        # 3. Assign p
        p = estimate_p(omega_a, omega_b, self._nerpa_cfg, f, g)
        return p


def estimate_p_on_condition(omega_a: List[Aminoacid], omega_b: List[ScoredAminoacid],
                            observations: List[Tuple[AlignedAminoacid, AlignedScoredAminoacid]],
                            condition: Callable[[Aminoacid, ScoredAminoacid], bool]) -> Dict[str, float]:
    alphabet: List[str] = [a.name for a in omega_a] + [b.name for b in omega_b]
    met = {a: 1 for a in alphabet}
    for s1, s2 in observations:
        if s1 == Gap() or s2 == Gap():
            continue
        if s1.name not in alphabet or s2.name not in alphabet:
            continue
        if condition(s1, s2):
            met[s1.name] += 1
            met[s2.name] += 1
    g: Dict[str, float] = defaultdict(float)
    for a, met_a in met.items():
        g[a] = met_a / sum([m for _, m in met.items()])
    return g


def _get_state_from_observation(symb1: AlignedAminoacid, symb2: AlignedScoredAminoacid):
    if symb1 != Gap() and symb2 != Gap():
        return 0
    elif symb2 == Gap():
        return 1
    else:
        return 2


def estimate_tau_mu(alignments: List[PairwiseAlignmentOutputWithLogs]) -> Tuple[np.ndarray, np.ndarray]:
    tau = np.ones((3, 3))
    mu = np.ones((3,))

    for alignment in alignments:
        states: List[int] = []
        for pos in range(len(alignment.aligned_sequence1)):
            symb1: AlignedAminoacid = alignment.aligned_sequence1.symbols[pos]
            symb2: AlignedScoredAminoacid = alignment.aligned_sequence2.symbols[pos]
            states.append(_get_state_from_observation(symb1, symb2))
        mu[states[0]] += 1
        for i in range(1, len(states)):
            tau[states[i - 1], states[i]] += 1

    tau = tau / tau.sum(axis=1)
    mu = tau.mean(axis=0)
    return tau, mu


def get_tau_mu_from_nerpa_cfg(nerpa_cfg: Dict) -> Tuple[np.ndarray, np.ndarray]:
    tau = np.ones((3, 3))
    tau[0, 1] = tau[1, 1] = tau[2, 1] = np.exp(nerpa_cfg['insertion'])
    tau[0, 2] = tau[1, 2] = tau[2, 2] = np.exp(nerpa_cfg['deletion'])
    tau = tau / tau.sum(axis=1)
    mu = tau.mean(axis=0)
    return tau, mu
