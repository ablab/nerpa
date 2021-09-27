import numpy as np
import os

from src.markov_probability_model.base.alphabet import Aminoacid, ScoredAminoacid, AlignedAminoacid, \
    AlignedScoredAminoacid, Gap
from src.markov_probability_model.parameters.utils import same_modifications, estimate_p_with_modifications
from src.markov_probability_model.pairwise_alignment.sequence_aligner import PairwiseAlignmentOutputWithLogs
from src.markov_probability_model.data_loader.data_loader import TwoSequenceListsData
from src.markov_probability_model.hmm.pairwise_alignment_hmm import PairwiseAlignmentHMMParameters
from src.markov_probability_model.parameters.utils import parse_nerpa_config, get_alphabets_from_data, estimate_qa_qb, \
    array_to_dict, get_names_alphabet, generate_p_score, generate_p_mods
from typing import List, Optional, Tuple, Dict
from collections import defaultdict


class MaxLikelihoodParametersEstimator:
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
        # 1. Estimate g(a, mod, meth) = P(a, mod, meth | match)
        # 2. Estimate f(a, mod1, meth1, mod2, meth2) = P(a, mod1, meth1, mod2, meth2 | mismatch)
        f, g = {}, {}
        for modification1 in ['@L', '@D']:
            for methylation1 in [True, False]:
                g[(modification1, methylation1)] = \
                    estimate_g(omega_a, omega_b, observations, modification1, methylation1)
                for modification2 in ['@L', '@D']:
                    for methylation2 in [True, False]:
                        f[(modification1, methylation1, modification2, methylation2)] = \
                            estimate_f(omega_a, omega_b, observations,
                                       modification1, methylation1, modification2, methylation2)
        # 3. Assign p
        p_score = generate_p_score(self._nerpa_cfg, self._data)
        p_mods = generate_p_mods(self._data)
        p = estimate_p_with_modifications(omega_a, omega_b, self._nerpa_cfg, f, g, p_score, p_mods)
        return p


def estimate_f(omega_a: List[Aminoacid], omega_b: List[ScoredAminoacid],
               observations: List[Tuple[AlignedAminoacid, AlignedScoredAminoacid]],
               modification1, methylation1, modification2, methylation2) -> Dict[str, float]:
    alphabet: List[str] = get_names_alphabet(omega_a, omega_b)
    met = {a: 1 for a in alphabet}
    div_term = len(alphabet)
    for s1, s2 in observations:
        if s1 == Gap() or s2 == Gap():
            continue
        if s1.name not in alphabet or s2.name not in alphabet:
            continue
        if s1.name != s2.name:
            continue
        if same_modifications(s1, modification1, methylation1) and same_modifications(s2, modification2, methylation2):
            met[s1.name] += 1
        div_term += 1
    f: Dict[str, float] = defaultdict(float)
    for a, met_a in met.items():
        f[a] = met_a / div_term
    return f


def estimate_g(omega_a: List[Aminoacid], omega_b: List[ScoredAminoacid],
               observations: List[Tuple[AlignedAminoacid, AlignedScoredAminoacid]],
               modification, methylation) -> Dict[str, float]:
    alphabet: List[str] = get_names_alphabet(omega_a, omega_b)
    met = {a: 1 for a in alphabet}
    div_term = len(alphabet)
    for s1, s2 in observations:
        if s1 == Gap() or s2 == Gap():
            continue
        if s1.name not in alphabet or s2.name not in alphabet:
            continue
        if s1.name == s2.name:
            continue
        if same_modifications(s1, modification, methylation):
            met[s1.name] += 1
        if same_modifications(s2, modification, methylation):
            met[s2.name] += 1
        div_term += 1
    g: Dict[str, float] = defaultdict(float)
    for a, met_a in met.items():
        g[a] = met_a / div_term
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
