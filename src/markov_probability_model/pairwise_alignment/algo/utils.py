import numpy as np

from src.markov_probability_model.base.alphabet import Gap
from src.markov_probability_model.base.sequence import AminoacidSequence, ScoredAminoacidSequence, AlignedAminoacid, \
    AlignedScoredAminoacid, AlignedAminoacidSequence, AlignedScoredAminoacidSequence
from src.markov_probability_model.base.utils import my_log, log_add_exp, my_exp
from src.markov_probability_model.hmm.pairwise_alignment_hmm import PairwiseAlignmentHmm, \
    PairwiseAlignmentHmmObservation
from typing import Tuple, Optional, List


def _get_symbols(seq1: AminoacidSequence, seq2: ScoredAminoacidSequence, i: int, j: int, di: int, dj: int) -> \
        Tuple[AlignedAminoacid, AlignedScoredAminoacid]:
    return Gap() if di == 0 else seq1.symbols[i - 1], Gap() if dj == 0 else seq2.symbols[j - 1]


def calculate_log_alpha(seq1: AminoacidSequence, seq2: ScoredAminoacidSequence,
                        hmm: PairwiseAlignmentHmm) -> np.ndarray:
    n, m = len(seq1), len(seq2)

    alpha = np.full((n + 1, m + 1, len(hmm.states)), -np.inf)
    d = {hmm.M: (1, 1), hmm.A: (1, 0), hmm.B: (0, 1)}
    for v in hmm.states:
        di, dj = d[v]
        symb1, symb2 = _get_symbols(seq1, seq2, di, dj, di, dj)
        alpha[di, dj, hmm.state_index(v)] = my_log(hmm.start_prob(v)) + my_log(
            hmm.observation_prob(v, PairwiseAlignmentHmmObservation(symb1, symb2)))
    for i in range(n + 1):
        for j in range(m + 1):
            for v in hmm.states:
                di, dj = d[v]
                if i < di or j < dj:
                    continue
                symb1, symb2 = _get_symbols(seq1, seq2, i, j, di, dj)
                s = []
                for prev in hmm.states:
                    s.append(alpha[i - di, j - dj, hmm.state_index(prev)] + my_log(hmm.transition_prob(prev, v)))

                alpha[i, j, hmm.state_index(v)] = log_add_exp([
                    alpha[i, j, hmm.state_index(v)],
                    my_log(hmm.observation_prob(v, PairwiseAlignmentHmmObservation(symb1, symb2))) + log_add_exp(s)])
    return alpha


def calculate_log_beta(seq1: AminoacidSequence, seq2: ScoredAminoacidSequence, hmm: PairwiseAlignmentHmm) -> np.ndarray:
    n, m = len(seq1), len(seq2)

    beta = np.full((n + 1, m + 1, len(hmm.states)), -np.inf)
    for v in hmm.states:
        beta[n, m, hmm.state_index(v)] = 0
    d = {hmm.M: (1, 1), hmm.A: (1, 0), hmm.B: (0, 1)}
    for i in range(n, -1, -1):
        for j in range(m, -1, -1):
            if i == 0 and j == 0:
                continue
            for v in hmm.states:
                for nxt in hmm.states:
                    di, dj = d[nxt]
                    if i + di > n or j + dj > m:
                        continue
                    symb1, symb2 = _get_symbols(seq1, seq2, i + di, j + dj, di, dj)
                    beta[i, j, hmm.state_index(v)] = log_add_exp([
                        beta[i, j, hmm.state_index(v)],
                        beta[i + di, j + dj, hmm.state_index(nxt)] +
                        my_log(hmm.observation_prob(nxt, PairwiseAlignmentHmmObservation(symb1, symb2))) +
                        my_log(hmm.transition_prob(v, nxt))])
    return beta


def calculate_log_marginal_prob(seq1: AminoacidSequence, seq2: ScoredAminoacidSequence, hmm: PairwiseAlignmentHmm,
                                alpha: np.ndarray, beta: np.ndarray) -> np.ndarray:
    n, m, k = len(seq1), len(seq2), len(hmm.states)
    sum_alpha = log_add_exp([alpha[n, m, x] for x in range(k)])
    marginal_prob = np.zeros((n + 1, m + 1, k))
    for i in range(n + 1):
        for j in range(m + 1):
            for s in range(k):
                marginal_prob[i, j, s] = alpha[i, j, s] + beta[i, j, s] - sum_alpha
    return marginal_prob


def calculate_marginal_prob(seq1: AminoacidSequence, seq2: ScoredAminoacidSequence, hmm: PairwiseAlignmentHmm,
                            alpha: np.ndarray, beta: np.ndarray) -> np.ndarray:
    log_marginal_prob = calculate_log_marginal_prob(seq1, seq2, hmm, alpha, beta)
    n, m, k = log_marginal_prob.shape
    marginal_prob = np.array(
        [[[my_exp(log_marginal_prob[i, j, p]) for p in range(k)] for j in range(m)] for i in range(n)])
    return marginal_prob


def _get_sequence_from_aligned(sequence):
    return [a for a in sequence.symbols if a != Gap()]


def log_marginal_prob_for_alignment(sequence1: AlignedAminoacidSequence, sequence2: AlignedScoredAminoacidSequence,
                                    hmm: PairwiseAlignmentHmm,
                                    marginal_prob: Optional[np.ndarray] = None) -> List[float]:
    if marginal_prob is None:
        seq1 = AminoacidSequence(sequence1.sequence_id, _get_sequence_from_aligned(sequence1))
        seq2 = ScoredAminoacidSequence(sequence2.sequence_id, _get_sequence_from_aligned(sequence2))
        alpha = calculate_log_alpha(seq1, seq2, hmm)
        beta = calculate_log_beta(seq1, seq2, hmm)
        marginal_prob = calculate_marginal_prob(seq1, seq2, hmm, alpha, beta)

    d = {hmm.M: (1, 1), hmm.A: (1, 0), hmm.B: (0, 1)}
    i, j = 0, 0
    probs: List[float] = []
    for k in range(len(sequence1)):
        symb1: AlignedAminoacid = sequence1.symbols[k]
        symb2: AlignedScoredAminoacid = sequence2.symbols[k]
        state = hmm.A if symb2 == Gap() else hmm.B if symb1 == Gap() else hmm.M
        i += d[state][0]
        j += d[state][1]
        probs.append(marginal_prob[i, j, hmm.state_index(state)])
    return probs
