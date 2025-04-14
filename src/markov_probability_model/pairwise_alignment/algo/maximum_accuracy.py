import numpy as np

from src.markov_probability_model.pairwise_alignment.sequence_aligner import ScoredPairwiseAlignmentOutput, \
    PairwiseAlignmentOutputWithLogs
from src.markov_probability_model.base.sequence import AminoacidSequence, ScoredAminoacidSequence, \
    AlignedAminoacidSequence, AlignedScoredAminoacidSequence
from src.markov_probability_model.hmm.pairwise_alignment_hmm import PairwiseAlignmentHmm
from src.markov_probability_model.pairwise_alignment.sequence_aligner import PairwiseSequenceAligner
from src.markov_probability_model.pairwise_alignment.score_augmentations import ScoreAugmentator
from src.markov_probability_model.pairwise_alignment.algo.utils import calculate_log_alpha, calculate_log_beta, \
    calculate_marginal_prob, log_marginal_prob_for_alignment
from src.markov_probability_model.base.alphabet import AlignedAminoacid, AlignedScoredAminoacid, Gap
from typing import List


class MaximumAccuracyOutput(PairwiseAlignmentOutputWithLogs, ScoredPairwiseAlignmentOutput):
    def __init__(self, aligned_sequence1: AlignedAminoacidSequence, aligned_sequence2: AlignedScoredAminoacidSequence,
                 mpd_score: float, logs: str):
        super(MaximumAccuracyOutput, self).__init__(aligned_sequence1, aligned_sequence2, logs)
        self.mpd_score = mpd_score

    def score(self):
        return self.mpd_score


class MaximumAccuracy(PairwiseSequenceAligner[MaximumAccuracyOutput]):
    def __init__(self, hmm: PairwiseAlignmentHmm, sa: ScoreAugmentator):
        self._hmm = hmm
        self._sa = sa

    def align(self, seq1: AminoacidSequence, seq2: ScoredAminoacidSequence) -> MaximumAccuracyOutput:
        hmm = self._hmm
        n, m, k = len(seq1.symbols), len(seq2.symbols), len(hmm.states)

        alpha = calculate_log_alpha(seq1, seq2, hmm)
        beta = calculate_log_beta(seq1, seq2, hmm)
        marginal_prob = calculate_marginal_prob(seq1, seq2, hmm, alpha, beta)

        d = {hmm.M: (1, 1), hmm.A: (1, 0), hmm.B: (0, 1)}
        max_accuracy = np.full((n + 1, m + 1, k), None)
        max_accuracy[0, 0] = 0
        prev_state = np.full((n + 1, m + 1, k), None)
        for i in range(n + 1):
            for j in range(m + 1):
                for v in range(k):
                    di, dj = d[hmm.states[v]]
                    if di > i or dj > j:
                        continue
                    for prev in range(k):
                        score = marginal_prob[i, j, v]
                        if max_accuracy[i - di, j - dj, prev] is None:
                            continue
                        new_score = score + max_accuracy[i - di, j - dj, prev]
                        if max_accuracy[i, j, v] is None or max_accuracy[i, j, v] < new_score:
                            max_accuracy[i, j, v] = new_score
                            prev_state[i, j, v] = prev

        def log_calculate_score(v, i, j):
            if v == 0:
                return 'P_M({}, {}) = {}'.format(i, j, marginal_prob[i, j, v])
            else:
                return ''

        i, j, v = n, m, 0
        max_score = max_accuracy[i, j, v]
        for tmp_v in range(k):
            if max_accuracy[i, j, tmp_v] > max_accuracy[i, j, v]:
                v = tmp_v
                max_score = max_accuracy[i, j, tmp_v]
        states = []
        logs = []
        while prev_state[i, j, v] is not None:
            states.append(v)
            logs.append(log_calculate_score(v, i, j))
            di, dj = d[hmm.states[v]]
            prev = prev_state[i, j, v]
            v, i, j = prev, i - di, j - dj
        logs.reverse()
        states.reverse()
        assert i == 0 and j == 0, 'Internal error'

        aligned1: List[AlignedAminoacid] = []
        aligned2: List[AlignedScoredAminoacid] = []
        for st in states:
            di, dj = d[hmm.states[st]]
            aligned1.append(Gap() if di == 0 else seq1.symbols[i])
            aligned2.append(Gap() if dj == 0 else seq2.symbols[j])
            i += di
            j += dj

        aligned_seq1: AlignedAminoacidSequence = AlignedAminoacidSequence(seq1.sequence_id, aligned1)
        aligned_seq2: AlignedScoredAminoacidSequence = AlignedScoredAminoacidSequence(seq2.sequence_id, aligned2)
        logs.append('Marginal probs:')
        logs.append(' '.join(map(str, log_marginal_prob_for_alignment(aligned_seq1, aligned_seq2, hmm, marginal_prob))))
        return MaximumAccuracyOutput(aligned_seq1, aligned_seq2,
                                     self._sa.recalculate_score(max_score, seq1, seq2, hmm.parameters),
                                     logs='\n'.join(logs))
