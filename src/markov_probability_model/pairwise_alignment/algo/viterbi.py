import numpy as np
import copy

from src.markov_probability_model.pairwise_alignment.sequence_aligner import ScoredPairwiseAlignmentOutput, \
    PairwiseAlignmentOutputWithLogs
from src.markov_probability_model.base.sequence import AminoacidSequence, ScoredAminoacidSequence, \
    AlignedAminoacidSequence, \
    AlignedScoredAminoacidSequence
from src.markov_probability_model.base.alphabet import Gap, Symbol, AlignedAminoacid, AlignedScoredAminoacid
from src.markov_probability_model.hmm.pairwise_alignment_hmm import PairwiseAlignmentHmm, \
    PairwiseAlignmentHmmObservation, PairwiseAlignmentHmmState
from src.markov_probability_model.pairwise_alignment.sequence_aligner import PairwiseSequenceAligner
from src.markov_probability_model.pairwise_alignment.score_augmentations import ScoreAugmentator
from src.markov_probability_model.base.utils import my_log
from typing import List


class ViterbiOutput(PairwiseAlignmentOutputWithLogs, ScoredPairwiseAlignmentOutput):
    def __init__(self, aligned_sequence1: AlignedAminoacidSequence, aligned_sequence2: AlignedScoredAminoacidSequence,
                 viterbi_score: float, log: str):
        super(ViterbiOutput, self).__init__(aligned_sequence1, aligned_sequence2, log)
        self.viterbi_score = viterbi_score

    def score(self):
        return self.viterbi_score


class Viterbi(PairwiseSequenceAligner[ViterbiOutput]):
    def __init__(self, hmm: PairwiseAlignmentHmm, sa: ScoreAugmentator):
        self._hmm = hmm
        self._sa = sa

    def align(self, seq1: AminoacidSequence, seq2: ScoredAminoacidSequence) -> ViterbiOutput:
        hmm = self._hmm
        n, m, k = len(seq1.symbols), len(seq2.symbols), len(hmm.states)

        f = np.full((n + 1, m + 1, k), None)
        prev_state = np.full((n + 1, m + 1, k), None)

        d = {hmm.M: (1, 1), hmm.A: (1, 0), hmm.B: (0, 1)}

        for st in hmm.states:
            i, j = d[st]
            s1 = seq1.symbols[0] if i == 1 else Gap()
            s2 = seq2.symbols[0] if j == 1 else Gap()
            f[i, j, hmm.state_index(st)] = \
                my_log(hmm.start_prob(st)) + my_log(hmm.observation_prob(st, PairwiseAlignmentHmmObservation(s1, s2)))

        for i in range(n + 1):
            for j in range(m + 1):
                for st in hmm.states:
                    di, dj = d[st]
                    if i < di or j < dj:
                        continue
                    for frm in hmm.states:
                        s1 = seq1.symbols[i - 1] if di == 1 else Gap()
                        s2 = seq2.symbols[j - 1] if dj == 1 else Gap()
                        sc = my_log(hmm.transition_prob(frm, st)) + my_log(
                            hmm.observation_prob(st, PairwiseAlignmentHmmObservation(s1, s2)))
                        frm_idx = hmm.state_index(frm)
                        st_idx: int = hmm.state_index(st)
                        if f[i - di, j - dj, frm_idx] is None:
                            continue
                        if f[i, j, st_idx] is None or \
                                f[i, j, st_idx] < f[i - di, j - dj, frm_idx] + sc:
                            f[i, j, st_idx] = f[i - di, j - dj, frm_idx] + sc
                            prev_state[i, j, st_idx] = copy.deepcopy(frm)

        i, j, best_st, best_score = n, m, hmm.M, -np.inf
        for st in hmm.states:
            if f[i, j, hmm.state_index(st)] > best_score:
                best_st, best_score = st, f[i, j, hmm.state_index(st)]
        states = [best_st]
        while prev_state[i, j, hmm.state_index(best_st)] is not None:
            prev_st = prev_state[i, j, hmm.state_index(best_st)]
            states.append(prev_st)
            i -= d[best_st][0]
            j -= d[best_st][1]
            best_st = prev_st

        states.reverse()
        assert i <= 1 and j <= 1, "Internal error"

        def _log_observation_score(state: PairwiseAlignmentHmmState, symb1: Symbol, symb2: Symbol):
            name = 'P' if state == hmm.M else 'Q_a' if state == hmm.A else 'Q_b'
            return 'log({}({}, {})) = log({}) = {}'.format(
                name, symb1, symb2, hmm.observation_prob(state, PairwiseAlignmentHmmObservation(symb1, symb2)),
                my_log(hmm.observation_prob(state, PairwiseAlignmentHmmObservation(symb1, symb2))))

        def _log_start_score(state: PairwiseAlignmentHmmState):
            return 'log(Mu({})) = log({}) = {}'.format(state, hmm.start_prob(state),
                                                       my_log(hmm.start_prob(state)))

        def _log_transition_score(prev: PairwiseAlignmentHmmState, nxt: PairwiseAlignmentHmmState):
            return 'log(Tau({} -> {})) = log({}) = {}'.format(prev, nxt, hmm.transition_prob(prev, nxt),
                                                              my_log(hmm.transition_prob(prev, nxt)))

        i, j = 0, 0
        aligned1: List[AlignedAminoacid] = []
        aligned2: List[AlignedScoredAminoacid] = []
        logs = []
        for prev_state, st in zip([None] + states[:-1], states):
            s1 = Gap() if d[st][0] == 0 else seq1.symbols[i]
            s2 = Gap() if d[st][1] == 0 else seq2.symbols[j]
            if prev_state is None:
                logs.append(_log_start_score(st))
            else:
                logs.append(_log_transition_score(prev_state, st))
            logs.append(_log_observation_score(st, s1, s2))
            aligned1.append(s1)
            aligned2.append(s2)
            i += d[st][0]
            j += d[st][1]

        return ViterbiOutput(AlignedAminoacidSequence(seq1.sequence_id, aligned1),
                             AlignedScoredAminoacidSequence(seq2.sequence_id, aligned2),
                             self._sa.recalculate_score(best_score, seq1, seq2, hmm.parameters),
                             log='\n'.join(logs))
