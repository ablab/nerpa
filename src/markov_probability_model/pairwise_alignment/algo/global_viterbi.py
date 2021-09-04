from src.markov_probability_model.pairwise_alignment.sequence_aligner import ScoredPairwiseAlignmentOutput, \
    PairwiseAlignmentOutputWithLogs
from src.markov_probability_model.base.sequence import AminoacidSequence, ScoredAminoacidSequence, \
    AlignedAminoacidSequence, \
    AlignedScoredAminoacidSequence
from src.markov_probability_model.hmm.pairwise_alignment_hmm import PairwiseAlignmentHmm
from src.markov_probability_model.pairwise_alignment.sequence_aligner import PairwiseSequenceAligner
from src.markov_probability_model.pairwise_alignment.score_augmentations import ScoreAugmentator
from src.markov_probability_model.pairwise_alignment.algo.utils import calculate_log_alpha
from src.markov_probability_model.pairwise_alignment.algo.viterbi import Viterbi, ViterbiOutput
from src.markov_probability_model.base.utils import log_add_exp


class GlobalViterbiOutput(PairwiseAlignmentOutputWithLogs, ScoredPairwiseAlignmentOutput):
    def __init__(self, aligned_sequence1: AlignedAminoacidSequence, aligned_sequence2: AlignedScoredAminoacidSequence,
                 global_viterbi_score: float, log: str = ''):
        super(GlobalViterbiOutput, self).__init__(aligned_sequence1, aligned_sequence2, log)
        self.global_viterbi_score = global_viterbi_score

    def score(self):
        return self.global_viterbi_score


class GlobalViterbi(PairwiseSequenceAligner[GlobalViterbiOutput]):
    def __init__(self, hmm: PairwiseAlignmentHmm, sa: ScoreAugmentator):
        self._hmm = hmm
        self._sa = sa
        self._viterbi = Viterbi(hmm, sa)

    def align(self, seq1: AminoacidSequence, seq2: ScoredAminoacidSequence) -> GlobalViterbiOutput:
        hmm = self._hmm
        n, m, k = len(seq1.symbols), len(seq2.symbols), len(hmm.states)

        alpha = calculate_log_alpha(seq1, seq2, hmm)
        score = log_add_exp([alpha[n, m, x] for x in range(k)])
        o: ViterbiOutput = self._viterbi.align(seq1, seq2)

        return GlobalViterbiOutput(o.aligned_sequence1, o.aligned_sequence2,
                                   self._sa.recalculate_score(score, seq1, seq2, hmm.parameters))
