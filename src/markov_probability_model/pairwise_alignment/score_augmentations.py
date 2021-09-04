import abc

from src.markov_probability_model.base.utils import my_log
from src.markov_probability_model.hmm.pairwise_alignment_hmm import PairwiseAlignmentHMMParameters
from src.markov_probability_model.base.sequence import AminoacidSequence, ScoredAminoacidSequence


class ScoreAugmentator(abc.ABC):
    @abc.abstractmethod
    def recalculate_score(self, score: float, seq1: AminoacidSequence, seq2: ScoredAminoacidSequence,
                          params: PairwiseAlignmentHMMParameters) -> float:
        pass


class IdentityScoreAugmentator(ScoreAugmentator):
    def recalculate_score(self, score: float, seq1: AminoacidSequence, seq2: ScoredAminoacidSequence,
                          params: PairwiseAlignmentHMMParameters) -> float:
        return score


class NullHypothesisScoreAugmentator(ScoreAugmentator):
    def recalculate_score(self, score: float, seq1: AminoacidSequence, seq2: ScoredAminoacidSequence,
                          params: PairwiseAlignmentHMMParameters) -> float:
        score -= sum(my_log(params.q_a[a]) for a in seq1.symbols) + sum(my_log(params.q_b[b]) for b in seq2.symbols)
        return score
