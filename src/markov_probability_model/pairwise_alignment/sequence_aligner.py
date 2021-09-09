import abc

from src.markov_probability_model.base.sequence import AlignedAminoacidSequence, AlignedScoredAminoacidSequence, \
    AminoacidSequence, ScoredAminoacidSequence
from typing import TypeVar, Generic


class PairwiseAlignmentOutput:
    def __init__(self, aligned_sequence1: AlignedAminoacidSequence,
                 aligned_sequence2: AlignedScoredAminoacidSequence):
        self.aligned_sequence1 = aligned_sequence1
        self.aligned_sequence2 = aligned_sequence2
        assert len(self.aligned_sequence1) == len(self.aligned_sequence2)

    def __len__(self):
        return len(self.aligned_sequence1)


class PairwiseAlignmentOutputWithLogs(PairwiseAlignmentOutput):
    def __init__(self, aligned_sequence1: AlignedAminoacidSequence,
                 aligned_sequence2: AlignedScoredAminoacidSequence, logs: str):
        super().__init__(aligned_sequence1, aligned_sequence2)
        self.logs = logs


class ScoredPairwiseAlignmentOutput(abc.ABC, PairwiseAlignmentOutput):
    @abc.abstractmethod
    def score(self):
        pass


O = TypeVar('O')


class PairwiseSequenceAligner(abc.ABC, Generic[O]):
    @abc.abstractmethod
    def align(self, seq1: AminoacidSequence, seq2: ScoredAminoacidSequence) -> O:
        pass
