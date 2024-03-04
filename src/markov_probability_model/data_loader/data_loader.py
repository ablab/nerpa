import abc

from src.markov_probability_model.base.sequence import AminoacidSequence, ScoredAminoacidSequence
from typing import List, Generic, TypeVar

D = TypeVar('D')


class DataLoader(abc.ABC, Generic[D]):
    @abc.abstractmethod
    def load_data(self) -> D:
        pass


class TwoSequenceListsData:
    def __init__(self, sequences1: List[AminoacidSequence], sequences2: List[ScoredAminoacidSequence]):
        self.sequences1 = sequences1
        self.sequences2 = sequences2


class PairwiseAlignmentDataLoader(DataLoader[TwoSequenceListsData], abc.ABC):
    pass
