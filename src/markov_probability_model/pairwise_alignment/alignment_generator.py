import abc

from typing import List, TypeVar, Generic
from src.markov_probability_model.pairwise_alignment.sequence_aligner import PairwiseAlignmentOutputWithLogs, \
    PairwiseSequenceAligner
from src.markov_probability_model.data_loader.data_loader import TwoSequenceListsData
from src.markov_probability_model.base.sequence import AminoacidSequence
from multiprocessing import Pool
from tqdm import tqdm


class AlignmentGenerator(abc.ABC):
    @abc.abstractmethod
    def generate_alignments(self) -> List[PairwiseAlignmentOutputWithLogs]:
        pass


O = TypeVar('O')


class AllPairsAlignmentGenerator(AlignmentGenerator, Generic[O]):
    def __init__(self, data: TwoSequenceListsData, sequence_aligner: PairwiseSequenceAligner[O],
                 pool_sz: int = 1):
        self._data = data
        self._sequence_aligner = sequence_aligner
        self._pool_sz = pool_sz

    def generate_alignments(self) -> List[O]:
        with Pool(self._pool_sz) as p:
            alignment_lists: List[List[O]] = list(
                tqdm(p.imap(self._generate_alignments_for_seq1, self._data.sequences1),
                     total=len(self._data.sequences1)))
        return [alignment for alignments in alignment_lists for alignment in alignments]

    def _generate_alignments_for_seq1(self, seq1: AminoacidSequence) -> List[O]:
        alignments: List[O] = []
        for seq2 in self._data.sequences2:
            alignments.append(self._sequence_aligner.align(seq1, seq2))
        return alignments
