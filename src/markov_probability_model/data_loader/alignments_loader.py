import abc

from typing import List, TypeVar
from src.markov_probability_model.pairwise_alignment.sequence_aligner import PairwiseAlignmentOutputWithLogs
from src.markov_probability_model.data_loader.data_loader import PairwiseAlignmentDataLoader, TwoSequenceListsData
from src.markov_probability_model.base.sequence import AminoacidSequence, ScoredAminoacidSequence
from src.markov_probability_model.base.alphabet import Gap, Aminoacid, ScoredAminoacid


class AlignmentsLoader(PairwiseAlignmentDataLoader):
    @abc.abstractmethod
    def load_alignments(self) -> List[PairwiseAlignmentOutputWithLogs]:
        pass

    def load_data(self) -> TwoSequenceListsData:
        alignments = self.load_alignments()
        seqs1: List[AminoacidSequence] = []
        seqs2: List[ScoredAminoacidSequence] = []
        for alignment in alignments:
            s1: List[Aminoacid] = []
            for s in alignment.aligned_sequence1.symbols:
                if s != Gap():
                    s1.append(s)
            seqs1.append(AminoacidSequence(alignment.aligned_sequence1.sequence_id, s1))
            s2: List[ScoredAminoacid] = []
            for s in alignment.aligned_sequence2.symbols:
                if s != Gap():
                    s2.append(s)
            seqs2.append(ScoredAminoacidSequence(alignment.aligned_sequence2.sequence_id, s2))
        return TwoSequenceListsData(seqs1, seqs2)
