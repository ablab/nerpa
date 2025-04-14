from src.markov_probability_model.base.base_sequence_id_resolver import SequenceId, BaseSequenceIdResolver, \
    SimpleBaseSequenceIdResolver
from src.markov_probability_model.base.alphabet import Aminoacid, ScoredAminoacid, AlignedAminoacid, \
    AlignedScoredAminoacid
from typing import List, Optional, Generic, TypeVar

S = TypeVar('S')


class Sequence(Generic[S]):
    def __init__(self, sequence_id: SequenceId, symbols: List[S],
                 base_seq_id_resolver: Optional[BaseSequenceIdResolver] = SimpleBaseSequenceIdResolver()):
        self.sequence_id = sequence_id
        self.symbols = symbols
        self._base_seq_id_resolver = base_seq_id_resolver

    @property
    def base_sequence_id(self):
        return self._base_seq_id_resolver.resolve(self.sequence_id)

    def __len__(self):
        return len(self.symbols)


class AminoacidSequence(Sequence[Aminoacid]):
    pass


class ScoredAminoacidSequence(Sequence[ScoredAminoacid]):
    pass


class AlignedAminoacidSequence(Sequence[AlignedAminoacid]):
    pass


class AlignedScoredAminoacidSequence(Sequence[AlignedScoredAminoacid]):
    pass
