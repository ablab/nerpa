import pickle
import os

from src.markov_probability_model.data_loader.data_loader import PairwiseAlignmentDataLoader, TwoSequenceListsData
from src.markov_probability_model.base.sequence import AminoacidSequence, ScoredAminoacidSequence, SequenceId
from src.markov_probability_model.base.alphabet import Aminoacid, ScoredAminoacid
from typing import List


class PickleDataLoader(PairwiseAlignmentDataLoader):
    def __init__(self, data_dir: str):
        self._data_dir = data_dir

    def load_data(self) -> TwoSequenceListsData:
        structures = self._load_sequences1()
        predictions = self._load_sequences2()
        return TwoSequenceListsData(structures, predictions)

    def _load_sequences1(self) -> List[AminoacidSequence]:
        structs = pickle.load(open(os.path.join(self._data_dir, 'structures.pickle'), 'rb'))
        structures: List[AminoacidSequence] = []
        for bgc, struct in structs.items():
            structures.append(AminoacidSequence(bgc, [Aminoacid(s.lower(), None, False) for s in struct]))
        return structures

    def _load_sequences2(self) -> List[ScoredAminoacidSequence]:
        pred = pickle.load(open(os.path.join(self._data_dir, 'pred_score.pickle'), 'rb'))
        predictions: List[ScoredAminoacidSequence] = []
        for bgc, prs in pred.items():
            for i, pr in enumerate(prs):
                predictions.append(ScoredAminoacidSequence(
                    SequenceId(f'{bgc}_{i}'), [ScoredAminoacid(s.lower(), None, False) for s in pr]))
        return predictions
