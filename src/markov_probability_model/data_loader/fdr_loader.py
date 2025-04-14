import abc
import pandas as pd
import os
import copy

from src.markov_probability_model.pairwise_alignment.fdr import FdrData, FdrParameters, FdrGenerator
from src.markov_probability_model.pairwise_alignment.sequence_aligner import ScoredPairwiseAlignmentOutput
from src.markov_probability_model.base.sequence import AlignedAminoacidSequence, AlignedScoredAminoacidSequence
from src.markov_probability_model.base.base_sequence_id_resolver import SequenceId
from typing import List, Dict


class FdrLoader(abc.ABC):
    @abc.abstractmethod
    def load_fdr(self) -> List[Dict[str, FdrData]]:
        pass


class FdrGeneratorFromReport(FdrLoader):
    def __init__(self, data_dir: str, fdr_parameters: List[FdrParameters]):
        self._data_dir = data_dir
        self._fdr_parameters = copy.deepcopy(fdr_parameters)
        for f in self._fdr_parameters:
            f.pairs_df_logpath = None
            f.best_pairs_df_logpath = None

    def load_fdr(self) -> List[Dict[str, FdrData]]:
        report = pd.read_csv(os.path.join(self._data_dir, 'report.csv'))
        alignments: List[ScoredPairwiseAlignmentOutput] = [
            self.NerpaFdrAlignment(score, seq1, seq2.split('/')[-1])
            for score, seq1, seq2 in zip(report['score'], report['mol id'], report['prediction id'])
        ]
        fdrs = FdrGenerator(alignments, self._fdr_parameters).generate_fdr()
        return [{'NERPA': f} for f in fdrs]

    class NerpaFdrAlignment(ScoredPairwiseAlignmentOutput):
        def __init__(self, score: float, sequence_id1: SequenceId, sequence_id2: SequenceId):
            super().__init__(AlignedAminoacidSequence(sequence_id1, []),
                             AlignedScoredAminoacidSequence(sequence_id2, []))
            self._score = score

        def score(self) -> float:
            return self._score


class CsvFdrLoader(FdrLoader):
    def __init__(self, data_dir: str, fdr_parameters: List[FdrParameters]):
        self._data_dir = data_dir
        self._fdr_parameters = fdr_parameters

    def load_fdr(self) -> List[Dict[str, FdrData]]:
        return [self._load_single_fdr(p) for p in self._fdr_parameters]

    def _load_single_fdr(self, p: FdrParameters) -> Dict[str, FdrData]:
        t = pd.read_csv(os.path.join(
            self._data_dir, f'FDR_top{p.topk}_{p.relative_to}_scores_with_garlic.csv'))
        return {'NERPA': FdrData(_fdr_row_to_array(t['FDR Nerpa'])),
                'GARLIC': FdrData(_fdr_row_to_array(t['FDR Garlic']))}


def _fdr_row_to_array(row, max_len=500) -> List[float]:
    a = []
    for c in row:
        if c == '-':
            break
        a.append(float(c))
    if len(a) > max_len:
        a = a[:max_len]
    return list(a)
