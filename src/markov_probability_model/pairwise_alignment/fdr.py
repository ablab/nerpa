import pandas as pd
import matplotlib.pyplot as plt
import os

from src.markov_probability_model.pairwise_alignment.sequence_aligner import ScoredPairwiseAlignmentOutput
from typing import List, NewType, Optional, Dict
from collections import defaultdict


class FdrParameters:
    def __init__(self, topk: int, relative_to: str,
                 pairs_df_logpath: Optional[str], best_pairs_df_logpath: Optional[str]):
        self.topk = topk
        self.relative_to = relative_to
        self.pairs_df_logpath = pairs_df_logpath
        self.best_pairs_df_logpath = best_pairs_df_logpath

    def __str__(self):
        return f'{self.relative_to}_top{self.topk}'


FdrData = NewType('FdrData', List[float])


class FdrGenerator:
    def __init__(self, alignments: List[ScoredPairwiseAlignmentOutput], fdr_parameters: List[FdrParameters]):
        self._fdr_parameters: List[FdrParameters] = fdr_parameters
        self._alignments: List[ScoredPairwiseAlignmentOutput] = sorted(alignments, key=lambda x: -x.score())

    def generate_fdr(self) -> List[FdrData]:
        return [self._generate_single_fdr(p) for p in self._fdr_parameters]

    def _generate_single_fdr(self, p: FdrParameters) -> FdrData:
        groups: Dict[List[ScoredPairwiseAlignmentOutput]] = defaultdict(list)
        for alignment in self._alignments:
            base_struct_id = alignment.aligned_sequence1.base_sequence_id if p.relative_to == 'mol' else \
                alignment.aligned_sequence2.base_sequence_id
            groups[base_struct_id].append(alignment)
        groups_list: List[List[ScoredPairwiseAlignmentOutput]] = list(groups.values())
        groups_list.sort(key=lambda x: -x[0].score())

        if p.pairs_df_logpath is not None:
            pd.DataFrame({
                'Seq1': [a.aligned_sequence1.sequence_id for g in groups_list for a in g],
                'Seq2': [a.aligned_sequence2.sequence_id for g in groups_list for a in g],
                'Score': [a.score() for g in groups_list for a in g],
                'Aligned1': [' '.join(map(str, a.aligned_sequence1.symbols)) for g in groups_list for a in g],
                'Aligned2': [' '.join(map(str, a.aligned_sequence2.symbols)) for g in groups_list for a in g],
            }).to_csv(p.pairs_df_logpath)

        correct: int = 0
        incorrect: int = 0
        fdr: List[float] = []
        result_alignments: List[ScoredPairwiseAlignmentOutput] = []
        for alignments in groups_list:
            cut_alignments: List[ScoredPairwiseAlignmentOutput] = alignments[:min(len(alignments), p.topk)]
            while len(alignments) > len(cut_alignments) and \
                    cut_alignments[-1].score() == alignments[len(cut_alignments)].score():
                cut_alignments.append(alignments[len(cut_alignments)])
            res: ScoredPairwiseAlignmentOutput = cut_alignments[0]
            for a in cut_alignments:
                if a.aligned_sequence1.base_sequence_id == a.aligned_sequence2.base_sequence_id:
                    res = a
                    break
            result_alignments.append(res)
            if res.aligned_sequence1.base_sequence_id == res.aligned_sequence2.base_sequence_id:
                correct += 1
            else:
                incorrect += 1
            fdr.append(incorrect / (correct + incorrect))

        if p.best_pairs_df_logpath is not None:
            pd.DataFrame({
                'Seq1': [a.aligned_sequence1.sequence_id for a in result_alignments],
                'Seq2': [a.aligned_sequence2.sequence_id for a in result_alignments],
                'Score': [a.score() for a in result_alignments],
                'Aligned1': [' '.join(map(str, a.aligned_sequence1.symbols)) for a in result_alignments],
                'Aligned2': [' '.join(map(str, a.aligned_sequence2.symbols)) for a in result_alignments],
            }).to_csv(p.best_pairs_df_logpath)

        return FdrData(fdr)


def plot_fdrs(fdr_parameters: List[FdrParameters], fdrs: List[Dict[str, FdrData]], save_dir: str):
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    for i, p in enumerate(fdr_parameters):
        plt.title(f'FDR_{str(p)}')
        legends: List[str] = []
        for name, others_fdr in fdrs[i].items():
            plt.plot(range(len(others_fdr)), others_fdr)
            legends.append(name)
        plt.legend(legends)
        plt.savefig(os.path.join(save_dir, f'FDR_{str(p)}.png'))
        plt.clf()
