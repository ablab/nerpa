import pandas

from src.markov_probability_model.data_loader.alignments_loader import AlignmentsLoader, PairwiseAlignmentOutputWithLogs
from src.markov_probability_model.base.alphabet import Aminoacid, ScoredAminoacid, Gap
from src.markov_probability_model.base.sequence import AlignedAminoacidSequence, AlignedScoredAminoacidSequence
from typing import List, Tuple, Optional


class MibigAlignmentsLoader(AlignmentsLoader):
    def __init__(self, mibig_filepath: str):
        self._mibig_filepath = mibig_filepath

    def load_alignments(self) -> List[PairwiseAlignmentOutputWithLogs]:
        mibig = pandas.read_csv(self._mibig_filepath)
        outp: List[PairwiseAlignmentOutputWithLogs] = []
        for bgc, alignment_data in mibig.groupby(mibig['BGC']):
            sequence, prediction = [], []
            for seq_symbol, base_seq_symbol, pred_symbols, e_domain, m_domain in zip(alignment_data['rBan AA-ID'],
                                                                                     alignment_data['rBan AA'],
                                                                                     alignment_data['PRED_TOP5'],
                                                                                     alignment_data['L-/D- (E domain)'],
                                                                                     alignment_data['M domain']):
                if pandas.isna(seq_symbol) or seq_symbol == '-':
                    s1 = '-'
                    sequence.append(Gap())
                else:
                    s1, modification1, methylation1 = _simplify_seq_symbol(seq_symbol, base_seq_symbol)
                    sequence.append(Aminoacid(s1, modification1, methylation1))

                if pandas.isna(pred_symbols) or pred_symbols == '-':
                    prediction.append(Gap())
                else:
                    s2 = _choose_pred_symbol(pred_symbols, s1)
                    modification2: str = None if pandas.isna(e_domain) or e_domain == '-' else '@' + e_domain
                    methylation2: bool = False if pandas.isna(m_domain) or m_domain == '-' else (
                                m_domain.lower() == 'true')
                    prediction.append(ScoredAminoacid(s2, modification2, methylation2))

            outp.append(PairwiseAlignmentOutputWithLogs(
                AlignedAminoacidSequence(bgc, sequence), AlignedScoredAminoacidSequence(bgc, prediction), logs=''))
        return outp


def _simplify_base_seq_symbol(s: str) -> str:
    return s.split('+')[0]


def _simplify_seq_symbol(s: str, base_s: str) -> Tuple[str, Optional[str], bool]:
    modification = None
    if '@L' in s:
        modification = '@L'
    elif '@D' in s:
        modification = '@D'
    return _simplify_base_seq_symbol(base_s), modification, ('NMe' in s)


def _choose_pred_symbol(pred_symbols, seq_symbol) -> str:
    if pandas.isna(pred_symbols):
        return '-'
    pred_symbols = pred_symbols.split(';')
    base_symbols, scores = zip(*[(s.split('(')[0], float(s.split('(')[1].split(')')[0])) for s in pred_symbols])
    max_score_base_symbols = [base_symbols[i] for i in range(len(base_symbols)) if scores[i] == scores[0]]
    if seq_symbol in max_score_base_symbols:
        return pred_symbols[max_score_base_symbols.index(seq_symbol)]
    return pred_symbols[0]
