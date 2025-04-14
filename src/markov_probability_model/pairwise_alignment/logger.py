import abc
import os

from src.markov_probability_model.pairwise_alignment.sequence_aligner import PairwiseAlignmentOutputWithLogs
from src.markov_probability_model.base.sequence import SequenceId
from src.markov_probability_model.base.alphabet import AlignedScoredAminoacid, Gap
from typing import TypeVar, Generic

O = TypeVar('O')


class AlignmentsLogger(Generic[O]):
    @abc.abstractmethod
    def log_alignment(self, alignment: O):
        pass


class HtmlAlignmentsLogger(AlignmentsLogger[PairwiseAlignmentOutputWithLogs]):
    def __init__(self, log_dir: str):
        self._log_dir = log_dir
        if not os.path.exists(log_dir):
            os.makedirs(log_dir)

    @staticmethod
    def _get_colour_by_score(score: float):
        if score > 90.0:
            return 'green'
        elif score > 60.0:
            return 'yellow'
        else:
            return 'red'

    def _get_string_for_output(self, s: AlignedScoredAminoacid) -> str:
        if s == Gap():
            return '<td>' + str(s) + '</td>'
        return '<td>' + str(s).split('(')[0] + \
               '<font style="background-color:%s;">%s<font></td>' % (
                   self._get_colour_by_score(s.score), s.score)

    def log_alignment(self, alignment: PairwiseAlignmentOutputWithLogs):
        id1: SequenceId = alignment.aligned_sequence1.sequence_id
        id2: SequenceId = alignment.aligned_sequence2.sequence_id
        logpath = os.path.join(self._log_dir, f'{id1}_{id2}.html')

        with open(logpath, 'w') as f:
            f.write('<html><body>')
            f.write(f'</br><h1>{id1} & {id2}</h1>')
            f.write('<table>')
            f.write('<tr>')
            for ind in range(len(alignment)):
                f.write(f'<td><font style="-color:black;">{str(alignment.aligned_sequence1.symbols[ind])}<font></td>')
            f.write('<tr/>')
            for ind in range(len(alignment)):
                f.write(self._get_string_for_output(alignment.aligned_sequence2.symbols[ind]))
            f.write('</tr>')
            f.write('</table>')
            f.write('</br>')
            f.write('<p> Logs: </p>' + ''.join('<p>&emsp;' + p + '</p>' for p in alignment.logs.split('\n')))
            f.write('</body></html>')
