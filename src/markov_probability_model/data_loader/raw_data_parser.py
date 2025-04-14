import os
import re
import sys

from src.markov_probability_model.data_loader.data_loader import PairwiseAlignmentDataLoader, TwoSequenceListsData
from src.markov_probability_model.base.sequence import AminoacidSequence, ScoredAminoacidSequence, SequenceId
from src.markov_probability_model.base.alphabet import Aminoacid, ScoredAminoacid, AminoacidAlphabet
from typing import List, Tuple, Optional


class RawDataParser(PairwiseAlignmentDataLoader):
    def __init__(self, data_dir: str, old_omega_a: AminoacidAlphabet):
        self._data_dir = data_dir
        self._old_omega_a = old_omega_a

    def load_data(self) -> TwoSequenceListsData:
        structures = self._load_sequences1()
        predictions = self._load_sequences2()
        return TwoSequenceListsData(structures, predictions)

    def _load_sequences1(self) -> List[AminoacidSequence]:
        structures: List[AminoacidSequence] = []
        with open(os.path.join(self._data_dir, 'structures.info'), 'r') as f:
            lines = [l.rstrip() for l in f.readlines() if len(l) > 0]
            for l in lines:
                bgc_name = SequenceId(l.split()[0])
                bgc_symbols = l.split()[1].split(';')[0].split(',')
                bgc_symbols = list(map(lambda x: _parse_symbol(x, self._old_omega_a), bgc_symbols))
                structures.append(
                    AminoacidSequence(bgc_name,
                                      [Aminoacid(s, modification, methylation) for s, modification, methylation in
                                       bgc_symbols]))
        return structures

    def _load_sequences2(self) -> List[ScoredAminoacidSequence]:
        pred: List[ScoredAminoacidSequence] = []
        path = os.path.join(self._data_dir, 'predictions')
        for subdir, _, files in os.walk(path):
            if subdir != path:
                continue
            for file in files:
                pred_name = SequenceId(file)
                with open(os.path.join(subdir, file), 'r') as f:
                    lines = f.readlines()
                pred_symbols = []
                for l in lines:
                    l_split = l.split()
                    base_name = l_split[2].split(';')[0]
                    modification = '@D' if 'd-' in l_split[1].lower() else '@L'
                    methylation = ('+mt' in l_split[1].lower())
                    pred_symbols.append(ScoredAminoacid(base_name, modification, methylation))
                pred.append(ScoredAminoacidSequence(pred_name, pred_symbols))
        return pred


def _parse_base_symbol(y: str, old_omega_a: AminoacidAlphabet) -> str:
    x = y.lower()
    for p in sorted(re.findall(r'[a-z]+', x), key=lambda a: -len(a)):
        if p in [a.name for a in old_omega_a] and p not in ['dhpg']:
            return p
    if x[-2:] in ['x0', 'x1', 'x2', 'x3', 'x4']:
        return 'none'
    if x[-4:] in ['mabu', 'adda', 'dhpg']:
        return x.lower()[-4:-1]
    if x[-4:] in ['bala', 'bphe', 'corn'] or x[-5:] in ['dhabu', 'bhend', 'dhcys']:
        return x.lower()[-3:]
    if x[-5:] in ['valol']:
        return 'val'
    res = sorted(re.findall(r'[a-z]+', x), key=lambda a: -len(a))[0]
    if res in ['put', 'end', 'hva', 'dbu', 'piz', 'hmp', 'hpr', 'dov', 'pha', 'suc', 'eta', 'hyv', 'dpb', 'ahp',
               'uda', 'cha', 'bht', 'dhb', 'cit', 'hty', 'tcl', 'met']:
        return res
    sys.stderr.write(f'No alphabet symbol for {x}, returning none\n')
    return 'none'


def _parse_symbol(y: str, old_omega_a: AminoacidAlphabet) -> Tuple[str, Optional[str], bool]:
    name = _parse_base_symbol(y, old_omega_a)
    modification: Optional[str] = None
    if '@L' in y:
        modification = '@L'
    elif '@D' in y:
        modification = '@D'
    methylation: bool = ('NMe' in y) or len(re.findall(r'[0-9]Me', y)) > 0
    return name, modification, methylation
