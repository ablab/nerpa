import abc
import re
from typing import NewType, List

SequenceId = NewType('SequenceId', str)
BaseSequenceId = NewType('BaseSequenceId', str)


class BaseSequenceIdResolver(abc.ABC):
    @abc.abstractmethod
    def resolve(self, sequence_id: SequenceId) -> SequenceId:
        pass


class SimpleBaseSequenceIdResolver(BaseSequenceIdResolver):
    def __init__(self):
        self._regexps: List[str] = [
            r'BGC[0-9]{7}',
            r'NPA[0-9]{6}',
            r'[A-Z]{3}[0-9]{5}_variant',
            r'[A-Z]{3}[0-9]{2}-[A-Z]{1}_variant',
            r'antimarin[0-9]{4}_[0-9]{4,5}_variant',
            r'streptomedb.[0-9]{2,4}_variant',
            r'mibig_[0-9]{3}_variant',
        ]

    def resolve(self, sequence_id: SequenceId) -> SequenceId:
        for seq_re in self._regexps:
            matches = re.findall(seq_re, str(sequence_id))
            if len(matches) > 0:
                res = matches[0]
                if res.endswith('_variant'):
                    return res[:-len('_variant')]
                return res
        raise IndexError(f'Cannot resolve base for {sequence_id}')
