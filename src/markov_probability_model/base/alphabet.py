import abc
from typing import List, NewType, TypeVar, Optional


class Symbol(abc.ABC):
    def __init__(self, name: str, modification: Optional[str], methylation: bool):
        self.name = name
        self.modification = modification
        self.methylation = methylation

    def __str__(self) -> str:
        res = ''
        if self.modification is not None:
            res += self.modification + '-'
        if self.methylation:
            res += 'NMe-'
        res += self.name
        return res

    def __eq__(self, other):
        return str(self) == str(other)

    def __hash__(self):
        return hash(str(self))

    def __repr__(self):
        return str(self)


class Gap(Symbol):
    def __init__(self):
        super().__init__('-', None, False)


class Aminoacid(Symbol):
    pass


class ScoredAminoacid(Symbol):
    def __init__(self, init_name: str, modification: Optional[str], methylation: bool):
        super().__init__(init_name.split('(')[0], modification, methylation)
        self._score_str: str = init_name.split('(')[1].split(')')[0]
        self.score: float = float(self._score_str)

    def __str__(self):
        return super().__str__() + '(' + self._score_str + ')'


AlignedAminoacid = TypeVar('AlignedAminoacid', Aminoacid, Gap)

AlignedScoredAminoacid = TypeVar('AlignedScoredAminoacid', ScoredAminoacid, Gap)

Alphabet = NewType('Alphabet', List[Symbol])

AminoacidAlphabet = NewType('AminoacidAlphabet', List[Aminoacid])

ScoredAminoacidAlphabet = NewType('ScoredAminoacidAlphabet', List[ScoredAminoacid])
