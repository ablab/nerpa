import numpy as np

from src.markov_probability_model.hmm.hmm import HMM, Prob
from src.markov_probability_model.base.alphabet import AminoacidAlphabet, ScoredAminoacidAlphabet, Symbol, Gap, \
    Aminoacid, ScoredAminoacid
from typing import Dict, NewType


class PairwiseAlignmentHMMParameters:
    def __init__(self,
                 omega_a: AminoacidAlphabet, omega_b: ScoredAminoacidAlphabet,
                 mu: np.ndarray, tau: np.ndarray,
                 p: Dict[Aminoacid, Dict[ScoredAminoacid, Prob]],
                 q_a: Dict[Aminoacid, Prob], q_b: Dict[ScoredAminoacid, Prob]):
        if mu.shape != (3,) or tau.shape != (3, 3):
            raise IndexError(f'Expected 3 HMM states for mu = {mu} and tau = {tau}')
        self.omega_a = omega_a
        self.omega_b = omega_b
        self.mu = mu
        self.tau = tau
        self.p = p
        self.q_a = q_a
        self.q_b = q_b

    def log_to(self, log_filepath):
        with open(log_filepath, 'w') as log:
            print('Omega_a: {}'.format(self.omega_a), file=log)
            print('', file=log)
            print('Omega_b: {}'.format(self.omega_b), file=log)
            print(file=log)
            print('Tau: {}'.format(self.tau), file=log)
            print('Mu: {}'.format(self.mu), file=log)
            # print(file=log)
            # print('P example:', file=log)
            # for a in [Aminoacid('arg', '@L', False)]:
            #     for b in [ScoredAminoacid('arg(100.0)', '@L', False),
            #               ScoredAminoacid('arg(80.0)', '@L', False),
            #               ScoredAminoacid('arg(70.0)', '@L', False),
            #               ScoredAminoacid('arg(60.0)', '@L', False),
            #               ScoredAminoacid('ala(100.0)', '@L', False),
            #               ScoredAminoacid('ala(80.0)', '@L', False),
            #               ScoredAminoacid('ala(70.0)', '@L', False),
            #               ScoredAminoacid('ala(60.0)', '@L', False)]:
            #         print(f'\tp[{a}][{b}] = {self.p[a][b]}', file=log)
            print(file=log)
            print('p: {}'.format(self.p), file=log)
            print(file=log)
            print('q_a: {}'.format(self.q_a), file=log)
            print('q_b: {}'.format(self.q_b), file=log)


PairwiseAlignmentHmmState = NewType('PairwiseAlignmentHmmState', str)


class PairwiseAlignmentHmmObservation:
    def __init__(self, first: Symbol, second: Symbol):
        self.first = first
        self.second = second

    def __hash__(self):
        return hash(str(self))

    def __str__(self):
        return str(self.first) + '-' + str(self.second)

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return str(self) == str(other)


class PairwiseAlignmentHmm(HMM[PairwiseAlignmentHmmState, PairwiseAlignmentHmmObservation]):
    def __init__(self, params: PairwiseAlignmentHMMParameters):
        self._params = params
        self.M = PairwiseAlignmentHmmState('M')
        self.A = PairwiseAlignmentHmmState('A')
        self.B = PairwiseAlignmentHmmState('B')
        states = [self.M, self.A, self.B]
        start_probs: Dict[PairwiseAlignmentHmmState, Prob] = dict(zip(states, list(params.mu)))
        transition_probs: Dict[PairwiseAlignmentHmmState, Dict[PairwiseAlignmentHmmState, Prob]] = {}
        for i, frm in enumerate(states):
            transition_probs[frm] = {}
            for j, to in enumerate(states):
                transition_probs[frm][to] = params.tau[i][j]
        observation_probs: Dict[PairwiseAlignmentHmmState, Dict[PairwiseAlignmentHmmObservation, Prob]] = {}
        for st in states:
            observation_probs[st] = {}
        for a in params.omega_a:
            observation_probs[self.A][PairwiseAlignmentHmmObservation(a, Gap())] = params.q_a[a]
            for b in params.omega_b:
                observation_probs[self.M][PairwiseAlignmentHmmObservation(a, b)] = params.p[a][b]
                observation_probs[self.B][PairwiseAlignmentHmmObservation(Gap(), b)] = params.q_b[b]
        super().__init__(start_probs, transition_probs, observation_probs)

    @property
    def parameters(self):
        return self._params
