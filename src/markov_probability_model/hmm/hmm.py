from typing import TypeVar, Generic, NewType, Dict, List

S = TypeVar('S')
O = TypeVar('O')

Prob = NewType('Prob', float)


class HMM(Generic[S, O]):
    def __init__(self, start_probs: Dict[S, Prob],
                 transition_probs: Dict[S, Dict[S, Prob]],
                 observation_probs: Dict[S, Dict[O, Prob]]):
        self._start_probs = start_probs
        self._transition_probs = transition_probs
        self._observation_probs = observation_probs

    @property
    def states(self) -> List[S]:
        return list(self._start_probs.keys())

    def observations(self, s: S) -> List[O]:
        return list(self._observation_probs[s].keys())

    def start_prob(self, s: S) -> Prob:
        return self._start_probs[s]

    def transition_prob(self, frm: S, to: S) -> Prob:
        return self._transition_probs[frm][to]

    def observation_prob(self, s: S, o: O) -> Prob:
        return self._observation_probs[s][o]

    def state_index(self, s: S) -> int:
        return self.states.index(s)
