import abc

from src.markov_probability_model.hmm.pairwise_alignment_hmm import PairwiseAlignmentHMMParameters


class ParametersCalculator(abc.ABC):
    @abc.abstractmethod
    def calculate_parameters(self) -> PairwiseAlignmentHMMParameters:
        pass
