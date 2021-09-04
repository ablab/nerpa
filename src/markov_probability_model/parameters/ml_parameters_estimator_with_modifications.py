import numpy as np

from src.markov_probability_model.parameters.ml_parameters_estimator import MaxLikelihoodParametersEstimator, \
    estimate_p_on_condition
from src.markov_probability_model.base.alphabet import Aminoacid, ScoredAminoacid
from src.markov_probability_model.parameters.utils import same_modifications_methylations, estimate_p_with_modifications
from typing import List


class MaxLikelihoodParametersEstimatorWithModifications(MaxLikelihoodParametersEstimator):
    def _estimate_p(self, omega_a: List[Aminoacid], omega_b: List[ScoredAminoacid]) -> np.ndarray:
        observations = [(alignment.aligned_sequence1.symbols[i],
                         alignment.aligned_sequence2.symbols[i])
                        for alignment in self._alignments
                        for i in range(len(alignment.aligned_sequence1))]
        # 1. Estimate g(a) = P(a | mismatch, mod1, meth1, mod2, meth2)
        # 2. Estimate f(a) = P(a | match,    mod1, meth1, mod2, meth2)
        f, g = {}, {}
        for modification1 in ['None', '@L', '@D']:
            for methylation1 in [True, False]:
                for modification2 in ['None', '@L', '@D']:
                    for methylation2 in [True, False]:
                        condition = lambda x, y: same_modifications_methylations(x, modification1, methylation1) and \
                                                 same_modifications_methylations(y, modification2, methylation2)
                        g[(modification1, methylation1, modification2, methylation2)] = \
                            estimate_p_on_condition(
                                omega_a, omega_b, observations,
                                condition=lambda x, y: x.name != y.name and condition(x, y))
                        f[(modification1, methylation1, modification2, methylation2)] = \
                            estimate_p_on_condition(
                                omega_a, omega_b, observations,
                                condition=lambda x, y: x.name == y.name and condition(x, y))
        # 3. Assign p
        p = estimate_p_with_modifications(omega_a, omega_b, self._nerpa_cfg, f, g)
        return p
