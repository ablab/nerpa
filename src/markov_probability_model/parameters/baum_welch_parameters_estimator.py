import numpy as np
import os

from src.markov_probability_model.parameters.parameters_calculator import ParametersCalculator
from src.markov_probability_model.data_loader.alignments_loader import PairwiseAlignmentOutputWithLogs
from src.markov_probability_model.data_loader.data_loader import TwoSequenceListsData
from src.markov_probability_model.base.alphabet import Gap, AminoacidAlphabet, ScoredAminoacidAlphabet, Aminoacid, \
    ScoredAminoacid
from src.markov_probability_model.base.sequence import AminoacidSequence, ScoredAminoacidSequence
from src.markov_probability_model.base.utils import my_log, my_exp, log_add_exp
from src.markov_probability_model.parameters.utils import estimate_p_with_modifications, parse_nerpa_config, \
    estimate_qa_qb, array_to_dict, \
    same_modifications_methylations
from src.markov_probability_model.pairwise_alignment.algo.utils import calculate_log_alpha, calculate_log_beta
from src.markov_probability_model.parameters.utils import get_alphabets_from_data
from src.markov_probability_model.hmm.pairwise_alignment_hmm import PairwiseAlignmentHmm, \
    PairwiseAlignmentHmmObservation, \
    PairwiseAlignmentHMMParameters
from typing import List, Optional, Dict, Callable, Tuple
from tqdm import tqdm
from multiprocessing import Pool

PSEUDOCOUNT = 1e-15


def _calculate_a(seqs1: List[AminoacidSequence], seqs2: List[ScoredAminoacidSequence],
                 gammas: List[np.ndarray], epss: List[np.ndarray], hmm: PairwiseAlignmentHmm):
    a = np.full((len(hmm.states), len(hmm.states)), 0.0)
    for i in range(len(hmm.states)):
        for j in range(len(hmm.states)):
            a[i, j] = my_exp(log_add_exp([eps[t, s, i, j] + gamma[t, s, i]
                                          for seq1, seq2, eps, gamma in zip(seqs1, seqs2, epss, gammas)
                                          for t in range(len(seq1) + 1) for s in range(len(seq2) + 1)]))
        a[i, :] /= a[i, :].sum()
    return a


def _calculate_gamma(seq1: AminoacidSequence, seq2: ScoredAminoacidSequence, hmm: PairwiseAlignmentHmm,
                     alpha: np.ndarray, beta: np.ndarray) -> np.ndarray:
    n, m = len(seq1), len(seq2)
    gamma = np.full((n + 1, m + 1, len(hmm.states)), -np.inf)
    sum_ = log_add_exp([alpha[n, m, hmm.state_index(i)] for i in hmm.states])
    for t in range(n + 1):
        for s in range(m + 1):
            for i in hmm.states:
                gamma[t, s, hmm.state_index(i)] = \
                    alpha[t, s, hmm.state_index(i)] + beta[t, s, hmm.state_index(i)] - sum_
    return gamma


def _calculate_eps(seq1: AminoacidSequence, seq2: ScoredAminoacidSequence, hmm: PairwiseAlignmentHmm,
                   alpha: np.ndarray, beta: np.ndarray) -> np.ndarray:
    d = {hmm.M: (1, 1), hmm.A: (1, 0), hmm.B: (0, 1)}
    n, m = len(seq1), len(seq2)
    eps = np.full((n + 1, m + 1, len(hmm.states), len(hmm.states)), -np.inf)
    for t in range(n + 1):
        for s in range(m + 1):
            for x in hmm.states:
                for y in hmm.states:
                    di, dj = d[y]
                    if t + di > n or s + dj > m:
                        continue
                    if beta[t, s, hmm.state_index(x)] == -np.inf:
                        continue
                    if alpha[t, s, hmm.state_index(x)] == -np.inf:
                        continue
                    symb1 = Gap() if di == 0 else seq1.symbols[t + di - 1]
                    symb2 = Gap() if dj == 0 else seq2.symbols[s + dj - 1]
                    eps[t, s, hmm.state_index(x), hmm.state_index(y)] = \
                        my_log(hmm.observation_prob(y, PairwiseAlignmentHmmObservation(symb1, symb2))) + my_log(
                            hmm.transition_prob(x, y)) + beta[t + di, s + dj, hmm.state_index(y)] - beta[
                            t, s, hmm.state_index(x)]
    return eps


def _calculate_pi(gammas: List[np.ndarray], hmm: PairwiseAlignmentHmm):
    d = {hmm.M: (1, 1), hmm.A: (1, 0), hmm.B: (0, 1)}
    pi = np.full((len(hmm.states),), 0.0)
    for v in hmm.states:
        di, dj = d[v]
        pi[hmm.state_index(v)] = my_exp(log_add_exp([
            gamma[di, dj, hmm.state_index(v)] for gamma in gammas]))
    for v in hmm.states:
        pi[hmm.state_index(v)] /= len(gammas)
    return pi


def _estimate_g(alphabet: List[str], seqs1: List[AminoacidSequence], seqs2: List[ScoredAminoacidSequence],
                gammas: List[np.ndarray], hmm: PairwiseAlignmentHmm,
                condition: Callable[[Aminoacid, ScoredAminoacid], bool]) -> Dict[str, float]:
    return {a: my_exp(
        log_add_exp([
                        gam[t, s, hmm.state_index(hmm.M)]
                        for seq1, seq2, gam in zip(seqs1, seqs2, gammas)
                        for t in range(1, len(seq1) + 1) for s in range(1, len(seq2) + 1)
                        if seq1.symbols[t - 1].name != seq2.symbols[s - 1].name and (
                    seq1.symbols[t - 1].name == a or seq2.symbols[s - 1].name == a) and condition(seq1.symbols[t - 1],
                                                                                                  seq2.symbols[s - 1])
                    ] + [my_log(PSEUDOCOUNT)]) -
        log_add_exp([
                        2 * gam[t, s, [hmm.state_index(hmm.M)]]
                        for seq1, seq2, gam in zip(seqs1, seqs2, gammas)
                        for t in range(1, len(seq1) + 1) for s in range(1, len(seq2) + 1)
                        if
                        seq1.symbols[t - 1].name != seq2.symbols[s - 1].name and condition(seq1.symbols[t - 1],
                                                                                           seq2.symbols[s - 1])
                    ] + [my_log(len(alphabet) * PSEUDOCOUNT)])
    ) for a in alphabet}


def _estimate_f(alphabet: List[str], seqs1: List[AminoacidSequence], seqs2: List[ScoredAminoacidSequence],
                gammas: List[np.ndarray], hmm: PairwiseAlignmentHmm,
                condition: Callable[[Aminoacid, ScoredAminoacid], bool]) -> Dict[str, float]:
    return {a: my_exp(
        log_add_exp([
                        gam[t, s, hmm.state_index(hmm.M)]
                        for seq1, seq2, gam in zip(seqs1, seqs2, gammas)
                        for t in range(1, len(seq1) + 1) for s in range(1, len(seq2) + 1)
                        if seq2.symbols[s - 1].name == a and seq1.symbols[t - 1].name == a and condition(
                seq1.symbols[t - 1],
                seq2.symbols[s - 1])
                    ] + [my_log(PSEUDOCOUNT)]) -
        log_add_exp([
                        gam[t, s, hmm.state_index(hmm.M)]
                        for seq1, seq2, gam in zip(seqs1, seqs2, gammas)
                        for t in range(1, len(seq1) + 1) for s in range(1, len(seq2) + 1)
                        if
                        seq2.symbols[s - 1].name == seq1.symbols[t - 1].name and condition(seq1.symbols[t - 1],
                                                                                           seq2.symbols[s - 1])
                    ] + [my_log(len(alphabet) * PSEUDOCOUNT)])
    ) for a in alphabet}


class BaumWelchParametersEstimator(ParametersCalculator):
    def __init__(self, alignments: List[PairwiseAlignmentOutputWithLogs],
                 data: TwoSequenceListsData,
                 nerpa_cfg_path: str,
                 init_params: PairwiseAlignmentHMMParameters,
                 n_iterations: int,
                 recalculate_transition_probs: bool = True,
                 log_dir: Optional[str] = None,
                 pool_sz: int = 1):
        self._alignments = alignments
        self._data = data
        self._nerpa_cfg = parse_nerpa_config(nerpa_cfg_path)
        self._cur_params = init_params
        self._hmm = PairwiseAlignmentHmm(self._cur_params)
        self._n_iterations = n_iterations
        self._cur_iter = 0
        self._recalculate_transition_probs = recalculate_transition_probs
        self._log_dir = log_dir
        self._omega_a, self._omega_b = get_alphabets_from_data(self._data, self._alignments)
        self._alphabet: List[str] = [a.name for a in self._omega_a] + [b.name for b in self._omega_b]
        self._pool_sz = pool_sz
        print('Will use Baum-Welch parameters estimator')

    def calculate_parameters(self) -> PairwiseAlignmentHMMParameters:
        for _ in range(self._n_iterations):
            print(f'Iteration {self._cur_iter}')
            new_params = self._calculate_parameters_on_iteration()
            d = _calculate_distance(new_params, self._cur_params)
            print(f' distance = {d}')
            self._cur_params = new_params
            self._hmm = PairwiseAlignmentHmm(self._cur_params)
        return self._cur_params

    def _calculate_alpha_beta_gamma_eps(self, alignment: PairwiseAlignmentOutputWithLogs) -> Tuple:
        seq1 = AminoacidSequence(alignment.aligned_sequence1.sequence_id,
                                 [i for i in alignment.aligned_sequence1.symbols if i != Gap()])
        seq2 = ScoredAminoacidSequence(alignment.aligned_sequence2.sequence_id,
                                       [i for i in alignment.aligned_sequence2.symbols if i != Gap()])
        alpha = calculate_log_alpha(seq1, seq2, self._hmm)
        beta = calculate_log_beta(seq1, seq2, self._hmm)
        gamma = _calculate_gamma(seq1, seq2, self._hmm, alpha, beta)
        eps = _calculate_eps(seq1, seq2, self._hmm, alpha, beta)
        return seq1, seq2, alpha, beta, gamma, eps

    def _calculate_g_f(self, case: Tuple[str, bool, str, bool, List, List, List[np.ndarray]]) -> Tuple:
        modification1, methylation1, modification2, methylation2, seqs1, seqs2, gammas = case
        condition = lambda x, y: same_modifications_methylations(x, modification1, methylation1) and \
                                 same_modifications_methylations(y, modification2, methylation2)
        g = _estimate_g(self._alphabet, seqs1, seqs2, gammas, self._hmm, condition)
        f = _estimate_f(self._alphabet, seqs1, seqs2, gammas, self._hmm, condition)
        return g, f

    def _calculate_parameters_on_iteration(self) -> PairwiseAlignmentHMMParameters:

        self._cur_iter += 1

        alignments: List[PairwiseAlignmentOutputWithLogs] = [a for a in self._alignments if
                                                             _check_observation(a, self._omega_a, self._omega_b)]
        with Pool(self._pool_sz) as p:
            prob_data: List[Tuple] = list(
                tqdm(p.imap(self._calculate_alpha_beta_gamma_eps, alignments), total=len(alignments)))

        seqs1, seqs2, alphas, betas, gammas, epss = map(list, zip(*prob_data))
        if self._recalculate_transition_probs:
            mu, tau = np.array(_calculate_pi(gammas, self._hmm)), np.array(
                _calculate_a(seqs1, seqs2, gammas, epss, self._hmm))
        else:
            mu, tau = self._cur_params.mu, self._cur_params.tau

        # 1. Estimate g(a) = P(a | mismatch, mod1, meth1, mod2, meth2)
        # 2. Estimate f(a) = P(a | match,    mod1, meth1, mod2, meth2)
        f, g = {}, {}
        cases = [(modification1, methylation1, modification2, methylation2, seqs1, seqs2, gammas)
                 for modification1 in ['@D', '@L', 'None']
                 for methylation1 in [False, True]
                 for modification2 in ['@D', '@L', 'None']
                 for methylation2 in [False, True]]
        with Pool(self._pool_sz) as p:
            prob_data: List[Tuple] = list(
                tqdm(p.imap(self._calculate_g_f, cases), total=len(cases)))
        for i, (modification1, methylation1, modification2, methylation2, _, _, _) in enumerate(cases):
            g[(modification1, methylation1, modification2, methylation2)] = prob_data[i][0]
            f[(modification1, methylation1, modification2, methylation2)] = prob_data[i][1]

        p = estimate_p_with_modifications(self._omega_a, self._omega_b, self._nerpa_cfg, f, g)
        q_a, q_b = estimate_qa_qb(p)
        p, q_a, q_b = array_to_dict(self._omega_a, self._omega_b, p, q_a, q_b)
        res = PairwiseAlignmentHMMParameters(self._omega_a, self._omega_b, mu=mu, tau=tau, p=p, q_a=q_a, q_b=q_b)
        if self._log_dir is not None:
            res.log_to(os.path.join(self._log_dir, f'baum_welch_parameters_iter{self._cur_iter}.txt'))
        return res


def _calculate_distance(p1: PairwiseAlignmentHMMParameters, p2: PairwiseAlignmentHMMParameters) -> float:
    p = np.array([[p1.p[a][b] - p2.p[a][b] for b in p1.omega_b] for a in p1.omega_a])
    q_a = np.array([p1.q_a[a] - p2.q_a[a] for a in p1.omega_a])
    q_b = np.array([p1.q_b[b] - p2.q_b[b] for b in p1.omega_b])
    return np.linalg.norm(p1.mu - p2.mu) + np.linalg.norm(p1.tau - p2.tau) + np.linalg.norm(p) + np.linalg.norm(
        q_a) + np.linalg.norm(q_b)


def _check_observation(alignment: PairwiseAlignmentOutputWithLogs, omega_a: AminoacidAlphabet,
                       omega_b: ScoredAminoacidAlphabet) -> bool:
    for s in alignment.aligned_sequence1.symbols:
        if s != Gap() and s.name not in [a.name for a in omega_a]:
            return False
    for s in alignment.aligned_sequence2.symbols:
        if s != Gap() and s.name not in [b.name for b in omega_b]:
            return False
    return True
