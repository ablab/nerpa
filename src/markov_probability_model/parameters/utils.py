import numpy as np

from src.markov_probability_model.base.alphabet import Aminoacid, ScoredAminoacid, Gap, AminoacidAlphabet, \
    ScoredAminoacidAlphabet, Symbol
from src.markov_probability_model.data_loader.data_loader import TwoSequenceListsData
from src.markov_probability_model.pairwise_alignment.sequence_aligner import PairwiseAlignmentOutput
from typing import Tuple, List, Dict
from collections import defaultdict


def get_alphabets_from_data(data: TwoSequenceListsData, alignments: List[PairwiseAlignmentOutput]) -> Tuple[
    AminoacidAlphabet, ScoredAminoacidAlphabet]:
    omega_a: List[Aminoacid] = [symb for sequence in data.sequences1 for symb in sequence.symbols if symb != Gap()] + \
                               [s for a in alignments for s in a.aligned_sequence1.symbols if s != Gap()]
    omega_b: List[ScoredAminoacid] = [symb for sequence in data.sequences2
                                      for symb in sequence.symbols if symb != Gap()] + \
                                     [s for a in alignments for s in a.aligned_sequence2.symbols if s != Gap()]
    return AminoacidAlphabet(list(set(omega_a))), ScoredAminoacidAlphabet(list(set(omega_b)))


def get_alphabets_from_alignments(data: List[PairwiseAlignmentOutput]) \
        -> Tuple[AminoacidAlphabet, ScoredAminoacidAlphabet]:
    omega_a: List[Aminoacid] = [symb for alignment in data for symb in alignment.aligned_sequence1.symbols
                                if symb != Gap()]
    omega_b: List[ScoredAminoacid] = [symb for alignment in data
                                      for symb in alignment.aligned_sequence2.symbols if symb != Gap()]
    return AminoacidAlphabet(list(set(omega_a))), ScoredAminoacidAlphabet(list(set(omega_b)))


def parse_nerpa_config(cfg_path: str) -> Dict:
    with open(cfg_path, 'r') as f:
        lines = f.readlines()
    cfg = {'insertion': float(lines[0].rstrip().split(' ')[1]), 'deletion': float(lines[1].rstrip().split(' ')[1]),
           'Scores': [float(l) for l in lines[2].rstrip().split(' ')[1:] if len(l) != 0],
           'ProbGenCorrect': [float(l) for l in lines[3].rstrip().split(' ')[1:] if len(l) != 0],
           'ProbGenIncorrect': [float(l) for l in lines[4].rstrip().split(' ')[1:] if len(l) != 0]}
    return cfg


def _estimate_p_match(a: Aminoacid, b: ScoredAminoacid, nerpa_cfg: Dict, f: Dict,
                      p_score: Dict, mods_prob: Dict) -> float:
    if a.modification is None:
        a_mods_prob = mods_prob.copy()
    else:
        a_mods_prob = defaultdict(int)
        a_mods_prob[a.modification] = 1.0
    if b.modification is None:
        b_mods_prob = mods_prob.copy()
    else:
        b_mods_prob = defaultdict(int)
        b_mods_prob[b.modification] = 1.0
    f_ = sum(f[(mod_a, a.methylation, mod_b, b.methylation)][a.name] * a_mods_prob[mod_a] * b_mods_prob[mod_b]
             for mod_a in ['@L', '@D'] for mod_b in ['@L', '@D'])
    return f_ * get_prob_gen(nerpa_cfg, b.score, p_score, correct=True)


def _estimate_p_miss(a: Aminoacid, b: ScoredAminoacid, nerpa_cfg: Dict, g: Dict,
                     p_score: Dict, mods_prob: Dict) -> float:
    if a.modification is not None:
        g_a = g[(a.modification, a.methylation)][a.name]
    else:
        g_a = sum(g[(mod, a.methylation)][a.name] * mods_prob[mod] for mod in ['@L', '@D'])
    if b.modification is not None:
        g_b = g[(b.modification, b.methylation)][b.name]
    else:
        g_b = sum(g[(mod, b.methylation)][b.name] * mods_prob[mod] for mod in ['@L', '@D'])
    return g_a * g_b * get_prob_gen(nerpa_cfg, b.score, p_score, correct=False)


def estimate_p_with_modifications(omega_a: List[Aminoacid], omega_b: List[ScoredAminoacid],
                                  nerpa_cfg: Dict, f: Dict, g: Dict, p_score: Dict, mods_prob: Dict):
    # Estimate p(a, b).
    # Model:
    #  p(@N-a, @M-a(score)) = f(a, N, M) * prob_gen_correct(score) * P(score)
    #  p(@N-a, @M-b(score)) = g(a, N) * g(b, M) * prob_gen_incorrect(score) * P(score)
    p = np.zeros((len(omega_a), len(omega_b)))
    for i, a in enumerate(omega_a):
        for j, b in enumerate(omega_b):
            if a.name == b.name:
                p[i, j] = _estimate_p_match(a, b, nerpa_cfg, f, p_score, mods_prob)
            else:
                p[i, j] = _estimate_p_miss(a, b, nerpa_cfg, g, p_score, mods_prob)
    p /= p.sum()
    return p


def estimate_qa_qb(p: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    return p.sum(axis=1), p.sum(axis=0)


def array_to_dict(omega_a: AminoacidAlphabet, omega_b: ScoredAminoacidAlphabet, p: np.ndarray, q_a: np.ndarray,
                  q_b: np.ndarray) -> Tuple[Dict, Dict, Dict]:
    p = {a: {b: p[i, j] for j, b in enumerate(omega_b)} for i, a in enumerate(omega_a)}
    q_a = {a: q_a[i] for i, a in enumerate(omega_a)}
    q_b = {b: q_b[i] for i, b in enumerate(omega_b)}
    return p, q_a, q_b


def _truncate_score_to_nerpa_cfg(config: Dict, score: float):
    for i in range(len(config['Scores'])):
        if config['Scores'][i] <= score:
            return config['Scores'][i]


def get_prob_gen(config: Dict, score: float, p_score: Dict, correct: bool):
    score = _truncate_score_to_nerpa_cfg(config, score)
    for i in range(len(config['Scores'])):
        if config['Scores'][i] == score:
            return np.exp(config['ProbGenCorrect' if correct else 'ProbGenIncorrect'][i]) * \
                   p_score[_truncate_score_to_nerpa_cfg(config, score)]


def same_modifications(x: Symbol, mod: str, methylation: bool) -> bool:
    return x.modification == mod and x.methylation == methylation


def get_names_alphabet(omega_a: AminoacidAlphabet, omega_b: ScoredAminoacidAlphabet) -> List[str]:
    return list(set([a.name for a in omega_a] + [b.name for b in omega_b]))


def generate_p_score(config, data: TwoSequenceListsData) -> Dict[int, int]:
    p_score = defaultdict(float)
    div_term = 0
    for seq in data.sequences2:
        for symb in seq.symbols:
            score = int(_truncate_score_to_nerpa_cfg(config, symb.score))
            p_score[score] += 1
            div_term += 1
    return {score: p / div_term for score, p in p_score.items()}


def generate_p_mods(data: TwoSequenceListsData) -> Dict[str, int]:
    p_mods = defaultdict(float)
    div_term = 0
    for seq in data.sequences1:
        for symb in seq.symbols:
            if symb.modification is not None:
                div_term += 1
                p_mods[symb.modification] += 1
    for seq in data.sequences2:
        for symb in seq.symbols:
            if symb.modification is not None:
                div_term += 1
                p_mods[symb.modification] += 1
    return {mod: p / div_term for mod, p in p_mods.items()}
