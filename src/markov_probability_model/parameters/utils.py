import numpy as np

from src.markov_probability_model.base.alphabet import Aminoacid, ScoredAminoacid, Gap, AminoacidAlphabet, \
    ScoredAminoacidAlphabet, Symbol
from src.markov_probability_model.data_loader.data_loader import TwoSequenceListsData
from src.markov_probability_model.pairwise_alignment.sequence_aligner import PairwiseAlignmentOutput
from typing import Tuple, List, Dict, Optional


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


def estimate_p(omega_a: List[Aminoacid], omega_b: List[ScoredAminoacid],
               nerpa_cfg: Dict, f: Dict[str, float], g: Dict[str, float]):
    # Estimate p(a, b).
    # Model:
    #  p(a, a(score)) = f(a) * prob_gen_correct(score)
    #  p(a, b(score)) = g(a) * g(b) * prob_gen_incorrect(score)
    p = np.zeros((len(omega_a), len(omega_b)))
    for i, a in enumerate(omega_a):
        for j, b in enumerate(omega_b):
            if a.name == b.name:
                p[i, j] = f[a.name] * get_prob_gen(nerpa_cfg, b.score, correct=True)
            else:
                p[i, j] = g[a.name] * g[b.name] * get_prob_gen(nerpa_cfg, b.score, correct=False)
    p /= p.sum()
    return p


def estimate_p_with_modifications(omega_a: List[Aminoacid], omega_b: List[ScoredAminoacid],
                                  nerpa_cfg: Dict, f: Dict, g: Dict):
    # Estimate p(a, b).
    # Model:
    #  p(a, a(score)) = f(a) * prob_gen_correct(score)
    #  p(a, b(score)) = g(a) * g(b) * prob_gen_incorrect(score)
    p = np.zeros((len(omega_a), len(omega_b)))
    for i, a in enumerate(omega_a):
        for j, b in enumerate(omega_b):
            idx = (a.modification if a.modification is not None else 'None', a.methylation,
                   b.modification if b.modification is not None else 'None', b.methylation)
            if a.name == b.name:
                p[i, j] = f[idx][a.name] * get_prob_gen(nerpa_cfg, b.score, correct=True)
            else:
                p[i, j] = g[idx][a.name] * g[idx][b.name] * get_prob_gen(nerpa_cfg, b.score, correct=False)
                p[i, j] = g[idx][a.name] * g[idx][b.name] * get_prob_gen(nerpa_cfg, b.score, correct=False)
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


def get_prob_gen(config, score, correct: bool):
    for i in range(len(config['Scores'])):
        if config['Scores'][i] <= score:
            return np.exp(config['ProbGenCorrect' if correct else 'ProbGenIncorrect'][i])


def same_modifications_methylations(x: Symbol, mod: Optional[str], methylation: bool) -> bool:
    return (x.modification == mod or mod == 'None') and x.methylation == methylation
