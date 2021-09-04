import os
import argparse

from src.markov_probability_model.data_loader.mibig_alignments_loader import MibigAlignmentsLoader, \
    PairwiseAlignmentOutputWithLogs
from src.markov_probability_model.data_loader.data_loader import TwoSequenceListsData
from src.markov_probability_model.data_loader.raw_data_parser import RawDataParser
from src.markov_probability_model.data_loader.fdr_loader import FdrGeneratorFromReport
from src.markov_probability_model.parameters.ml_parameters_estimator_with_modifications import \
    MaxLikelihoodParametersEstimatorWithModifications
from src.markov_probability_model.parameters.baum_welch_parameters_estimator import BaumWelchParametersEstimator
from src.markov_probability_model.parameters.utils import get_alphabets_from_alignments
from src.markov_probability_model.hmm.pairwise_alignment_hmm import PairwiseAlignmentHmm, PairwiseAlignmentHMMParameters
from src.markov_probability_model.pairwise_alignment.algo.viterbi import Viterbi, ViterbiOutput
from src.markov_probability_model.pairwise_alignment.algo.maximum_posterior_decoding import \
    MaximumPosteriorDecodingOutput, MaximumPosteriorDecoding
from src.markov_probability_model.pairwise_alignment.algo.maximum_accuracy import MaximumAccuracyOutput, MaximumAccuracy
from src.markov_probability_model.pairwise_alignment.algo.global_viterbi import GlobalViterbi, GlobalViterbiOutput
from src.markov_probability_model.pairwise_alignment.score_augmentations import NullHypothesisScoreAugmentator, \
    IdentityScoreAugmentator
from src.markov_probability_model.pairwise_alignment.alignment_generator import AllPairsAlignmentGenerator
from src.markov_probability_model.pairwise_alignment.fdr import FdrGenerator, FdrParameters, plot_fdrs, FdrData
from typing import List, Dict


def generate_alignments_for_parameters(algo: List[str], pool_sz: int, topk: List[int],
                                       data: TwoSequenceListsData, parameters: PairwiseAlignmentHMMParameters,
                                       all_fdr: List[Dict[str, FdrData]], res_fdr_folder: str):
    hmm = PairwiseAlignmentHmm(parameters)
    assert len(algo) > 0
    for algo in algo:
        print(f' Generating {algo} alignments...')
        score_augmentator = {'viterbi': NullHypothesisScoreAugmentator,
                             'global_viterbi': NullHypothesisScoreAugmentator,
                             'maximum_posterior_decoding': IdentityScoreAugmentator,
                             'maximum_accuracy': IdentityScoreAugmentator}[algo]()
        sequence_aligner = {'viterbi': Viterbi, 'global_viterbi': GlobalViterbi,
                            'maximum_posterior_decoding': MaximumPosteriorDecoding,
                            'maximum_accuracy': MaximumAccuracy}[algo](hmm, score_augmentator)
        alignment_type = {'viterbi': ViterbiOutput, 'global_viterbi': GlobalViterbiOutput,
                          'maximum_posterior_decoding': MaximumPosteriorDecodingOutput,
                          'maximum_accuracy': MaximumAccuracyOutput}[algo]
        alignments = AllPairsAlignmentGenerator[alignment_type](
            data, sequence_aligner, pool_sz=pool_sz).generate_alignments()

        print(f' Generating {algo} FDR...')
        fdr_parameters: List[FdrParameters] = [
            FdrParameters(k, rt,
                          os.path.join(res_fdr_folder, f'pairs_{algo}_{rt}.csv'),
                          os.path.join(res_fdr_folder, f'best_pairs_{algo}_{rt}_top{k}.csv'),
                          ) for k in topk for rt in ['mol', 'genome']
        ]
        algo_fdrs: List[FdrData] = FdrGenerator(alignments, fdr_parameters).generate_fdr()
        for i, fdr in enumerate(algo_fdrs):
            all_fdr[i][algo] = fdr

    print(' Saving FDRs...')
    plot_fdrs(fdr_parameters, all_fdr, save_dir=res_fdr_folder)


def run(data_dir: str, prob_gen_filepath: str,
        results_dir: str, mibig_path: str, pool_sz: int, algo: List[str],
        use_bw: bool, bw_iters: int, topk: List[int]):
    print('Starting alignments generation using Hidden Markov Model...')
    print(f'Generating results for {data_dir}')

    res_parameters_folder = os.path.join(results_dir, 'parameters')
    res_fdr_folder = os.path.join(results_dir, 'fdr')
    for folder in [res_parameters_folder, res_fdr_folder]:
        if not os.path.exists(folder):
            os.makedirs(folder)

    print(' Loading Mibig alignments...')
    ground_truth_alignments: List[PairwiseAlignmentOutputWithLogs] = \
        MibigAlignmentsLoader(mibig_path).load_alignments()

    print(' Loading data...')
    data: TwoSequenceListsData = RawDataParser(
        data_dir, get_alphabets_from_alignments(ground_truth_alignments)[0]).load_data()

    print(' Estimating parameters...')
    parameters: PairwiseAlignmentHMMParameters = \
        MaxLikelihoodParametersEstimatorWithModifications(ground_truth_alignments, data, prob_gen_filepath,
                                                          log_dir=res_parameters_folder).calculate_parameters()

    fdr_parameters: List[FdrParameters] = [
        FdrParameters(
            topk,
            rt,
            os.path.join(res_fdr_folder, f'pairs_v_{rt}.csv'),
            os.path.join(res_fdr_folder, f'best_pairs_v_{rt}_top{topk}.csv'),
        ) for topk in topk for rt in ['mol', 'genome']
    ]
    all_fdr: List[Dict[str, FdrData]] = FdrGeneratorFromReport(data_dir, fdr_parameters).load_fdr()

    if use_bw:
        baum_welch = BaumWelchParametersEstimator(ground_truth_alignments, data, prob_gen_filepath,
                                                  init_params=parameters, n_iterations=1,
                                                  recalculate_transition_probs=True,
                                                  log_dir=res_parameters_folder,
                                                  pool_sz=pool_sz)
        for it in range(bw_iters):
            parameters: PairwiseAlignmentHMMParameters = baum_welch.calculate_parameters()
            cur_res_fdr_folder = os.path.join(res_fdr_folder, f'iter{it}')
            if not os.path.exists(cur_res_fdr_folder):
                os.makedirs(cur_res_fdr_folder)
            generate_alignments_for_parameters(algo, pool_sz, topk, data, parameters, all_fdr,
                                               cur_res_fdr_folder)
    else:
        generate_alignments_for_parameters(algo, pool_sz, topk, data, parameters, all_fdr, res_fdr_folder)

    print(' Done!')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_dir', type=str, default='data/sequences/res_small_with_modifications')
    parser.add_argument('--results_dir', type=str, default='results/res_small')
    parser.add_argument('--mibig_path', type=str, default='data/sequences/mibig.csv')
    parser.add_argument('--pool_sz', type=int, default=8)

    parser.add_argument('--algo', nargs='*', action='append',
                        default=['viterbi', 'global_viterbi', 'maximum_accuracy', 'maximum_posterior_decoding'],
                        choices=['viterbi', 'global_viterbi', 'maximum_accuracy', 'maximum_posterior_decoding'])
    parser.add_argument('--use_bw', type=bool, default=True)
    parser.add_argument('--bw_iters', type=int, default=1000)

    parser.add_argument('--topk', type=list, default=[1, 3, 5, 10])

    args = parser.parse_args()

    run(args.data_dir, 'data/parameters/prob_gen.cfg', args.results_dir, args.mibig_path, args.pool_sz,
        args.algo, args.use_bw, args.bw_iters, args.topk)
