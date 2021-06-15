#!/usr/bin/env python3
import argparse

import sys
import os
import datetime
import shutil
import logging
import numpy as np
import pandas as pd
from collections import defaultdict

PATH_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(PATH_REPO)

import nerpa
import logger

class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


class ScoringModel(object):
    def __init__(self, path_to_config):
        self.df_aminoacids = None
        self.df_monomers   = None
        self.pred_thresh = None
        self.path_config = path_to_config
        self._init(path_to_config)

    def _init(self, path_to_config):
        self.df_aminoacids = pd.read_csv(os.path.join(path_to_config, 'aminoacids.tsv'), index_col=None, sep='\t')
        self.df_monomers = pd.read_csv(os.path.join(path_to_config, 'monomersLogP.tsv'), index_col=None, sep='\t')
        with open(os.path.join(path_to_config, 'prob_gen.cfg')) as f:
            for line in f:
                if line.startswith('Score'):
                    self.pred_thresh = sorted(map(float, line.strip().split()[1:]), reverse=True)
                    break
            else:
                raise

    def _logp_insertion(self, df):
        mask = df['Monomer_AA'].isna()
        logp = np.log(mask.sum() / df.shape[0])
        return logp

    def _logp_deletion(self, df):
        mask = df['Pred_TopAA'].isna()
        logp = np.log(mask.sum() / df.shape[0])
        return logp

    def _logp_probgen_df(self, df, thresholds):
        raw_scores = df['Pred_TopAAScore'].copy().astype(float)
        bins = np.array(sorted(thresholds, reverse=True))
        ranks = pd.Series(np.digitize(raw_scores, bins=bins, right=False), index=df.index)
        ranks[ranks >= bins.shape[0]] = bins.shape[0] - 1

        mask_match = ~df['Pred_TopAA'].isna() & ~df['Monomer_AA'].isna()
        mask_correct = ~df['Matched_AA'].isna() & mask_match

        counts_correct = ranks.loc[mask_correct].value_counts()
        counts = ranks.loc[mask_match].value_counts()

        logp_correct = np.log(counts_correct) - np.log(counts)
        logp_incorrect = np.log(counts - counts_correct) - np.log(counts)

        df_probgen = pd.concat([pd.Series(dtype=float), logp_correct, logp_incorrect], axis=1)
        df_probgen.columns = ['Scores', 'ProbGenCorrect', 'ProbGenIncorrect']
        df_probgen['Scores'] = bins
        df_probgen.fillna(0, inplace=True)

        return df_probgen.transpose()

    def _compute_logp_from_masks(self, mask_pT, mask_nT):
        """
        input: two binary vectors
            (TODO: arrays with values from {True, False, np.nan})
        output: array with log probs
        """
        counts = np.zeros((2,2))
        counts[0,0] = (~mask_pT & ~mask_nT).sum()
        counts[1,0] = (mask_pT & ~mask_nT).sum()
        counts[0,1] = (~mask_pT & mask_nT).sum()
        counts[1,1] = (mask_pT & mask_nT).sum()

        logp = np.log(1+counts) + np.log(mask_pT.shape[0])
        # PF_NT, PF_NF
        logp[0,:] -= np.log(1 + (~mask_pT).sum())
        # PT_NT, PT_NF
        logp[1,:] -= np.log(1 + mask_pT.sum())
        # PT_NF, PF_NF
        logp[:,0] -= np.log(1 + (~mask_nT).sum())
        # PT_NT, PF_NT
        logp[:,1] -= np.log(1 + mask_nT.sum())
        # PF_NF, PF_NT, PT_NF, PT_NT -> PT_NT, PT_NF, PF_NT, PF_NF
        return logp.flatten()[::-1]

    def _logprobs_mt(self, df):
        mask = ~df['Pred_MT'].isna() & ~df['Monomer_MT'].isna()
        mask_pT = df.loc[mask, 'Pred_MT'].astype(bool)
        mask_nT = df.loc[mask, 'Monomer_MT'].astype(bool)
        #     mask_pF = df['Prediction_Modifications'].str.contains('MT') & ~df['Prediction_Modifications'].isna()
        #     mask_nF = df['Monomer_Modifications'].str.contains('MT') & ~df['Monomer_Modifications'].isna()
        return self._compute_logp_from_masks(mask_pT, mask_nT)

    def _logprobs_d(self, df):
        mask = ~df['Pred_D_conf'].isna() & ~df['Monomer_D_conf'].isna()
        mask_pT = df.loc[mask, 'Pred_D_conf'].astype(bool)
        mask_nT = df.loc[mask, 'Monomer_D_conf'].astype(bool)
        return self._compute_logp_from_masks(mask_pT, mask_nT)

    def _logp_monomers(self, df, monomers):
        counts_dict = df['Monomer_AA'].value_counts().to_dict()
        #     display(counts_dict)
        counts = np.array([counts_dict.get(x, 1) for x in monomers])
        log_freqs = np.log(counts) - np.log(counts.sum())
        return log_freqs

    # def _logp_predictions(self, df, monomers):
    #     counts_dict = df['Matched_Residue'].value_counts().to_dict()
    #     counts = np.array([counts_dict.get(x, 1) for x in monomers])
    #     log_freqs = np.log(counts) - np.log(counts.sum())
    #     return log_freqs

    def generate_modification_scores_str(self, df):
        mask_indel = df['Pred_TopAA'].isna() | df['rBan AA-ID'].isna()
        res = 'Name	predT-nrpT	predT-nrpF	predF-nrpT	predF-nrpF\n'
        res += 'MT ' + ' '.join(f'{x:.3f}' for x in self._logprobs_mt(df)) + '\n'
        res += '@D ' + ' '.join(f'{x:.3f}' for x in self._logprobs_d(df))
        return res

    def generate_prob_gen_str(self, df):
        mask_indel = df['Pred_TopAA'].isna() | df['rBan AA-ID'].isna()
        res  = f'insertion {self._logp_insertion(df)}\n'
        res += f'deletion {self._logp_deletion(df)}\n'
        df_probgen = self._logp_probgen_df(df[~mask_indel], self.pred_thresh)
        # print(df_probgen)
        res += df_probgen.to_string(header=False)
        return res

    # def generate_aa_df(self, df):
    #     df_res = self.df_aminoacids.copy()
    #     monomers = df_res['NameId'].to_list()
    #     logp = self._logp_predictions(df, monomers)
    #     df_res['NRPSPred2_LogP'] = logp
    #     return df_res

    def generate_monomers_df(self, df):
        df_res = self.df_monomers.copy()
        monomers = df_res['NameID'].to_list()
        logp = self._logp_monomers(df, monomers)
        df_res['LogP'] = logp
        return df_res

    def generate_config(self, path_out, df):
        """ generates new config files in the specified directory

        :param path_out: path to new config to be put
        :param df: training dataset
        :return:
        """

        shutil.copytree(self.path_config, path_out)
        with open(os.path.join(path_out, 'modifications.tsv'), 'w') as f:
            f.write(self.generate_modification_scores_str(df))
        with open(os.path.join(path_out, 'prob_gen.cfg'), 'w') as f:
            f.write(self.generate_prob_gen_str(df))
        df_monomers = self.generate_monomers_df(df)
        df_monomers.to_csv(os.path.join(path_out, 'monomersLogP.tsv'), index=False, sep='\t', na_rep='none')
        # df_aminoacids = self.generate_aa_df(df)
        # df_aminoacids.to_csv(os.path.join(path_out, 'aminoacids.tsv'), index=False, sep='\t')


def parse_details_csv(path):
    df = pd.read_csv(path, index_col=False, na_values=['-', 'NA', 'none', ''])
    mask_to_drop = df['PRED_TOP5'].isna() & df['rBan AA-ID'].isna()
    df.drop(df[mask_to_drop].index, inplace=True)

    mask = ~df['rBan AA'].isna()
    df.loc[mask, 'Monomer_AA'] = df.loc[mask, 'rBan AA'].map(lambda x: x.split('+', 1)[0])
    df.loc[mask, 'Monomer_MT'] = df.loc[mask, 'rBan AA'].map(lambda x: '+' in x and 'MT' in x.split('+',1)[1])

    def parse_pred_tops(s):
        chunks = s.split(';')
        def chunk_to_tup(c):
            aa, rest = c.split('(')
            score = float(rest[:-1])
            return aa, score
        return sorted([chunk_to_tup(c) for c in chunks], key=lambda x: x[1], reverse=True)

    def get_max_preds(preds):
        def chunk_to_tup(c):
            aa, rest = c.split('(')
            score = float(rest[:-1])
            return aa, score
        parsed = [chunk_to_tup(c) for c in preds.split(';')]
        max_score = max(x[1] for x in parsed)
        return max_score, [aa for aa, s in parsed if s >= max_score]

    def get_matched_aa(row):
        s, aa = get_max_preds(row['PRED_TOP5'])
        if row['Monomer_AA'] in aa:
            return row['Monomer_AA'], s
        return np.nan, np.nan

    mask = ~df['rBan STRUCT_CONFIGURATION'].isna()
    df.loc[mask, 'Monomer_D_conf'] = df.loc[mask, 'rBan STRUCT_CONFIGURATION'].map(lambda x: x == 'D')

    mask = ~df['PRED_TOP5'].isna()
    df.loc[mask, 'Pred_D_conf'] = df.loc[mask, 'L-/D-'].map(lambda x: x == 'D')
    df.loc[mask, 'Pred_TopAA'] = df.loc[mask, 'PRED_TOP5'].map(lambda s: parse_pred_tops(s)[0][0])
    df.loc[mask, 'Pred_TopAAScore'] = df.loc[mask, 'PRED_TOP5'].map(lambda s: parse_pred_tops(s)[0][1])

    mask = ~df['PRED_TOP5'].isna() & ~df['Monomer_AA'].isna()
    df.loc[mask, 'Matched_AA'] = df.loc[mask].apply(lambda row: get_matched_aa(row)[0], axis=1)
    df.loc[mask, 'Matched_AAScore'] = df.loc[mask].apply(lambda row: get_matched_aa(row)[1], axis=1)

    cols = ['BGC', 'CONTIG', 'ORF', 'A-ID', 'rBan VERTEX', 'rBan AA-ID',
            'Pred_TopAA', 'Pred_TopAAScore', 'Pred_D_conf', 'MT',
            'Matched_AA', 'Matched_AAScore',
            'Monomer_AA', 'Monomer_D_conf', 'Monomer_MT'
            ]

    df_training = df[cols].copy()
    df_training.rename(columns={'MT': 'Pred_MT', 'BGC': 'Id'}, inplace=True)
    return df_training


def parse_args(log):
    parser = argparse.ArgumentParser()
    parser.add_argument("--training-data", dest="train_set",
                        help="csv file with true alignments", type=str, required=True)
    parser.add_argument("--num-folds", dest="n_folds", default=5, type=int, help="number of CV folds", action="store")
    parser.add_argument("--seed", dest="seed", default=0, type=int, help="rng seed", action="store")
    parser.add_argument("--output_dir", "-o", help="output dir", type=str)
    args, unknown = parser.parse_known_args()

    old_sys_argv = sys.argv
    sys.argv = ['nerpa.py'] + unknown + ['--output_dir', args.output_dir]
    nerpa_args = nerpa.parse_args(log)
    sys.argv = old_sys_argv

    if not args.output_dir:  # 'output dir was not specified with -o option'
        args.output_dir = os.path.join(
            PATH_REPO, 'results_cv',
            f"results_{datetime.datetime.now().strftime('%Y_%m_%d_%H_%M_%S')}")

    return args, nerpa_args


def k_fold_split(df, k=5, seed=0):
    rng = np.random.default_rng(seed)
    accessions = df['Id'].unique().copy()
    rng.shuffle(accessions)
    size = accessions.shape[0] // k
    folds = []
    for i in range(k - 1):
        folds.append(accessions[i * size: (i + 1) * size])
    folds.append(accessions[(k - 1) * size:])
    masks = [df['Id'].isin(fold) for fold in folds]
    return folds, masks


def run_on_fold(args, log, path_to_fold):
    out_path = os.path.join(path_to_fold, 'run/')
    log.info('\n============================ [fold Nerpa run] ============================')
    log.info(f'Results: {out_path}')

    args_dict = {**vars(args)}
    input_keys = ['smiles', 'smiles_tsv', 'rban_output', 'antismash_out', 'antismash']
    for key in input_keys:
        args_dict[key] = None
    args_dict['output_dir'] = out_path
    # args_dict['structures'] = os.path.join(path_to_fold, 'structures.info')
    # # args_dict['predictions'] = [x.strip().split()[0] for x in open(os.path.join(path_to_fold, 'predictions.info'))]
    # args_dict['predictions'] = [os.path.join(path_to_fold, 'predictions.info')]
    args_dict['configs_dir'] = os.path.join(path_to_fold, 'configs')
    log.info('Running arguments:')
    for k,v in args_dict.items():
        log.info(f'\t{k} : {v}')

    # log_nerpa = logger.NerpaLogger()
    cwd = os.getcwd()
    logger_nerpa = logger.NerpaLogger()
    nerpa.run(AttrDict(args_dict), logger_nerpa)
    logger_nerpa._logger.handlers.clear()
    os.chdir(cwd)
    log.info('\nDone.')


def main():
    log = logging.getLogger("nerpa_cv")
    log.setLevel(logging.DEBUG)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.DEBUG)
    log.addHandler(console_handler)

    args, nerpa_args = parse_args(log)
    os.makedirs(args.output_dir)
    file_handler = logging.FileHandler(os.path.join(args.output_dir, 'nerpa_cv.log'))
    file_handler.setLevel(logging.DEBUG)
    log.addHandler(file_handler)

    df = parse_details_csv(args.train_set)
    df.to_csv(os.path.join(args.output_dir, 'df_training.tsv'), index=False, sep='\t')
    conf_gen = ScoringModel(os.path.join(PATH_REPO, 'configs'))

    all_out_path = os.path.join(args.output_dir, 'all_data')
    conf_gen.generate_config(os.path.join(all_out_path, 'configs'), df)
    run_on_fold(nerpa_args, log, all_out_path)

    log.info(f'Splitting into {args.n_folds} folds. Using seed={args.seed}.')
    folds, masks = k_fold_split(df, k=args.n_folds, seed=args.seed)
    for i, fold in enumerate(folds):
        log.info(f'Fold #{i:02d} :')
        log.info(f'\tsize: {len(fold)}')
        log.info(f'\tcontent: {fold}')

    for i, (fold, mask) in enumerate(zip(folds, masks)):
        fold_out_path = os.path.join(args.output_dir, f'fold-{i:02d}')
        os.makedirs(fold_out_path, exist_ok=True)
        with open(os.path.join(fold_out_path, 'holdout.ids'), 'w') as f:
            f.write('\n'.join(fold))
        df[~mask].to_csv(os.path.join(fold_out_path, f'matches_training_{i:02d}.tsv'), index=False, sep='\t', na_rep='NaN')
        conf_gen.generate_config(os.path.join(fold_out_path, 'configs'), df[~mask])
        run_on_fold(nerpa_args, log, fold_out_path)


    results = []
    for i, (fold, mask) in enumerate(zip(folds, masks)):
        fold_out_path = os.path.join(args.output_dir, f'fold-{i:02d}')
        df_fold = pd.read_csv(os.path.join(fold_out_path, 'run', 'report.csv'))
        mask = df_fold['prediction id'].map(lambda x: x.rsplit('/', 1)[1].split('_',1)[0] in fold)
        df_fold.loc[mask, 'Fold_id'] = i
        results.append(df_fold.loc[mask])

    df_results = pd.concat(results, axis=0)
    df_results.to_csv(os.path.join(args.output_dir, 'report_cv.csv'), index=False)


if __name__ == '__main__':
    main()