#!/usr/bin/env python3
"""
Rank amino acid positions by ESM zero-shot scores and select top k mutants
per position. Runs ESM DMS once and writes two CSVs:

  1. position_improvability.csv      — all positions, ranked by
                                       (n_above_wt desc, mean_logratio desc)
  2. top_positions_by_improvability.csv — top k mutants for each of the top N
                                          positions, with per-model means

Position columns:
  pos, wt_aa, n_above_wt, mean_logratio, max_logratio, model_1_mean, ...

Mutant columns:
  pos, wt_aa, mutations, rank_within_pos, average_model_logratio,
  total_model_pass, n_above_wt, mean_logratio

Usage (run from scripts/):
python rank_positions_by_improvability.py \
    --wt-file ../data/dgoA/Q6BF16_V85A.fasta \
    --top-n 10 \
    --top-k 5 \
    --excluded-positions 14,37,85,154,126
"""

import argparse
import os
import pandas as pd
from Bio import SeqIO

from multievolve import zero_shot_esm_dms


def parse_args():
    parser = argparse.ArgumentParser(
        description='Rank positions by ESM improvability and select top k mutants'
    )
    parser.add_argument('--wt-file', required=True,
                        help='Path to wildtype FASTA file')
    parser.add_argument('--top-n', type=int, default=10,
                        help='Number of top positions to output (default: 10)')
    parser.add_argument('--top-k', type=int, default=5,
                        help='Number of top mutants per position (default: 5)')
    parser.add_argument('--excluded-positions', required=False,
                        help='Comma-separated positions to exclude (1-indexed)')
    parser.add_argument('--out-dir', required=False,
                        help='Output directory (default: same dir as --wt-file)')
    return parser.parse_args()


def main():
    args = parse_args()

    excluded = (
        [int(p) for p in args.excluded_positions.split(',')]
        if args.excluded_positions else []
    )
    out_dir = args.out_dir or os.path.dirname(os.path.abspath(args.wt_file))
    os.makedirs(out_dir, exist_ok=True)

    wt_seq = str(SeqIO.read(args.wt_file, 'fasta').seq)
    print(f'Sequence length: {len(wt_seq)} residues')
    print(f'Excluded positions: {excluded}')
    print('\nRunning ESM zero-shot DMS...')
    df = zero_shot_esm_dms(wt_seq)

    df['pos'] = df['mutations'].apply(lambda x: int(x[1:-1]))
    df['wt_aa'] = df['mutations'].apply(lambda x: x[0])
    df = df[~df['pos'].isin(excluded)]

    model_logratio_cols = [c for c in df.columns if c.startswith('model_') and c.endswith('_logratio')]
    n_models = len(model_logratio_cols)

    # --- CSV 1: all positions ranked ---
    pos_stats = (
        df.groupby(['pos', 'wt_aa'], sort=False)
        .agg(
            n_above_wt=('average_model_logratio', lambda s: (s > 0).sum()),
            mean_logratio=('average_model_logratio', 'mean'),
            max_logratio=('average_model_logratio', 'max'),
            **{f'model_{i+1}_mean': (f'model_{i+1}_logratio', 'mean') for i in range(n_models)}
        )
        .reset_index()
        .sort_values(['n_above_wt', 'mean_logratio'], ascending=False)
    )

    pos_out = os.path.join(out_dir, 'position_improvability.csv')
    pos_stats.to_csv(pos_out, index=False)
    print(f'\nAll positions ranked -> {pos_out}')

    top_positions = pos_stats.head(args.top_n)
    print(f'\nTop {args.top_n} positions:')
    print(f'{"pos":>5}  {"wt":>3}  {"n_above_wt":>10}  {"mean_logratio":>13}  {"max_logratio":>12}')
    print('-' * 52)
    for _, row in top_positions.iterrows():
        print(f'{int(row["pos"]):>5}  {row["wt_aa"]:>3}  {int(row["n_above_wt"]):>10}  '
              f'{row["mean_logratio"]:>13.4f}  {row["max_logratio"]:>12.4f}')

    # --- CSV 2: top k mutants per top N positions ---
    top_muts = df[df['pos'].isin(top_positions['pos'])].copy()
    top_muts['rank_within_pos'] = (
        top_muts.groupby('pos')['average_model_logratio']
        .rank(ascending=False, method='first')
        .astype(int)
    )
    top_muts = top_muts[top_muts['rank_within_pos'] <= args.top_k]
    top_muts = top_muts.merge(
        top_positions[['pos', 'n_above_wt', 'mean_logratio']], on='pos', how='left'
    )
    out_cols = ['pos', 'wt_aa', 'mutations', 'rank_within_pos',
                'average_model_logratio', 'total_model_pass', 'n_above_wt', 'mean_logratio']
    top_muts = (
        top_muts[out_cols]
        .sort_values(['mean_logratio', 'rank_within_pos'], ascending=[False, True])
    )

    muts_out = os.path.join(out_dir, 'top_positions_by_improvability.csv')
    top_muts.to_csv(muts_out, index=False)
    print(f'\nTop {args.top_k} mutants per position -> {muts_out}')


if __name__ == '__main__':
    main()
