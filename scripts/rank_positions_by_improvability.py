#!/usr/bin/env python3
"""
Rank amino acid positions by ESM zero-shot scores and select top k mutants
per position. Runs ESM DMS once and writes two CSVs:

  1. position_improvability.csv      — all positions, ranked by
                                       (n_above_wt desc, mean_logratio desc)
  2. top_positions_by_improvability.csv — top k mutants for each of the top N
                                          positions, with per-model means

If --pdb-files and --chain-id are provided, ESM-IF (structure-based) logratios
are averaged across PDB files and combined with ESM sequence scores into a single
combined_logratio, which is used for all ranking and aggregation.

Position columns:
  pos, wt_aa, n_above_wt, mean_logratio, max_logratio, model_1_mean, ...
  (+ esm_if_mean, combined_mean if ESM-IF is used)

Mutant columns:
  pos, wt_aa, mutations, rank_within_pos, combined_logratio,
  average_model_logratio, esm_if_logratio (if used),
  n_above_wt, mean_logratio

Usage (run from scripts/):
python rank_positions_by_improvability.py \
    --wt-file ../data/dgoA/Q6BF16_V85A.fasta \
    --pdb-files ../data/dgoA/Q6BF16_V85C_model_0.cif \
    --chain-id A \
    --top-n 10 \
    --top-k 5 \
    --excluded-positions 14,37,85,154,126
    --out-dir ../data/
"""

import argparse
import os
import pandas as pd
from Bio import SeqIO

from multievolve import zero_shot_esm_dms, zero_shot_esm_if_dms


def parse_args():
    parser = argparse.ArgumentParser(
        description='Rank positions by ESM improvability and select top k mutants'
    )
    parser.add_argument('--wt-file', required=True,
                        help='Path to wildtype FASTA file')
    parser.add_argument('--pdb-files', required=False,
                        help='Comma-separated list of PDB/CIF structure files (enables ESM-IF scoring)')
    parser.add_argument('--chain-id', required=False, default='A',
                        help='Chain ID to use for ESM-IF (default: A)')
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

    pdb_files = (
        [f.strip() for f in args.pdb_files.split(',')]
        if args.pdb_files else []
    )

    wt_seq = str(SeqIO.read(args.wt_file, 'fasta').seq)
    print(f'Sequence length: {len(wt_seq)} residues')
    print(f'Excluded positions: {excluded}')
    print('\nRunning ESM zero-shot DMS...')

    df = zero_shot_esm_dms(wt_seq)

    # Extract position and wildtype amino acid from 'mutations' column
    df['pos'] = df['mutations'].apply(lambda x: int(x[1:-1]))
    df['wt_aa'] = df['mutations'].apply(lambda x: x[0])
    df = df[~df['pos'].isin(excluded)]

    # Extract model logratio columns 
    model_logratio_cols = [c for c in df.columns if c.startswith('model_') and c.endswith('_logratio')]
    n_models = len(model_logratio_cols)

    # --- Optional ESM-IF scoring ---
    '''
    I haven't tested use_esm_if.
    The score of ESM-IF is not scaled the same as ESM sequence-bsed score in the original code. So I normalized both scores to z-scores before averaging them into a combined score. 
    '''
    use_esm_if = bool(pdb_files)
    if use_esm_if:
        esm_if_dfs = []
        for pdb_file in pdb_files:
            print(f'Running ESM-IF zero-shot DMS on {pdb_file}...')
            esm_if_df = zero_shot_esm_if_dms(wt_seq, pdb_file, chain_id=args.chain_id)
            esm_if_dfs.append(esm_if_df[['mutations', 'logratio']])

        # Average ESM-IF logratio across all PDB files
        esm_if_merged = esm_if_dfs[0].rename(columns={'logratio': 'logratio_pdb0'})
        for j, esm_if_df in enumerate(esm_if_dfs[1:], start=1):
            esm_if_merged = esm_if_merged.merge(
                esm_if_df.rename(columns={'logratio': f'logratio_pdb{j}'}),
                on='mutations', how='outer'
            )
        pdb_cols = [c for c in esm_if_merged.columns if c.startswith('logratio_pdb')]
        esm_if_merged['esm_if_logratio'] = esm_if_merged[pdb_cols].mean(axis=1)
        esm_if_merged = esm_if_merged[['mutations', 'esm_if_logratio']]

        # Merge into main df and compute combined score
        df = df.merge(esm_if_merged, on='mutations', how='left')
        # Z-score normalize each signal globally before combining so that
        # neither score dominates due to scale differences
        df['esm_z'] = (
            (df['average_model_logratio'] - df['average_model_logratio'].mean())
            / df['average_model_logratio'].std()
        )
        df['esm_if_z'] = (
            (df['esm_if_logratio'] - df['esm_if_logratio'].mean())
            / df['esm_if_logratio'].std()
        )
        df['combined_logratio'] = (df['esm_z'] + df['esm_if_z']) / 2
    else:
        df['combined_logratio'] = df['average_model_logratio']

    # --- CSV 1: all positions ranked ---
    extra_aggs = {f'model_{i+1}_mean': (f'model_{i+1}_logratio', 'mean') for i in range(n_models)}

    if use_esm_if:
        extra_aggs['esm_if_mean'] = ('esm_if_logratio', 'mean')

    pos_stats = (
        df.groupby(['pos', 'wt_aa'], sort=False)
        .agg(
            n_above_wt=('combined_logratio', lambda s: (s > 0).sum()),
            mean_logratio=('combined_logratio', 'mean'),
            max_logratio=('combined_logratio', 'max'),
            **extra_aggs
        )
        .reset_index()
        .sort_values(['n_above_wt', 'mean_logratio'], ascending=False)
    )

    # output the csv with all positions ranked by improvability
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
        top_muts.groupby('pos')['combined_logratio']
        .rank(ascending=False, method='first')
        .astype(int)
    )
    top_muts = top_muts[top_muts['rank_within_pos'] <= args.top_k]
    top_muts = top_muts.merge(
        top_positions[['pos', 'n_above_wt', 'mean_logratio']], on='pos', how='left'
    )

    out_cols = ['pos', 'wt_aa', 'mutations', 'rank_within_pos',
                'combined_logratio', 'average_model_logratio', 
                'n_above_wt', 'mean_logratio']
    
    if use_esm_if:
        out_cols.insert(out_cols.index('average_model_logratio') + 1, 'esm_if_logratio')
        
    top_muts = (
        top_muts[out_cols]
        .sort_values(['mean_logratio', 'rank_within_pos'], ascending=[False, True])
    )

    muts_out = os.path.join(out_dir, 'top_positions_by_improvability.csv')
    top_muts.to_csv(muts_out, index=False)
    print(f'\nTop {args.top_k} mutants per position -> {muts_out}')


if __name__ == '__main__':
    main()
