#!/usr/bin/env python3
"""
Score amino acid positions by the average ESM log-likelihood ratio across
all possible substitutions and all ESM models, then output the top N positions.

For each position the score is:
    position_score = mean(average_model_logratio over all 19 substitutions)

This answers: "at which positions does the model think mutations are most
likely to be beneficial, on average across all possible substitutions?"

Usage (run from scripts/):
cd scripts
python zeroshot_top_positions_esm.py \
    --wt-file ../data/dgoA/Q6BF16_V85A.fasta \
    --top-n 10 \
    --excluded-positions 14,37,85,154,126
"""

import argparse
import os

import pandas as pd
from Bio import SeqIO

from multievolve import zero_shot_esm_dms


def parse_args():
    parser = argparse.ArgumentParser(
        description='Rank positions by average ESM mutation score'
    )
    parser.add_argument(
        '--wt-file',
        required=True,
        help='Path to the wildtype FASTA file'
    )
    parser.add_argument(
        '--top-n',
        type=int,
        default=10,
        help='Number of top positions to output (default: 10)'
    )
    parser.add_argument(
        '--excluded-positions',
        required=False,
        help='Comma-separated positions to exclude (1-indexed)'
    )
    parser.add_argument(
        '--output',
        required=False,
        help='Output CSV path (default: zeroshot_top_positions_esm.csv next to wt-file)'
    )
    return parser.parse_args()


def main():
    args = parse_args()

    excluded = (
        [int(p) for p in args.excluded_positions.split(',')]
        if args.excluded_positions else []
    )

    wt_seq = str(SeqIO.read(args.wt_file, 'fasta').seq)
    print(f'Sequence length: {len(wt_seq)} residues')
    print('Running ESM zero-shot DMS...')

    dms_df = zero_shot_esm_dms(wt_seq)

    # Parse position from mutation string e.g. "A109F" -> 109
    dms_df['pos'] = dms_df['mutations'].apply(lambda x: int(x[1:-1]))
    dms_df['wt_aa'] = dms_df['mutations'].apply(lambda x: x[0])

    # Drop excluded positions
    if excluded:
        dms_df = dms_df[~dms_df['pos'].isin(excluded)]

    # Per-model logratio columns returned by zero_shot_esm_dms
    model_logratio_cols = [c for c in dms_df.columns if c.startswith('model_') and c.endswith('_logratio')]
    n_models = len(model_logratio_cols)

    # Aggregate per position: mean across all substitutions for each metric
    pos_df = (
        dms_df.groupby(['pos', 'wt_aa'])
        .agg(
            position_score=('average_model_logratio', 'mean'),
            max_logratio=('average_model_logratio', 'max'),
            n_positive=('average_model_logratio', lambda x: (x > 0).sum()),
            **{f'model_{i+1}_mean': (f'model_{i+1}_logratio', 'mean') for i in range(n_models)}
        )
        .reset_index()
        .sort_values('position_score', ascending=False)
    )

    top = pos_df.head(args.top_n).reset_index(drop=True)

    print(f'\nTop {args.top_n} positions by average ESM mutation score:')
    print(top[['pos', 'wt_aa', 'position_score', 'max_logratio', 'n_positive']].to_string(index=False))

    out_path = args.output or os.path.join(
        os.path.dirname(args.wt_file), 'zeroshot_top_positions_esm.csv'
    )
    top.to_csv(out_path, index=False)
    print(f'\nWrote {len(top)} positions -> {out_path}')


if __name__ == '__main__':
    main()
