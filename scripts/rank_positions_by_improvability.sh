conda activate multievolve
cd /homes/xlian/dgoA/MULTI-evolve/scripts

python rank_positions_by_improvability.py \
    --wt-file ../data/dgoA/Q6BF16_V85A.fasta \
    --top-n 10 \
    --top-k 5 \
    --excluded-positions 14,37,85,154,126