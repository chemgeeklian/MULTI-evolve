conda activate multievolve
cd /homes/xlian/dgoA/MULTI-evolve/scripts

FASTA_DIR='/homes/xlian/dgoA/NF006600_DgoA_no_KDPG/selected_seq/singleseqfastas/'
OUTPUT_BASE='../data/dgoA/'

for fasta in "$FASTA_DIR"*.fasta; do
    name=$(basename "$fasta" .fasta)
    out_dir="${OUTPUT_BASE}${name}/"
    mkdir -p "$out_dir"

    python rank_positions_by_improvability.py \
        --wt-file "$fasta" \
        --top-n 10 \
        --top-k 5 \
        --excluded-positions 14,37,85,154,126 \
        --out-dir "$out_dir"
done