conda activate multievolve

plm_zeroshot_ensemble.py \
--wt-file /homes/xlian/dgoA/MULTI-evolve/data/dgoA/Q6BF16_V85A.fasta \
--pdb-files /nfs/ml_lab/projects/ml_lab/xlian/dgoa_sim/boltz_results_Q6BF16_V85A/predictions/Q6BF16_V85A/Q6BF16_V85A_model_0.pdb \
--chain-id A \
--variants 24 \
--normalizing-method aa_substitution_type \
--excluded-positions 14,37,85,154,126