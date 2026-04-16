Here's exactly what the script does, start to finish. It doesn't use anything in `multievolve/proposers/` — zero-shot has its own simple proposal logic living in these two files:

## Input → Output

**Input:** a wildtype FASTA + one or more PDB/CIF structures + `--variants N`
**Output:** `plm_zeroshot_ensemble_nominated_mutations.csv` listing up to `4 × N` nominated single-point mutations.

## Step 1 — Enumerate every possible single mutation

In `zero_shot_esm_dms()` (and `zero_shot_esm_if_dms()`), the "proposer" is just a double for-loop: for every position in the WT sequence, for every amino acid except the WT residue, construct the mutation string `"{wt}{pos}{mt}"`. That's it — no heuristics, just exhaustive DMS over the sequence.

````python path=multievolve/utils/zeroshot_utils.py mode=EXCERPT
for i, residue in enumerate(wt_seq):
    for aa in amino_acids:
        if wt_seq[i] == aa:
            continue
        mutations.append(wt_seq[i] + str(i + 1) + aa)
````

## Step 2 — Score every mutation with two complementary PLM families

### (a) ESM sequence ensemble — `zero_shot_esm_dms`

Loads **6 models**: ESM-1v seeds 1-5 + ESM-2 3B. For each model:

1. Feed WT sequence through the transformer → get logit tensor over (position, vocab).
2. For each mutation `W{pos}M`, compute the **log-likelihood ratio**
   `score = log P(M | context) − log P(W | context)` at that position.
3. This is the `wt-marginals` strategy (the default). A `masked-marginals` alternative mask-fills each position instead.

Each mutation ends up with 6 scores. The script stores them as `model_1_logratio … model_6_logratio`, plus:
- `average_model_logratio` — mean across the 6 models
- `total_model_pass` — how many of the 6 voted `logratio > 0` (i.e. predicted beneficial)

### (b) ESM-IF structure model — `zero_shot_esm_if_dms`

Uses `esm_if1_gvp4_t16_142M_UR50`, which conditions on 3D coordinates (this is why you need the PDB/CIF). Same idea:
- Pass coords + WT sequence through the inverse-folding model.
- Score each mutation as `logit[M] − logit[W]` at that position.

If you pass multiple PDB files (e.g., different conformations), the script runs this per-structure and averages them into `average_model_logratio`.

## Step 3 — Rank, then sample ≤1 mutation per position

Back in `plm_zeroshot_ensemble.py`, there are **four ranking variants** produced, each a different way of ordering the same DMS table:

| Nomination set | How it's ranked |
|---|---|
| `muts_esm` | ESM table sorted by `(total_model_pass desc, average_model_logratio desc)` — favors mutations where all 6 ESM models agree |
| `muts_esm_if` | ESM-IF table sorted by `average_model_logratio desc` |
| `muts_esm_z` | ESM table grouped by `--normalizing-method` (e.g. `aa_substitution_type` = `"A-V"`, `"A-L"`, …), converted to within-group **z-scores**, then sorted |
| `muts_esm_if_z` | Same z-score normalization, but applied to ESM-IF scores |

The z-scored variants exist to correct for the fact that raw log-ratios are systematically biased per substitution class (some "A→X" changes just have higher baseline probability than others). Within each group, a z-score measures *"how unusually good is this mutation compared to other mutations of the same type?"*.

Then `sample_mutations(df, variants, excluded_positions)` walks the sorted list top-down and **takes at most one mutation per residue position**, skipping any `--excluded-positions`, until it has `variants` rows:

````python path=scripts/plm_zeroshot_ensemble.py mode=EXCERPT
for index, row in df.iterrows():
    if row['pos'] not in pos:
        muts.append(row.to_frame().T)
        pos.append(row['pos'])
    if len(muts) == total_muts:
        break
````

This one-per-position rule is the key design choice: it prevents the output from being dominated by several substitutions at the same "hot" position, spreading nominations across the protein.

## Step 4 — Merge the four sets into one CSV

`merge_mutation_dfs` does an outer-join on `mutations` and writes a CSV with one row per nominated mutation and four 0/1 columns (`esm_sampled`, `esm_if_sampled`, `esm_z_sampled`, `esm_if_z_sampled`) indicating which of the four methods nominated it.

Mutations picked by multiple methods → higher confidence. This is the "ensemble" in the name — not an ensemble of *scores*, but an ensemble of *ranking criteria*, each casting a vote.

## Diagram

```
WT FASTA ─┬──► zero_shot_esm_dms ──► DMS table with 6 logratios
          │                          (avg logratio, total_model_pass)
          │
PDB/CIF ──┴──► zero_shot_esm_if_dms ─► DMS table with IF logratio
                                       (averaged across structures)

                    │
                    ▼
     Rank 4 different ways:                     Take top N
     ┌────────────────────┐                     with ≤1 per
     │ 1. ESM raw         │ ──► sample_mutations ──► position,
     │ 2. ESM-IF raw      │                          excluding
     │ 3. ESM z-score     │                          --excluded-
     │ 4. ESM-IF z-score  │                          positions
     └────────────────────┘
                    │
                    ▼
     Outer-join → plm_zeroshot_ensemble_nominated_mutations.csv
```

## TL;DR

The script "proposes" mutations by:
1. Enumerating **all** single-point mutations (this is the entire "proposal" step — no clever search).
2. Scoring each with an **ESM sequence ensemble** and **ESM-IF structure model** (log-likelihood ratio vs. WT).
3. Ranking the same DMS table **four different ways** (raw vs. substitution-type-normalized z-score, for each of ESM / ESM-IF).
4. Greedy-picking the top `--variants` mutations from each ranking under a **≤1-per-position** constraint and merging into a single CSV where you can see which methods voted for each mutation.

The `multievolve/proposers/` file you have open is **not used** here — that's for the supervised multi-mutant workflow.
