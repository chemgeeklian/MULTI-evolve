## (1) Which `-er` modules matter for zero-shot prediction

The repo has four `-er` module families:

| Module | Purpose | Matters for zero-shot? |
|---|---|---|
| **`featurizers/`** | Turn sequences → numerical features | ✅ **Yes** (the zero-shot scorers live here / in `utils/zeroshot_utils.py`) |
| **`proposers/`** | Generate candidate mutants to score | ✅ **Yes** — but only implicitly. The script `plm_zeroshot_ensemble.py` does proposing inline (enumerates all single mutants); it does **not** import a `Proposer` class. The `BaseProposer` in `proposers/base_proposers.py` is an alternative, richer entry point. |
| **`predictors/`** (regressors) | Supervised regression models trained on labeled data | ❌ **No** — these need training labels. Skip entirely for zero-shot. |
| **`splitters/`** | Split labeled datasets into train/val/test | ❌ **No** — irrelevant without training data. |

So: **`featurizers` (core) + optionally `proposers`**. Ignore `predictors` and `splitters`.

## (2) What they do, and the order of execution in `plm_zeroshot_ensemble.py`

For the specific script you're running, there is **no explicit Featurizer/Proposer class** — it calls the underlying utility functions directly. Here's the actual order:

```
scripts/plm_zeroshot_ensemble.py  (main)
  │
  ├─ 1. Parse args, read wt FASTA                        → wt_seq
  │
  ├─ 2. PROPOSE (inline, implicit)
  │     Inside zero_shot_esm_dms / zero_shot_esm_if_dms:
  │     enumerate all single-residue substitutions      → ~L × 19 mutants
  │     (equivalent to a deep_mutational_scan proposer)
  │
  ├─ 3. FEATURIZE / SCORE
  │     a. zero_shot_esm_dms(wt_seq)                    → 6 ESM LM logratios
  │             └── multievolve/utils/zeroshot_utils.py L13
  │     b. for each pdb: zero_shot_esm_if_dms(...)      → ESM-IF logratio
  │             └── multievolve/utils/zeroshot_utils.py L167
  │
  ├─ 4. AGGREGATE
  │     • average logratios across ESM-IF structures    → average_model_logratio
  │     • for ESM ensemble: count total_model_pass,
  │       then sort by (total_model_pass, avg logratio)
  │
  ├─ 5. NORMALIZE (optional z-scoring)
  │     calculate_z_scores(df, normalizing_method='aa_substitution_type')
  │     → z-score each mutation within its AA-substitution type bucket
  │
  ├─ 6. SELECT / RANK (four parallel sampling strategies)
  │     sample_mutations(df, 24, excluded_positions) picks top-N mutants
  │     with one-per-position + exclusion filter, for each of:
  │        • esm (raw logratio rank)
  │        • esm_if (raw logratio rank)
  │        • esm_z (z-scored rank)
  │        • esm_if_z (z-scored rank)
  │
  └─ 7. MERGE & WRITE
        merge_mutation_dfs(...) unions the four sets, marks origin
        → plm_zeroshot_ensemble_nominated_mutations.csv
```

### Mapping to the `-er` abstractions

If you wanted to express the same pipeline in the repo's class-based style, it would look like:

| Step | Class-based equivalent | File |
|---|---|---|
| Propose all single mutants | `BaseProposer` (with `mutation_pool=None` → auto deep-mutational scan) | `multievolve/proposers/base_proposers.py` L27 |
| Score with ESM ensemble | `ZeroshotEsmFeaturizer` | `multievolve/featurizers/zeroshot_featurizers.py` (`ZeroshotBaseFeaturizer` at L7, ESM subclass below) |
| Score with ESM-IF | `ZeroshotEsmIfFeaturizer` | `multievolve/featurizers/zeroshot_featurizers.py` L576 |
| Aggregate/rank/filter | (done inline in `plm_zeroshot_ensemble.py`) | — |

**But note:** the script `plm_zeroshot_ensemble.py` deliberately bypasses the Featurizer classes and calls `zero_shot_esm_dms` / `zero_shot_esm_if_dms` directly from `multievolve/utils/zeroshot_utils.py`. So for your run, the **only files that actually execute** on the zero-shot path are:

1. `scripts/plm_zeroshot_ensemble.py` — orchestration
2. `multievolve/utils/zeroshot_utils.py` — the two scoring functions (`zero_shot_esm_dms` at L13, `zero_shot_esm_if_dms` at L167)
3. `multievolve/utils/other_utils.py` — provides `AAs`, `read_msa`, etc. (used by the ESM functions)

The `featurizers/` and `proposers/` trees are **not imported** by your CLI invocation. They're relevant only if you extend the pipeline into supervised learning later (where you'd then also need `predictors/` and `splitters/`).
