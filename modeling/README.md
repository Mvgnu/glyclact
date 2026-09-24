# Modeling

This directory contains the CarveMe/FBA workflow.

| Script | Responsibility |
|---|---|
| `prepare_candidates.py` | stage genomes declared in the manual panel and write a checksummed manifest |
| `reconstruct_models.py` | build CarveMe models |
| `curate_models.py` | repair known exchange/transport representation problems |
| `audit_models.py` | report model feasibility, gaps, and redox diagnostics |
| `medium_diagnostics.py` | exchange and single-dropout diagnostics |

`medium_diagnostics.py` does **not** infer a biological minimal medium. A CarveMe model gapfilled on M9 can report synthesis routes absent from the organism, and single dropouts are not combinatorial minimization. Treat every model-derived growth rate, supplement requirement, secretion capacity, and product capacity as a diagnostic until experimentally measured.

All modeling scripts default to `candidate_genomes/manifest_pool_v5.tsv`; pass another
manifest explicitly for a historical cohort. This prevents the 89-model directory from
silently becoming the current panel.
