# Results index

Result files span multiple historical cohorts. A filename under `pool_v5/` does not by itself mean “the final 16”.

| Scope | Files |
|---|---|
| current manual panel: 14 producers + 2 held | `panel/pool_v5_panel.yaml`, `candidate_genomes/manifest_pool_v5.tsv` |
| 36-candidate screening cohort | `pool_v5/screening_candidates.tsv`, `pool_v5/virulence_screen.tsv`, `pool_v5/lactate_consumption.tsv` |
| 36-candidate ANI run used during pool-v5 review | `pool_v5/identity/strain_identity.tsv`, `pool_v5/identity/ani_all_hits.tsv` |
| earlier 24-genome ANI output | `identity/strain_identity.tsv`, `identity/ani_all_hits.tsv` (a rerun with current defaults will replace this scope) |
| current 16-genome functional summary | `pool_v5/guild_diversity.tsv` |
| 17-member exploratory draft (included the now-dropped supporter P3_C11) | `pool_v5/capacity_table*.tsv` |
| 89 reconstructed models | `pool_v5/model_gate.tsv`, `pool_v5/minimal_medium.tsv` |
| current interpreted decision | `pool_v5/POOL_V5.md` |

Historical tables are retained as evidence, not silently relabelled as final-panel output. The old `lactate_consumption.tsv` also contains a `coconsumes_*` column from an invalid objective; ignore that column. The repaired script no longer emits it.
