# Staged genome inputs

Genome directories and local staging manifests are ignored from version control.
The manifests record provenance and checksums when the external source library is
available; this README describes their expected scopes, not files distributed here.

| Manifest | Scope |
|---|---|
| `manifest_pool_v5.tsv` | **current** manual panel: 14 producers + 2 held |
| `manifest.tsv` | historical 25-genome modeling cohort |
| `manifest_d_guild.tsv` | archived 45-genome D-lactate guild |

New current-panel work should default to `manifest_pool_v5.tsv`. Refresh it from `panel/pool_v5_panel.yaml` with `modeling/prepare_candidates.py`.
