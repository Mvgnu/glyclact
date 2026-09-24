# Identity

`identify_strains.py` searches ANI references and records the closest hits. Its direct species call requires both high ANI and adequate aligned fraction; weaker hits remain closest-reference evidence for manual interpretation.

Default input: `candidate_genomes/manifest_pool_v5.tsv`.

The current checked-in `results/identity/` tables predate that default and cover 24
genomes; the 36-candidate pool-v5 run is retained under `results/pool_v5/identity/`.
Scope is indexed in `results/README.md` and should not be inferred from the directory name.

Identity is not a pool generator. The six current producers without direct species-level ANI support remain conservatively labelled in the panel YAML rather than inheriting a species name from a low-ANI nearest neighbour.
