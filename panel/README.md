# Current panel

`pool_v5_panel.yaml` is the manually curated, authoritative pool definition.

It contains 14 producers and 2 `held_lar_operon` strains. The held strains require chiral-HPLC evidence before promotion. No script automatically adds or removes candidates, and virulence screens are reviewed manually.

After an intentional YAML change, refresh the staged manifest with:

```bash
python modeling/prepare_candidates.py \
  --config panel/pool_v5_panel.yaml \
  --manifest candidate_genomes/manifest_pool_v5.tsv
```

The command requires the configured strain-library source volume.
