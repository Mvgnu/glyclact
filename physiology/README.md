# Physiology and diversity

These scripts answer bounded diagnostic questions:

| Script | Responsibility |
|---|---|
| `lactate_capacity.py` | independent product-capacity envelopes at a growth constraint |
| `lactate_consumption.py` | lactate sole-carbon capacity with oxygen, and without oxygen/alternative acceptors |
| `guild_diversity.py` | KEGG gene-content distance and functional trait summaries |
| `overflow_markers.py` | marker scan for an explicitly supplied genome list |
| `detox_screen.py` | crude-glycerol tolerance/detox marker summary |

Competing-product columns from `lactate_capacity.py` are separate optimizations, not products secreted in one shared flux solution.

The oxygen-free consumption condition contains no nitrate or fumarate. Do not use it to predict product retention in alternative-acceptor media. Likewise, gene-content clustering is a way to inspect functional redundancy; it is not a taxonomic or inclusion/exclusion gate.

The diagnostics default to `candidate_genomes/manifest_pool_v5.tsv`, so a plain run covers
the current 14 producers plus 2 held strains rather than every historical staged genome.
