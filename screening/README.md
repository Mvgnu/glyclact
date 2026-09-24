# Screening

These scripts create evidence tables. None of them edits or regenerates the current panel.

| Script | Responsibility |
|---|---|
| `screen_genomes.py` | broad marker screen across the genome library |
| `markers.py` | marker definitions used by the broad screen |
| `build_ldh_references.py` | build the target-plus-decoy LDH reference panel |
| `assign_folds.py` | HMM-based fold typing of references and candidates |
| `adjudicate_ldh.py` | four-part LDH adjudication: fold, best class, class margin, and RBH |
| `virulence_screen.py` | report VFDB/AMR homology hits for manual review |

The LDH reference builder removes ambiguous accessions before it counts the final classes.

Virulence output is intentionally not integrated into a pool selector. Inspect the hit identities, coverage, locus completeness, and biological context manually. Automated tiers are navigation aids only and include incorrect/context-dependent calls.
