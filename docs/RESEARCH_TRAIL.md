# Literature to a candidate community panel

## Research question

Which combinations of isolates are worth testing for net L-lactate accumulation
from glycerol, and what evidence would justify testing one combination rather
than another? The aim is to generate testable biological hypotheses, not to
equate a gene annotation or a model optimum with a demonstrated phenotype.

## Read the decision process in this order

| Stage | Record | What changed or was challenged |
|---|---|---|
| Biological motivation | [Why lactate](../results/lactate_biology/WHY_LACTATE.md), [project brief](PROJECT_BRIEF.md) | Product requirements, carbon accounting and the evidence needed before claiming a community benefit |
| Literature synthesis | [Literature cross-reference](../results/literature/LITERATURE_CROSSREF.md), [pool-v3 literature review](../results/literature/LITERATURE_POOL_V3.md) | Published mechanisms and counterexamples used to challenge candidate roles and process assumptions |
| Broad genome screen | [Screen-v2 findings](../results/screening/SCREEN_V2_FINDINGS.md), [screening code](../screening/README.md) | Pathway coverage, production versus consumption, and limits of annotation-derived scores |
| Initial mechanistic panel | [Pool v2](../results/pools/POOL_V2.md) | Producer/supporter and substrate-sharing hypotheses; taxonomic reassessment affected candidate interpretation |
| Wider candidate set | [Pool v3](../results/pools/POOL_V3.md) | Broader taxonomic and functional coverage; conflicts between literature, annotations and model predictions |
| Product-purity constraint | [Pool v4](../results/pools/POOL_V4.md) | Excluding D-lactate capacity changed the available partners and their proposed functions |
| Alternative direction | [Archived D-lactate work](../archive/README.md) | A different product route was considered and abandoned; its annotation checks remained useful |
| Sequence-level review | [LDH adjudication](../results/pool_v5/LDH_ADJUDICATION.md), [racemase audit](../results/pool_v5/RACEMASE_AUDIT.md) | Re-examination of enzyme identity and possible isomer conversion, rather than accepting annotation names alone |
| Revised selection | [Pool v5](../results/pool_v5/POOL_V5.md), [manual panel](../panel/pool_v5_panel.yaml) | Retained candidates, explicit holds and reversals of exclusions that depended on unreliable model behaviour |
| Proposed validation | [Experimental plan v5](../results/planning/EXPERIMENT_PLAN_V5.md) | Experiments proposed to distinguish biological hypotheses and resolve remaining uncertainty |

The original committed record covers pools v2–v4. Later computational review and
reorganisation extend it to pool v5. Pool v1 is not fully preserved: the pool-v2
record explicitly says the earlier screen results were replaced. The documents
are successive research states, not a complete version-by-version event log.

## What makes the revisions informative

- **A plausible annotation was not always enough.** Sequence and domain evidence
  were used to revisit whether a candidate enzyme supported the intended role.
- **A nearest reference was not automatically a species identification.** The
  final record preserves unresolved identities instead of upgrading weak matches.
- **An exclusion could be reversed.** Model transport, gap-filling and redox
  problems were distinguished from evidence that the actual organism was unsuitable.
- **Uncertainty could lead to a hold rather than a forced yes/no decision.**
  The current panel records two racemase-risk holds requiring experimental resolution.
- **Selection remained a scientific judgment.** Screening scripts provide
  evidence; they do not automatically regenerate the authoritative candidate panel.

These are records of reasoning and computational checks. They do not establish
that the selected communities outperform monocultures, that the intended product
is optically pure, or that genome content alone predicts community function.

## How to interpret the computational results

The [results index](../results/README.md) identifies the different cohort sizes.
Do not treat every table in `pool_v5/` as the final panel. Metabolic-model tables
are retained as decision-support diagnostics, not trained community predictors
or experimental growth measurements. Historical confident wording should be
read alongside the later corrections and the current panel's explicit cautions.

The current computational panel is 14 producer candidates plus two held strains.
Here, “producer” is a candidate classification from the reviewed evidence, not
confirmation from a culture assay. Safety-related homology screens are also
evidence for review, not clearance to handle a strain.

## What is deliberately outside this record

Raw strain-library genomes and annotations, reconstructed model artifacts,
simulation outputs, wet-lab measurements and assay validation, physical strain
mapping, Echo assembly/run logs, and adaptive community-acquisition runs stay
outside this repository's shareable scope. No wet-lab result is needed to follow
the literature-to-shortlist decision process above.
