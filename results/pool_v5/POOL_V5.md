# Pool v5: current adjudicated L-lactate panel

## Decision

The current panel is **14 producers plus 2 strains held for chiral-HPLC resolution of a lactate-racemase risk**. There is no non-producing supporter.

The authoritative definition is `panel/pool_v5_panel.yaml`. It is manually curated and is not regenerated from a Boolean selector. `candidate_genomes/manifest_pool_v5.tsv` is the matching 16-genome staged-input manifest.

## Why the decision is manual

The evidence streams do not all operate at the same confidence level:

| Evidence | What is reliable | How it is used |
|---|---|---|
| LDH adjudication | fold, target-plus-decoy best hit, class margin, and RBH | hard purity evidence |
| racemase audit | LarA-domain and `larBCDE` architecture | hold ambiguous complete systems for chiral HPLC |
| ANI | species calls only with adequate identity and aligned fraction | conservative label/relatedness evidence |
| virulence/AMR homology | underlying hits and complete loci after inspection | manual biosafety judgment; automated tiers are not gates |
| CarveMe/FBA | route feasibility and reconstruction defects | diagnostic only |
| gene-content diversity | functional similarity and unusual accessory traits | manual redundancy/diversity aid |

## Current producers

| Genome | Conservative identity | Main reason retained | Main caution |
|---|---|---|---|
| `P2_C2` | *Oerskovia merdavium* | clean LDH/racemase evidence; distinct Actinomycetota physiology | model redox leak |
| `P3_D2` | *Microbacterium* sp.; closest *M. forte* | fast modeled glycerol growth; distinct guild member | aerobic lactate-use capacity; species unresolved |
| `P2_C10` | *Microbacterium* sp. | distinct producer architecture | strongest modeled aerobic lactate-use signal; species unresolved |
| `P3_G8` | *Rhodococcus erythropolis* | strongest crude-glycerol tolerance profile | oxygen/acceptor-linked lactate oxidation risk |
| `P3_C3` | *Lactococcus lactis* | three near-identical characterized L-LDH orthologs; no D-LDH | auxotrophy/model mismatch and `lctO` |
| `P1_B11` | *Bacillus velezensis* | clean sequence evidence and robust biosynthesis | manually reviewed virulence-family hit |
| `P1_F7` | *B. altitudinis* group | faster member of a close within-group comparison | modeled aerobic lactate use; direct species call incomplete |
| `P3_F11` | *Bacillus altitudinis* | clean slower comparator; no modeled lactate sole-carbon growth | near-redundant with `P1_F7` |
| `P3_A9` | *Paenibacillus taichungensis* | robust biosynthesis; deliberate transport-model challenge | glycerol transport is not represented reliably |
| `P3_A8` | *Paenibacillus amylolyticus* | distinct Paenibacillus producer | marginal model redox leak |
| `P3_C5` | *Paenibacillus odorifer* | oxidative glycerol route (`gldA` + `dhaK`) | model lacks a defensible glycerol transport route |
| `P3_G9` | *Paenibacillus* sp.; closest *P. lupini* | deliberate test of a model no-growth false negative | species unresolved; model behavior unreliable |
| `P2_D9` | Bacillaceae sp.; closest *N. tiangongensis* | strong detox/accessory profile | low-ANI identity and model redox leak |
| `P2_B6` | *Pseudomonas* sp. | only Pseudomonadota producer; major phylogenetic spread | species unresolved; defined supplements predicted |

Full strain-level rationales remain in the YAML so the decision and its exceptions stay together.

## Held, not included in the 14-member search

- `P1_B3` — *Priestia aryabhattai*
- `P1_B8` — *Priestia megaterium*

Both carry a `larBCDE` maturation system plus an ambiguous Lar_N-domain protein. That is not proof of racemase activity, but it is exactly the purity failure mode sequence-only exclusion could miss. A monoculture chiral-HPLC measurement decides promotion.

## Where the screens agree

All 16 current/held genomes have the intended L-LDH/D-LDH sequence evidence after fold-aware target/decoy adjudication. The final panel does not depend on eggNOG naming alone.

The broad marker screen, dedicated racemase audit, and manual panel agree that the two Priestia strains need a purity hold. Gene-content analysis also supports carrying broad functional spread while treating the `P1_F7`/`P3_F11` pair as an intentional close comparison.

## Where they disagree, and the resolution

- **Annotation versus LDH sequence adjudication:** automatic family labels disagree with the fold-aware screen. Keep the sequence adjudication.
- **KO-only versus domain/operon racemase screening:** KO-only screening undercounts plausible systems. Keep the domain/operon audit and experimental hold.
- **Nearest ANI hit versus taxonomic label:** low-ANI nearest neighbours are not species identifications. Keep conservative “sp.” labels.
- **Automated virulence tier versus manual locus reasoning:** generic environmental homologs and partial systems create incorrect calls. Keep the hits; keep the manual decisions; do not integrate the tier into selection code.
- **FBA ranking versus biological plausibility:** M9-gapfilled models previously excluded plausible producers for transport, growth, or redox artifacts. Keep FBA as diagnostics and reverse those automated exclusions.
- **Phylogenetic dereplication versus functional diversity:** pairwise ANI is useful evidence but is not needed as a second automatic pruning layer here. Use guild gene-content diversity for functional redundancy and retain deliberate close comparators where experimentally informative.

## Process condition is unresolved

The 35/36 oxygen-free result in `lactate_consumption.tsv` used no oxygen **and no alternative respiratory acceptor**. It does not show that nitrate or fumarate will protect product.

Quinone/cytochrome-linked L-lactate oxidation can feed anaerobic respiration. Nitrate or fumarate may improve redox balance while simultaneously reopening lactate consumption. Acetate is also a plausible sink, but transfer from the cited organism to this panel is not established. Measure net L-lactate retention in each exact proposed medium before fixing the process condition.

## Consolidation outcome

Retained and repaired:

- LDH reference/fold/adjudication scripts
- conservative ANI hit reporting
- manual virulence evidence generation
- model reconstruction/curation/audit diagnostics
- lactate capacity/consumption diagnostics
- guild diversity and crude-feed marker summaries

Retired:

- automatic pool-v5 selector
- standalone ANI dereplication candidate selector
- obsolete hard-coded community designer
- unused LDH `anchors.json`
- the invalid lactate co-consumption optimization

Historical result tables are indexed by cohort in `results/README.md`; they are evidence, not current-panel definitions.
