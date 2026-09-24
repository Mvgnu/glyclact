# Pool v4 — the pure L-lactate pool

Optical purity treated as a hard specification rather than a soft penalty. PLA-grade
L-lactate needs ≥98–99% optical purity, so a member carrying a **D-LDH** can take a batch
off-spec no matter how much total lactate the community makes. Consumption is a different
class of failure: a member that eats L-lactate costs titre but leaves the product on-spec.

So the filter is: **no D-LDH anywhere** (hard), **no lactate oxidase** (preferred).

Counts recomputed directly from the eggNOG annotations in `results/purity_scan.tsv`, not
from the screen's score columns — those disagreed on L-LDH presence in **19 of 260**
genomes, and the whole argument rests on the counts being right.

## The pool

**Producers — all L-only, none consumes lactate**

| Genome | Species | µ glycerol | Lactate mol/mol | e⁻ recovery |
|---|---|---:|---:|---:|
| `P1_B3` | *Priestia aryabhattai* | 0.117 | 1.00 | 1.0000 |
| `P1_B11` | *Bacillus velezensis* | 0.129 | 1.00 | 1.0000 |
| `P3_F11` | *Bacillus altitudinis* | 0.120 | 1.00 | 1.0000 |

**Supporters — six that clear purity, growth and the redox gate**

| Genome | Species | Clean? | Prov | AA exp | Cat | DHA | µ glycerol |
|---|---|:--|---:|---:|---:|:--|---:|
| `P3_B12` | *Lactococcus lactis* | **clean** | 11 | 1 | 0 | — | 0.179 |
| `P2_E7` | *Chryseobacterium pennae* | **clean** | 18 | 1 | 2 | periplasmic | 0.604 |
| `P2_G8` | *Chishuiella changwenlii* | **clean** | 15 | 0 | 1 | periplasmic | 0.403 |
| `P3_C11` | *Gluconobacter frateurii* | safe, eats | 15 | 0 | 1 | PQQ ×4 | 0.127 |
| `P3_G8` | *Rhodococcus erythropolis* | safe, eats | 20 | **6** | 3 | — | 0.217 |
| `P3_G12` | *Brevibacterium aurantiacum* | safe, eats | 18 | 3 | 2 | — | 0.211 |

**Nine members: 3 producers, 3 clean supporters, 3 purity-safe consumers.** Zero D-LDH
anywhere in the pool.

Three independent DHA-donor mechanisms survive the filter — PQQ quinoproteins in
`P3_C11`, periplasmic glycerol dehydrogenases in `P2_E7` and `P2_G8` — so the substrate-
division axis, the only one with published lactate co-culture support, is intact.

## What the purity constraint costs

Four of the original eight supporters carry a D-LDH and had to go:

| Dropped | Species | AA exporters | Why |
|---|---|---:|---|
| `P2_G9` | unplaced *Leclercia*-like | **9** | D-LDH |
| `P3_F8` | *Ewingella americana* | **7** | D-LDH |
| `P2_D5` | *Pantoea pleuroti* | 6 | D-LDH |
| `P2_D12` | *Rahnella inusitata* | 6 | D-LDH |

Those were the pool's four best amino-acid providers. Total annotated AA-exporter count
across the supporter bench falls from **36 to 11**, and the peak from 9 to 6 (`P3_G8`).

**This is the real price of the purity constraint, and it is not small.** The H5
precursor-supply hypothesis — that partners raise lactate by relieving biosynthetic burden,
the mechanism worth up to 38% of production cost as a yeast-extract substitute — is the
hypothesis most weakened. It can still be tested against the yeast-extract medium arm, but
the strain-based version of it is now much thinner. The D-LDH risk sits precisely in the
Enterobacteriaceae, which is where the amino-acid provisioning lives.

**All four are recoverable with two deletions each**, and they are the most tractable hosts
in the pool — see `GENOMIC_MODIFICATION_OPPORTUNITIES.md` §5. Each carries a single `ldhA`
(D-LDH, EC 1.1.1.28) and a single `lldD` (L-lactate oxidase, EC 1.1.2.3), both dispensable,
no `dld` and no redundancy. Deleting the pair makes each strain purity-safe *and* stops it
consuming the product. So the choice is not "purity or provisioning" permanently — it is
"purity or provisioning **without strain work**", which is a different and much better
decision to be facing.

## Two proposed replacements failed the gate

Both were picked on annotation and did not survive modelling:

- **`P2_E9`** (*Leucobacter*-like) — chosen as the second-best purity-safe AA provider at 5
  exporters. **Cannot grow on glycerol at all**: µ = 0.0000 aerobic and anaerobic, and only
  0.0639 on glucose, which is low enough to suggest a heavily gapped reconstruction. Its
  ANI is also unplaced (83.9% to *L. muris*, 68% aligned). Rejected on all three counts.
- **`P3_A8`** (*Paenibacillus amylolyticus*, 95.4%) — the cleanest option on paper: 18
  pathways, two L-LDH, no D-LDH, no oxidase, only needs biotin. **Fails the redox gate at
  1.0161.** That is marginal — the tolerance is 1.015 — unlike `P2_C2`'s 1.0714, so it is
  the one worth trying to repair rather than discard. Excluded for now.

## Alternates on the bench

| Genome | Species | Why held back |
|---|---|---|
| `P3_G9` | *Paenibacillus* | Best profile in the whole qualifying set — 19 pathways, catalase 3, **zero requirements**, clean. **Confirmed cannot grow on glycerol** (µ 0.0000; 1.19 on glucose). Only usable as an obligate cross-feeder. |
| `P3_A8` | *P. amylolyticus* | Marginal redox-gate failure (1.0161). Repairable. |
| `P2_C2` | *Oerskovia merdavium* | The only clean qualifier with both `glpK` and `glpF`, 19 pathways — but redox recovery 1.0714 and the leak does not localise. Strain probably fine, model is not. |

## Design consequence

Nine members, with at least one producer present:

- producer subsets (≥1 of 3): **7**
- supporter subsets (6): **64**
- **448 communities** per medium condition

At 6 supporters the 64 subsets still cover every main effect, all 15 pairs, all 20 triples
and beyond — a complete factorial, so every interaction term is estimable rather than
aliased. That property was lost at 8 supporters, where the pairwise block had to be
fractional.

## Open question this pool is built to answer

The landscape currently running (1792 communities on the v3 pool, glycerol-minimal and
yeast-extract arms) already shows the producer-only result: **total biomass is conserved at
~0.24–0.26 regardless of how many producers are present**, and **0 of 7 producer
combinations beat the best monoculture** (`P1_B11`, 0.258). Producers split the glycerol
rather than adding to it.

So the whole case for a community rests on whether a *supporter* raises producer biomass
above 0.258. Pool v4 is the version of that question where the answer, if positive, gives a
saleable product.
