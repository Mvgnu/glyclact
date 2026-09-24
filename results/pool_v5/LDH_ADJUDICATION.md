# LDH calls, re-decided by sequence

Pool v4's purity gate rests entirely on per-genome L-LDH and D-LDH counts, and those
counts came from eggNOG EC/KO/name fields. Two tables in this repo disagreed about them
in 19 of 260 genomes. Agreement would not have been reassuring either, because both
tables derive from the same annotation — they could only have been wrong together.

The specific worry raised against the screen was that D-LDH calls are annotation
artifacts. `ldhD` belongs to the D-isomer-specific 2-hydroxyacid dehydrogenase family
and shares its fold, PFAM architecture (PF00389 + PF02826) and His/Glu/Arg catalytic
triad with `serA`, `ghrA`/`ghrB`, `pdxB`, glycerate dehydrogenase, NAD-formate
dehydrogenase and D-2-hydroxyisocaproate dehydrogenase — several of which never touch
pyruvate in vivo. An automatic annotator handed one of those will often say
"D-lactate dehydrogenase".

That worry is correct in principle. **It does not hold for this library**, and something
else does.

## Method

All 1,687,322 proteins from all 384 genomes were searched against a Swiss-Prot panel
(`references/ldh_panel/`) that holds the targets and their decoys *together*, so a
family member is out-competed by its own family rather than merely failing a threshold.
Four lines of evidence per protein, recorded separately in
`results/screening/ldh_adjudication.tsv`:

| Evidence | What it rules out |
|---|---|
| Fold, by hmmsearch, on **both** sides | wrong-family references and wrong-family candidates |
| Best-hit class against targets + decoys | family members that are not orthologs |
| Bitscore margin over the best decoy | borderline assignments |
| Reciprocal best hit within the genome | "closest thing I have" masquerading as orthology |

The catalytic triad explains why annotation-only calls are ambiguous within this family;
it was not implemented as an independent sequence gate and is not claimed as one here.

## The reference panel had the disease it was built to diagnose

The first pass produced 186 genomes that gained an L-LDH. The cause was one reference:
**Swiss-Prot Q07251** (*Cupriavidus necator*, gene `ldh`) carries EC 1.1.1.27 but its
only Pfam domain is **PF02615 `Ldh_2`** — the LDH2/MDH2 family, a different fold from
the canonical L-LDH. It recruited 141 Enterobacteriaceae `ybiC` proteins as false
L-LDH.

Fold-typing every panel member by HMM found exactly one such contradiction in 1198
sequences, and it was that one. With the fold gate in place the disagreement with
annotation fell from 211 genomes to 59.

The lesson generalises past this project: a curated reference set assembled by EC number
inherits every EC-number error in the source, and EC number is exactly the field under
suspicion. Fold is measured, so it is the only part of this that does not depend on
somebody's label being right.

## Result: D-LDH calls were already right

| | genomes |
|---|---:|
| Agree exactly between annotation and sequence | **367 / 384** |
| Annotation said D-LDH, sequence says none (rescued) | 2 |
| Annotation said none, sequence found D-LDH | 4 |

The 4 newly-disqualified are `P1_E2`, `P3_E1`, `CEent1` and `GCF_002966775.1` — all four
have **empty or near-empty eggNOG annotation files**, so the annotator saw no genes at
all rather than seeing them wrongly. 17 of the 384 genomes are in this state, and the
old screen scored them as carrying nothing.

The two rescued genomes (`MYb382`, `P2_D8`) are genuine annotation artifacts of exactly
the predicted kind. Neither is a pool candidate.

**The four Enterobacteriaceae supporters pool v4 dropped carry real `ldhA`** — 80–91%
identity to characterized references, reciprocal-best-hit confirmed, correct fold. They
stay dropped. The purity constraint's cost to amino-acid provisioning documented in
`POOL_V4.md` is real and is not recoverable by re-reading the annotation.

## Result: L-LDH was under-counted

147 genomes carry a sequence-confirmed L-LDH against 132 by annotation. The gains are
concentrated in the unannotated genomes and in strains whose L-LDH sat below the
annotator's confidence threshold.

## *Lactococcus lactis* was never excluded for D-LDH

The suggestion to re-check a "4–2 LDH split" in *L. lactis* does not correspond to
anything in the data. `P3_B12` has **zero D-LDH** by every method and always did — both
tables agree. What it has is:

| Protein | Best characterized match | Identity | RBH |
|---|---|---:|:--:|
| `FGFGBCHD_01174` | `ldh1`, *L. lactis* IL1403 | **100.0%** | yes |
| `FGFGBCHD_00648` | `ldh3`, *L. lactis* IL1403 | 99.7% | yes |
| `FGFGBCHD_00439` | `ldh2`, *L. lactis* IL1403 | 99.7% | yes |
| `FGFGBCHD_00507` | `ldh`, *Fervidobacterium nodosum* | 35.8% | no — loses to `ldh2` |

Three L-LDH orthologs at essentially perfect identity to the canonical reference strain,
plus one divergent paralog that the old count included and this one does not (4 → 3).

Where *L. lactis* was excluded is the **D-lactate guild screen**, which required
`makes_D_lactate`. Failing that gate is the correct result for a homofermentative
L-lactate producer, and it is irrelevant now that D-lactate is off the table. In pool v4
*L. lactis* is present — as a *supporter*.

That is the finding worth acting on. `P3_B12` is not a marginal producer: it has the
cleanest L-LDH complement in the library and is the most industrially characterized
homolactic organism there is.

### Do not use its in-silico growth rate as an argument

`P3_B12`'s model reports µ = 0.179 on glycerol-minimal with **zero required
supplements** — faster than all three pool v4 producers. That number is an artifact and
should not be quoted.

*L. lactis* is a textbook multi-auxotroph: branched-chain amino acids, methionine,
histidine, arginine, glutamate and several vitamins. The model says otherwise for two
compounding reasons. The reconstruction was gapfilled on **M9**, so CarveMe was required
to make it grow on a minimal medium and added what that took; `minimal_medium.tsv` shows
it synthesising **pantothenate and pyridoxal with no genomic evidence at all**, and
biotin, folate and lipoate by transport only. Upstream of that, the annotation was
already wrong in the same direction — `SCREEN_V2_FINDINGS.md` flagged that `P3_B12`
scores as *providing* branched-chain amino acids, "contradicting known *L. lactis*
Ile/Leu/Val auxotrophy".

So the promotion to producer rests on the LDH evidence, which is sequence-level and
solid. The growth and medium-cost picture is the opposite of what the model says, and it
is the real argument against this strain: **`P3_B12` has the best L-LDH in the library
and the worst medium cost**, because it will need amino-acid supplementation that the
Bacillaceae producers do not.

### One further correction against it

`P3_B12` carries a **lactate oxidase** — `FGFGBCHD_01264`, 48.1% identity to
characterized `lctO` of *Streptococcus iniae*, e = 4.6e-113, RBH-confirmed. Pool v4
records it as having none, because the annotation did not call it. Pool v4's own filter
is "no D-LDH (hard), **no lactate oxidase (preferred)**", so this matters to the choice
even though it does not disqualify.

`lctO` is O₂-dependent and cannot turn over anaerobically, so under fermentative
operation it is idle. Under the aeration these strains need to grow on glycerol, it is
not. Together with the `noxE` and `poxL` already documented in `SCREEN_V2_FINDINGS.md`,
`P3_B12` has three separate O₂-linked routes competing with lactate formation.
