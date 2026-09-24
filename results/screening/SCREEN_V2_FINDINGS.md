# Screen v2: what v1 missed, and what the rescan changed

Rescan of all 384 annotated genomes with `screening/screen_genomes.py` v2 (marker sets in
`screening/markers.py`). Every claim below was checked against the annotation files, not
inferred from the code.

## Oversights in v1

### DHA export was invisible, and the donor score was inverted

v1 scored `10 * int(sld_count >= 1 or gld >= 1)`, treating periplasmic `sldAB`
and cytoplasmic `gldA` as equivalent. They are opposite: SldAB is membrane-bound
with its active site in the periplasm, so DHA appears outside the cell with no
transporter; GldA makes DHA in the cytoplasm, where it is stranded. v2 reports
`dha_export_route` as `periplasmic_oxidation`, `cytoplasmic_via_facilitator`,
`cytoplasmic_no_exporter` or `none`.

v1 then added `4 * int(not dhaklm_complete)` — a reward for *lacking* DHA
consumption — but only tested the three-subunit PTS route. `P3_C11` and `P3_H7`
carry a standalone ATP-dependent `dhaK` and collected that bonus anyway. The same
strict cutoff excluded receivers: `P3_A10` has `dhaK` and was never credited.
v2 recognises both kinase families and separates gene-name evidence
(`dha_kinase_confident`) from description-only matches, because description hits
are near-universal — scoring them made the donor role unassignable to anything.

`glpF` was recorded and never scored, though aquaglyceroporins pass DHA as well
as glycerol and are the most plausible cytoplasmic DHA efflux route.

### No lactate export marker at all

This is what caused the GEM surprise in `MODEL_GATE.md`:

| | `lctP` annotated | GEM lactate exchange |
|---|---|---|
| `P3_A3`, `P3_F11` | yes | `L_LACt2r` + `EX_lac__L_e` |
| `P3_A10`, `P3_B12` | no | none |

CarveMe mirrored the annotation faithfully. The gap was real and the screen
simply did not look. v2 emits `lactate_export_call` = `annotated` /
`unannotated_gap` / `no_ldh`.

### Aerobic pyruvate competition had zero coverage

The process is aerobic, which makes this the yield-determining axis:

| Marker | Found in | Effect under O₂ |
|---|---|---|
| `pdhABCD`/`aceE`/`lpdA` | all six panel members | pyruvate → acetyl-CoA → TCA → CO₂, not lactate |
| `noxE` | `P3_B12` only | oxidises NADH with O₂, competing with LDH; shifts *L. lactis* off homolactic metabolism |
| `poxL`/`poxB` | `P3_A10`, `P3_B12` | pyruvate → acetyl-P + H₂O₂ |

v1 penalised `poxB` but not the Gram-positive name `poxL`, so `P3_B12` escaped a
penalty it should have taken.

### A community mechanism better supported than the DHA shuttle

`P3_B12` has **no catalase** (only `ahpC`/`ahpF`/`sodA`/`tpx`) while carrying both
H₂O₂ generators, and is **heme-dependent** and **lipoate salvage-only**. Every
other panel member has `katA`. "Catalase-negative lactate producer rescued by a
catalase-positive partner under aeration" is a named, measurable community
function of the kind the README demands, and it falls out of the data.

The heme dependency compounds it: catalase is a haemoprotein, so the heme
auxotrophy is *why* the strain is catalase-negative, and supplying heme is a
medium lever that changes its respiratory behaviour.

The lipoate result points the other way and is useful: lipoate is the PDH
cofactor, and `P3_A3`, `P3_B12`, `P3_F11` carry only the salvage gene `lplA`. A
**lipoate-free minimal medium suppresses PDH and pushes pyruvate toward lactate**.

### Auxotrophy was averaged, and covered 7 of 20 amino acids

v1's `9*aa_fraction + 9*vitamin_fraction` gave a strain complete in six of seven
pathways and empty in the seventh a high supporter score, when it is an auxotroph
needing rescue. v2 reports `requires` and `provides` as separate lists and never
averages them.

Cysteine, serine, glycine, alanine, aspartate, asparagine, proline,
phenylalanine, tyrosine and the whole shikimate pathway were absent from v1.
Added, along with cobalamin, respiratory quinone, lipoate and heme.

## Corrections made during the rescan

Three of these were caught only by checking outputs against the annotations, and
each would have distorted the result:

1. **`cat` is not catalase.** It matched "Cation efflux family" in `P3_B12` and
   chloramphenicol acetyltransferase in `P3_F11`, turning a catalase-negative
   strain into a false positive. Removed.
2. **`mecA` is not methicillin resistance here.** It matched the
   Bacillus/Lactococcus ClpC protease adaptor in three genomes. Gene-name
   virulence screening collides badly (`esxA`, `spa`, `plc`, `inlA`, `sak` are the
   same problem), so the column is now `biosafety_review_required` plus
   `virulence_name_matches` for review — not a determination. A real assessment
   needs a curated VF database (VFDB/abricate), which is not installed.
3. **Requirements were over-called in two ways that inflate the medium:**
   - *Respiratory quinone*: `P1_B1`/`P3_D6` were called menaquinone auxotrophs
     while carrying 8–9 `ubi*` genes — they use ubiquinone. Now satisfied by any
     of the men/mqn/ubi routes.
   - *Cobalamin*: absence of B12 synthesis is not a B12 requirement. Every genome
     carries `metE`, the B12-independent methionine synthase. Requiring B12 only
     when a B12-dependent enzyme lacks its substitute drops the count from most of
     the pool to **24 of 384**.
4. **Lactate consumers were ranking as producers.** `lctP`/`lldP` are catabolic
   *uptake* permeases, and v2 initially credited them as export while penalising
   consumption too weakly. `P1_B1`, `P3_D6`, `MYb115` and 18 others have **zero
   L-LDH, one D-LDH and 1–4 consumption genes** — they are lactate catabolisers.
   A `lactate_direction` call now excludes `consumer_favoured` strains from the
   producer role, cutting it from 50 genomes to 29.

## Rescan outcome

Requirement frequency across 384 genomes — the real minimal-medium concerns:

| Requirement | Genomes | |
|---|---:|---:|
| biotin (B7) | 97 | 25% |
| thiamine (B1) | 83 | 22% |
| cysteine / sulfur assimilation | 77 | 20% |
| pyridoxal (B6) | 70 | 18% |
| glutamate/glutamine | 57 | 15% |
| niacin (B3) | 50 | 13% |
| heme | 45 | 12% |

Assembly quality from the previously-unused QUAST output: 272 good, 48
acceptable, 21 fragmented, 43 without QUAST.

### L-capable producers, prototrophic tier

Zero annotation-predicted requirements, L-lactate capable, producer-favoured:

| Genome | Genus | % | Score | Drain (PDH/POX/NOX) | Export | Assembly | Biosafety |
|---|---|---:|---:|---|---|---|---|
| `P1_B11` | Bacillus | 99.38 | 39.7 | 4/0/0 | annotated | good | review |
| `P1_B3` | Aneurinibacillus | 100.00 | 38.3 | 4/0/1 | annotated | good | clear |
| `P1_B8` | Aneurinibacillus | 59.43 | 38.3 | 4/0/1 | annotated | good | clear |
| `P3_A10` | Paenibacillus | 100.00 | 35.45 | 4/1/0 | **unannotated** | acceptable | review |
| `P3_B5` | Bacillaceae | 100.00 | 34.55 | 4/0/0 | annotated | good | review |

`P1_B3` is the cleanest candidate in the pool: prototrophic, L-only,
producer-favoured, annotated export, good assembly, and not in a risk genus.
`P1_B8` matches it but its genus call is only 59% — ANI required.

One requirement (biotin): `P1_F7` Raoultella (39.7, genus call 57%), `P3_F11`
Bacillus (39.7), `P3_C5` Paenibacillus (35.3, export gap).

### DHA donors with the periplasmic route

| Genome | Genus | % | Score | PQQ | sldAB | Requires |
|---|---|---:|---:|---:|---:|---:|
| `MYb596` | unresolved | – | 24 | 4 | 1 | 5 |
| `MYb592` | unresolved | – | 21 | 4 | 1 | 3 |
| `MYb595` | unresolved | – | 21 | 4 | 1 | 3 |
| `P3_C11` | Gluconobacter | 100.00 | 21 | 4 | 1 | 4 |
| `P3_C9` | Gluconobacter | 99.86 | 21 | 4 | 1 | 3 |
| `P3_H7` | Gluconobacter | 100.00 | 21 | 4 | 1 | 4 |

`P3_C9` matters because `P3_C11` and `P3_H7` carve to **identical** models, so the
panel has no independent Gluconobacter. The `MYb59x` trio carries the same
periplasmic signature but has unresolved taxonomy.

The Gluconobacters reach 12/17 cobalamin genes — far the highest in the panel, but
below the 80% completeness bar, so the screen does not credit B12 provision. Worth
testing rather than assuming.

### The original panel under v2

| Genome | Genus | Producer | Donor | Receiver | Supporter | Isomer | Roles |
|---|---|---:|---:|---:|---:|---|---|
| `P3_A3` | Staphylococcus | 43.8 | -3 | 42.5 | 19.8 | L+D | producer, receiver |
| `P3_F11` | Bacillus | 39.7 | 0 | 41.8 | 27.0 | L | producer, receiver, supporter |
| `P3_A10` | Paenibacillus | 35.45 | -3 | 34.5 | 29.8 | L | producer, receiver, supporter |
| `P3_B12` | Lactococcus | 26.4 | -3 | 27.5 | 11.0 | L | **none** |
| `P3_C11` | Gluconobacter | 0.45 | 21 | 10.8 | 16.0 | none | donor |
| `P3_H7` | Gluconobacter | 0.45 | 21 | 10.8 | 16.0 | none | donor |

**`P3_B12` loses every role.** Under aeration it combines `noxE`, `poxL`, no
catalase, heme dependence, lipoate salvage-only, eight requirements and an
unannotated lactate export. The original panel's centrepiece is the wrong choice
for an aerobic process — a reversal that only became visible once the screen
looked at oxygen.

## Known limits

- Requirements are annotation priors. `P3_B12` is scored as providing
  branched-chain amino acids, contradicting known *L. lactis* Ile/Leu/Val
  auxotrophy. The authoritative call must come from GEM minimal-medium and
  dropout analysis.
- Prototrophy is not export. A supporter must overproduce and secrete; annotated
  exporters are sparse and clade-biased, so provision needs GEM or bench support.
- `P2_G9`'s genus reads as `BIGb0383`, a strain ID leaked into the
  diamond/MEGAN taxonomy file. Genus calls below ~90% support (`P1_B8` 59%,
  `P1_F7` 57%, `P1_B4` 64%) need ANI.
