# Pool v2 — ten strains for a presence/absence community landscape

Selected from the screen v2 rescan of all 384 annotated genomes. Ten strains give
1023 non-empty subsets, or 375 communities at sizes 2–4. Copied into
`candidate_genomes/` with sha256 provenance in `manifest.tsv`; declared in
`candidate_panel.yaml`.

The v1 scan results have been discarded and replaced — they were produced by a
screen that could not see lactate export, aerobic pyruvate competition, or the
difference between periplasmic and cytoplasmic DHA oxidation.

## The pool

| Genome | Genus | % | Role | Producer | Donor | Receiver | Supporter | Isomer | Direction | Needs | Provides | Catalase |
|---|---|---:|---|---:|---:|---:|---:|---|---|---:|---:|---:|
| `P1_B3` | Aneurinibacillus | 100.00 | keystone producer/receiver/supporter | 38.3 | 5 | 40.4 | 41.8 | L | producer | 0 | 18 | 3 |
| `P1_B11` | Bacillus | 99.38 | producer, lowest aerobic drain | 39.7 | 8 | 41.8 | 40.8 | L | producer | 0 | 18 | 3 |
| `P3_F11` | Bacillus | 100.00 | producer, biotin-dependent | 39.7 | 0 | 41.8 | 39.6 | L | producer | 1 | 18 | 1 |
| `P3_A10` | Paenibacillus | 100.00 | producer, export-gap control | 35.45 | -3 | 34.5 | 42.4 | L | producer | 0 | 19 | 1 |
| `P3_G9` | Paenibacillus | 99.96 | non-consuming supporter | 33.3 | -3 | 34.0 | 42.4 | L | producer | 0 | 19 | 3 |
| `P3_C11` | Gluconobacter | 100.00 | DHA donor, B12 source | 0.45 | 21 | 10.8 | 36.4 | none | consumer | 3 | 15 | 1 |
| `P3_C9` | Gluconobacter | 99.86 | independent DHA donor | 0.45 | 21 | 10.8 | 32.4 | none | consumer | 3 | 15 | 0 |
| `P3_B12` | Lactococcus | 100.00 | rescue test case | 26.4 | -3 | 27.5 | 13.0 | L | producer | 8 | 11 | 0 |
| `P2_D12` | Rahnella | 99.90 | rich supporter, eats lactate | 22.95 | 5 | 28.5 | 52.2 | L+D | ambiguous | 0 | 22 | 3 |
| `P1_G5` | Ochrobactrum | 100.00 | lactate sink, negative control | 3.65 | 0 | 14.0 | 45.8 | none | consumer | 0 | 18 | 2 |

## Species identity (ANI) — and why it changed two decisions

Full table in `results/strain_identity.tsv`. Each assembly was placed against one
RefSeq genome per species in its genus, escalating to family when no species matched.

| Genome | Screen genus | ANI species | ANI | Aligned | Verdict |
|---|---|---|---:|---:|---|
| `P1_B3` | Aneurinibacillus | ***Priestia aryabhattai*** | 99.47 | 0.94 | **genus was wrong** |
| `P1_B11` | Bacillus | *Bacillus velezensis* | 98.81 | 0.87 | confirmed |
| `P3_F11` | Bacillus | *Bacillus altitudinis* | 98.17 | 0.95 | confirmed |
| `P3_A10` | Paenibacillus | *Paenibacillus taichungensis* | 98.54 | 0.94 | confirmed |
| `P3_C11` | Gluconobacter | *Gluconobacter frateurii* | 96.73 | 0.90 | confirmed |
| `P3_C9` | Gluconobacter | *Gluconobacter cerinus* | 96.45 | 0.90 | confirmed |
| `P3_B12` | Lactococcus | *Lactococcus lactis* | 98.23 | 0.83 | confirmed |
| `P2_D12` | Rahnella | *Rahnella inusitata* | 98.83 | 0.95 | confirmed |
| `P1_G5` | Ochrobactrum | ***Brucella pseudogrignonensis*** | 99.56 | 0.94 | **reclassified** |

The two Gluconobacters resolve to **different species**, independently confirming
that the 573-reaction difference between their models is real signal rather than
reconstruction noise.

### Both surprises were post-2020 reclassifications

The annotation's reference taxonomy predates two genus splits, and neither could be
caught by looking at gene content:

- **`P1_B3` is *Priestia aryabhattai*** (formerly *Bacillus aryabhattai*; *Priestia*
  was split from *Bacillus* in 2020). The MEGAN lineage was self-contradictory —
  order Bacillales 99.96%, family Bacillaceae 99.98%, genus "Aneurinibacillus 100%",
  though Aneurinibacillus is not in Bacillaceae. That 100% is normalised *within the
  genus row*, so a single sparse assignment reads as unanimous. Corroborating
  evidence: GC 37.86% against 43–46.5% for every Aneurinibacillus reference, and no
  ANI hit against 9 Aneurinibacillus, 60 Paenibacillaceae or all 132 *Bacillus*
  species — because *Priestia* is no longer in *Bacillus*.
- **`P1_G5` is *Brucella pseudogrignonensis*** (formerly *Ochrobactrum
  pseudogrignonense*; *Ochrobactrum* was folded into *Brucella* in 2020) — see the
  biosafety note below.

`screen_genomes.py` now cross-checks each genus against its parent family and flags
`verify_by_ani_family_mismatch`. On this pool the flag has perfect recall: all six
`clean` calls were ANI-confirmed and both wrong calls were caught, with one false
alarm (`P2_D12`, flagged only because MEGAN still places *Rahnella* in
Enterobacteriaceae rather than Yersiniaceae). The flag deliberately does **not**
overrule the genus — doing so would have discarded two correct calls.

## Biosafety: `P1_G5` must be reconsidered

`P1_G5` is a ***Brucella*** at 99.56% ANI with 94% alignment — strain-level, not
borderline. *B. pseudogrignonensis* is a documented opportunistic human pathogen
with bacteraemia case reports.

The screen scored it `biosafety_review_required = 0`, because `Ochrobactrum` was not
in the risk-genus list and the taxonomy predates the reclassification. A genus-name
lookup against a stale reference cleared a *Brucella*.

Its role — "lactate sink negative control" — is the most replaceable in the pool;
any strain with lactate consumption genes and no L-LDH would serve. **Recommend
swapping it** unless there is a specific reason to keep this isolate.

## In-pool requirement coverage

This is the structure that makes the landscape mechanistic rather than arbitrary:
every dependency below is switched on or off by adding or removing one strain.

| Requirement | Needed by | Suppliers in pool |
|---|---|---|
| biotin B7 | `P3_F11`, `P3_B12` | 8 members |
| heme | `P3_B12` | 9 members |
| **lipoate** | `P3_B12` | **none — medium only** |
| methionine | `P3_C11`, `P3_C9` | `P3_A10`, `P3_B12` |
| niacin B3 | `P3_B12` | 9 members |
| pantothenate B5 | `P3_C11`, `P3_C9`, `P3_B12` | 6 members |
| pyridoxal B6 | `P3_B12` | 9 members |
| respiratory quinone | `P3_C11`, `P3_C9` | 7 members |
| **sulfate assimilation** | `P3_B12` | **`P2_D12` only** |
| thiamine B1 | `P3_B12` | 9 members |

### Three sharp trade-offs fall out of this

1. **Lipoate is the only unmet requirement, and supplying it is double-edged.**
   `P3_B12` needs exogenous lipoate (it has only the salvage gene `lplA`). But
   lipoate is the pyruvate dehydrogenase cofactor, and PDH is the main aerobic
   drain away from lactate in every other member. A lipoate-free medium suppresses
   PDH and favours lactate, while excluding `P3_B12` entirely. Directly testable as
   a medium axis crossed with community composition.
2. **Rescuing the Lactococcus costs product.** `P3_B12`'s sulfate-assimilation
   requirement is covered only by `P2_D12`, which consumes lactate. Any community
   that makes `P3_B12` viable also contains a product sink.
3. **The DHA shuttle can be closed nutritionally.** Both Gluconobacter donors need
   methionine, pantothenate and a respiratory quinone; `P3_A10` supplies all three
   and consumes no lactate. So donor + `P3_A10` is a self-sufficient pair, and the
   donors return B12 (complete aerobic `cob*` route) plus DHA as carbon.

## Aerobic peroxide handling

`P3_B12` is the only member that generates H₂O₂ (`poxL`, `noxE`) while having no
catalase — eight members do have catalase. That makes catalase rescue a clean
presence/absence contrast rather than a confound.

| Genome | Catalase | Other peroxidase | POX | NOX | Exposed |
|---|---:|---:|---:|---:|---|
| `P1_B3` | 3 | 3 | 0 | 1 | no |
| `P1_B11` | 3 | 4 | 0 | 0 | no |
| `P3_F11` | 1 | 2 | 0 | 0 | no |
| `P3_A10` | 1 | 2 | 1 | 0 | no |
| `P3_G9` | 3 | 2 | 0 | 1 | no |
| `P3_C11` | 1 | 3 | 0 | 0 | no |
| `P3_C9` | 0 | 3 | 0 | 0 | no |
| `P3_B12` | 0 | 4 | 1 | 1 | **yes** |
| `P2_D12` | 3 | 4 | 1 | 0 | no |
| `P1_G5` | 2 | 3 | 0 | 0 | no |

## Why two strains were dropped from the old panel

- **`P3_A3`** — highest producer score in the whole pool (43.8) but 99.78%
  *Staphylococcus*, nearest reference *S. aureus* USA300. Not an industrially
  usable host. Model moved to `models/excluded/`.
- **`P3_H7`** — carves to a reaction set identical to `P3_C11` (2484 reactions,
  zero symmetric difference), so it is simulation cost with no information.
  `P3_C9` replaces it as the independent Gluconobacter.

## Status and next steps

Reconstructions: `P3_A10`, `P3_B12`, `P3_C11`, `P3_F11` already exist; `P1_B3`,
`P1_B11`, `P3_G9`, `P3_C9`, `P2_D12`, `P1_G5` are being built with CarveMe.

1. Re-run `audit_models.py` over the full ten — the model gate in `MODEL_GATE.md`
   describes the superseded six-member panel.
2. GEM minimal-medium and dropout analysis, which supersedes the annotation-based
   requirement table above. The annotation prior is already known to disagree with
   biology in at least one place: `P3_B12` scores as providing branched-chain amino
   acids, contradicting documented *L. lactis* Ile/Leu/Val auxotrophy.
3. Add product labels to the COMETS batch summarizer — it currently reports only
   biomass, survival and diversity, so it cannot measure glycerol consumed,
   lactate titre, DHA or yield.
4. Run the presence/absence landscape crossed with the medium axes that matter:
   lipoate present/absent, methionine present/absent, biotin present/absent.
