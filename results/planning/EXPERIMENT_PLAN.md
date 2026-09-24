# Experiment plan — community × medium landscape for glycerol → lactate

Two questions in one design:

1. **Can a model predict optimal community composition for lactate from glycerol?** —
   answered honestly only if a mechanistic baseline is fitted alongside the ML model on
   the same held-out data.
2. **Does it find a community with good lactate production?** — answered by Bayesian
   optimisation over the same feature space.

The design is a fixed producer plus a presence/absence vector over partners, crossed with
a binary medium vector. Ridge on main effects and pairwise interactions, then expected
improvement.

---

## 1. Feature vector

```
y  =  lactate (g/L, and mol/mol glycerol consumed)

x  =  [ s1 … s8 | m1 m2 m3 m4 ]        12 binary features
```

**Anchor (always present): `P1_B3` — *Priestia aryabhattai*.**
Highest carbon yield (0.850), most robust to under-aeration (13.54 capacity at
O₂/glycerol 0.1 versus 3.91 and 1.52 for the others), and — decisively for interpretation
— **prototrophic**. An auxotrophic anchor would confound "partner supports" with "partner
rescues", and `P3_F11` needs biotin. Fixing the anchor also removes the trivially
predictable effect (LDH presence/absence explains most of the variance and is already
known from the genomes), so the model spends its capacity on interactions, which is the
part nobody can predict.

**Partners `s1…s8`** — chosen to span one functional axis each, so a main effect is
interpretable:

| | Strain | Species | Axis it tests |
|---|---|---|---|
| `s1` | `P3_C11` | *Gluconobacter frateurii* | **substrate division** — glycerol → DHA, which the anchor can take up. The only mechanism with published lactate co-culture support |
| `s2` | `P2_A5` | *Lactiplantibacillus plantarum* | second producer; keystone organism in the published LAB pairs; FBP-independent LDH |
| `s3` | `P2_G9` | unplaced Enterobacteriaceae | max amino-acid provision (9 exporters), **no POX, no NOX** |
| `s4` | `P2_E7` | *Chryseobacterium pennae* | provisioner that **consumes no lactate**; proteolytic, distinct phylum |
| `s5` | `P3_F8` | *Ewingella americana* | max provision (23 pathways) + best oxidative shield (cat 3 / perox 4 / SOD 2) |
| `s6` | `P2_D12` | *Rahnella inusitata* | rich provisioner that **does** consume lactate — the cost-of-consumption contrast against `s4` |
| `s7` | `P3_B12` | *Lactococcus lactis* | LAB specialist, catalase-negative and heme-dependent — the rescue test |
| `s8` | `P3_G8` | *Rhodococcus erythropolis* | robustness / distant lineage; documented consortium partner |

Held out of the design: `P2_D5` (functional near-duplicate of `s5`/`s6`) and `P3_F11`
(biotin-dependent — better used as a second anchor later than as a partner).

**Medium `m1…m4`** — each binary, each testing a named hypothesis:

| | Variable | Hypothesis it tests |
|---|---|---|
| `m1` | **glucose co-feed**, ~5% of glycerol carbon | H3: raises FBP to activate LDH. *And* represses pyruvate oxidase, the route by which `P2_A5` burns its own lactate. One lever, two mechanisms |
| `m2` | **complex N** (yeast extract 4 g/L or casamino acids) | H5 benchmark: this is the thing the supporters must beat on cost. Also tests whether minimal medium drives the parasitic interactions the LAB study reported |
| `m3` | **lipoate** (2 µM) | The PDH fork: rescues `P2_A5`/`P3_B12`, but enables pyruvate dehydrogenase pool-wide |
| `m4` | **pH 6.0 vs 7.0** (buffered) | H3 (LDH works without FBP at low pH) against H4 (low pH induces `alsSD` → acetoin) |

**Aeration is a block, not a feature.** It is the single most important variable
stoichiometrically (target OTR ≈ 0.5 × glycerol uptake), but it cannot be varied
well-by-well in a plate. Run the whole design at **two shaking speeds / fill volumes** and
treat it as a blocking factor. Note the known conflict: `P2_A5` peaks at O₂/glycerol 0.6
and is harmed above it, while the *Bacillus* producers plateau from 0.5 — so one setpoint
may not suit both, and the block will show it.

---

## 2. Sizing

12 binary features → 2¹² = 4096 full factorial, which is out. Ridge on **main effects +
pairwise interactions** needs 12 + 66 = **78 coefficients**, which sets the scale.

| Round | Unique conditions | Purpose |
|---|---:|---|
| **1 — structured** | **96** | identify all pairwise terms |
| **2 — BO** | 60 | 48 expected-improvement + **12 random (honest test set)** |
| **3 — BO + validation** | 48 | 36 EI + predicted optimum in replicate |
| **Total** | **~204** | ≈ 5–7 plates at n=2, or 8–10 at n=3 |

**Round 1 composition (96):**

- anchor alone × all 4 media singly + baseline = 5
- anchor + each partner, baseline medium = 8
- **anchor + all partner pairs, baseline medium = 28** ← this is what identifies the
  interaction terms directly, no acquisition function needed
- anchor + each partner × each medium variable (fractional, D-optimal) = 24
- space-filling higher order (3–5 partners, mixed media) = 27
- on-plate controls = 4

**Why the 28 pairs matter.** With all singles and all pairs measured, every pairwise
coefficient is estimable directly. You then *predict* the 3–5 member communities and test
whether the predictions hold. That is the design-versus-ML question answered as a
measurement, not an argument — and it is affordable.

---

## 3. Measurement

### Tier 1 — every well

| Analyte | Method | Notes |
|---|---|---|
| **L-lactate and D-lactate, separately** | enzymatic, NAD⁺/LDH, A340 | Megazyme K-DLATE or equivalent; D-GPT traps pyruvate to pull the equilibrium. **Isomer separation is not optional** — the pool spans L-only, D-only and DL strains, and optical purity is a product spec |
| **glycerol consumed** | enzymatic (glycerokinase/PK/LDH), A340 | Megazyme K-GCROL. Needed for **yield**; titre alone is not the objective you want |
| biomass | OD600 (or bioluminescence, as in the previous work) | total only — says nothing about composition |
| **final pH** | plate reader / probe | cheap, and the strongest single predictor of a failed well |

All four are 340 nm plate-reader assays, so an enzymatic panel is genuinely the right
answer here. At ~200 conditions × replicates the kit cost adds up — making the reagents
in-house (LDH, NAD⁺, hydrazine/glycine buffer) is substantially cheaper at this scale.

### Tier 2 — every well, cheap addition

| **acetoin** | Voges–Proskauer colorimetric (α-naphthol/KOH, A560) | Acetoin is the competing pyruvate sink the models under-represent and the literature says dominates *Bacillus* aerobic overflow. A classic plate-readable assay — do not skip it |

### Tier 3 — ~15% subset (~30 wells, spanning the response range)

**HPLC, Aminex HPX-87H or equivalent.** One run gives lactate, acetate, ethanol, formate,
succinate, glycerol, 2,3-butanediol and 1,3-propanediol together. Two jobs: cross-validate
the enzymatic assays, and **close the carbon balance**. 2,3-BD has no convenient plate
assay and is an 8.5 mmol/gDW/h capacity in two of the producers, so it needs this tier.

### Tier 4 — top ~20 and bottom ~20 performers only

**Community composition** — 16S amplicon or strain-specific qPCR. Expensive per sample and
*not needed for the ridge/BO model*, whose features are the inoculum vector. It is needed
to explain **why** the winners won, and to fit a gLV if you want one later.

---

## 4. Practical failure modes to design against

1. **pH crash.** Lactate acidifies; in unbuffered microplates growth stops early and every
   well converges to the same low number. **Buffer properly** (100 mM MOPS or phosphate) or
   add CaCO₃. This is the single most common way plate-scale lactate work fails.
2. **Growth-rate mismatch.** The enterobacteria grow at µ ≈ 0.55 h⁻¹ against ≈ 0.12 for the
   *Bacillus* producers. Fast growers will dominate regardless of design. Inoculate at
   **equal OD**, record it as a fixed choice, and expect composition ≠ inoculum.
3. **Plate effects.** Randomise well positions, put a **fixed reference community on every
   plate** for plate-to-plate normalisation, and use edge wells for blanks only.
4. **Endpoint ambiguity.** Communities finish at different times. Sample **two timepoints**
   (e.g. 24 h and 48 h) rather than one, or the model learns growth rate rather than yield.
5. **Pure glycerol first.** Crude glycerol adds methanol, salts and soap; add it as a later
   block once the landscape is known.

---

## 5. Modelling

**Objective.** Use **lactate per glycerol consumed** (yield), not titre, as the primary
response. Titre rewards whatever grew most; yield is the process-relevant quantity and is
what every model number in this project is expressed in.

**Models to fit, on identical data and identical held-out folds:**

| Model | What it represents |
|---|---|
| **M0 — mechanistic baseline** | predicted lactate = anchor capacity, adjusted only by known members' capacity and consumption genes. **Zero fitted parameters.** |
| **M1 — main effects ridge** | 12 coefficients |
| **M2 — ridge + pairwise** | 78 coefficients; the primary model |
| **M3 — GP / BO surrogate** | for expected improvement |

**The comparison M0 versus M2 on held-out data is the answer to the design-versus-ML
question.** If M0 is competitive, deliberate design was sufficient and that is a real,
publishable negative result. If M2 beats it substantially, the gap is the size of the
non-derivable interaction landscape — a number nobody currently has.

**Reporting.** Hold out ~20% of round 1, never used for training or acquisition, and
report Spearman and precision@k on it. Do **not** report BO's best-found as evidence the
model learned anything — those are different claims, and conflating them is the main way
this kind of study goes wrong.

**Acquisition.** Expected improvement on M3, but constrain the space: exclude
partner-only communities without the anchor (no anchor → no design), and cap community
size at 5 unless a prediction says otherwise.

---

## 6. Optional second block — worth the cost

Repeat round 1 with **`P2_A5` (*L. plantarum*) as anchor** instead of `P1_B3`.

Two returns for one extra plate set: it directly tests the FBP-independence hypothesis
(H3) — an FBP-independent LDH should be *less* penalised on glycerol than the
FBP-dependent *Bacillus* enzymes — and it tests whether the learned interaction structure
transfers between anchors, which is the strongest available check on whether the model
learned biology or fitted one organism's quirks.
