# Pool v3 — a wider panel, and what the models say about it

The design brief changed on a fair objection: community behaviour stops tracking genome
content as members are added, so a pairwise design reasoned from annotation is the kind
that decays fastest. Breadth beats optimality. This expands the panel from 9 to 19
evaluated strains and cuts to **12** — with the cut models moved to `models/excluded/` so
`models/` now contains exactly the pool.

Roles below are hypotheses from annotation and FBA. Nothing here is a phenotype claim.

## The pool

| # | Genome | Species | Role | µ on glycerol | Lactate mol/mol |
|---:|---|---|---|---:|---:|
| 1 | `P1_B3` | *Priestia aryabhattai* | producer, prototroph | 0.117 | **1.00** |
| 2 | `P1_B11` | *Bacillus velezensis* | producer, fastest of the three | 0.129 | **1.00** |
| 3 | `P3_F11` | *Bacillus altitudinis* | producer, biotin-dependent (rescue test) | 0.120 | **1.00** |
| 4 | `P2_A5` | *Lactiplantibacillus plantarum* | producer, FBP-independent LDH test | 0.189 | **1.00** ¹ |
| 5 | `P3_C11` | *Gluconobacter frateurii* | DHA donor, sole aerobic B12 route | 0.127 | — ² |
| 6 | `P3_B12` | *Lactococcus lactis* | lactate specialist, heme+catalase rescue test | 0.179 | 0.00 ³ |
| 7 | `P2_G9` | unplaced Enterobacteriaceae | supporter: 9 AA exporters, **no POX, no NOX** | 0.551 | — |
| 8 | `P3_F8` | *Ewingella americana* | supporter: 23 pathways, best oxidative shield | 0.566 | — |
| 9 | `P2_D5` | *Pantoea pleuroti* | supporter + second DHA donor (5 pqq) | 0.561 | — |
| 10 | `P2_D12` | *Rahnella inusitata* | supporter, 22 pathways, 6 AA exporters | 0.564 | — |
| 11 | `P3_G8` | *Rhodococcus erythropolis* | breadth: Actinobacteria, 6 AA exporters | 0.217 ⁴ | — ² |
| 12 | `P2_E7` | *Chryseobacterium pennae* | breadth: Bacteroidetes, DHA donor, **eats no lactate** | 0.604 | — |

¹ after curation.  ² LDH assignment not genuine or ambiguous.  ³ has L-LDH but no
annotated transporter; the export gap is left in place deliberately.  ⁴ after the
oxygen-exchange fix.

Four phyla, eleven genera, four confirmed producers, four high-provision supporters,
three independent DHA donors. 4095 non-empty subsets; 781 communities at sizes 2–4.

## Who makes which isomer, and who eats lactate

This was in `genome_scores.tsv` from the start but never surfaced in a results document,
which was an omission — optical purity is a product specification, not a detail.

| Genome | Species | L-LDH | D-LDH | Isomer | `lldD`/`dld` | Eats lactate |
|---|---|---:|---:|:--|:--|:--|
| `P1_B3` | *Priestia aryabhattai* | 1 | 0 | **L only** | 0 / 0 | no |
| `P1_B11` | *Bacillus velezensis* | 1 | 0 | **L only** | 0 / 0 | no |
| `P3_F11` | *Bacillus altitudinis* | 1 | 0 | **L only** | 0 / 0 | no |
| `P3_B12` | *Lactococcus lactis* | 4 | 0 | **L only** | 0 / 0 | no |
| `P2_A5` | *Lactiplantibacillus plantarum* | 6 | 2 | **DL** | 0 / 0 | yes, via reverse LDH + POX |
| `P2_G9` | *Leclercia*-like | 1 | 1 | DL | 1 / 0 | yes |
| `P3_F8` | *Ewingella americana* | 1 | 1 | DL | 1 / 0 | yes |
| `P2_D12` | *Rahnella inusitata* | 2 | 1 | DL | 1 / 0 | yes |
| `P2_D5` | *Pantoea pleuroti* | 0 | 1 | **D only** | 1 / 0 | yes |
| `P3_C11` | *Gluconobacter frateurii* | 0 | 0 | none | 1 / 1 | yes |
| `P3_G8` | *Rhodococcus erythropolis* | 0 | 0 | none | 2 / 1 | yes |
| `P2_E7` | *Chryseobacterium pennae* | 0 | 0 | none | 0 / 0 | **no** |

**Two consequences that were not being carried forward.**

**All three *Bacillus*-type producers are L-only** — one L-LDH, no D-LDH. That is optically
pure L-lactate, which is what polylactide requires. **`P2_A5` is DL** (six L-LDH plus two
D-LDH), which is the textbook *L. plantarum* phenotype and commercially a much weaker
product. So the FBP-independence that makes `P2_A5` scientifically attractive comes bundled
with an optical-purity problem, and the two pull against each other in the anchor choice.

**`P2_A5` has no `lldD` or `dld` at all.** Its aerobic lactate consumption runs through
*reverse* NAD-LDH plus pyruvate oxidase, not a respiratory oxidase — so the edit that stops
it eating its own product is `poxB`/`poxF`, not `lldD`. `P2_E7` is the only supporter in the
pool that consumes no lactate by either route.

## Which producers are real

Any strain reporting 1.00 mol/mol had to clear a gene-provenance check, because that is
how the two *Gluconobacter* strains previously entered the producer table: CarveMe had
built their `LDH_L` from `lldD` (EC 1.1.2.3), the respiratory quinone-coupled L-lactate
**oxidase** — a lactate-*consuming* enzyme credited as production capacity.

| Genome | `LDH_L` built from | Verdict |
|---|---|---|
| `P1_B3`, `P1_B11`, `P3_F11` | `ldh`, EC 1.1.1.27, "LDH family" | **genuine** |
| `P2_A5` | 4 genes, 2 unambiguous `ldh` EC 1.1.1.27 | **genuine** |
| `P3_C11`, `P3_C9` | `lldD`, EC 1.1.2.3 | artifact |
| `P1_G3` | unnamed "FMN-dependent dehydrogenase", no EC | artifact |
| `P3_G8` | dual EC 1.1.1.27/1.1.1.37, no "LDH family" | ambiguous, probably MDH |

**`P2_A5` is the addition that matters.** *L. plantarum* LDH is the documented
**FBP-independent** bacterial L-LDH, where the *Bacillus* and *Lactococcus* enzymes need
fructose-1,6-bisphosphate to activate. Glycerol enters glycolysis below the aldolase step,
so FBP is a gluconeogenic trickle — hypothesis H3 in `OVERFLOW_TRIGGERS.md` predicts that
should specifically cripple FBP-dependent producers on glycerol. `P2_A5` against
`P1_B3`/`P1_B11`/`P3_F11` on the same medium is a direct test of that, and it brings
decades of process literature that the three *Bacillus*-type producers entirely lack.

Two caveats: it carries **`nox`** (H₂O-forming NADH oxidase, which competes with LDH for
NADH aerobically) and **five pyruvate-oxidase copies** — the first strain in the pool with
either. And **eleven auxotrophies**, which makes it the strongest beneficiary of the
supporters and the natural partner for the precursor-supply experiment.

---

# The redox gate — what it is, and why it replaced the workarounds

This section answers a fair challenge: closing sulfide and proton uptake to keep ethanol
under 1.167 mol/mol looks like committing to an arbitrary ceiling by fiat.

## The ceiling is not arbitrary

Available electrons are conserved. For a compound CcHhOoNnSsPp with net charge z,
referenced to CO₂/H₂O/NH₃/SO₄²⁻/H₃PO₄:

```
gamma = 4c + h - 2o - 3n + 6s + 5p - z
```

This is a conservation law, not a fitted constant, and it reproduces every textbook value
when run against the models' own metabolite formulas:

| | formula | charge | γ |
|---|---|---:|---:|
| glycerol | C3H8O3 | 0 | **14** |
| ethanol | C2H6O | 0 | **12** |
| L-lactate | C3H5O3 | −1 | **12** |
| acetate | C2H3O2 | −1 | 8 |
| CO₂, H₂O | | 0 | 0 |
| O₂ | O2 | 0 | **−4** |
| NH₄⁺, SO₄²⁻ | | | 0 |
| **H⁺** | H | **+1** | **0** |

So 10 glycerol carry 140 available electrons, and 140/12 = **11.67 ethanol** is the
absolute maximum regardless of pathway, thermodynamics or enzyme complement. The 1.167
figure is arithmetic.

**The last row is the decisive one.** A proton carries zero available electrons. In a
model whose reactions all conserve charge, proton uptake *cannot* supply reducing power,
and closing it would change nothing at all. That closing it changes results is therefore
not a modelling choice — it is a direct measurement of how badly the charge balance is
broken. These reconstructions carry **232–555 charge-unbalanced reactions each**.

## So were the closures legitimate?

They are two different things, and only one of them is a workaround.

- **Sulfide is a medium correction, not a patch.** `h2s` was in the mineral base as a
  sulfur source, but it is an electron donor worth 8 electrons per molecule. Our process
  feeds sulfate, not sulfide. Removing it makes the medium match the intended experiment,
  and it happens to put `P3_F8` and `P1_G3` back on exactly 1.167. Sulfate still supplies
  sulfur, and no model lost growth.
- **Closing proton uptake is a genuine workaround**, and it is now labelled as one in the
  code. Its justification is the γ = 0 argument above: it should be a no-op, so anything
  it changes was spurious. Proton *export* stays open, so fermentation products are still
  secreted with their protons.

## But closing channels was the wrong primary answer

The right answer is to **measure the electron budget and reject models that fail it**,
rather than to keep hunting for the exchange that has to be shut. That is now implemented
in `audit_models.py` as `electron_recovery`.

**How it works.** Maximise a reduced product, then sum electrons over the exchange
reactions only:

```
recovery = sum(gamma x secretion) / sum(gamma x uptake)
```

Above 1.0 is impossible. Below 1.0 is normal — the balance went into biomass.

Three design points that matter:

1. **It uses only exchange metabolites** — glycerol, lactate, ethanol, O₂, water. Those
   have clean, well-formed formulas, so the measurement is independent of the internal
   annotation errors described below. This is why it is trustworthy where per-reaction
   scanning is not.
2. **It probes several products, not one.** The first version tested ethanol only, and
   five models have no ethanol exchange — they came back blank, i.e. silently untested.
   It now tries ethanol, L- and D-lactate, acetate and succinate, and keeps the **worst**
   case: one leaking route is enough to disqualify a model, and averaging would let a
   clean probe mask a broken one. `P3_A3` was caught only by this change.
3. **Tolerance is 1.5%**, absorbing rounding and the handful of metabolites (metals,
   R-groups) whose formulas cannot be scored. Skipping those is safe in this direction —
   a skipped *uptake* can only make measured recovery too low, never too high.

## What it found

| Genome | recovery | |
|---|---:|---|
| pool members | **1.0000** | except the two below |
| `P2_E7` | 1.0077 | marginal, under tolerance |
| `P3_A3` | **1.0909** | **fails** |
| `P1_A3` (cut) | **1.0244** | **fails** |
| `P2_C2` (cut) | **1.0714** | **fails** |

**This independently flagged three models, two of which had already been cut on unrelated
grounds** — which is the useful property of a gate that measures something real.

## Can the leaks be fixed rather than gated?

Sometimes. Localising them by scanning per-reaction electron balance works only
partially, and understanding why is worth recording:

- **`P1_A3` localises cleanly.** `SUCD2` is written
  `fad + 2 h + mqn7 + succ --> fadh2 + fum + mql7` — a succinate dehydrogenase that
  reduces **two** electron carriers while oxidising one substrate. Succinate → fumarate
  releases 2 electrons; this reaction banks them twice. Blocking `SUCD2` restores recovery
  to exactly 1.0000. It is present in `P3_A3`, `P1_A3` and `P3_A10`.
- **`P3_A3` does not.** Blocking its worst offender makes recovery *worse* (1.0909 →
  1.1077) because the solver reroutes. Chasing reactions one at a time is the wrong tool.
- **`P2_C2` does not localise at all**: its per-reaction imbalances sum to zero, so the
  leak runs through a metabolite whose formula cannot be scored.

**A separate, systematic annotation error underlies much of this.** `g3p_c` is written
`C3H5O6P` with charge **0**, while its exact isomer `dhap_c` is `C3H5O6P` with charge
**−2** — the BiGG value. `g3p_c` is mis-charged by +2 in *every* model in the panel. That
makes `TPI` appear to destroy 2 electrons and `GAPD` appear to create 2. They cancel along
a linear pathway, which is why well-formed models still measure 1.0000, but any branch
point around G3P breaks the cancellation and makes per-reaction diagnostics unreliable.
Charges do not enter the stoichiometric matrix, so this affects *accounting*, not FBA
results — and it does not affect the exchange-level test at all.

## The lesson for the gate

The existing free-lunch test asks "can carbon or biomass appear from nothing". Every model
in this panel passed it while three were manufacturing **electrons**. A mass balance
cannot see that, because carbon balanced perfectly in every case — `P2_C2`'s leak produced
12.5 ethanol + 5 CO₂ from 10 glycerol, which is 30 carbon in and 30 carbon out, and 150
electrons out against 140 in. **Mass and redox are independent conservation laws and a
gate needs both.**

---

## Model amendments

`curate_models.py` applies only amendments that are (a) genuine BiGG universe reactions,
(b) backed by named genes in that strain's own annotation, and (c) logged with evidence in
`results/model_curation.tsv`. **3 of 19 models were touched.**

| Genome | Added | Justification |
|---|---|---|
| `P2_A5` | `L_LACt3`, `EX_lac__L_e` | `AIJAFFKC_00135` = `lctP`, KO K03303, annotated **"L-lactate permease"**; plus two genuine `ldh` (EC 1.1.1.27, "LDH family"). Model had only a D-lactate exchange. Lactate yield 0.00 → **1.00** |
| `P3_G8` | `O2tex`, `EX_o2_e` | had `o2_c`, periplasmic transport and cytochrome oxidases but **no oxygen exchange** — an obligate aerobe that could not respire. µ 0.085 → **0.217** |
| `P2_C2` | `O2tex`, `EX_o2_e` | same defect |

**On name collisions.** The lactate rule requires KO **and** gene name **and** EC, because
K00016 is shared with malate dehydrogenase in eggNOG — `P1_A3`'s two K00016 hits are both
`mdh`, and a looser rule would have manufactured a producer out of the sink control. The
`P2_A5` evidence was checked by hand afterwards: `lctP`/K03303/"L-lactate permease" is
unambiguous.

Deliberately **not** amended: `P3_A10` and `P3_B12` make L-lactate but have no annotated
transporter, so their export gaps remain live test cases. `P1_F12` likewise.

### A bug in the first version of this curation, and why it mattered

The added exchanges were built as plain boundary `Reaction` objects. `model.exchanges` is
resolved from the **SBO term**, not from a reaction's shape, so those reactions read as
`boundary=True` but `in model.exchanges=False`. Every routine in this project that
iterates `model.exchanges` — to close uptake, or to sum an electron budget — skipped them
silently. The symptom was `P2_A5` reporting an electron recovery of **0.0000**: the
lactate it had just been given the ability to secrete was invisible to the accounting.
Fixed by using `model.add_boundary`, which sets the term, and by reading the model's own
extracellular compartment name (`C_e`, not `e`) rather than assuming it.

---

## Lipoate: the one requirement unmet in-pool

Every requirement is covered in-pool except lipoate, and it is worth being precise about
why that is a designed trade-off rather than a gap.

**What lipoate does.** Lipoic acid is the covalently-bound swinging-arm cofactor of the E2
subunit of **pyruvate dehydrogenase** (and of 2-oxoglutarate dehydrogenase and the glycine
cleavage system). Without it, PDH cannot function.

**Who needs it.** `P2_A5` and `P3_B12` — both lactic acid bacteria — are annotated as
**salvage-only**: they can take lipoate up and attach it, but cannot synthesise it. This
is the normal LAB phenotype. No pool member has the complete de novo route, so nothing
in-pool can supply it.

**Why that is useful rather than a problem.** PDH is the main aerobic drain away from
lactate: it pulls pyruvate into acetyl-CoA and the TCA cycle, which is exactly the carbon
we want going to LDH instead. So the medium choice is a genuine fork:

- **Lipoate supplied** → both LAB can grow, *and* PDH is enabled in all ten other members,
  increasing the drain away from lactate everywhere.
- **Lipoate withheld** → PDH is suppressed pool-wide, which favours lactate, but `P2_A5`
  and `P3_B12` cannot participate.

**This is the nutritional route to suppressing PDH**, and it is why PDH deletion does not
appear as a strain-engineering target: withholding the cofactor achieves the same thing,
applies to every member at once, and is reversible mid-run — which a deletion is not. It
belongs in the medium sweep as a deliberate variable, crossed with the presence or
absence of the two LAB.

One other single point of failure worth noting: **methionine**, required by `P3_C11`, is
supplied in-pool only by `P3_B12`. Any community containing the Gluconobacter but not the
Lactococcus needs methionine in the medium.

---

## What was cut

Seven strains moved to `models/excluded/`, which now documents each with its reason. In
brief: `P3_A3` (redox gate failure plus a pending AMR screen), `P1_A3` and `P2_C2` (redox
gate failures), `P1_F12` (cannot grow on glycerol at all), `P3_C9` and `P1_G3` (functional
duplicates whose producer status was an artifact), `P3_A10` (export-gap test case,
lower priority now). All have finished models and species calls, so re-adding any is free.

### The one worth revisiting: `P3_A3` / *Mammaliicoccus sciuri*

It was previously excluded as "99.78% Staphylococcus, nearest reference *S. aureus*
USA300". **That was wrong**, and the correction is instructive: with the reference set
escalated to Staphylococcaceae, ANI places it at **99.50% to *Mammaliicoccus sciuri***
with 96% aligned, and it matches *nothing* among 68 *Staphylococcus* species. *M. sciuri*
left the genus in 2020 — the third post-2020 reclassification to matter in this project.
"Nearest reference *S. aureus*" only ever meant "least-bad hit in a reference set that did
not contain the right organism".

It is excluded again, for two different reasons:

1. **Its model fails the redox gate** (recovery 1.0909). Its reported 1.00 mol/mol lactate
   yield collapses to **0.676** once oxygen is held to the stoichiometric 0.5 per
   glycerol — the headline number was bought with electrons the model does not have. An
   earlier draft of this document called it "the standout producer"; that was wrong.
2. ***M. sciuri* is a documented `mecA` reservoir.** Direct inspection found only name
   collisions — the `mecA` hit is the ClpC adaptor, `femA`/`femB` are native peptidoglycan
   genes present in methicillin-*sensitive* staphylococci, the `msrA` hits are methionine
   sulfoxide reductase, and there is no `blaZ`, PBP2a, PVL, TSST or coagulase. **That is
   not a clearance.** eggNOG name matching is the wrong instrument; this needs
   AMRFinderPlus or CARD, and the decision is the supervisor's or biosafety officer's.

The strain's genome is otherwise genuinely clean — one fermentative `ldh`, no `lldD`/`dld`,
and no ethanol or acetate exchange in the model at all. If the AMR screen clears it and
the leak is repaired, it is a strong bench candidate.

---

# The identification rescan — limits dropped, and what it changed

The 80-reference cap and the size/GC prefilter were memory guards, not taxonomy. Both are
now unnecessary: references are sketched, indexed and queried **in batches**, each index
discarded before the next is built, so peak memory is set by `--batch-size` (default 40)
rather than by how well-sequenced a clade happens to be. **A 436-species family search now
peaks at 0.17 GB.** Both guards remain available behind `--prefilter` and
`--max-references` for a quick pass.

The GC window was the more dangerous of the two. A cap fails loudly — you can count what
was dropped. A size/GC window assumes the correct species resembles the query in exactly
the properties that vary within a genus, which is the assumption a reclassified or novel
organism breaks.

## Every call held under exhaustive search

Eleven pool members were re-searched with no cap and no filter. **Not one identification
moved** — including `P2_E7` against all **138** *Chryseobacterium* species, `P1_B11` and
`P3_F11` against all **132** *Bacillus*, and `P3_A3` against all 68 *Staphylococcus* plus
all 114 Staphylococcaceae.

**The two unplaced members stayed unplaced, and that now means something.** `P2_G9` sits
at 81.9% against all **203** Enterobacteriaceae and `P1_G3` at 84.6% against all 26
*Erwinia* plus all 90 Erwiniaceae. Previously this rested on a search that skipped 123 of
203 species; it is now a complete search, so these are genuinely unrepresented organisms
rather than an artifact of where the cap fell.

## `P1_B3` was mislabelled, and the reason generalises

Its exhaustive escalation returned **zero hits across all 436 Paenibacillaceae species**.
Its own MEGAN output explains why, by contradicting itself:

```
Order:   Bacillales 99.96%    vs  Aneurinibacillales 0.024%
Family:  Bacillaceae 99.98%   vs  Aneurinibacillaceae 0.024%
Genus:   Aneurinibacillus 100%
```

*Aneurinibacillus* belongs to Aneurinibacillaceae. The genus row says one thing; the family
and order rows say the opposite at ~100%.

**The genus row is the unreliable one.** MEGAN percentages are normalised *within each rank
row*, so a row holding a single name always reads 100% however few reads support it.
"Aneurinibacillus, 100% genus support" never meant 100% confidence — it meant "of the
genus-level assignments made, all went to Aneurinibacillus". That label had been in the
panel since pool v2. ANI settles it: nothing in Paenibacillaceae, against ***Priestia
aryabhattai* at 99.47%** in Bacillaceae — a **fourth** post-2020 reclassification to matter
here, after *Priestia* ex-*Bacillus*, *Brucella* ex-*Ochrobactrum* and *Mammaliicoccus*
ex-*Staphylococcus*.

`identify_strains.py` now reads MEGAN's family row directly and escalates there first when
it disagrees with the family derived from the genus.

## Robustness fixes the long runs forced

- **A truncated download killed a run.** `http.client.IncompleteRead` is not an
  `HTTPError`, so it slipped past the retry and ended a 436-genome search on its ~300th
  download. Now retried, with the partial file deleted rather than cached as complete.
- **Failed downloads silently shrank coverage.** A search could report "40/40" having
  indexed 37. It now counts and names what it skipped — otherwise "no species match"
  quietly reverts to meaning "we did not look", which is the failure the cap removal was
  meant to end.
- **A sequencing bug in the new escalation**: the `break` on a successful family match
  fired *before* the hits were adopted, which would have discarded the very match it broke
  for. Caught before the run.

---

## Caveats carried forward

- **Two members are taxonomically unplaced**, now on complete searches: `P2_G9` at 81.9%
  against all 203 Enterobacteriaceae, `P1_G3` at 84.6% against all 90 Erwiniaceae. Their
  roles here are annotation-based and do not depend on the species call.
- **`P2_E7` sits at 1.0077 electron recovery** — under tolerance, but it is the pool member
  closest to the line.
- **The screen's biosafety flag is advisory only.** It also flags `P1_B11`, a long-standing
  pool member, and most matches look like the familiar name-collision problem (`hlb` in
  *Bacillus*, `vanB` in *Rahnella*). Every pool member carries zero virulence name matches.
