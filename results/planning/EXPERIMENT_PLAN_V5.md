# Experiment plan v5 — random communities, then EI

Updates `EXPERIMENT_PLAN.md`; does not replace it. **Sections 3 (measurement tiers), 4
(failure modes) and 5 (modelling) of that document stand** — the enzymatic L/D panel, the
acetoin assay, the HPLC subset and the composition tier are all still right.

Three things changed from the first draft of this plan, and each removed something rather
than adding it.

---

## The rule that decides what goes in the search space

**A factor belongs in the acquisition space only if you would accept it in the final
process.**

This is not a stylistic preference. Expected improvement optimises whatever it is given,
and it optimises it *correctly*. Put yeast extract in the feature matrix and EI will find
its optimum with yeast extract in the well — because yeast extract genuinely helps. The
run then succeeds on its own terms and hands back a community that only works in a medium
the project exists to avoid. The optimiser did its job; the design asked the wrong
question.

So each candidate factor is sorted by one question — *would we ship this?*

| Factor | In the search space? | Why |
|---|---|---|
| the 14 members | **yes** | this is the thing being designed |
| N source, buffering | **only if block 0 says they interact** | genuinely interchangeable axes — see below |
| acetate | **no — evaluate then fix in block 0** | plausible redox mechanism, but the sign and dose must be measured in this panel/feedstock |
| ~~fumarate~~ | **drop** | adds carbon *and* its product succinate is already on the analyte panel |
| ~~crude vs pure glycerol~~ | **not a factor** | crude is the only substrate. Pure glycerol is a solved problem and reproducing it gains nothing |
| **yeast extract** | **no — reference arm only** | this is the cost we are trying to remove. A fixed ceiling to beat, never a variable to optimise |
| **dissolved O₂** | **no — fixed condition** | not settable reproducibly in microplates; see below |

### Which medium axes are BO-shaped, and which are not

The criterion is not "medium in or out". Plenty of published work puts medium in BO and is
right to — but for **genuinely interchangeable elements**, alternatives with no strong prior
ordering. Not for known improvements where you already know the sign.

Sorted that way:

| Axis | Shape | Why |
|---|---|---|
| **acetate** | directed validation | ¹³C tracing establishes a mechanism in another organism; transfer, dose, and net product retention still need block-0 measurement |
| **nitrogen source** | **interchangeable** | crude glycerol is nitrogen-free, so N *must* be added regardless. Ammonium sulfate, ammonium chloride, and urea are genuine alternatives. Nitrate belongs here only for strains shown to assimilate it; respiratory `narGH` alone does not establish that |
| **buffering / pH strategy** | **interchangeable** | lactate crashes the pH and something must hold it. CaCO₃, phosphate, base titration, Ca(OH)₂ — genuinely different, large effect, and the downstream consequences differ (CaCO₃ precipitates calcium lactate, which changes the whole recovery step) |
| **manganese** | marginally interchangeable | Mn matters for LDH activity in LAB and crude glycerol's metal content varies. Cheap either way, low stakes |

### The real question is whether medium and community interact

If the best medium is roughly the same whichever community is in the well, medium and
membership can be optimised **sequentially** at almost no loss — settle medium, fix it,
search members. If the medium optimum *flips* with composition, sequential optimisation
finds a good community for the wrong medium.

That is testable rather than arguable, and block 0 is where it gets tested.

**Decision rule:**

- **Ranking stable across communities** → decouple. Fix the medium, run BO over the
  14-member vector alone. Cheapest path, and it keeps the space at 2¹⁴.
- **Ranking flips** → the interaction is real. Fold N source and buffering into the feature
  vector, knowing you are paying for something that exists.

Default to sequential. Every binary medium dimension roughly doubles the space — two of them
takes 16,384 to ~65,000 on the same well budget, which thins coverage badly. Block 0 is what
earns the right to deviate.

---

## Oxygen is a fixed condition, not a factor

The previous draft made O₂ the primary factor. That was wrong twice over.

**Practically**, oxygen in shaken microplates is not controllable to a level worth
treating as a design variable — the transfer rate depends on fill volume, shaking orbit,
plate geometry and evaporation, none of which stay constant across a screen. A factor you
cannot set reproducibly is noise wearing a label.

**Scientifically**, glycerol → pyruvate produces more reducing equivalents than pyruvate →
lactate consumes, so a sink is needed. Low oxygen is a plausible operating point, but it
does not by itself settle the medium. The existing model result called “anaerobic” removed
oxygen **and every alternative respiratory acceptor**.

Nitrate and fumarate are therefore tradeoffs, not automatic fixes. `narGH` and `frdAB`
support respiratory electron transfer, but quinone/cytochrome-linked L-lactate oxidation
can use the same respiratory chain. Adding an acceptor may improve redox balance while
also reopening product consumption. Block 0 must measure net L-lactate retention in the
exact acceptor-containing medium; gene presence alone does not decide the sign.

### Acetate is a strong candidate sink, not a settled answer

*Enterococcus faecalis* QU 11,
fed-batch on **30 g/L glycerol + 10 g/L acetic acid**, gives **55.3 g/L optically pure
L-lactic acid at 0.991 mol/mol glycerol**
([Kobayashi et al., J Biosci Bioeng 2015](https://pubmed.ncbi.nlm.nih.gov/26168904/)).
¹³C₂-acetate traced by GC-MS shows the mechanism directly: acetate is reduced to
**ethanol**, which reoxidises the excess NADH that glycerol → lactate generates. Total
carbon yield 1.06 C-mol product per C-mol substrate — the acetate carbon leaves as ethanol
rather than CO₂, hence "without carbon loss".

**That is pure glycerol on a rich medium, not crude.** mMRS: 10 g/L peptone, 8 g/L beef
extract, 4 g/L yeast extract, Tween 80. And *E. faecalis* is an opportunistic pathogen and
VRE reservoir. **The mechanism transfers; the process does not.** An earlier version of
this section claimed a crude-glycerol result — that came from a search summary and is not
in the paper.

That is the same redox problem this plan is built around, solved with a cheaper and better
reagent than either salt:

- **Carbon-neutral in the accounting sense** — the added carbon exits as a second saleable
  product rather than confounding the balance the way fumarate→succinate does.
- **Cheap**, and acetate is a fermentation product that multiple pool members can make,
  so it may partly self-supply in a community.
- Works **aerobically** as well: a companion paper reports lactate fermentation as the main
  aerobic route in *E. faecalis* on high glycerol
  ([Appl Microbiol Biotechnol 2018](https://pubmed.ncbi.nlm.nih.gov/30232536/)), which
  loosens the oxygen constraint rather than tightening it.

**Which pool members can run it.** The route is acetate → acetyl-CoA → acetaldehyde →
ethanol, needing `ackA`/`pta` (or `acs`) plus the bifunctional `adhE`. `adhE` is present in
**7 of 14**: `P3_G8`, `P2_C2`, `P3_A8`, `P2_D9`, `P3_A9`, `P3_C3`, `P3_C5`. Four of those
also carry `pta` and `ackA` for the activation step: **`P3_C3`, `P2_C2`, `P2_D9`, `P3_G8`**.

**`P3_C3` (*L. lactis*) is the one to watch.** It has `adhE`, `pta`, `ackA`, the most `ldh`
copies in the pool, and it is the GRAS/QPS organism that *E. faecalis* — an opportunistic
pathogen and VRE reservoir — is not. If the acetate trick transfers, it transfers to the
member that was already the best producer on paper.

**Design consequence: acetate is the first sink to test in the medium factorial.** The
mechanism is plausible, but transfer to these strains and this crude feedstock is not
established. Settle dose, net lactate retention, and by-products in block 0 before fixing
it for the community search.

### Nitrate and fumarate — neither is present in crude glycerol, and they are not equivalent

Crude glycerol carries no nitrate — transesterification uses NaOH/KOH neutralised with HCl
or H₂SO₄, so the ash is chloride and sulfate and **no nitrogen enters the process at all**.
It carries no fumarate either; the organic non-glycerol fraction is soap, FAME, free fatty
acids, mono/diglycerides and methanol. Both are deliberate supplements.

**Nitrate is not a clean lever by default.**

- Contributes no carbon, so it does not touch the mol-lactate / mol-glycerol denominator.
- Cheap at fertiliser grade.
- It doubles as the nitrogen source only where assimilatory nitrate use is demonstrated;
  respiratory `narGH` is not evidence of nitrogen assimilation.
- `narGH` is O₂-repressed and anaerobically induced, which is the regime being run.
- The same respiratory chain can accept electrons from quinone-linked lactate oxidation,
  so nitrate may increase rather than suppress product re-consumption.
- Watch for nitrite accumulation, which is toxic and inhibitory; and if denitrification
  runs to completion the nitrogen leaves as N₂O/N₂ rather than reaching biomass.

**Fumarate is the messy one, despite covering more members.**

- **It is a carbon source.** Adding it adds non-glycerol carbon and confounds the yield.
- **Its reduction product is succinate**, which is already on the competing-product analyte
  panel. Dosing in fumarate means dosing in the precursor of something being measured as a
  competing product, and the carbon balance stops closing cleanly.
- More expensive than nitrate.

**So nitrate enters block 0 only as an explicitly labelled acceptor condition**, paired
with a nitrogen source unless nitrate assimilation has been shown. It advances only if
the exact medium improves net L-lactate retention. Fumarate remains dropped from the main
factorial; if tested later it needs separate carbon and succinate accounting.

### There is no co-product route to fall back on

The classic way to balance glycerol fermentation is to send part of the feed to 1,3-PDO —
that is how *Klebsiella* and *Clostridium* processes close their redox. **No member of this
pool can do it: zero 1,3-propanediol genes across all 14 producers.** So the excess NADH has to leave
via an external acceptor (O₂, nitrate) or via the acetate → ethanol route above, which is
the same trade made deliberately rather than by accident,
which costs lactate yield directly. That is a real constraint on the design and it is not
negotiable by strain choice within this pool.

**Fix the condition: sealed or low-headspace plates, static or gently shaken, one regime
throughout.** Record the actual condition; do not pretend it is a variable.

---

## No monoculture screening block

The previous draft opened with 17 monocultures as a gate into the community rounds. That
is the wrong experiment, and the reason is in this project's own reference case.

**LS30 showed the consortium outperforming its members on lag phase and production rate.**
If monoculture performance predicted community contribution, there would be no consortium
literature to cite. A strain that is slow or feeble alone can be the one that shortens lag,
detoxifies an impurity, or supplies a cofactor — and those are precisely the strains a
monoculture gate deletes. Screening on monocultures selects *against* the mechanism that
makes a community worth building.

It also cannot be rescued by arguing the gate is cheap. A cheap filter that removes the
right answer is worse than no filter.

**So: go straight into communities.** Random sampling first, then the kernel, then EI.

### The one thing monocultures are still needed for

Attribution of optical purity for `P1_B3` and `P1_B8`. If a community containing `P1_B3`
produces D-lactate, that is not attributable — the `lar` question needs the strain alone.

That is **two wells, not a block**. Run `P1_B3` and `P1_B8` as monocultures alongside round
1 purely for the chiral HPLC readout, and let them into the pool from round 2 if clean.
Everything else goes in from the start.

---

## The design

**Feature vector: 14 binary members. Nothing else.** 2¹⁴ = 16,384 candidate points.
`P1_B3` and `P1_B8` join from round 2 if their chiral HPLC clears, taking it to 2¹⁶.

The medium is **fixed** at whatever block 0 selects, so it is a constant and not a
dimension. That halves the space against the earlier draft and puts every acquired sample
on the question that is actually unknown.

| Round | Plates | Unique conditions | n | Design wells | Composition |
|---|---:|---:|---:|---:|---|
| **0 — medium × community interaction** | 1 | 36 | 2 | 72 | see below |
| **1 — random** | 2 | 88 | 2 | 176 | uniform over community size 2–14. The 2 monoculture purity wells sit in the control block |
| **2 — EI** | 1 | 44 | 2 | 88 | 33 expected improvement + **11 uniform random** |
| **3 — EI** | 1 | 44 | 2 | 88 | 33 EI + 11 random |
| **4 — confirm** | 1 | 22 | 4 | 88 | top predictions and best observed, in replicate, ± yeast-extract reference |
| | **6** | | | **512** | plus 48 control wells = **560 wells** |

Round sizes are **88 design wells + 8 control wells = 96** per plate. An earlier version of
this table counted design wells only and so allocated no room for blanks, the glycerol-free
control, or the plate-normalising reference community — which would have pushed the campaign to
ten plates once controls were added. Sizing the rounds to 88 keeps it at six. The cost is
trimming round 1 from 96 conditions to 88 and rounds 2–4 by four and two, which changes
nothing that matters.

**Bench execution is in `protocols/SOP — Community Assembly and Crude-Glycerol Fermentation
Screen.md`** — medium, Echo picklist and volumes, plate layout, O₂ regime, harvest, and the
validation set that must pass before block 0.

### Block 0 — one plate that decides the shape of everything after it

**Three deliberately different communities**, chosen to span the pool rather than to win:

- Actinomycetota-heavy — `P2_C2` + `P3_D2` + `P3_G8`
- Bacillales-heavy — `P1_B11` + `P3_F11` + `P3_A9`
- mixed across groups — `P2_C2` + `P3_C3` + `P1_B11` + `P2_B6`

Crossed with **± acetate**, two assimilable nitrogen sources (ammonium salt and urea), and
a separately labelled nitrate-acceptor arm that retains an assimilable N source, at one
buffering condition. Roughly 36 unique conditions at n = 2.

It answers three things at once:

1. **Does a proposed sink improve net L-lactate retention** — acetate and nitrate are read
   directly rather than assumed to have the same sign.
2. **Which assimilable N source** works best without confusing respiration with assimilation.
3. **Does the medium ranking depend on the community** — the decision rule above. Same
   ordering across all three → decouple and keep the search at 2¹⁴. Ordering flips → fold
   the medium axes in.

**And it is the cheapest possible failure detector.** If none of the three communities makes
measurable L-lactate on crude glycerol under any medium condition, that is one plate spent
rather than five hundred wells — and it is a real result, not a wasted run. A defined panel
assembled from genomes against a product specification, which then does not produce on the
target feedstock, says something specific about the limits of genome-driven design that an
enrichment study cannot say. It would be a worse thesis and a perfectly publishable one.

**Why 96 random to start.** A GP over 14 binary dimensions needs an initial design that
spans it; 96 points with sizes drawn uniformly over 2–14 gives coverage of both sparse and
dense communities, which is where the interaction structure lives. Sampling sizes uniformly
rather than sampling members independently at p = 0.5 matters — independent sampling
concentrates almost everything near size 7 and never tests the extremes.

**Why 12 random points in every EI round.** They are the only honest held-out test of
whether the surrogate generalises. Without them a BO run cannot distinguish a good model
from a lucky one, and the project's stated first question — *can a model predict optimal
community composition* — becomes unanswerable.

**Objective:** L-lactate in g/L and in mol/mol glycerol consumed. **Optical purity is a
hard constraint in the acquisition, not a term in the objective** — a well below the 98%
spec is a failure whatever its titre, and averaging it into a score lets EI trade purity
for yield.

---

## The reference arms, outside the loop

Run alongside, never in the acquisition:

- **Yeast extract**, on round 1's best few and on the round 4 confirmations. This is the
  ceiling. The in-silico value of relieving biosynthetic burden is +18.8%; the observed
  value on at least one strain here is ~6-fold. Knowing the gap between the defined-medium
  optimum and the YE ceiling tells you how much is left on the table — without letting EI
  chase it.
- **B12**, as a plain supplement on a few round-1 communities. Eleven of the panel carry
  B12-dependent enzymes and, with `P3_C11` dropped, no member makes it by a route that
  runs under O₂ limitation — the aerobic `cob*` route left with the *Gluconobacter*, and
  the anaerobic `cbi*` route sits in `P1_B3`/`P1_B8`, which are held. So B12 may simply be
  a required medium component rather than a cross-feed. Cheap to test, and it is a defined
  supplement rather than a rich one.

---

## The substrate is crude glycerol. There is no pure-glycerol arm.

Pure glycerol → lactate is a solved problem: ~85% of theoretical yield in an engineered
monoculture (Mazumdar 2010). Running it here would reproduce a published result with worse
strains. **Crude glycerol is the reason the project exists**, so it is the base medium from
round 1, not a confirmation block.

This also settles the question the rule raises about feedstock, which was the wrong
question to ask: optimising on pure and confirming on crude has exactly the yeast-extract
structure — EI returns a community optimal in a condition nobody intends to use, and it
would look best precisely because nothing asked it to tolerate anything.

**The one thing to lock down before round 1: a single batch of crude glycerol large enough
for the whole campaign.** If the feedstock changes between rounds the surrogate is fitting
a moving target, and that failure is invisible in the data — it looks like a noisy response
surface rather than a shifted one. Characterise the batch (methanol, Na/K, free fatty
acids, soap, pH, water) and record it with the results.

### What the pool already tells you about surviving it

`detox_screen.py` scored all 384 genomes on the crude-glycerol stressors, and the ranking
of the existing pool is in `POOL_V5.md`. Three things carry into the design:

- **`P3_G8` (*Rhodococcus erythropolis*) scores 86 against 34 for the next best** — 49 LCFA
  genes, 23 halotolerance. Its lactate oxidases create an oxygen-linked risk, and any
  quinone-linked route must be retested with nitrate or another acceptor. Its crude-feed
  robustness and product-consumption risk must therefore be measured together.
- **`P3_C3` (*L. lactis*) scores lowest at 10.** Best LDH complement in the library,
  halotolerance 3, against a feedstock carrying 50–200 g/L Na/K. Expect it to fail in the
  crude arm and do not read that as the community failing.
- **No genome in the library oxidises methanol** — none of 384. Methanol has to be handled
  upstream, which is standard practice since biodiesel plants recover it, but the batch
  spec needs checking rather than assuming.

## Practical risks

- **Strain-specific media.** Many of the 14 producers need at least one predicted supplement and `P3_C3`
  (*L. lactis*) needs eight. A single shared minimal medium silently selects against the
  auxotrophs and reproduces the models' own bias. Either supplement to the union of
  requirements — accepting that this is a defined, cheap addition and not yeast extract —
  or record that the minimal arm is testing a smaller pool than it appears to.
- **Inoculum ratio.** The feature vector is presence/absence, so **total** inoculum is held
  constant and split equally among the members present — not a fixed dose per member. A fixed
  per-member dose makes total density scale with community size, and the surrogate then learns
  "more members is better" as an artifact of inoculum density rather than of community
  structure. Record the actual delivered volumes from the Echo transfer log; a community that
  fails because one member was under-inoculated is a measurement error the model will learn as
  biology. Note also that equal OD is not equal cell number across this panel — a *Rhodococcus*
  and a *Lactococcus* differ severalfold — which is a documented bias, not a fixable one.

- **Two timepoints, superseded.** `EXPERIMENT_PLAN.md` §4.4 calls for sampling at 24 h and
  48 h so the model learns yield rather than growth rate. That still holds as an *objective*,
  but the O₂ regime makes destructive double-sampling impossible: breaking the seal at 24 h
  destroys the one condition being held constant, and duplicate plates double the campaign.
  **Kinetic OD600 read through the seal replaces the second timepoint** — harvest once, after
  every well has been flat for ≥12 h, with the trace as the evidence that growth had finished.
  See the protocol, Part H2.
- **pH.** Lactate accumulation crashes unbuffered wells, and a crashed well reads as a
  failed community. Buffer, and measure final pH in every well — it is the cheapest and
  strongest single predictor of a dead well.
- **Composition at the end.** Tier 4 of `EXPERIMENT_PLAN.md` — 16S or strain-specific qPCR
  on the top and bottom performers. Not needed to fit the surrogate, whose features are the
  inoculum. Needed to explain why the winners won, and required if a gLV is fitted later.
