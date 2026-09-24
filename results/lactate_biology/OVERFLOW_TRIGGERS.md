# What actually triggers the pyruvate spill — hypotheses and evidence

`WHY_LACTATE.md` established *why* lactate exists (NAD⁺ regeneration) and that plain FBA
cannot predict it. The open question it left is the useful one: **glycolytic flux has to
exceed PDH and TCA capacity for the spill to happen — so what conditions and what partners
create that excess?** Knocking out PDH or the TCA cycle is not available to us (lethal, and
we are not engineering strains), so everything here is about *conditions*.

This is a literature pass cross-referenced against our own genomes and models. Six
hypotheses, ranked by leverage. Three of them are directly actionable at the bench; two
are levers the models are structurally blind to; one is a liability specific to glycerol
that nothing in the project had flagged.

It also corrects three earlier results — see **Corrections** at the end. One of them
reverses a conclusion.

---

## H1 — Oxygen supply at the stoichiometric ratio. **Strongest, quantitative, testable.**

**Claim.** Lactate output is governed by O₂ *rate* relative to glycerol uptake, with a
sharp optimum at **0.5 mol O₂ per mol glycerol**. Below it, carbon leaks to ethanol and
lactate capacity collapses. Above it, respiration takes the carbon to CO₂.

**Our evidence.** Sweeping the O₂ cap as a ratio to glycerol (20 mmol glycerol, growth held
at 90% of max, lactate capacity in mmol/gDW/h):

| O₂ / glycerol | `P1_B3` | `P1_B11` | `P3_F11` | µ (`P1_B11`) |
|---:|---:|---:|---:|---:|
| 1.00 | 17.00 | 16.71 | 16.68 | 0.257 |
| 0.60 | 17.00 | 16.71 | 16.68 | 0.257 |
| **0.50** | **17.04** | **16.74** | 16.17 | 0.255 |
| 0.40 | 16.76 | 15.87 | 13.26 | 0.230 |
| 0.25 | 15.20 | 9.89 | 7.64 | 0.191 |
| 0.10 | 13.54 | 3.91 | 1.85 | 0.153 |

The plateau ends exactly at the stoichiometric ratio. There is nothing to gain above 0.6
and a steep, near-linear loss below 0.5. `P1_B3` is markedly more robust to under-aeration
than the other two (13.54 versus 3.91 and 1.85 at ratio 0.1) — a robustness ranking that
did not exist before, and one that matters in a real vessel with dissolved-oxygen
gradients.

**Why 0.5 exactly, mechanistically.** All three producers route glycerol *only* through
`GLYK` → `G3PD4` (menaquinone-linked) in the model solution — `GLYCDx`, the NAD-linked
route, carries zero flux. That means glycerol's surplus electrons land on the **quinone
pool**, not on NADH:

```
glycerol + ATP → G3P                      (GlpK)
G3P + MK      → DHAP + MKH2               (G3PD4 / GlpD-type)
DHAP → pyruvate                           +1 NADH, +2 ATP
pyruvate + NADH → lactate                 −1 NADH
MKH2 + ½ O2 → MK + H2O                    ← the only O2 demand
```

The NAD pool balances itself through LDH; the *only* thing oxygen is needed for is
reoxidising the quinone. That is 0.5 O₂ per glycerol, and the solver consumes exactly
5.000 for 10 glycerol. Supply exactly that much and there is no oxygen left for NADH
dehydrogenase — so glycolytic NADH has nowhere to go except LDH. Supply more and
respiration outcompetes.

**Literature.** This is precisely the trick used in the best published glycerol→lactate
process. [Mazumdar et al., AEM
2010](https://journals.asm.org/doi/10.1128/aem.00664-10) reached **85% of theoretical
yield** (0.83 g/g, >99.9% chiral purity) by running *E. coli* **microaerobically** and
overexpressing `glpK`–`glpD` — the respiratory route — which alone gave a **~5-fold**
improvement over the fermentative GldA route, "because GlpD donates electrons directly to
quinones." **Our Bacillus strains carry that route natively** — `glpF`, `glpK` and `glpD`
in all nine genomes, zero anaerobic `glpB/glpC` — so the single most important engineering
step in that paper is something we get for free.

Independent confirmation that the optimum is real and two-sided:

- *B. coagulans* chemostat: aeration **halved** lactic acid (19.5 → 11 g/L) while raising
  biomass 3.9 → 4.6 g/L ([Enzyme Microb
  Technol](https://www.sciencedirect.com/science/article/abs/pii/S0141022998000982)).
- Micro-aerated lactate fermentation: moderate aeration gave **+35.5% lactate over strict
  anaerobiosis**, and **excessive aeration decreased it** ([Water
  Research](https://www.sciencedirect.com/science/article/abs/pii/S030147972600633X);
  [related](https://pubmed.ncbi.nlm.nih.gov/41412395/)).

**What to do.** Treat OTR as the primary controlled variable, targeted at
**0.5 × the measured glycerol uptake rate**. Use respiratory quotient as the online
signal — it is the standard readout for exactly this regime, since DO alone carries little
information once oxygen is limiting ([Bioeng
Transl](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7447875/)).

---

## H2 — Rex derepression through the NADH/NAD⁺ ratio. **Same knob, second mechanism.**

**Claim.** `ldh` is not constitutive. It is transcriptionally gated by Rex, a repressor
whose DNA binding is released as NADH/NAD⁺ rises. O₂ limitation raises that ratio, so H1
and H2 are the same physical lever acting through two independent mechanisms — which is
why the effect should be larger than either predicts alone.

**Our evidence.** `rex` (K03710) is present in **all nine** genomes (1–4 copies).
`cydAB` — the high-affinity cytochrome *bd* oxidase, co-regulated with `ldh` in the same
Rex regulon — is present in all nine as well. So the regulatory machinery is intact
throughout the pool.

**Literature.**

- *B. subtilis* Rex directly regulates `cydABCD` and the NADH-linked fermentative
  `lctP-ldh` operon ([Ravcheev et al., J
  Bacteriol](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3294762/)). Rex binding depends
  on the NADH/NAD⁺ ratio ([McLaughlin et al., Mol
  Cell](https://www.cell.com/molecular-cell/fulltext/S1097-2765(10)00368-0)).
- `lctE` induction under fermentative conditions is `fnr`- and `resDE`-dependent, and
  `lctE` mutants grow poorly anaerobically — NADH reoxidation by LDH is load-bearing
  ([Nakano & Zuber, J Bacteriol
  182:3072](https://jb.asm.org/content/182/11/3072); [cytochrome *bd* / LDH
  co-expression](https://www.microbiologyresearch.org/content/journal/micro/10.1099/mic.0.28124-0)).
- `rex` deletion in *B. cereus* changes carbon flux through the NADH-recycling lactate
  pathway, in opposite directions aerobically and anaerobically ([Laouami et
  al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4162614/)).
- The dependence is close to linear and quantified in *L. lactis*: NADH/NAD⁺ from 0.711 to
  0.383 took lactate from 21.2 to 9.9 mM while diacetyl rose 1.07 → 4.16 mM ([PLOS
  ONE](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0036296)).

**What to do.** Measure intracellular NADH/NAD⁺ alongside the OTR series. It converts H1
from "we set an oxygen level" into a mechanistic dose-response, and it is the variable that
would transfer to any other strain. One favourable note from the earlier screen: **no
H₂O-forming NADH oxidase (`nox`, K00359/K17869) in any of the nine genomes** — nothing
competes with LDH to dump NADH without making product.

---

## H3 — FBP activation of LDH: the glycerol-specific liability. **Nothing had flagged this.**

**Claim.** Firmicute L-LDH is allosterically activated by fructose-1,6-bisphosphate. FBP
tracks glycolytic flux — but **glycerol enters glycolysis at DHAP, *below* the aldolase
step**, so on glycerol FBP is a gluconeogenic trickle rather than a glycolytic pool. The
allosteric switch that turns LDH on may simply never be thrown. This is a mechanism by
which glycerol is a structurally *worse* lactate substrate than sugar, independent of
redox.

**Our evidence.** In the maximum-lactate flux solution, aldolase (`FBA`) carries only
**0.36–1.12** against a glycerol uptake of 10 and a pyruvate kinase flux of **8.7–9.4**.
So 4–12% of carbon touches the FBP node while ~90% goes straight from DHAP to pyruvate.
Whatever FBP pool exists is small and gluconeogenic. Our three producers are all
Firmicutes, the clade where FBP-dependent LDH is the norm.

**Literature.**

- Bacterial NAD-dependent L-LDH is "activated specifically by low concentrations of
  fructose-1,6-diphosphate"; activation **lowers the apparent Km for both pyruvate and
  NADH**, and the enzyme is inhibited by Pi, ATP and high NADH ([Garvie, Microbiol
  Rev](https://pubmed.ncbi.nlm.nih.gov/6997721/); [LDH from *B.
  subtilis*](https://pubmed.ncbi.nlm.nih.gov/236452/)).
- The structural basis is worked out in *B. stearothermophilus* — the FBP-bound R state
  versus the T state ([Iwata et al., J Mol
  Biol](https://www.sciencedirect.com/science/article/abs/pii/S0022283683711228)).
- FBP sensitivity varies by species: *L. lactis* LDH is strongly FBP-activated,
  *L. plantarum* LDH is **not regulated by FBP at all** ([JBC
  review](https://www.jbc.org/article/S0021-9258(20)45563-2/fulltext)). So an
  FBP-independent LDH is a real, selectable property.
- FBP is the canonical glycolytic-flux signalling metabolite ([Nat Commun
  2017](https://www.nature.com/articles/s41467-017-01019-z)), and it drops >5-fold when
  carbon is diverted away from upper glycolysis ([Nat Commun
  2021](https://www.nature.com/articles/s41467-021-25142-0)). One source states the link
  directly: high FBP activates **both PYK and LDH** and directs flux toward lactate.
- **pH partially substitutes for FBP.** FBP is absolutely required around pH 8, but
  measurable activity persists without it at **pH 6.0–7.2** — and *L. casei* LDH, which
  strictly requires FBP at neutral pH, becomes active under acidic conditions through
  homotropic activation by pyruvate itself ([Protein Eng Des
  Sel](https://academic.oup.com/peds/article/15/1/35/1528652)).

**What to do — three independent levers, cheapest first.**

1. **Run at pH 6.0–6.5**, not 7. Free, and it is the condition under which the enzyme works
   without its activator. (But see the conflict flagged in H4.)
2. **Co-feed a small amount of sugar** — glucose at a few percent of the glycerol carbon,
   enough to fill the FBP pool without meaningfully diluting the glycerol feedstock or the
   carbon yield. This is a hypothesis, not an established result: co-feeding is a
   documented way to control the ATP/NADH/acetyl-CoA balance ([Sci
   Rep](https://www.nature.com/articles/s41598-020-76761-4)), but I found nothing testing
   sugar specifically as an FBP trigger for lactate on glycerol. If it works it is cheap
   and immediately transferable; it is the highest-information experiment in this document.

   > **Strengthened by the pool v3 literature pass.** Pyruvate oxidase — the enzyme that
   > burns lactate to acetate aerobically in *L. plantarum* — is **repressed by glucose**
   > ([AEM 72](https://journals.asm.org/doi/10.1128/aem.00659-06)). So the same co-feed
   > raises FBP *and* suppresses the largest competing drain in `P2_A5`, by two independent
   > mechanisms. See `LITERATURE_POOL_V3.md`.
3. **Treat FBP-independence as a strain-selection criterion** if the pool is ever
   re-screened.

**Nothing in our models can see any of this.** FBA has no allostery. The capacity numbers
are therefore upper bounds that implicitly assume LDH is fully active.

---

## H4 — Acetoin and 2,3-butanediol are the competing sink, and it is the wrong one to have

**Claim.** In *Bacillus*, the dominant aerobic overflow route is not lactate — it is
acetate, acetoin and 2,3-butanediol. Acetolactate synthase consumes **two** pyruvate and
regenerates **no** NAD⁺, so it competes for pyruvate without doing the redox job lactate
does. If the spill happens but goes to acetoin, we get nothing.

**Our evidence.** `alsS` **and** `alsD` are present in all three Bacillus producers
(`P1_B3`, `P1_B11`, `P3_F11`), and every model in the pool carries `ACLS` plus an
`EX_actn__R_e` exchange. **This was missing from the earlier capacity analysis** — acetoin
and 2,3-BD were never in the product list, so the competition was never scored. It is in
now. In silico the competition looks mild (acetoin capacity ≤2.0 against lactate ~17), but
that is not reassuring, because the models under-represent the route: the literature says
otherwise, loudly.

**Literature.**

- Aerobic growth of *B. subtilis* on excess glucose produces **acetate, acetoin and
  2,3-butanediol** ([Microb Cell
  Fact](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3728045/)); the switch is
  CcpA-dependent, activated at high glucose ([J
  Bacteriol](https://journals.asm.org/doi/10.1128/jb.181.22.6996-7004.1999)).
- ALS is "active at a high pyruvate concentration and a low pH", and `alsSD` is induced by
  acetate accumulation, low pH and aerobic stationary phase via AlsR ([Frädrich et al., J
  Bacteriol](https://journals.asm.org/doi/10.1128/jb.06425-11); [factors in acetoin
  synthesis](https://pmc.ncbi.nlm.nih.gov/articles/PMC3964831/)). Acetoin is pH-neutral,
  and its role is thought to be *avoiding* acidification — i.e. it is actively favoured in
  the pH range H3 wants.
- **The decisive number**: inactivating `alsS` in *B. subtilis* raised L-lactate yield by
  **15%**, and that strain reached **95% of theoretical yield** and 105 g/L ([Romero et
  al., Microb Cell Fact 8:23](https://pmc.ncbi.nlm.nih.gov/articles/PMC2680810/)).

**A genuine conflict to resolve at the bench.** H3 wants lower pH (LDH works without FBP);
H4 says low pH plus accumulating acetate is exactly what induces `alsSD`. These pull
opposite ways. Best guess: pH 6.0–6.5 with tight base control so acetate does not
accumulate, and acetoin measured every time so we know which way it went. This is a real
tension, not a detail — it is the reason to measure acetoin in the very first experiment.

**Model caveat, self-caught.** Blocking `ACLS` in silico is *lethal* for `P1_B11` and
`P3_F11`. That is an artifact: CarveMe merged the catabolic AlsS with the anabolic IlvBH
acetolactate synthase, so deleting it also deletes branched-chain amino acid biosynthesis.
Real *B. subtilis* `alsS` mutants are viable and are the basis of the 95%-yield result
above. **Do not read the in-silico lethality as a prediction.**

---

## H5 — Precursor supply, reframed with a price tag

**Claim.** Supplying amino acids and nucleotides lets the producer stop diverting glycerol
carbon into building blocks, so more reaches pyruvate. This is the supporter hypothesis
from `WHY_LACTATE.md`; the literature both supports it and prices it.

**Literature.**

- The 95%-of-theoretical *B. subtilis* result required **mineral medium plus 8 g/L yeast
  extract**. Base mineral medium alone worked but gave less; full rich medium was *not*
  needed ([Romero et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC2680810/)). So the
  effect is real and saturates early.
- Individual amino acids matter: cysteine supplementation gave 94.1 g/L at 0.90 g/g yield
  in *Sporolactobacillus* ([Sci Rep](https://www.nature.com/articles/s41598-024-61247-4)).
- **Yeast extract accounts for up to 38% of total lactic-acid production cost**, which is
  why the field keeps trying to replace it ([Fermentation
  6:97](https://www.mdpi.com/2311-5637/6/4/97)).

**Why this is the strongest economic argument for a community in the whole project.** The
value of a partner is not rescuing an auxotrophy — every member is prototrophic — it is
**delivering the yeast-extract effect without buying yeast extract.** That is a cost line
worth up to ~38%, and it is the first thing in this project that gives the community form a
concrete commercial reason to exist rather than a scientific one.

`P2_D12` (*Rahnella inusitata*) secretes ten amino acids **including cysteine**, the one
with the best published effect. Its drawback stands unchanged — it consumes lactate — so
the experiment is a three-way: producer alone, producer + `P2_D12`, producer + 8 g/L yeast
extract. The third arm is the benchmark the partner has to beat.

---

## H6 — Export thermodynamics, and why titre and pH interact

**Claim.** Lactate efflux is electrogenic H⁺-symport that *generates* proton-motive force
at low external lactate — but the driving force collapses as titre rises and pH falls, at
which point export becomes an energy cost.

**Literature.** In *S. cremoris*, the H⁺/lactate stoichiometry falls from 1.4 at 1 mM to
**0.7 at 50 mM** at pH 5.5, and the lactate gradient collapses from ~100 mV to ~40 mV
([Otto et al.](https://journals.asm.org/doi/10.1128/jb.149.2.733-738.1982); [energy
recycling by lactate
efflux](https://journals.asm.org/doi/10.1128/jb.162.1.383-390.1985)). At high product
concentration and low external pH the energy cost of export can equal or exceed the energy
yield of making the product.

**Why it matters here.** Lactate export has now been the binding constraint in **six**
separate places in this project (`P3_A10`, `P3_B12`, `P1_A3`, `P2_D12`, the ethanol
recycle, and now the pH/titre window). It reinforces the standing recommendation to add
`L_LACt3` to the models that make lactate but cannot export it, and it adds a process
constraint: **pH-controlled fed-batch with base addition, or in-situ product removal** —
not free-acid accumulation.

---

## Corrections to earlier results

**1. Nitrate in the mineral base invalidated every oxygen conclusion.** `MINERAL_BASE`
opened `EX_no3_e` at −1000 as a nitrogen source. Nitrate is also a respiratory electron
acceptor, *Bacillus* respires it, and every model in the pool carries a nitrite efflux
exchange — the giveaway. With nitrate open, growth was flat against oxygen all the way down
to O₂/glycerol = 0.02, because the cell simply switched acceptor. Ammonium already supplies
nitrogen, so nothing is lost by closing it.

**This reverses a conclusion.** `CAPACITY_RANKING.md` said *"Oxygen is a growth constraint
here, not a routing lever — aeration should be set for growth, not used to push carbon
toward lactate; that lever does not exist in these models."* That is wrong. With nitrate
closed, oxygen is the cleanest lever in the analysis, with a sharp knee at the
stoichiometric ratio (H1).

Fixed in `physiology/lactate_capacity.py`, `modeling/audit_models.py` and the now-renamed
`modeling/medium_diagnostics.py`; nitrate and the
other alternative acceptors are now closed explicitly. `capacity_table.tsv` is regenerated
(567 rows, oxygen now swept as a ratio to glycerol). The **strain ranking and carbon yields
are unchanged** — `P1_B3` 0.850, `P1_B11` 0.835, `P3_F11` 0.834 — so the lab-ready
recommendation stands. Only the oxygen conclusion changes. The auxotrophy and
minimal-medium results are aerobic and unaffected in substance, though those two scripts
have not been re-run.

**2. The two *Gluconobacter* strains are not lactate producers.** `P3_C11` (0.833) and
`P3_C9` (0.806) appeared in the producer table. Neither genome has **any** fermentative
L-LDH — zero K00016 hits. CarveMe built their `LDH_L` from `lldD` (EC 1.1.2.3, K00101),
the *respiratory, quinone-coupled* L-lactate **oxidase** that runs lactate → pyruvate, and
their `LDH_D` reactions from `dld` (EC 1.1.2.4). Those are lactate-consuming enzymes
credited as production capacity. Their lactate numbers should be struck; they are DHA
donors and potential lactate **sinks**, which is also the more useful thing to know about
them. The three Bacillus producers are unaffected — each has a genuine `ldh` (e.g.
`P1_B3`/`BBFJKGGH_04870`, EC 1.1.1.27), and `P1_A3`'s two K00016 hits are both `mdh`,
confirming the earlier "no L-LDH" call for the sink control.

**3. Acetoin and 2,3-butanediol were missing from the product list.** Added to
`physiology/lactate_capacity.py` (H4).

**4. A framing refinement in `WHY_LACTATE.md`.** That document said glycerol → lactate
leaves "+1 NADH surplus per glycerol". The electron count is right, the carrier is not: via
the `GlpK`/`G3PD4` route the models actually use, the surplus sits on the **quinone pool**
and the NAD pool balances through LDH. The 0.5 mol O₂ per mol glycerol figure is unchanged
either way, but the distinction is the whole mechanism behind H1 — it is why oxygen
limitation forces NADH to LDH rather than merely starving the cell.

---

## A risk flag worth stating plainly

Searching the literature for lactate fermentation by *Priestia aryabhattai*, *Bacillus
velezensis* or *Bacillus altitudinis* returns **nothing**. All three are documented as
biocontrol and plant-growth-promoting organisms, not fermentation organisms. Our capacity
ranking is built entirely on genome-derived models with no published physiology behind it
for these species — unlike *B. subtilis*, *B. coagulans* and *L. lactis*, which have
decades of process literature.

That is not a reason to drop them; it is the reason the first experiment is a measurement,
not an optimisation. But it argues for **carrying one strain with a track record as a
positive control** — a *B. coagulans* or a `lctE`-intact *B. subtilis* — so that if all
three pool producers make no lactate, we can tell "wrong strains" from "wrong conditions".

---

## Ranked experiments

1. **OTR series at O₂/glycerol = 0.1, 0.25, 0.5, 0.6, 1.0** on `P1_B3`, `P1_B11`, `P3_F11`.
   Measure lactate, **ethanol**, **acetoin/2,3-BD**, acetate, glycerol consumed, biomass.
   This tests H1 directly, tests H4 in the same run, and produces the measured product
   ratios that turn every model number in this project from capacity into prediction.
2. **pH 6.0 versus 7.0** crossed with the two best O₂ points. Tests H3 and H4 against each
   other where they conflict. Cheap.
3. **Glucose co-feed at ~5% of glycerol carbon**, at the best O₂/pH point. Tests H3's FBP
   hypothesis. Highest information per unit effort if it works.
4. **Three-way nutrition arm**: minimal / + `P2_D12` / + 8 g/L yeast extract. Tests H5 and
   prices the partner against the thing it has to replace.
5. **NADH/NAD⁺ measurement** on the OTR series, if the assay is available. Converts H1+H2
   into a mechanism that transfers to other strains.
