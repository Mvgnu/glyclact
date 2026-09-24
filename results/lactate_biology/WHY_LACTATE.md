# What is lactate actually for, and what does that mean for the design?

Prompted by the right question: lactate is not a growth requirement, so under what
conditions is it made at all — and does that change which supporters we want?

## Lactate is a redox valve, not a product the cell wants

Lactate dehydrogenase exists to regenerate NAD⁺. Glycolysis oxidises
glyceraldehyde-3-phosphate and produces NADH; if NAD⁺ is not returned, glycolysis
stops. The cell has two ways to return it — respiration, which needs an electron
acceptor, or fermentation, which dumps the electrons onto pyruvate. **Lactate is made
when NADH production outruns the capacity to reoxidise it any other way.** Nothing
about lactate is wanted for its own sake.

That makes it conditional by nature, exactly as suspected. The conditions are worth
being precise about, because for *glycerol* they are not the familiar ones.

## Glycerol is more reduced than sugar, and that changes the arithmetic

Degree of reduction per molecule: glycerol 14, lactate 12, ethanol 12, 1,3-propanediol
16, CO₂ 0.

Glucose → 2 lactate is perfectly redox balanced — the 2 NADH from glycolysis are
exactly consumed by 2 LDH steps. That is why homolactic fermentation of sugar works
without oxygen.

**Glycerol → lactate is not balanced.** Glycerol must first be *oxidised* to DHAP, then
glycolysis makes an NADH, and LDH consumes only one. Net **two surplus electrons per
glycerol**, which have to go somewhere.

> **Refinement (see `OVERFLOW_TRIGGERS.md`).** This originally read "+1 NADH surplus". The
> electron count is right but the carrier is not. All three producers route glycerol
> through `GLYK` → `G3PD4`, the menaquinone-linked G3P dehydrogenase, so the surplus lands
> on the **quinone pool** and the NAD pool balances itself through LDH. The 0.5 mol O₂ per
> mol glycerol below is unchanged either way — but the distinction is the whole mechanism:
> it is why supplying *exactly* that much oxygen leaves none for NADH dehydrogenase and
> forces glycolytic NADH into LDH.

Confirmed in the models — maximising lactate from 10 mmol glycerol:

| O₂ available | Lactate | O₂ consumed | Co-products |
|---:|---:|---:|---|
| 20 | **10.000** | **5.000** | water only |
| 5 | **10.000** | **5.000** | water only |
| 1 | 7.810 | 1.000 | methanol 6.57 |
| 0 | 3.333 | 0 | **ethanol 8.33**, CO₂ 3.33 |

Ten glycerol leave ten surplus NADH; ten NADH need five O₂ (each O₂ accepts four
electrons). The model consumes **exactly 5.000**.

Across the producers, O₂ per lactate:

| Genome | O₂ / lactate |
|---|---:|
| `P1_B11`, `P3_F11`, `P3_C11` | **0.500** (theoretical) |
| `P3_C9` | 0.464 |
| `P1_B3` | 0.446 |

## Consequence 1: aerobic operation is obligatory, not a preference

Glycerol → lactate at 1.0 mol/mol **requires ~0.5 mol O₂ per mol glycerol**. Without
oxygen the ceiling collapses from 1.0 to 0.33 mol/mol and ethanol takes over. The
decision to run aerobically was correct for a harder reason than convenience — the
stoichiometry does not permit anything else.

It also gives a concrete process target: **oxygen transfer rate ≈ 0.5 × glycerol
uptake rate.** Below it, carbon leaks to ethanol; there is no benefit above it.

## Consequence 2: the window, and why our models cannot find it

Too little O₂ and lactate is impossible. Plenty of O₂ and the cell can simply respire
pyruvate all the way to CO₂ for far more ATP — which is what our simulations do, and
why they report almost no lactate.

Real cells make lactate in between, through **overflow metabolism**: when glycolytic
flux exceeds the capacity of pyruvate dehydrogenase and the TCA cycle, pyruvate
accumulates and spills to lactate even with oxygen present. The driver is *enzyme
capacity*, not oxygen availability.

Plain FBA has no enzyme capacity constraints. It assumes any flux is achievable if
stoichiometry allows, so it sees no reason to overflow — which is precisely why
lactate came out degenerate (min 0, max ~8) in `PRODUCT_PREDICTION_LIMITS.md`. **The
degeneracy is not a quirk of these reconstructions; it is the absence of the mechanism
that causes lactate production in the first place.**

Predicting it properly needs an enzyme-constrained model (ecFBA / GECKO-style, or
proteome-allocation), or measured flux ratios imposed as bounds.

## Consequence 3: this changes what a supporter is for

The supporter question was framed as auxotrophy rescue — who supplies biotin, who
supplies methionine. On this analysis that is the less important axis.

If lactate is driven by flux overflow, then the useful partner is one that **raises
glycolytic flux relative to biosynthetic demand**. Supplying amino acids and
nucleotides means the producer diverts less glycerol carbon into building blocks, so
more passes through glycolysis to pyruvate, so more spills to lactate. This is the
established reason rich media raise lactate yields in lactic acid bacteria
fermentations.

Our capacity analysis concluded supplements do not improve carbon yield — and that
conclusion is sound *within FBA*, but it cannot see this effect at all, for the same
reason it cannot see overflow. The two results are not in conflict; they are answers
to different questions.

So the supporter hypothesis worth testing is not "who rescues an auxotrophy" but
**"who supplies biosynthetic precursors and thereby pushes carbon through glycolysis"**
— and the pool already contains a strong candidate in `P2_D12` (*Rahnella inusitata*),
which secretes alanine, aspartate, cysteine, glutamate, glycine, serine, threonine,
tryptophan, phenylalanine and methionine. Its drawback remains that it also consumes
lactate, so the test has to weigh precursor supply against product loss.

## What to do with this

**Followed up in `OVERFLOW_TRIGGERS.md`**, which cross-references this against the
published literature and adds five further mechanisms — Rex derepression, FBP activation
of LDH, acetoin competition, precursor supply and export thermodynamics — plus a
correction: contrary to the earlier capacity table, **oxygen *is* a routing lever**, with a
sharp knee at exactly the 0.5 ratio derived here. That table had nitrate open, which let
the models respire their way around every oxygen cap.

1. **Set aeration deliberately** to ~0.5 mol O₂ per mol glycerol consumed, and treat it
   as a controlled variable rather than "aerobic". This is the single most important
   process parameter that follows from the stoichiometry.
2. **Measure the lactate/ethanol split against O₂ transfer rate** in the first
   experiments. That is the axis the models say matters and cannot themselves resolve.
3. **Test precursor supply as a lactate lever**, not just as auxotrophy rescue —
   e.g. producer alone versus producer plus `P2_D12` versus producer plus casamino
   acids, measuring lactate per glycerol in each.
4. **Treat the supplement conclusion as scope-limited.** "Supplements do not improve
   carbon yield" holds for capacity in FBA; it says nothing about overflow, which is
   the mechanism that actually makes lactate.
