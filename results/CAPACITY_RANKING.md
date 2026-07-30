# Lactate capacity ranking — what to take to the bench

From `lactate_capacity.py`; full table in `results/capacity_table.tsv` (756 conditions:
12 pool members × 3 glycerol levels × 7 oxygen ratios × 3 supplement sets).

**Capacity, not titre.** dFBA cannot predict lactate for these models — at maximum growth
the minimum lactate is zero and the maximum ~8, so any titre a simulation reports is
solver-arbitrary (`PRODUCT_PREDICTION_LIMITS.md`). Capacity — the maximum lactate reachable
while still growing at 90% of maximum — has a unique answer and a clear meaning: **the best
this strain could do on this medium.** It is an upper bound; where a real culture sits
inside it is what the bench decides.

**Oxygen is swept as a ratio to glycerol**, not as an absolute rate, because the ratio is
what the stoichiometry fixes: glycerol → lactate leaves two surplus electrons per glycerol,
so **0.5 mol O₂ per mol glycerol** is the exact respiratory demand (`OVERFLOW_TRIGGERS.md`).

## The producers

Four of twelve have genuine lactate capacity. Two further strains report capacity that is
an artifact of a mis-assigned LDH and are excluded from the ranking — see below.

On glycerol minimal medium, 20 mmol glycerol, O₂/glycerol = 1.0:

| Strain | Species | Capacity | Carbon yield | Ethanol | Acetoin | 2,3-BD |
|---|---|---:|---:|---:|---:|---:|
| `P1_B3` | *Priestia aryabhattai* | 17.00 | **0.850** | 17.32 | 2.00 | 8.54 |
| `P1_B11` | *Bacillus velezensis* | 16.71 | 0.835 | 17.17 | 2.00 | 8.46 |
| `P3_F11` | *Bacillus altitudinis* | 16.66 | 0.833 | **4.98** | 2.00 | **0.00** |
| `P2_A5` | *Lactiplantibacillus plantarum* | 8.40 ¹ | 0.363 | **0.00** | 1.43 | **0.00** |

¹ at its own optimum of O₂/glycerol = 0.6; 7.27 at 1.0.

**Excluded from the ranking**: `P3_C11` (16.66) and `P3_G8` (17.03) both report high
capacity, but neither has a genuine fermentative L-LDH. CarveMe built `P3_C11`'s `LDH_L`
from `lldD` (EC 1.1.2.3), the respiratory L-lactate **oxidase**, and `P3_G8`'s from a gene
carrying both LDH and MDH EC numbers without the "LDH family" descriptor. Reporting them
would credit lactate *consumption* as production.

## The trade-off is now three-way, and acetoin/2,3-BD is what changed it

The earlier ranking scored only lactate, ethanol, formate, acetate and succinate. Adding
**acetoin and 2,3-butanediol** — the textbook *Bacillus* overflow sinks, absent from the
product list until `OVERFLOW_TRIGGERS.md` flagged them — separates strains that looked
nearly identical:

- **`P1_B3` — highest yield and the most robust.** 0.850 mol/mol, and it holds 13.54
  capacity even at O₂/glycerol = 0.1, against 3.91 for `P1_B11` and 1.52 for `P3_F11`.
  That robustness matters in a real vessel with dissolved-oxygen gradients. Its cost is
  the largest competing capacity in the pool: 17.3 ethanol **plus** 8.5 2,3-BD.
- **`P3_F11` — the selectivity pick.** Carbon yield within 0.02 of the leader, but ethanol
  capacity 4.98 against ~17, and **2,3-butanediol capacity of exactly zero**. Its total
  competing-product capacity is by far the lowest of the three *Bacillus*-type producers.
  The earlier recommendation flagged its low ethanol; the 2,3-BD result strengthens it
  considerably. Its weakness is fragility to under-aeration (1.52 at ratio 0.1).
- **`P2_A5` — less than half the yield, but structurally the cleanest.** 0.363 mol/mol,
  yet **zero ethanol and zero 2,3-BD capacity** — the only genuine producer with no
  reduced co-product route at all. It also carries the pool's only real fermentation
  pedigree, and it is the test of the FBP-independence hypothesis.

### `P2_A5` behaves oppositely to the *Bacillus* strains on oxygen

The three *Bacillus*-type producers plateau: capacity is flat from O₂/glycerol = 0.5
upward and falls only below it. `P2_A5` **peaks at 0.6 and is actively harmed by excess
oxygen** — 8.40 at ratio 0.6 against 4.00 at 2.0.

| O₂/glycerol | 2.0 | 1.0 | 0.6 | 0.5 | 0.4 | 0.25 | 0.1 |
|---|---:|---:|---:|---:|---:|---:|---:|
| `P1_B3` | 17.00 | 17.00 | 17.00 | 17.04 | 16.76 | 15.20 | 13.54 |
| `P1_B11` | 16.71 | 16.71 | 16.71 | 16.74 | 15.87 | 9.89 | 3.91 |
| `P3_F11` | 16.66 | 16.66 | 16.66 | 16.17 | 13.26 | 7.64 | 1.52 |
| `P2_A5` | 4.00 | 7.27 | **8.40** | 8.00 | 7.40 | 4.00 | 2.80 |

This is consistent with lactic acid bacterium physiology and it means the two groups want
**different aeration setpoints** — which is a real constraint on co-culturing them, and
one that would not have appeared without sweeping oxygen as a ratio.

## Oxygen is a routing lever, with a knee at the stoichiometric ratio

The previous version of this document said oxygen was "a growth constraint, not a routing
lever". That was an artifact of nitrate being open in the mineral base, letting the models
respire their way around every oxygen cap. With nitrate closed, the plateau ends exactly at
O₂/glycerol ≈ 0.5–0.6 and capacity falls steeply below it. Nothing is gained above 0.6.

Two further electron-donor artifacts of the same family — sulfide, and unbounded proton
uptake against several hundred charge-unbalanced reactions — were found and removed, and
the gate now measures electron recovery directly rather than relying on the closures
(`POOL_V3.md`).

## Supplements still do not improve carbon efficiency

Re-confirmed on the corrected medium, at O₂/glycerol = 1.0:

| Strain | Minimal | + vitamins | + vitamins & AA | Supplement carbon share |
|---|---:|---:|---:|---:|
| `P1_B3` | **0.850** | 0.850 | 0.730 | 21.8% |
| `P1_B11` | **0.835** | 0.650 | 0.697 | 28.5% |
| `P3_F11` | **0.833** | 0.833 | 0.797 | 14.5% |
| `P2_A5` | 0.363 | 0.363 | **0.364** | 26.2% |

Measured per carbon *actually consumed*, supplements make yield worse or leave it
unchanged, because the amino acids contribute 15–29% of the carbon taken up. `P1_B11`'s
vitamin row (0.650) sits below both its neighbours and is more likely an alternative
optimum than a real effect — treat it as noise.

**For a glycerol-valorisation process this remains the desirable answer: the cheapest
medium is also the most carbon-efficient.** But the scope limit from `OVERFLOW_TRIGGERS.md`
stands — this is capacity in FBA, and it says nothing about overflow, which is the
mechanism that actually makes lactate. The precursor-supply hypothesis (H5) is untestable
here by construction.

## Lab-ready recommendation

**Medium:** glycerol minimal — mineral salts plus glycerol. Best or equal-best carbon
yield for every producer, and the cheapest. Lipoate is the one deliberate variable
(`POOL_V3.md`): withholding it suppresses pyruvate dehydrogenase pool-wide but excludes
both lactic acid bacteria.

**Aeration:** target OTR ≈ **0.5 × the measured glycerol uptake rate** for the *Bacillus*
producers; **0.6** for `P2_A5`, which is harmed above it.

**Strains, by what you want to optimise:**

1. `P3_F11` — **selectivity.** Yield within 0.02 of the leader, lowest competing capacity
   by a wide margin. Needs aeration held near setpoint.
2. `P1_B3` — **yield and robustness.** Highest carbon yield, tolerates poor aeration best.
   Largest competing capacity, so the product spectrum matters most here.
3. `P2_A5` — **the informative one.** Half the yield, but zero reduced co-products, real
   process literature, and the FBP-independence test.

**What to measure** — unchanged, and it is what converts capacity into prediction:
glycerol consumed; L- versus D-lactate separately; **ethanol, acetoin and 2,3-butanediol**
(the models say up to 17 ethanol and 8.5 2,3-BD sit alongside the lactate, and which one a
strain picks is exactly what the model cannot tell us); acetate and formate; glyceric acid
for `P3_C11`; growth rate and final biomass.
