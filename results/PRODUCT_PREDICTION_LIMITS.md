# These GEMs cannot predict lactate titre — and what to measure instead

Found while validating the new COMETS product labels end to end. The labels work; the
quantity they were built to report does not exist in these models.

## What the smoke run showed

Two species (`P1_B3` + `P1_A3`) on glycerol minimal, 200 cycles, aerobic:

| | Value |
|---|---:|
| glycerol consumed | 5.91 mmol |
| lactate produced | **0.008 mmol** |
| carbon into tracked products | **0.13%** |

Essentially all carbon went to biomass and CO₂.

## It is not the objective style, the oxygen, or the medium

Each candidate explanation was tested and ruled out:

| Hypothesis | Test | Result |
|---|---|---|
| `MAX_OBJECTIVE_MIN_TOTAL` suppresses secretion | reran with `MAXIMIZE_OBJECTIVE_FLUX` | lactate **0.0000** — worse, not better |
| Oxygen supply too high | swept O₂ = 100, 10, 2, 0.5, 0 mmol | lactate 0.0025–0.0083 at every level; at O₂=0 no growth at all |
| Models cannot make lactate | growth-maximising FBA in cobra | secretes **8.23 mmol** (0.82 mol/mol) |

That last row is the clue: cobra and COMETS solve the *same* model and disagree
completely, which means the answer is not determined by the model.

## The actual cause: lactate is degenerate

At maximum growth, holding growth at its optimum and then minimising and maximising
lactate:

| Genome | µ max | lactate MIN | lactate MAX |
|---|---:|---:|---:|
| `P1_B3` | 0.1170 | **0.0000** | 8.229 |
| `P1_B11` | 0.1287 | **0.0000** | 8.059 |
| `P3_F11` | 0.1201 | **0.0000** | 8.103 |
| `P3_C11` | 0.1274 | **0.0000** | 7.796 |

**Lactate secretion carries no growth cost or benefit.** Anywhere from zero to ~8
mmol is equally optimal, so which value a run reports is an arbitrary property of the
solver's vertex choice, not a prediction. cobra picks a vertex with lactate; COMETS
picks one without.

### Redox limitation forces *a* fermentation product, but not lactate

Capping the O₂ uptake *rate* (not supply) does force overflow — but the cell chooses
ethanol:

| O₂ cap | µ max | lactate MIN | ethanol MIN |
|---:|---:|---:|---:|
| 20 | 0.1170 | 0.0000 | 0.0000 |
| 1.0 | 0.0540 | 0.0000 | 0.0000 |
| 0.5 | 0.0463 | 0.0000 | **2.852** |
| 0.25 | 0.0424 | 0.0000 | **4.698** |

And blocking the competing sinks does not rescue it: with ethanol, acetate and
formate all closed, lactate minimum stays 0 in all three producers, and for `P1_B11`
and `P3_F11` the lactate *maximum* collapses to ~0 as the models route carbon
elsewhere.

**Conclusion: lactate is never obligate in any of these models, under any oxygen
regime, with or without competing sinks blocked.**

## Correction

I previously wrote, from the production envelope, that "a growth-maximising dFBA will
therefore still make lactate, which is what COMETS needs to be informative." That was
wrong. The envelope gives the *upper* bound of the optimal face; it says nothing about
which point on that face a simulation lands on, and the lower bound is zero
throughout. The envelope was the right calculation for capacity and the wrong one for
prediction.

This is the standard FBA limitation for products that do not affect growth, not a
defect specific to these reconstructions.

## What the landscape can and cannot deliver

**Determinate, worth measuring:**

- growth, biomass, coexistence and survival — the existing 25 labels
- glycerol consumed
- **ethanol** under O₂-rate limitation, which is genuinely forced
- **lactate *capacity***: maximum lactate at the growth rate the simulation actually
  achieved. Deterministic, non-degenerate, and interpretable as "the best this
  community could do in this medium".

**Not determinate — do not report as a prediction:**

- lactate titre, yield or selectivity from the dFBA trajectory

## Recommended change

Keep the media-log labels — they are correct and needed for glycerol consumption and
the forced products. Add a per-simulation **capacity** metric computed post hoc:
constrain each member to its realised growth rate, maximise lactate, and record the
result alongside the titre. Report titre and capacity as separate columns so the
degeneracy is visible rather than hidden.

Two options for making titre itself meaningful, both real modelling commitments:

1. **Constrain secretion experimentally.** Measure single-strain product ratios on
   glycerol and impose them as bounds. This is the reality-gate step the brief already
   calls for, and it is the only route that makes the titre a measurement rather than
   an assumption.
2. **Break the degeneracy in silico** — flux sampling over the optimal face to get a
   distribution instead of a point, or minimisation of metabolic adjustment against a
   reference. Both report uncertainty honestly; neither invents information.

Option 1 is the one the project's own gates already require.
