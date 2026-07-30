# Recycling ethanol carbon back to lactate — assessment

Idea: the producers carry large ethanol capacity alongside their lactate capacity
(~17 mmol/gDW/h in four of five). If that carbon is lost to ethanol, a partner that
eats ethanol could return it to lactate.

**Verdict: metabolically sound, blocked by the same lactate-export gap as everything
else, and better implemented as a resting-cell second stage than as a co-culture.**

## The route works, at the efficiency theory predicts

Ethanol to lactate must run acetyl-CoA → glyoxylate shunt → malate → pyruvate,
because pyruvate dehydrogenase is effectively irreversible. Stoichiometry:

```
2 ethanol (4 C) → 2 acetyl-CoA → succinate → malate → pyruvate (3 C) + CO2 (1 C)
```

so **75% carbon efficiency** on the recycle step. Measured in the models, ignoring
export: `P1_A3` and `P2_D12` both reach exactly **5.000 mmol cytosolic lactate from
10 mmol ethanol** — 15 C from 20 C, matching the prediction.

Thermodynamically it is favourable: ethanol is more reduced per carbon than lactate
(degree of reduction 6 vs 4), so the conversion releases electrons and yields ATP
aerobically. Nothing blocks it energetically.

## What blocks it: lactate export, again

| Strain | Grows on ethanol | Cytosolic lactate from ethanol | Export route |
|---|---:|---:|---|
| `P1_A3` | 0.098 | 5.000 | **none** — `EX_lac__L_e` orphaned from the cytosol |
| `P2_D12` | 0.352 | 5.000 | `L_LACtex` only; inner membrane is import-only |
| `P3_A10` | 0.045 | 1.818 | **none** — no exchange at all |

Adding the CarveMe universe exporter `L_LACt3`
(`h_e + lac__L_c --> h_c + lac__L_e`, already present in three pool members) completes
the route: `P1_A3` goes from 0.000 to 0.500 mmol lactate at 90% of max growth.

This is now the **fifth** place lactate export has been the binding constraint —
`P3_A10`, `P3_B12`, `P1_A3`, `P2_D12`, and now the ethanol recycle. Adding a lactate
exporter is the single highest-leverage curation available in this project, and unlike
the vitamin-export case it is well supported: `lctP` is annotated in several strains,
every real lactate producer exports lactate, and `L_LACt3` is a genuine universe
reaction rather than an invention.

## The catch: a growing digester makes biomass, not lactate

| Operating point | Lactate from 10 mmol ethanol | Carbon efficiency |
|---|---:|---:|
| 90% of maximum growth | 0.500 | 7.5% |
| zero growth (resting cells) | 5.000 | **75%** |

Ethanol is an energy-rich growth substrate, so a growing partner spends it on biomass
and CO₂. The recycle only pays at near-zero growth.

**Implication: this is a two-stage process, not a co-culture.** Grow the producer on
glycerol, then convert accumulated ethanol with resting cells — a standard industrial
pattern, and one that matches where the carbon efficiency actually is.

## There is no ethanol specialist in the library

Screening all 384 genomes for the digester profile — isocitrate lyase + malate
synthase + lactate transport + NAD-dependent L-LDH — returns **82 candidates**. Every
single one also carries `glpK`.

So there is no strain that eats ethanol but not glycerol. Any digester added to a
co-culture would compete for the primary substrate, and the three pool members that
grow on ethanol all grow *better* on glycerol (`P1_A3` 0.213 vs 0.098; `P2_D12` 0.640
vs 0.352; `P3_A10` 0.182 vs 0.045), so they would preferentially take glycerol. That
independently points away from the co-culture form of the idea.

## The most interesting lead: the producer could recycle its own ethanol

`P1_B3` — the top producer — has the complete gene set itself: isocitrate lyase,
malate synthase, malic enzyme, a lactate transporter and eleven aldehyde-dehydrogenase
family proteins. A single self-recycling strain would be more elegant than any
community.

Its model cannot grow on ethanol, and the cause is **not** transport: ethanol,
acetaldehyde and acetate are all reachable in the cytosol at full flux. The blocked
step is downstream and not yet identified — an earlier reading that acetyl-CoA was
unreachable was an artefact of the probe draining the CoA pool, so that conclusion does
not stand. `ACS`, `PPA`, `PTAr` and `ACKr` are all present with working pyrophosphate
sinks.

Worth diagnosing properly if the self-recycling route is attractive, since it would
remove the need for a partner entirely.

## Recommendation

1. **Add the lactate exporter** (`L_LACt3`) to the models that make lactate but cannot
   export it. It unblocks five separate findings and is the best-supported curation in
   the project.
2. **Treat ethanol recycling as a downstream resting-cell stage**, not a community
   member — that is where the 75% carbon efficiency lives.
3. **Measure ethanol first.** The whole premise rests on the producers actually making
   ethanol, which is a capacity, not a prediction. If the bench shows they make little
   ethanol on glycerol, the recycle is unnecessary; if they make a lot, it is worth the
   second stage. This is one more reason single-strain product ratios are the next
   experiment.
