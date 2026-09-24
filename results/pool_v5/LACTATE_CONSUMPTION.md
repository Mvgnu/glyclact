# Lactate-consumption diagnostic

## Scope

`lactate_consumption.tsv` contains **36 screening candidates**, not the final panel. It asks whether each CarveMe model can grow with L-lactate as sole carbon source after opening model lactate transport.

The historical table has a `coconsumes_*` column from an invalid objective formulation. Ignore that column. The repaired `physiology/lactate_consumption.py` no longer computes or emits it.

## What the two valid conditions mean

- **Aerobic:** oxygen exchange is available.
- **Oxygen-free:** oxygen is closed and no nitrate, fumarate, or other alternative respiratory acceptor is added.

The second condition is not a general prediction for “low-oxygen production medium.” It is specifically a no-terminal-acceptor condition.

## Results

Thirty of 36 models have aerobic L-lactate sole-carbon growth capacity. The six without that capacity are `P1_B11`, `P1_B8`, `P2_B6`, `P3_F11`, `P2_C2`, and `P1_B3`. Four are current producers and two are held for the separate racemase question.

Thirty-five of 36 lack L-lactate sole-carbon growth in the oxygen-free/no-alternative-acceptor condition. The exception, `P3_G9`, also has an unreliable model-level glycerol no-growth result, so it is not strong biological evidence.

The qualitative aerobic risk is worth following up, especially for `P3_G8` and `P3_C3`, which have sequence-level lactate-oxidation evidence. The quantitative growth ratios are uncalibrated capacities from M9-gapfilled models and are not selection gates.

## Correct interpretation

The result supports only this bounded statement:

> Removing oxygen and all alternative respiratory acceptors closes modeled respiratory L-lactate use in nearly all screened models.

It does **not** establish that nitrate or fumarate prevents product re-consumption. Quinone/cytochrome-linked L-lactate oxidation can feed anaerobic respiration when a compatible terminal acceptor is available. An acceptor may improve glycerol redox balance while also reopening lactate oxidation.

Therefore nitrate, fumarate, acetate, and oxygen transfer must be evaluated in the exact intended medium by measuring net L-lactate retention over time. This is a process-design measurement, not an automated strain exclusion rule.
