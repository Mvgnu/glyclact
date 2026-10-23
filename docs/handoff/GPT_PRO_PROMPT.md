# Independent lactate screen audit prompt

You are auditing my master's experiment on lactate production by bacterial
monocultures and synthetic communities in crude-glycerol medium. Use the connected
repository https://github.com/Mvgnu/glyclact. I want a rigorous,
data-driven explanation of why these measurements look inconsistent and what
the smallest useful next experiment is. Do not just repeat the existing reports.

First read `docs/handoff/START_HERE.md` and the two current SOPs in `protocols/`.
The operator-confirmed routine is: 24 h nutrient-rich preculture; OD-based
dilution WITHOUT WASHING; transfer into crude-glycerol medium; 48 h incubation;
spin bacteria down with centrifuge set to 0 °C; 5 µL supernatant + 10 µL
crude-glycerol medium; 15 µL lactate assay mix; 1 h protected from light;
15 µL 1 M acetic acid stop; absorbance read. There is no extra cleaning.
The recorded protocol wavelength is 490 nm, but the actual reader setting needs
confirmation. Do not invent measured oxygenation, pH, endpoint OD or timed steps.

Old SOPs, proposals, generated layouts and empty templates are NOT execution
records. Use them for recipes/design intent, while prioritising raw exports,
full Echo transfer reports and explicit operator corrections. If evidence
conflicts, retain the conflict and state what depends on resolving it.

Please work through the evidence with code and provide:

1. **A run-by-run evidence map.** Identify culture versus assay plates, conditions,
   raw files, standards/blanks, technical culture replicates, analytical aliquots,
   repeat reader runs, read delays, and plate-ID/date ambiguity. Do not count
   multiple reads or analytical aliquots as independent biological replication.

2. **Independent assembly reconstruction.** Check full reports against planned
   picklists, using actual transferred volume, member identity and dose. Separate
   failed/cancelled run records from successful ones. Exceptions are not a second
   delivery. Compare repeatability within genuinely composition-and-dose-matched
   groups, and show which wells remain usable or ambiguous. Reported delivery is
   not proof that a species survived or grew. Do not explain all poor repeatability
   by assembly failures if clean groups still disagree.

3. **An assay and arithmetic audit.** Recalculate standards and blanks, raw and
   justified blank-corrected signals, analytical-pair differences, repeated-read
   changes and culture-to-culture variation. Inspect spatial patterns and compare
   ordinary versus SPR reads without assuming what SPR means. Audit reagent stock
   identities, units, enzyme activity, recipe arithmetic, incubation/stop timing,
   wavelength assumptions, and the distinction between pre-reagent and final-well
   lactate concentrations. Never fit a persuasive concentration scale to flat
   standards or use plate-median subtraction as proof of assay validity.

4. **Explain the apparent standards/spike contradiction.** For October 2 the
   actual spike was 7.5 µL of 12 mM lactate + 5 µL medium + 2.5 µL supernatant,
   versus ordinary 5 µL supernatant + 10 µL medium. Both received reagent and
   acid. Quantify the changed sample/matrix fractions. Ask what the increased
   signal actually establishes; do not call it validated matched spike recovery.
   Distinguish stock pH, reaction pH before acid, and stopped-mixture pH.

5. **Rank explanations rather than listing everything.** Separate demonstrated
   issues, plausible hypotheses and weak speculation. For each leading explanation,
   give the evidence for/against, the exact files/wells supporting it, what remains
   unidentified, and a test that could falsify it. Consider analytical handling,
   matrix/reaction interference, unwashed preculture carryover, actual delivered
   dose, culture growth variability, oxygen/evaporation, and lactate production
   versus consumption or D/L specificity—but assess their quantitative plausibility.
   Do not assume a single cause or that an unmeasured variable was definitely wrong.

6. **What survives biologically?** Are any communities consistently above their
   own constituent monocultures across separately grown wells? Distinguish stable
   effects from outlier-driven means, compare medians and individual replicates,
   and account for the larger total starting biomass of multispecies cultures.
   Genomic capacity is not demonstrated flux. Identify data suitable for raw-signal
   exploratory modelling versus validated lactate modelling; do not train yet.

7. **A minimal decision-making experiment.** We already planned a small October 21
   panel of 41, 200, original-stock 196, and their three-member community, with six
   cultures each and six medium blanks. Results are not supplied yet. Recommend
   the fewest informative controls/checks while retaining the current volumes
   where possible. State exactly what each possible result would change. Do not
   prescribe another four-plate acquisition campaign or a long wish list. Include
   how a free KU Leuven Echo-MS feasibility pilot could independently test the
   samples, remembering it measures total D+L, the current assay targets L, and
   they requested at least 50 µL clear supernatant per submitted well.

Use primary literature and manufacturer documentation for mechanistic claims,
especially the original assay paper https://doi.org/10.1371/journal.pone.0271818.
Check the paper's scope rather than assuming a mammalian-cell assay validates
crude bacterial media or a 384-well miniaturisation. Clearly separate conclusions
from my data from literature-based hypotheses. Do not wander into unrelated
pathogen-screening or taxonomy work unless essential to the assay question.

Deliver a plain-English bottom line first, followed by key diagnostic plots,
a ranked explanation table, a per-run retain/hold decision and a short practical
next-step plan. Cite repository-relative paths and exact wells for important
claims. Include reproducible code and preserve the raw files. If attachments
cannot be read, say what is missing rather than pretending you analysed them.
