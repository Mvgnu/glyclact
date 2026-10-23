# 23 September 2026 — one-plate analytical diagnostic round

## Purpose and scope

Assemble one small culture panel now. Measure growth, then run its samples,
standards and analytical controls together in one later assay session using one
fresh assay mix. No standards-only session today. One culture plate and one
separate assay plate are required; dilution/staging plates are handling plates,
not additional experimental culture plates.

This round asks whether low readings come from poor growth, low endpoint
L-lactate accumulation, sample interference or measurement variability. It does
not select EI/TM/random/max-diversity candidates, fit models, update acquisition
done lists or automatically approve observations for training.

Historical OD blocks, designs, measurements and models remain unchanged.
Precultures grow in nutrient medium and are transferred into crude-glycerol
screen medium, as in earlier rounds. Pellets were noticed in some previous wells,
but endpoint OD and saved untreated supernatants are unavailable for those runs.

## 1. Today's editable OD block and export

Enter fresh **1:10 OD600 readings** in
`OD600_READING_AT_1_TO_10_20260923` at the end of `adaptive_run.py`.
The calculation multiplies each reading by ten. `None` means pending; do not
copy the September 16 readings or use zero for missing data.

Run the new final section of `experiment_setup.ipynb`, or run
`analytical_round_20260923.py`. Do not run all historical notebook cells.
Running `adaptive_run.py` directly still invokes the historical September 2 run.

Each execution creates a new timestamped snapshot under
`experiment_setup_output/20260923_analytical_one_plate/`.
With missing OD, it exports the layouts and recording templates but **no Echo
picklist**. Once all ODs are entered, import only the resulting
`plate-1_echo_picklist_20260923_analytical.csv` into Echo. Other CSVs are maps,
audits and recording sheets, not instrument picklists.

## 2. Source plate and unchanged inoculation

Source plate **396**, one row only, **50 µL in each A1–A12**:

| Source well | Species | Source well | Species |
|---|---|---|---|
| A1 | 256 | A7 | 158 |
| A2 | 75 | A8 | 145 |
| A3 | 269 | A9 | 200 |
| A4 | 294 | A10 | 196 |
| A5 | 275 | A11 | 291 |
| A6 | 41 | A12 | 264 |

296 is excluded. The generated source plan gives culture and diluent volumes
**per 50 µL source well**, targeting working OD600 1.20 where dilution is possible.
Weaker cultures remain undiluted. Use the same source-suspension preparation and
diluent as the preceding run; record any deviation.

Keep **0.100 OD600·µL per member**, rounded to the nearest **25 nL**.
No new minimum transfer volume, no 15 µL destination cap and no automatic dose
switch. Mixtures therefore receive more total biomass than monocultures; this
round does not isolate a biomass-independent synergy effect.

The source-demand output sums all ejections from each physical well and reports
liquid remaining from 50 µL. No acoustic dead-volume guarantee is inferred from
that subtraction: inspect the Echo source survey. Save the full transfer report
and any warnings, even if the run appears successful. No refill is planned.

## 3. Culture plate 1

**50 µL screen medium per occupied culture/control well, plus inoculum for the
cultures.** Do not reduce the medium to compensate for inoculum. Do not change
volumes to accommodate the analytical checks.

| IDs | Conditions | Independently grown wells each |
|---|---|---|
| M01–M12 | 256, 75, 269, 294, 275, 41, 158, 145, 200, 196, 291, 264 | 6 |
| C01 | 41 + 200 + 196 | 6 |
| C02 | 256 + 269 + 294 + 41 + 145 + 200 | 6 |
| C03 | 294 + 275 | 6 |
| BLANK | Uninoculated screen medium, handled and incubated alongside cultures | 6 |

Total **90 inoculated + 6 control wells**. All are interior: rows **B, D, F, H,
J, L**, columns **2–17**. Each row contains all 15 conditions and one control.
Condition order rotates across rows to avoid fixing a species/community to one
column. Use the exported map, not a guessed species-to-column order.

Keep the preceding run's incubation temperature, shaking and seal, and its
48-hour endpoint. Record actual settings and timestamps in `session_record.json`.
48 h is the retained comparison endpoint, not a claim that all species produce
lactate specifically in stationary phase. Unassigned/edge wells are not
automatically blank controls; record any liquid added for evaporation management.

## 4. Harvest and growth measurement

1. Read **OD600 before centrifugation/clearing**, using consistent mixing and
   the six incubated medium controls. Keep the raw reader export and well IDs.
2. Record the read time, harvest time and any visible abnormalities in
   `plate-1_endpoint_OD600_and_harvest.csv`. Missing OD stays blank, not zero.
3. Clear cultures using the existing procedure. Keep each culture well separate.
   Do not pool the six independently grown wells.
4. Record recovered supernatant volume, especially for the wells used for extra
   analytical checks. Preserve a reserve if available. The generated sample budget
   excludes pellet retention, evaporation, pipette and staging-plate dead volume.

## 5. One combined assay session

Prepare assay reagents once, when the culture samples are ready, using the
existing high-LDH recipe in the assay SOP. Plan a reader slot for the prompt read.
Normal analytical wells remain **15 µL sample/standard + 15 µL assay mix**, then
**1 h at room temperature in the dark**, then **15 µL 1 M acetic acid**.
Final volume remains **45 µL**, read at **A490** using the existing method.

The operator confirmed the actual preceding procedure on 23 September:
**5 µL supernatant + 10 µL fresh crude-glycerol screen medium**, giving 15 µL
sample aliquot. Keep this unchanged. The diluent is NOT water. This dilutes
culture-derived changes threefold but does not dilute the common medium matrix
threefold: the aliquot is one-third spent medium and two-thirds fresh medium.
The generator refuses plans whose sample demand exceeds original culture volume.

### Primary assay duplicates and VIAFLO 384 head

Each of the 96 harvested samples contributes two **15 µL aliquots of the same
prepared sample**. First aliquot has the same well coordinate as the culture;
second is one row below at the same column:

| Source/dilution row | Assay aliquot 1 | Assay aliquot 2 |
|---|---|---|
| B | B | C |
| D | D | E |
| F | F | G |
| H | H | I |
| J | J | K |
| L | L | M |

The second dispense therefore uses a **+1-row destination offset**, not a new
culture source. Keep source mapping fixed. Every dispense is still 15 µL; do not
split one 15 µL dispense into two 7.5 µL wells. Use the usual contamination-safe
tip practice and do not return tips that contacted reaction mixture to samples.

This is a **liquid-well map**, not an automated instrument program. Only rows
B/D/F/H/J/L and columns 2–17 carry samples in the primary dilution/staging map.
The intervening positions must not carry other samples that would be transferred
into the duplicate wells. With a 384 head, check the active tip pattern, unused
channels and physical clearance of the one-row shift before liquid transfer.
Do not assume a fully populated 384-tip footprint can be shifted beyond a plate
edge safely. If that motion is not supported with the loaded tips, retain the
same well map/volumes and transfer the pairs with a compatible tip pattern or
multichannel method. Record the actual method.

The 12 primary assay rows occupy columns 2–17: **192 assay wells**.
These are 2 analytical aliquots from each culture, not 12 independent cultures
per condition. Do not count analytical duplicates as independent replicates.

### Extra checks, prepared from the same harvested samples

For **culture replicate 1 of each of the 15 biological conditions**, preselected
without looking at its signal:

- Two lactate-spiked analytical aliquots, in adjacent rows **N/O**, columns 2–16.
- Two sample-background aliquots, in adjacent rows in columns **18–22**. Replace
  the 15 µL assay mix with 15 µL plain 0.2 M Tris pH 8.2, then add the same stop
  solution. These report optical background, not every possible chemical
  interference. Do not dispense full reaction mix into these positions.

Spiked samples retain the **same sample fraction and total 15 µL sample aliquot**
as unspiked samples. The planned analytical spike is **+3 mM L-lactate before
adding assay reagents**, made using a 12 mM working standard **in the same
full-strength fresh crude-glycerol screen medium**, not the water-standard stock.
Keep medium components and solvent vehicle matched in the fortified and plain
medium preparations. At the confirmed 1:3 sample dilution, per aliquot this is:

| Aliquot | Original supernatant | 12 mM lactate in fresh screen medium | Plain fresh screen medium | Total |
|---|---|---|---|---|
| Unspiked or sample background | 5 µL | 0 µL | 10 µL | 15 µL |
| Spiked | 5 µL | 3.75 µL | 6.25 µL | 15 µL |

Prepare enough homogeneous working sample for each pair and the actual handling
losses. Do not add spike on top of a completed 15 µL aliquot. Match the diluent
vehicle between treatments. Both treatments contain 5 µL spent plus 10 µL
fresh-medium vehicle. Replacing 3.75 µL medium with an aqueous spike would change
the matrix and confound recovery, so do not use that earlier provisional recipe.
+3 mM here becomes +1 mM in the final 45 µL reaction;
retain the pre-reagent concentration convention for standards and recovery.

At 1:3, ordinary wells require **10 µL original supernatant** for their two assay
aliquots. Selected replicate-1 wells require **30 µL original supernatant** for
all six aliquots, before losses. The original 50 µL culture fill is unchanged.
Confirm recovered volume and instrument dead volume before preparing aliquots;
if insufficient, record omitted extra checks rather than silently reducing
dispense volumes or inventing results.

### Standards and assay controls in the same session

Use **0, 0.375, 0.75, 1.5, 3, 6 and 12 mM** before assay reagents, three wells
per level, in water and fresh screen-medium matrix: **42 wells**. The screen
standards use **full-strength fresh screen medium**, not one-third-strength medium.
Maintain medium composition throughout the dilution series using a matched
medium blank as diluent. Record actual preparation. Fresh-medium standards
approximate the base medium matrix; they do not reproduce bacterial changes in
the spent-medium component. That is why sample-specific spike recovery is needed.
This corrects the earlier provisional medium-fraction assumption without changing
culture conditions or sample dispense volumes.

The map intersperses both standard matrices within the interior. Neither matrix
is restricted to one plate edge. Do not use empty wells or a plate median as a
substitute for the explicitly mapped blanks.

Total assay occupancy: **192 primary + 30 spiked + 30 sample-background + 42
standards = 294 wells** on one 384-well assay plate. Keep standards and supplemental
checks out of the primary offset-transfer path until those transfers are done.
Use the assay map to stage/add reaction mix versus plain Tris correctly; a
blanket full-head reaction-mix dispense would invalidate the background controls.
Reagent quantities must include the actual dispensing method's reservoir/staging
dead volume, not just 15 µL multiplied by occupied wells.

Read promptly after stop, recording time. A later reread of the **same stopped
plate** can check drift without preparing new assay reagents; record delay,
temperature and light exposure. It is optional and does not replace the prompt
read. Compare each well to itself and standards to themselves.

## 6. Analysis to perform after readout

Preserve raw OD600 and absorbance exports. Join by culture and assay well maps.
Do not infer success from a smooth standard curve alone, exclude discordant wells
without a recorded reason, or replace low/negative corrected values with zero.

1. Plot standard replicate scatter and response, especially near the blank.
   Establish the useful measurement range from this batch, not an imported
   universal detection threshold.
2. Compare identical-liquid analytical duplicates before averaging them. Inspect
   row/column offsets using known identical samples and standards.
3. Compare spiked and matched unspiked aliquots. Recovery is the measured
   increment divided by the known added increment, on the same dilution scale.
   Investigate sample-specific suppression or enhancement. No universal pass/fail
   percentage is imposed here.
4. Inspect sample backgrounds for the selected culture wells. Do not apply one
   culture replicate's background to all six as if independently measured.
5. Only after checking analytical behaviour, compare the six culture wells per
   condition alongside their endpoint OD. OD is a growth indicator, not proof
   of lactate production or a validated cross-species biomass conversion.
6. If cultures grew, repeated measurements agree and spikes recover but the
   unspiked samples stay near blank, conclude low **endpoint L-lactate signal**
   within the validated range, not necessarily no lactate production at any time.

No models are retrained and no communities become acquisition-done merely by
appearing on this diagnostic plate. Decide data eligibility after reviewing the
assembly report and these analytical results.

## Sources and rationale

- Existing local assay SOP: `../SOP — Colorimetric L-Lactate Assay for Bacterial Supernatants.md`.
- Existing culture SOP: `../SOP — Community Assembly and Crude-Glycerol Fermentation Screen.md`.
- Assay method and documented pyruvate interference:
  https://doi.org/10.1371/journal.pone.0271818
- Candidate selection: September 19 unspiked triplets were 0.118/0.115/0.123 AU
  for C01 and 0.124/0.132/0.124 AU for C02. C03 had complete recorded transfers
  and high readings on September 4, but was near background on September 19.
  These are reasons to investigate, not confirmed high-lactate producers.
