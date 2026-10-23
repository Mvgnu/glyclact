# September 19 readout: 20-hour delay after acid stop

Four reader files copied unchanged from the USB folder. Plate IDs in each CSV match the September 16 designs. Analysis uses all readings; no outlier deletion, model fitting, or historical label changes. The user reports approximately 20 hours between acid stop and reading because the machine was occupied. Storage/light/temperature during the delay were not specified. The CSV reader temperatures are metadata, not storage temperatures.

## Main result

The validation communities do **not** reproduce their earlier advantage over their constituent monocultures. Their technical wells have small absolute differences, but are near the two interior controls (0.085 and 0.088). Close low readings do not establish that the earlier strong biological signal reproduced or that the assay is fully validated.

| members | wells | raw_readings | median |
|---|---|---|---|
| 294 + 196 + 291 + 264 | E3 + I18 + N11 | 0.072 + 0.079 + 0.063 | 0.0720 |
| 294 + 291 + 264 | C10 + H3 + L18 | 0.079 + 0.077 + 0.090 | 0.0790 |
| 294 + 275 | F11 + K4 + O19 | 0.073 + 0.076 + 0.106 | 0.0760 |
| 294 | F5 + J20 + O13 | 0.069 + 0.098 + 0.072 | 0.0720 |
| 275 | D18 + I11 + N4 | 0.103 + 0.092 + 0.076 | 0.0920 |
| 196 | B22 + G15 + L8 | 0.086 + 0.079 + 0.071 | 0.0790 |
| 291 | F9 + K2 + O17 | 0.076 + 0.083 + 0.094 | 0.0830 |
| 264 | C9 + H2 + L17 | 0.101 + 0.070 + 0.088 | 0.0880 |

Same-plate comparison against each community's own constituent monocultures:

| community | community_median | highest_constituent_median | median_difference |
|---|---|---|---|
| 294 + 196 + 291 + 264 | 0.0720 | 0.0880 | -0.0160 |
| 294 + 291 + 264 | 0.0790 | 0.0880 | -0.0090 |
| 294 + 275 | 0.0760 | 0.0920 | -0.0160 |

The four-member candidate's two previously fully transferred wells read 0.573 and 0.777 on September 4. Its third historical well lacked 196 and read 0.529. The repeated 294+275 pair previously read 0.543, 0.602 and 0.484 with full reported transfers. Today's values are much smaller, but batches and readout timing differ; this does not identify the delay, growth or assay preparation as the cause. The operator confirms correct assembly with no Echo warnings this time. A full transfer export was not supplied, so that execution status is operator-reported, not independently log-audited.

## Edge versus interior references

The operator identifies A1–F1 as crude-medium standards and A24–F24 as water standards. BOTH six-well series are excluded from the edge reference, leaving 64 edge wells. The two designed interior no-bacteria wells are F16 and K9. Exact blank media/reagent equivalence was not explicitly confirmed. The whole-plate median includes all 384 wells, matching the previous modeling convention. The first quick comparison included right-hand standards and is superseded by these corrected results.

| plate | plate_median | edge_nonstandard_median | control_F16 | control_K9 | interior_control_median |
|---|---|---|---|---|---|
| 1 | 0.1040 | 0.0930 | 0.0770 | 0.0940 | 0.0855 |
| 2 | 0.0900 | 0.0875 | 0.0850 | 0.0880 | 0.0865 |
| 3 | 0.0895 | 0.1075 | 0.0750 | 0.0980 | 0.0865 |
| 4 | 0.0870 | 0.0985 | 0.1000 | 0.0710 | 0.0855 |

After excluding both standard series, the edge reference is higher than the two-control median by 0.0075, 0.0010, 0.0210 and 0.0130 AU on plates 1–4, respectively. This supports investigating an edge/interior offset, especially on plates 3 and 4, but two interior controls per plate do not establish a precise spatial background surface. F16/K9 disagree substantially on plates 1, 3 and 4, so their mean/median is not automatically a validated blank either. On plate 2 the edge and interior-control medians differ by only 0.001 AU.

Whole-plate median subtraction is centering, **not media-blank subtraction**: it contains the cultures' signals and can remove a real broad shift in community performance. Edge subtraction can instead impose an edge-specific offset. Keep raw readings and show both as sensitivity analyses. Do not pick the method that makes a favored result look best.

Subtracting one common number from every well on plate 2 cannot change which community exceeds which monoculture: (community − baseline) − (mono − baseline) = community − mono. It also cannot improve within-plate replicate agreement. Therefore an edge-versus-median choice does not rescue the validation result.

## Repeatability across the full experiment

ICC(1,1), dimensionless: 1 is perfect agreement; near zero indicates weak separation of communities relative to within-triplet scatter. Negative estimates indicate greater within-triplet than between-triplet variation. The range is maximum minus minimum raw absorbance per triplet, then median across triplets. Validation repeats are excluded only from the new-acquisition ICC, not from the all-condition ICC.

| plate | path | icc_all | icc_new | median_triplicate_range |
|---|---|---|---|---|
| 1 | bayesian_ei | 0.2982 | 0.2982 | 0.0205 |
| 2 | max_diversity | 0.0454 | 0.0367 | 0.0240 |
| 3 | bayesian_posterior_mean | -0.0188 | -0.0188 | 0.0330 |
| 4 | random | 0.1447 | 0.1447 | 0.0410 |

For reference, September 11 intended-triplet ICCs were approximately 0.016, 0.073, 0.037 and 0.040 for plates 1–4; that run had documented composition errors. Current agreement is better for EI and random by this metric, but not uniformly better. Small absolute scatter near the blank is different from reliably distinguishing high- and low-performing communities. Different candidate sets and signal ranges also affect ICC.

## Acquisition summaries, excluding validation repeats

| plate | path | conditions | median_raw | median_plate_centered | median_edge_centered | best_median |
|---|---|---|---|---|---|---|
| 1 | bayesian_ei | 102 | 0.1055 | 0.0015 | 0.0125 | 0.1310 |
| 2 | max_diversity | 94 | 0.0895 | -0.0005 | 0.0020 | 0.1380 |
| 3 | bayesian_posterior_mean | 102 | 0.0860 | -0.0035 | -0.0215 | 0.1160 |
| 4 | random | 102 | 0.0805 | -0.0065 | -0.0180 | 0.1570 |

One physical plate represents each strategy, so differences are confounded with plate effects. Cross-plate ranks depend on the reference choice. These summaries are not proof that one acquisition strategy is superior. The pre-recorded predictions are compared with the new measurements in prospective_rank_checks.csv; no models were refit to these data.

## What the delayed read can establish

Both confirmed standard series show an overall concentration-associated response, especially the crude-medium series. This supports the narrower statement that the assay still responds to lactate after the delay. It does not establish precise concentrations in the cultures. The following descriptive straight-line fits assume the previous 12 mM top, twofold dilution series (12, 6, 3, 1.5, 0.75, 0.375 mM); the user confirmed matrix and positions but did not explicitly restate concentrations. No culture concentration conversion is made.

| plate | matrix | slope_AU_per_nominal_mM | intercept | r_squared | highest_minus_lowest_standard | n_dilution_order_reversals |
|---|---|---|---|---|---|---|
| 1 | crude | 0.0109 | 0.0863 | 0.9960 | 0.1280 | 0 |
| 1 | water | 0.0055 | 0.1110 | 0.6573 | 0.0800 | 1 |
| 2 | crude | 0.0085 | 0.0946 | 0.9593 | 0.1030 | 0 |
| 2 | water | 0.0035 | 0.1082 | 0.4462 | 0.0390 | 1 |
| 3 | crude | 0.0089 | 0.0987 | 0.9326 | 0.0980 | 1 |
| 3 | water | 0.0040 | 0.1263 | 0.4320 | 0.0790 | 2 |
| 4 | crude | 0.0092 | 0.0784 | 0.9854 | 0.1030 | 0 |
| 4 | water | 0.0054 | 0.0989 | 0.6706 | 0.0600 | 1 |

Water standards have dilution-order reversals on every plate, including plate 2's highest standard below its next standard. Matrix and position are confounded because crude is always left and water always right. Neither slope differences nor baseline differences alone establish a biochemical matrix effect. There are no within-plate standard replicates at each dose.

A single delayed measurement cannot show that the stop solution preserved signal or rankings for 20 hours: a matched prompt-read comparison would be required. Nor can it prove that delay caused the loss of the earlier community signal.

Current conclusions: the community-over-monoculture claim did not reproduce in this read; there is a plausible edge-reference offset; some within-plate variation is measurable, but broad repeatability remains mixed. Preserve this run as a delayed-read batch rather than silently treating its absolute signal as interchangeable with earlier batches.

## Files and checks

validation_repeats.png shows individual validation wells and the interior-control band. plate_heatmaps.png shows all four plates on one scale. standards.png shows both standard series. CSVs retain every culture well, all 384 positions, raw and differently centered values, condition roles, and summaries. manifest.json records input hashes, plate metadata, and the 20-hour-delay note. Raw inputs, runbook, saved designs and fitted models remain unchanged.
