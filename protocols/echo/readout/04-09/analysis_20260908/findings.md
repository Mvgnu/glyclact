# Community readouts from 4 September 2026

Analysis date: 8 September 2026. Four destination plates, 102 communities per plate, three technical replicates per community, and two measurement modes of the same wells. The eight reader exports contain 3,072 readings, of which 2,448 belong to the 1,224 inoculated wells. The repeated measurement modes are not additional biological replicates.

The strongest candidate supported by all three wells in both modes is **294 + 196 + 291 + 264**, on plate 2, wells B23, G16, L9. The dataset contains useful community-associated variation, especially on plates 2 and 3. The previous GP, however, did not successfully rank the newly tested communities.

## Mapping and processing

Plate numbers were taken from each CSV's ID1 metadata and mapped using the renamed Echo picklist filenames. Source wells were resolved through `source_plate_map_and_dilution_20260902.csv`. For every destination well, the species reconstructed from the picklist agreed with the community-design table.

| Physical plate | Strategy | Picklist filename |
|---|---|---|
| 1 | Expected improvement | echo_picklist_20260902_bayesian_ei-plate-1.csv |
| 2 | Repeat / max-diversity backfill | echo_picklist_20260902_maxdiv_repeat-plate-2.csv |
| 3 | Posterior mean | echo_picklist_20260902_posterior_mean-plate-3.csv |
| 4 | Random | echo_picklist_20260902_random-plate-4.csv |

All 408 current communities have three assigned wells. Plate 2 includes 54 exact historical compositions and 48 backfills. No current community includes 296.

Each numeric 16 × 24 matrix was located by its numeric content rather than a fixed line offset: these exports include a temperature header absent from the earlier file. Rankings use triplicate medians, retain every reading, and do not delete apparent outliers. The complete 76-well edge ring was used to calculate a descriptive baseline independently for each plate and mode. F16 and K9 are the two uninoculated interior wells.

No concentration standards were available. All numbers below are assay absorbance signals, not calibrated lactate concentrations. The numerical edge offsets are valid descriptive comparisons to the measured empty wells; interpreting them as biochemical blank correction requires those wells to contain the same media and assay reagents without bacteria. Their physical contents were queried and had not yet been confirmed when this report was prepared. No result here requires a concentration conversion.

## Strong candidates

Replicate values are raw readings, listed in the same order as the wells. Larger values are treated as higher assay signal, consistent with the preceding run's interpretation.

| Plate | Community ID | Members | Wells | SPR readings | Single readings |
|---|---|---|---|---|---|
| 2 | MAXDIV_BACKFILL_022 | 294, 196, 291, 264 | B23, G16, L9 | 1.119, 1.120, 1.103 | 0.529, 0.573, 0.777 |
| 3 | BAYESIAN_POSTERIOR_MEAN_013 | 269, 158, 145, 291, 264 | B7, F22, K15 | 1.034, 0.949, 1.003 | 0.480, 0.454, 0.501 |
| 2 | C050 | 269, 158, 145, 196, 291, 264 | D19, I12, N5 | 1.050, 0.877, 1.018 | 0.531, 0.343, 0.602 |
| 2 | C019 | 294, 275 | B3, F18, K11 | 0.997, 0.833, 0.917 | 0.543, 0.602, 0.484 |

The four-member plate-2 candidate exceeds the highest reading of every constituent monoculture in all three wells under both modes. Its lowest community reading is 1.103 in SPR and 0.529 in the single mode. The highest constituent monoculture reading is respectively 0.891 and 0.258. The six-member C050 also exceeds all constituent monoculture readings under both modes. These same-plate comparisons do not depend on a blank subtraction.

The five-member plate-3 candidate and C050 differ only by the presence of 196, making the shared five-member composition interesting for follow-up. They were measured on different plates, so their difference cannot identify the effect of 196.

The pair 294 + 275 gives a compact additional candidate: its single-mode readings are all above every reading of either monoculture, and its SPR median is also above both monoculture medians. Higher starting total biomass in mixtures remains a confound; these observations alone do not establish synergy or a specific interaction mechanism.

Some random communities have still higher medians but less consistent replicates. RANDOM_010 (256, 75, 269, 294, 275, 145, 196; B22/G15/L8) has SPR 1.147/1.113/0.809 and single 0.859/0.724/0.348. RANDOM_050 (256, 269, 294, 275, 41, 158, 264; E22/J15/O8) has SPR 1.113/1.133/0.591 and single 0.680/0.817/0.080. Both are candidates, but the latter includes a third well near the measured empty baseline.

## Four strategies

The table reports the median across the 102 community medians, after subtracting each plate's edge median. This is descriptive: one physical plate represents each strategy.

| Plate / strategy | SPR typical signal | Single typical signal | SPR best community median | Single best community median |
|---|---:|---:|---:|---:|
| 1 / EI | 0.1000 | 0.1315 | 0.3375 | 0.5555 |
| 2 / Repeat and backfill | 0.0435 | 0.0675 | 0.4440 | 0.4780 |
| 3 / Posterior mean | 0.0565 | 0.0835 | 0.3355 | 0.5990 |
| 4 / Random | 0.0700 | 0.0815 | 0.4575 | 0.6235 |

EI has the highest typical observed signal in both modes. Random finds some of the highest community medians. Neither establishes a winning sampling strategy: plate effects are confounded with strategy, no common inoculated reference was repeated across the four plates, and deduplication gave EI first access to overlapping candidates before posterior mean and random. The three new-strategy batches therefore do not constitute an independent comparison of unmodified acquisition policies.

The forecast-versus-observation Spearman correlations within the three entirely new batches are:

| Strategy | SPR | Single |
|---|---:|---:|
| EI | −0.050 | −0.131 |
| Posterior mean | −0.080 | −0.130 |
| Random | 0.158 | 0.146 |

Pooling within-plate ranks over those 306 new communities gives approximately 0.010 for SPR and −0.039 for the single mode. These are prospective evaluations of the stored predictions; no model was refitted. The posterior-mean batch's forecasts themselves occupy a very narrow interval, about 0.019804–0.019954, consistent with weakly differentiated predictions.

For the 54 exact repeated compositions, old-versus-new Spearman correlation is 0.184 with SPR and 0.298 with the new single readings. The corresponding Pearson correlations are 0.085 and 0.103. The original landscape has not transferred strongly across runs. This can reflect biology, preparation, measurement conditions, or their combination; the dataset does not identify a single cause.

## Measurement modes and position

The CSVs contain test names, plate IDs, times, one chromatic channel, one cycle and temperature, but no wavelength, scan geometry, number of sampled positions or other detailed reader settings. The plate-1 SPR-labelled filename has `Testname: MAGNUS_ABS` internally; it is assigned to SPR using the filename and paired export context. That discrepancy needs the saved reader method to resolve. The letters SPR alone do not identify a verified instrument operation.

| Plate | Edge median SPR / single | Interior control F16 SPR / single | Interior control K9 SPR / single | Community-rank agreement between modes | Raw triplicate ICC SPR / single |
|---|---|---|---|---:|---|
| 1 | 0.6685 / 0.0865 | 0.706 / 0.083 | 0.627 / 0.085 | 0.659 | 0.092 / 0.042 |
| 2 | 0.6750 / 0.0950 | 0.677 / 0.085 | 0.659 / 0.093 | 0.787 | 0.446 / 0.408 |
| 3 | 0.6675 / 0.1010 | 0.689 / 0.105 | 0.637 / 0.094 | 0.859 | 0.495 / 0.359 |
| 4 | 0.6555 / 0.1005 | 0.631 / 0.099 | 0.580 / 0.084 | 0.907 | 0.159 / 0.081 |

ICC here is one-way single-measure intraclass correlation over the triplicates; higher values indicate more between-community variation relative to within-community scatter. It is strongest on plates 2 and 3. SPR improves this variance-based metric on all four plates, although rank-based triplicate agreement is not uniformly improved. Its high background also makes raw CVs look artificially small: CV alone would be a misleading basis for preferring SPR.

SPR edge readings have a top-to-bottom gradient: top-row medians are 0.691–0.751, while bottom-row medians are 0.640–0.652. In plate 4, the three interior replicate blocks have raw medians 0.7305, 0.7775, 0.6645 in SPR and 0.1555, 0.2875, 0.1190 in the single mode. Every block contains one replicate of each community, so its composition set is identical: this difference is an experimental/position effect, not different community selection.

As a sensitivity check, removing only additive block-median offsets preserves most primary rankings (Spearman 0.868–0.998) and 7–10 of each top ten. It does not make weak or discrepant replicates disappear. A single reading of 2.043 at plate 1 N15 corresponds to 0.743 in SPR; it should not be treated as a standout community on its own.

## Previously discussed candidates

- C077 without 296, EI plate B19/G12/L5: SPR median 0.738, rank 71/102; single median 0.206, rank 59/102. It is not among this run's strongest candidates. This is a changed composition, so it is not a repeat test of C077 itself.
- C039, plate 2 F10/K3/O18: SPR median 0.696, rank 77/102; single median 0.183, rank 37/102. Modest signal, not a leading result.
- C027, plate 2 B21/G14/L7: SPR median 0.799, rank 24/102; single median 0.220, rank 23/102. Some signal, with one much higher well.
- 294 alone, plate 2 E9/J2/N17: SPR median 0.701; single median 0.196, compared with 0.309 in the previous single-mode experiment. Its earlier top monoculture result did not reproduce at the same magnitude.

## Source files

All paths below are relative to `glycerol_lactate_community/protocols/echo/`.

- Reader data: all eight CSVs in `readout/04-09/`, timestamps 215539, 220017, 220631, 221030, 221716, 222253, 222759 and 223454 on 2026-09-04.
- Plate identity and transfers: four renamed `echo_picklist_20260902_*-plate-*.csv` files in `experiment_setup_output/20260902_adaptive_12strains/`.
- Composition and replicate assignment: `community_design_all_plates_20260902.csv` in that same output directory.
- Source identities: `source_plate_map_and_dilution_20260902.csv` in that same output directory.
- Prospective predictions: `candidate_predictions_20260902.csv` in that same output directory.
- Historical observations: `training_data_20260902.csv` in that same output directory, derived from the 2026-08-30 rerun.

Original reader exports, picklists, design files and the runbook were not modified by this analysis. The two accompanying PNGs show the distribution and replicate comparisons, and the spatial patterns across all eight readouts.

## Follow-up: subtracting each plate's median

At the user's request, a second analysis centers each readout on its own median across **all 384 wells**, including the 76 edge wells and the two uninoculated interior wells. The formula is `community triplicate median − median of the 384 raw readings on that plate in that measurement mode`. A sensitivity comparison also uses the median of the 306 inoculated wells alone. This is a relative-to-plate score, not blank correction or a concentration estimate.

| Plate | SPR median, all 384 | Single median, all 384 | SPR median, inoculated 306 | Single median, inoculated 306 |
|---|---:|---:|---:|---:|
| 1 / EI | 0.7535 | 0.1825 | 0.7685 | 0.2115 |
| 2 / Repeat/backfill | 0.7130 | 0.1490 | 0.7200 | 0.1635 |
| 3 / Posterior mean | 0.7110 | 0.1465 | 0.7220 | 0.1665 |
| 4 / Random | 0.6960 | 0.1225 | 0.7250 | 0.1740 |

| Candidate | SPR, minus edge median | SPR, minus plate median | Single, minus edge median | Single, minus plate median |
|---|---:|---:|---:|---:|
| P2: 294 + 196 + 291 + 264 | 0.4440 | 0.4060 | 0.4780 | 0.4240 |
| P3: 269 + 158 + 145 + 291 + 264 | 0.3355 | 0.2920 | 0.3790 | 0.3335 |
| P2: 269 + 158 + 145 + 196 + 291 + 264 | 0.3430 | 0.3050 | 0.4360 | 0.3820 |
| P2: 294 + 275 | 0.2420 | 0.2040 | 0.4480 | 0.3940 |

All within-plate community rankings are exactly preserved. Subtracting a constant also preserves replicate differences, standard deviations, within-plate Pearson/Spearman correlations, and same-plate differences from constituent monocultures. It does not improve biological reproducibility or resolve spatial gradients. The weak prospective GP rank correlations and weak old-versus-new rank correlations are unchanged.

Across plates, 9 of the previous top 10 remain in the top 10 under each measurement mode after this change. The strongest reliable candidate remains 294 + 196 + 291 + 264. The median-centered results still contain the inconsistent high random candidates; the subtraction does not resolve their low third replicates.

The median across the 102 community medians now becomes:

| Plate | SPR, minus all-well median | Single, minus all-well median |
|---|---:|---:|
| 1 / EI | 0.0150 | 0.0355 |
| 2 / Repeat/backfill | 0.0055 | 0.0135 |
| 3 / Posterior mean | 0.0130 | 0.0380 |
| 4 / Random | 0.0295 | 0.0595 |

These residual medians should not be used to rank strategy success. Centering each differently selected plate removes its typical response, which can contain genuine biological differences as well as technical offsets. All-well medians also include the control wells. The median of triplicate medians need not equal the median of individual inoculated wells: the two aggregation orders are not interchangeable.

The alternative figures are `community_results_overview_plate_median.png` and `plate_position_comparison_plate_median.png`. The earlier edge-centered figures are retained for comparison.
