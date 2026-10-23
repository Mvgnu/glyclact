# Echo assembly audit, 14 September 2026

## Main finding

The later exception file reports 2,179 zero-volume transfers. Every row matches the final September continuation picklists on source well, destination well and requested volume. Under the explicit reconstruction assumptions below, 188 of 408 intended triplicates split into different inoculum compositions or doses. The data therefore cannot be averaged under their intended labels without reviewing assembly.

**This is a provisional reconstruction, not verified composition or approved training data.** The user may have supplied the wrong attempt. Internal Echo timestamps are 12 September 00:24–00:38, later than the associated reader timestamps on 11 September 23:30–23:43. The protocol name and filename suffix point to 9 September, but names cannot settle execution identity. Full run/transfer reports and retries are needed.

## What each supplied file establishes

| Setup / readout | Available evidence | Decision now |
|---|---|---|
| 28 Aug / 30 Aug | All 1,971 planned transfers are listed as Not processed, zero volume, with no timestamps or destination identity | Execution unresolved. Do not label every culture empty, delete the whole experiment as failed, or call it clean from this file. |
| 2 Sep / 4 Sep | No matching execution log supplied | Keep intended labels historical, but do not certify them for clean retraining. |
| 9 Sep design / 11 Sep readout | 2,179 timed zero-volume exceptions match plates 2–4. Association to measured plates remains provisional. | Reconstruct alternatives and hold all candidate training labels. |

The August file is an exact 1,971/1,971 match to the complete original picklist. Its Not processed state describes the supplied protocol/source entry. It does not distinguish a never-run or unused entry, an aborted attempt, a later successful retry, or the final physical run. A generic zero in that context is not equivalent to a timed instrument failure with source-survey metadata.

All 2,179 later exceptions explicitly report Actual Volume = 0 and MM0202001: Current fluid height below minimum transfer level. There are no positive-actual-volume exceptions or duplicate transfer keys in either supplied file. These later entries are evidence of zero delivery for the recorded attempt, not proof that no subsequent attempt delivered the species.

## Conditional reconstruction of the September continuation

Assume the matched log is the final attempt for these measured plates, no unrecorded retries or manual replacements occurred, and every unlisted picklist transfer completed. Keep source identities and per-species requested doses from the archived run. Remove only the logged zero-volume transfers. No readout value influences these classifications.

| plate | path | reported_zero_transfers | affected_wells | unchanged_triplicates | same_reduced_triplicates | split_triplicates | all_empty_triplicates | empty_wells |
|---|---|---|---|---|---|---|---|---|
| 1 | bayesian_ei | 0 | 0 | 102 | 0 | 0 | 0 | 0 |
| 2 | max_diversity | 71 | 71 | 50 | 2 | 50 | 0 | 0 |
| 3 | bayesian_posterior_mean | 966 | 299 | 0 | 2 | 100 | 0 | 7 |
| 4 | random | 1142 | 304 | 0 | 61 | 38 | 3 | 11 |

This leaves 152 intended triplicates (456 wells) with no reported exception, plus 94 unflagged wells in partially affected triplicates. 65 triplicates would share the same reduced, nonempty composition. 3 would be completely empty in all three wells. In total 18 wells would receive no bacteria under this scenario.

Plate 1 has no rows in this exception file. That is absence of reported exceptions, not a successful-transfer report. Plates 3 and 4 are affected throughout almost all culture wells. Apparent replication under the original design labels is particularly misleading there.

## Does matching composition fix repeatability?

The calculations below provisionally associate the log with the 11 September readouts. ICC(1,1) uses the three original technical wells of each design, after restricting by execution category. It is dimensionless, with 1 indicating perfect agreement. Subsets contain different communities and signal ranges, so changes in ICC are descriptive, not a causal effect estimate. No low/high absorbance values or statistical outliers are deleted.

| plate | path | subset | community_count | icc_1_1 | median_triplicate_range |
|---|---|---|---|---|---|
| 1 | bayesian_ei | all_intended_triplicates | 102 | 0.016 | 0.0295 |
| 1 | bayesian_ei | no_reported_exception_triplicates | 102 | 0.016 | 0.0295 |
| 1 | bayesian_ei | same_reduced_nonempty_triplicates | 0 | n.a. | n.a. |
| 1 | bayesian_ei | different_provisional_compositions | 0 | n.a. | n.a. |
| 2 | max_diversity | all_intended_triplicates | 102 | 0.073 | 0.0200 |
| 2 | max_diversity | no_reported_exception_triplicates | 50 | 0.064 | 0.0185 |
| 2 | max_diversity | same_reduced_nonempty_triplicates | 2 | -0.400 | 0.0275 |
| 2 | max_diversity | different_provisional_compositions | 50 | 0.069 | 0.0200 |
| 3 | bayesian_posterior_mean | all_intended_triplicates | 102 | 0.037 | 0.0310 |
| 3 | bayesian_posterior_mean | no_reported_exception_triplicates | 0 | n.a. | n.a. |
| 3 | bayesian_posterior_mean | same_reduced_nonempty_triplicates | 2 | -0.431 | 0.0335 |
| 3 | bayesian_posterior_mean | different_provisional_compositions | 100 | 0.036 | 0.0310 |
| 4 | random | all_intended_triplicates | 102 | 0.040 | 0.0180 |
| 4 | random | no_reported_exception_triplicates | 0 | n.a. | n.a. |
| 4 | random | same_reduced_nonempty_triplicates | 61 | 0.075 | 0.0190 |
| 4 | random | different_provisional_compositions | 38 | 0.010 | 0.0175 |

The EI plate still has the same low ICC when all triplicates have no reported exception. The exception file therefore cannot explain all disagreement. Unlogged delivery errors, growth variation, position effects and assay variation remain possible. Conversely, split triplicates cannot be used as evidence of biological repeatability for their original intended composition.

Pairwise differences within the original triplets, divided by whether reconstructed inoculum signatures agree:

| plate | same_provisional_signature | pairs | median_absolute_difference |
|---|---|---|---|
| 1 | True | 306 | 0.018000000000000002 |
| 2 | False | 106 | 0.012999999999999998 |
| 2 | True | 200 | 0.012999999999999998 |
| 3 | False | 265 | 0.018000000000000002 |
| 3 | True | 41 | 0.024999999999999994 |
| 4 | False | 81 | 0.009999999999999995 |
| 4 | True | 225 | 0.010999999999999996 |

## What could be retained for retraining

1. **Original-composition candidates:** the held export lists the 152 triplicates with no reported exception, including source evidence and individual reader values. Retain only after final-attempt identity and successful transfer coverage resolve. Isolated unflagged wells are available in the full well table with their original replicate identities.
2. **Relabelled candidates:** regroup wells by inferred species AND the complete per-species dose signature, within the same physical plate and run. The held regrouped export retains every original design and well. It includes the number of wells, scatter and individual readings. Empty reconstructed inocula are excluded from community-training candidates.
   The provisional regrouping creates 404 nonempty plate/composition/dose groups, of which 205 have at least three wells. These are not that many independent biological replicates. Replicate counts of one or two remain explicit rather than being padded or treated as triplicates.
3. **Unresolved:** retain raw data but do not fit composition-response pairs from the August Not processed file, the unlogged September 2 run, ambiguous attempts, or unverified reduced communities.
4. **Model paths:** after confirmation, the EI model should receive audited shared max-diversity data and audited EI-path observations. The posterior-mean model should receive the corresponding shared and own-path data. Keep random outcomes out of both training paths, as in the experiment design. Do not move a relabelled random observation into an EI model merely because its inferred composition resembles an EI candidate.

**No model has been retrained. `training_approved.csv` contains headers and zero observations. All held candidates have `training_approved=False`.** This prevents the provisional reconstruction being silently used as ground truth.

The original plate-median-centered response is retained as a traceable reference in the held exports. It is not recomputed from the retained subset. Filtering and regrouping can change the composition distribution and the background represented by the old plate center. Decide the training response normalization explicitly once valid wells and controls are established.

## Needed next

Obtain full transfer reports and run summaries for the 28 August, 2 September and 9 September setups. Include every attempt or retry, source/destination identities or barcodes, requested and actual transfer volumes, timestamps, survey records if available, and the saved protocol/run file. A final actual-volume record can resolve logged failures that were subsequently retried. Do not sum duplicated exports as if they were independent deliveries.

The full report for the later attempt should also explain why its internal timestamp follows the reader run. Exact picklist overlap distinguishes the design but cannot uniquely identify an execution of that design.

## Files

- `well_reconstruction.csv`: all 2,754 culture-well measurements across the three rounds, intended and provisional labels, flags and raw readout.
- `transfer_evidence.csv`: one planned source-to-destination transfer per row with original log line references.
- `triplicate_audit.csv`: 918 original triplicate groups with individual well values and inferred alternatives.
- `provisional_unchanged_triplicate_wells_HELD.csv`: candidate unchanged triplicates, not approved for fitting.
- `provisional_regrouped_observations_HELD.csv`: conditional nonempty composition/dose groups, not approved for fitting.
- `repeatability_comparison.csv`: all original round/plate ICCs and the conditional latest-round subsets.
- `manifest.json`: assumptions and SHA-256 hashes of every used input.

Raw source files, historical experiment outputs, runbook, slides and saved models remain unchanged.
