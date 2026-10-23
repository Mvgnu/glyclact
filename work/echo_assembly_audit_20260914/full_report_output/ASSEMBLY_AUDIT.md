# Echo assembly audit — full reports, 14 September 2026

This supersedes the exception-only analysis in `../output`. Raw files and previous analysis are preserved. The new folder was copied byte-for-byte from the USB drive (15 files, including four macOS metadata sidecars).

## The August reports are different run records

Run **1787930779** (header 28 August 17:31:02) contains 1,971 **Not processed** entries without transfer timestamps. Run **1787930946** (header 17:33:49) contains all **1,971 transfers at the full requested volume**, blank error status, and timed records. Both match the original picklist exactly. The first report has an exception file because every entry is unprocessed; the second report has no error entries. We cannot tell whether the first was cancelled, a pre-execution export, or another unused run record. We do not add the two reports together as deliveries.

**All 102 August triplicates have complete, intended transfer coverage in the successful report.** That removes the apparent blanket failure suggested by the exception file alone. It does not guarantee growth or assay quality.

## September: full transfer accounting

All 7,035 September 2 and all 7,044 September 9 planned transfers match their respective full reports exactly on plate, source well, destination well and requested volume. Their exception files are exact subsets of the full reports; they are not additional transfer attempts. Every reported transfer is either zero or its full requested volume: no partial positive transfers appear.

For September 2, 2,306 transfers report zero. Of those, 396 carry `MM0202006: Problem calc. eo fc...`; 1,910 report fluid height below the minimum. We use the explicit Actual Volume = 0, without guessing what the calculation error means physically. September 9 has 2,179 zero transfers, all minimum-height errors.

The September 9 file is strongly matched to the correct design, and its header says **10 September 00:18**. Its detail timestamps nevertheless say **12 September**, after the reader run on 11 September. The inconsistency is inside the full export, not just its filename. This makes a clock/export problem plausible, but does not establish it. Its readout association remains held pending confirmation. Blank plate barcodes do not independently resolve physical identity.

## What the triplicates received, according to the reports

“Unchanged” means every planned transfer in all three wells is recorded at its requested dose. “Same reduced” means all three received the same smaller nonempty species-and-dose set. “Split” means the three recorded inocula differ. “Empty” means no bacterial transfer was recorded, not that the well contained no medium.

| run | plate | path | reported_zero_transfers | unchanged_triplicates | same_reduced_triplicates | split_triplicates | all_empty_triplicates | empty_wells |
|---|---|---|---|---|---|---|---|---|
| 20260828 | 1 | max_diversity | 0 | 102 | 0 | 0 | 0 | 0 |
| 20260902 | 1 | bayesian_ei | 150 | 52 | 50 | 0 | 0 | 0 |
| 20260902 | 2 | max_diversity | 196 | 35 | 24 | 42 | 1 | 4 |
| 20260902 | 3 | bayesian_posterior_mean | 546 | 7 | 10 | 85 | 0 | 0 |
| 20260902 | 4 | random | 1414 | 0 | 4 | 93 | 5 | 41 |
| 20260909 | 1 | bayesian_ei | 0 | 102 | 0 | 0 | 0 | 0 |
| 20260909 | 2 | max_diversity | 71 | 50 | 2 | 50 | 0 | 0 |
| 20260909 | 3 | bayesian_posterior_mean | 966 | 0 | 2 | 100 | 0 | 7 |
| 20260909 | 4 | random | 1142 | 0 | 61 | 38 | 3 | 11 |

On the September 2 EI plate, all 150 failed transfers were strain **269**: it was missing in all three wells of 50 intended communities. These triplicates remain matched to each other, but need smaller-composition labels. The other 52 EI triplicates have full intended transfer coverage. In particular, within-triplet composition mismatch cannot explain the low EI repeatability in that round.

These reconstructed inocula assume the reported source identities match the loaded strains and there were no additional unrecorded transfers or manual rescue. Positive reported volume is stronger evidence than absence from an exception file, but is not a direct measurement of viable species in the destination.

## Repeatability after assembly filtering

ICC(1,1) is a dimensionless measure of agreement between the three original wells: 1 is perfect agreement; near 0 means within-triplet scatter dominates the between-community differences. Negative estimates mean within-triplet variation exceeds between-triplet variation. Very small subsets cannot support reliable ICC estimates. The range column is median (maximum minus minimum) **raw absorbance**, not lactate concentration. Centering each plate by one constant does not change either statistic.

| run | plate | subset | community_count | icc_1_1 | median_triplicate_range |
|---|---|---|---|---|---|
| 20260828 | 1 | all_intended_triplicates | 102 | 0.089 | 0.0330 |
| 20260828 | 1 | no_reported_exception_triplicates | 102 | 0.089 | 0.0330 |
| 20260828 | 1 | same_reduced_nonempty_triplicates | 0 | n/a | n/a |
| 20260902 | 1 | all_intended_triplicates | 102 | 0.042 | 0.1660 |
| 20260902 | 1 | no_reported_exception_triplicates | 52 | 0.004 | 0.1860 |
| 20260902 | 1 | same_reduced_nonempty_triplicates | 50 | 0.105 | 0.1515 |
| 20260902 | 2 | all_intended_triplicates | 102 | 0.408 | 0.0950 |
| 20260902 | 2 | no_reported_exception_triplicates | 35 | 0.387 | 0.1110 |
| 20260902 | 2 | same_reduced_nonempty_triplicates | 24 | 0.285 | 0.0815 |
| 20260902 | 3 | all_intended_triplicates | 102 | 0.359 | 0.0840 |
| 20260902 | 3 | no_reported_exception_triplicates | 7 | 0.743 | 0.0710 |
| 20260902 | 3 | same_reduced_nonempty_triplicates | 10 | 0.586 | 0.0820 |
| 20260902 | 4 | all_intended_triplicates | 102 | 0.081 | 0.2480 |
| 20260902 | 4 | no_reported_exception_triplicates | 0 | n/a | n/a |
| 20260902 | 4 | same_reduced_nonempty_triplicates | 4 | -0.021 | 0.1595 |
| 20260909 | 1 | all_intended_triplicates | 102 | 0.016 | 0.0295 |
| 20260909 | 1 | no_reported_exception_triplicates | 102 | 0.016 | 0.0295 |
| 20260909 | 1 | same_reduced_nonempty_triplicates | 0 | n/a | n/a |
| 20260909 | 2 | all_intended_triplicates | 102 | 0.073 | 0.0200 |
| 20260909 | 2 | no_reported_exception_triplicates | 50 | 0.064 | 0.0185 |
| 20260909 | 2 | same_reduced_nonempty_triplicates | 2 | -0.400 | 0.0275 |
| 20260909 | 3 | all_intended_triplicates | 102 | 0.037 | 0.0310 |
| 20260909 | 3 | no_reported_exception_triplicates | 0 | n/a | n/a |
| 20260909 | 3 | same_reduced_nonempty_triplicates | 2 | -0.431 | 0.0335 |
| 20260909 | 4 | all_intended_triplicates | 102 | 0.040 | 0.0180 |
| 20260909 | 4 | no_reported_exception_triplicates | 0 | n/a | n/a |
| 20260909 | 4 | same_reduced_nonempty_triplicates | 61 | 0.075 | 0.0190 |

Compare subsets cautiously: filtering changes which communities and signal range are represented. These are descriptive comparisons, not a causal estimate of the effect of assembly failure. Low agreement persists in the September 9 EI plate even though every transfer is reported successful (ICC 0.016). Assembly errors therefore do not explain all of the poor repeatability. Conversely, split triplicates are not repeated measurements of their intended community and should not be averaged under that label.

## What to retain for a future fit

- **Original labels:** `reported_intact_triplicate_wells.csv` contains 348 report-intact original triplicates (1044 wells), with the three individual readouts and provenance. September 9 rows retain an explicit timing/identity hold. Other report-intact wells from partially broken triplicates remain in the complete well table; do not fabricate triplicates for them.
- **Changed labels:** `reported_regrouped_candidates.csv` groups nonempty recorded species AND per-species dose, within each run and plate. It retains every original well/design ID, the individual readings, and the real replicate count. There are 961 groups; 510 have at least three wells. Repeated wells are technical replicates, not independent biological experiments. Reduced inocula also have reduced total biomass; presence-only modeling should not silently treat dose as irrelevant.
- **Empty inocula:** keep as distinct assay/control observations, not as their intended communities.
- **Training paths:** keep the shared max-diversity plus own EI or posterior-mean lineage. A relabelled random observation stays random. Do not reassign acquisition history based on its reconstructed membership.
- **Responses:** preserve raw absorbance and the historical plate-median-centered response. Assembly filtering does not validate the lactate calibration; agree on background handling and response noise before fitting. Do not filter on whether a result looks good.

**No models were retrained. All candidate exports retain `training_approved=False`.** This is a review gate, not a claim that the data are worthless. The August run now has positive transfer evidence; the two later runs can be separated into report-intact, relabelled and split groups instead of throwing everything away or training on the wrong labels.

## Traceability

`transfer_evidence.csv` links every planned transfer to a full-report file and original CSV line. `well_reconstruction.csv` links every culture well to its intended and reported inoculum and reader value. `triplicate_audit.csv` lists each original triplet and its individual reconstructions. `manifest.json` records input SHA-256 hashes and run headers. Historical designs, source maps, readouts, runbook, slides and fitted models are unchanged.
