# Assembly round 16 September 2026

Picklists ready: A/B/C/D feed plates 1/2/3/4 respectively, no refills. Minimum calculated source remainder 14.000 uL. Inspect the Echo source survey; this calculation is not a guarantee of acoustic transferability.

Source plate 396: A1–A12 → plate 1 EI; B1–B12 → plate 2 max diversity; C1–C12 → plate 3 posterior mean; D1–D12 → plate 4 random. Each well starts with 50 µL. No source well is shared between destination plates and no refill is assumed.

Columns 1–12: 256, 75, 269, 294, 275, 41, 158, 145, 200, 196, 291, 264. Species 296 is excluded. Each species needs 200 µL working suspension in total, split across its four wells. Preparation totals exclude pipetting losses.

OD inputs are readings at 1:10 dilution, multiplied by ten. Dilute toward OD600 1.2; weaker cultures remain undiluted. Dose remains 0.100 OD600·µL per species, rounded to 25 nL. No minimum transfer bound, 15 µL rule or automatic dose switch. Destination prefill is 50 µL medium plus inoculum. Review per-source remainder and the instrument survey before execution.

Four plates, 102 conditions each, three technical wells per condition, no culture wells on the edge. Plate 2 contains 94 new max-diversity conditions and eight validation repeats (24 wells), tagged separately from acquisition picks. Plates 1, 3 and 4 each contain 102 new conditions. The two remaining interior wells and the edge receive no bacterial transfer. No new standards are assigned by these picklists.

Validation communities: 294+196+291+264; 294+291+264; 294+275. Monocultures: 294, 275, 196, 291, 264. All are repeated on plate 2 and use source row B. Triplicates are dispersed across the interior like the other conditions; plate-2_validation_wells.csv identifies their exact wells. Exclude the labelled validation conditions from summaries of newly acquired max-diversity communities.

The original September 4 four-member candidate had complete transfers in G16/L9, but B23 lacked 196. The 196 monoculture also had one zero-inoculum well (B22); its other two wells were complete. The 294+275 pair and both constituent monocultures had full recorded transfers. The accompanying prior-assembly CSV retains those individual wells and raw readouts. The repeat tests the community-versus-constituent signal comparison; unchanged per-member dosing means mixtures still receive more total starting biomass than monocultures, so this alone is not a synergy test.

Corrected history is reconstructed from all three rounds' full reports. Zero-volume species are removed; measured wells are grouped only with the same recorded inoculum and dose. Empty inocula remain in the history well table but are excluded from fitting and the done set. Replicate counts of one or two remain explicit. No filtering by response or apparent repeatability. Ordinary absorbance minus the original full-plate median remains the response, not calibrated lactate concentration.

EI and posterior-mean models use their own corrected path plus shared max-diversity data. Random outcomes never enter either fit. The existing presence-only GP is retained; biomass and growth variation remain unmodeled. Intended compositions never actually measured are eligible again unless another path measured them.

The September 9 report header is September 10, while detail records say September 12. Use the complete design-matched report as the run record for this requested corrected retraining; retain the inconsistent timestamps unchanged. No additional unrecorded transfers assumed.

Historical files are unchanged. This folder contains the exact fitted data, delivered-composition done registry, model parameters, designs, source map, and input/code hashes. Missing ODs mean no Echo picklist CSVs are emitted. Enter ODs in the appended 20260916 block in adaptive_run.py, then run acquisition_round_20260916.py to create a new snapshot.
