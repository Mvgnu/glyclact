# Readout review: 11 September 2026

Analysis: 13 September. Original files and acquisition models remain unchanged.

## Interpretation

Most culture-well readings are close to the edge reference. Replicate agreement is weak. Reported Echo assembly warnings mean intended composition is not verified. Do not treat a low reading as a confirmed low-producing community or train the next model on this round without resolving assembly status.

The user confirmed standards in column 1: A1 = 12 mM, B1 = 6 mM, C1 = 3 mM, halving each row. All column-1 wells are excluded from the descriptive edge reference. The top standards respond on every plate, but lower levels overlap the background and spatial variation. No true zero-standard well was specified, and standard-versus-sample dilution equivalence is unconfirmed. No concentration conversion, arbitrary outlier deletion, or failed-well imputation is performed.

| Plate | Culture median | Edge reference | ICC(1,1) | Median replicate range |
|---|---:|---:|---:|---:|
| 1 | 0.1070 | 0.1070 | 0.016 | 0.0295 |
| 2 | 0.1010 | 0.1005 | 0.073 | 0.0200 |
| 3 | 0.0980 | 0.0970 | 0.037 | 0.0310 |
| 4 | 0.0940 | 0.0985 | 0.040 | 0.0180 |

ICC compares between-community variation with within-triplicate variation. Negative values are retained.

## Source loading audit

The historical statement that all wells retain at least 26 µL is correct **per destination plate with a new 50 µL fill**. It is not true for all four plates from one initial fill. Seven sources require more than 50 µL across all four plates. Nine exceed the 35 µL usable volume associated with the previous 15 µL reserve. The user confirmed there were no refills. The planning check therefore used the wrong loading assumption for the actual run. The cumulative demand explains why the run could not complete from one fill, even if all wells initially contained exactly 50 µL. No instrument report has been found in the scoped project files, so the exact failed transfers and affected compositions remain unresolved. This does not exclude additional assay or instrument effects.

| Source | Species | Demand across four plates (µL) | Minimum remainder per plate after 50 µL fill |
|---|---:|---:|---:|
| A1 | 256 | 55.125 | 29.300 |
| A2 | 75 | 28.575 | 37.850 |
| A3 | 269 | 55.125 | 28.625 |
| A4 | 294 | 56.475 | 31.775 |
| A5 | 275 | 45.900 | 34.700 |
| A6 | 41 | 26.325 | 38.075 |
| A7 | 158 | 31.275 | 37.850 |
| A8 | 145 | 38.475 | 36.500 |
| A9 | 200 | 56.025 | 32.225 |
| A10 | 196 | 70.800 | 26.000 |
| A11 | 291 | 56.100 | 34.100 |
| A12 | 264 | 61.875 | 27.050 |

## Earlier evidence worth retaining

The September 4 four-member candidate 294 + 196 + 291 + 264 has ordinary readings [0.529, 0.573, 0.777]. The highest constituent-monoculture reading on the same plate is 0.258. This is candidate evidence, not proof of synergy: starting total biomass increases with richness and no independent validation is established.

## Needed to resolve this round

Confirm blank identities, standard/sample assay dilutions and plate execution order. Match Echo survey/transfer exceptions to each requested transfer before deciding which compositions were assembled. Inspect source plate type and liquid-class settings against the actual labware. Do not infer successful assembly from the commanded picklist alone.

## Sources

- `protocols/echo/experiment_setup_output/community_design_20260828_13strains.csv` (SHA-256 `30215be73988f295ee5626d14198613f30ceef44569f4de4b8cf334b9bd161d7`)
- `protocols/echo/experiment_setup_output/20260902_adaptive_12strains/community_design_all_plates_20260902.csv` (SHA-256 `7635751f4c0ca78468d3de37c932c5e8d007e030c0c5c6b99f0289a51524c0e6`)
- `protocols/echo/experiment_setup_output/20260902_adaptive_12strains/source_plate_map_and_dilution_20260902.csv` (SHA-256 `d02738c720c815cc0fbdc5f39ba2e583bd9c8509ac8b667abd5d46348ba8346b`)
- `protocols/echo/readout/Karl_20260830_201548_MAGNUS_ABS-Lactate-assay-rerun.csv` (SHA-256 `8e0499d6013c6f1cd6e738da23d49dbc1d13cdb08e0e8dfc67ed596473968496`)
- `protocols/echo/experiment_setup_output/20260902_adaptive_12strains/echo_picklist_20260902_bayesian_ei-plate-1.csv` (SHA-256 `40f06b2e38ba5348aa5963e1f3f1eebe20d2eff8123ab5a1065738e85854c9f4`)
- `protocols/echo/readout/04-09/Karl_20260904_220017_MAGNUS_ABS-Lactate-assay-rerun.csv` (SHA-256 `7048deb2cdc236eafb6b35938d32e57089ff4f350551b22804b3dec3ece66b99`)
- `protocols/echo/experiment_setup_output/20260902_adaptive_12strains/echo_picklist_20260902_maxdiv_repeat-plate-2.csv` (SHA-256 `5180684b3cb10d176a19b698fc6416010bf9039ca96865a9edfbd499618bc92d`)
- `protocols/echo/readout/04-09/Karl_20260904_221030_MAGNUS_ABS-Lactate-assay-rerun.csv` (SHA-256 `a12e0bb48e18e601dccea7be30f21024905749a523d34b6805d2c12cada2fd09`)
- `protocols/echo/experiment_setup_output/20260902_adaptive_12strains/echo_picklist_20260902_posterior_mean-plate-3.csv` (SHA-256 `81a996bd0aa5808dbbab7808456313021c318966115d4115ac3596f24e6fdd59`)
- `protocols/echo/readout/04-09/Karl_20260904_221716_MAGNUS_ABS-Lactate-assay-rerun.csv` (SHA-256 `5c6e77ff4af77ec2f7e21f8c5bf07bfd35b6c3eceb5dbceb259127f78ecb65d1`)
- `protocols/echo/experiment_setup_output/20260902_adaptive_12strains/echo_picklist_20260902_random-plate-4.csv` (SHA-256 `c3222b26cf534b9870db24bc3495a2224e3d8a2c511bc38ee7dd37611fd6df7e`)
- `protocols/echo/readout/04-09/Karl_20260904_222759_MAGNUS_ABS-Lactate-assay-rerun.csv` (SHA-256 `39b1340e2a30a59daec71dcd3b4588cedc84296c499a0994e1543eef0aabc812`)
- `protocols/echo/experiment_setup_output/20260909_paths/20260909T220918_350704Z/12strains/design.csv` (SHA-256 `f9138e7ae5036babac4a3b79ac9dba4b4e4eee282e85be2e48a962cb30923ab8`)
- `protocols/echo/experiment_setup_output/20260909_paths/20260909T220918_350704Z/12strains/candidates.csv` (SHA-256 `d2516bbd06c930df2014372d94d4fc027a8b9d34bf7e29b003eb5fe9a5e67cb8`)
- `protocols/echo/experiment_setup_output/20260909_paths/20260909T220918_350704Z/12strains/source_demand.csv` (SHA-256 `db98dc5707da36bb62e9c72a3e85b660c253a9985a4211d4bea76a3338dd5fa3`)
- `protocols/echo/experiment_setup_output/20260909_paths/20260909T220918_350704Z/12strains/source_plan.csv` (SHA-256 `9f4d3c188ca086ae3797edff3ed07536e40e6ad25a98f052776f173b37b5f130`)
- `protocols/echo/readout/11-09/Karl_20260911_233338_MAGNUS_ABS-Lactate-assay-rerun.csv` (SHA-256 `0efd8d1f19f43e73a6df95f1f95f3047afec1c9918a0cbdef0cf402fb002af31`)
- `protocols/echo/readout/11-09/Karl_20260911_233743_MAGNUS_ABS-Lactate-assay-rerun.csv` (SHA-256 `e4eb0b1dc996d9c0ed5b6b3d81ab2f9cf8a73ba0f8cb51a7a7c78292fbe992db`)
- `protocols/echo/readout/11-09/Karl_20260911_234146_MAGNUS_ABS-Lactate-assay-rerun.csv` (SHA-256 `30c4f942c7b630e3b63115115feb1f7ab40966322f51e0095125ec4d4519821a`)
- `protocols/echo/readout/11-09/Karl_20260911_234548_MAGNUS_ABS-Lactate-assay-rerun.csv` (SHA-256 `3713bc30450ecdc32f2b6692f994ffd717c6e8ba1996f12ceea96de09b579651`)
- `results/screening/genome_scores.tsv` (SHA-256 `9c7bfd078e46b2aff74af100c2a470599865d4b1cff0482d0dec82d25b7390c8`)
- `results/screening/ldh_adjudication.tsv` (SHA-256 `672a4e783c9d1a42a50a645f2e7d802f1278b47612dfccd00acc12a90df84a15`)
