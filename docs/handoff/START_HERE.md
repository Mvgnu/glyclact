# Lactate experiment analysis handoff

Prepared **23 October 2026** for an independent analysis of inconsistent
community-screen and lactate-assay results. The question is what the existing
data support, what likely went wrong, and what minimal experiment would separate
the remaining explanations. Do not assume the previous analyses are correct.

## Actual workflow takes precedence over old plans

Magnus clarified the routine on October 23:

1. Precultures grow in nutrient-rich medium for **24 h**.
2. Measure OD and dilute to the run's OD-based input dose. **No cell washing**.
3. Transfer into crude-glycerol screen medium; incubate **48 h**.
4. Spin bacteria down with the centrifuge set to **0 °C**.
5. Ordinary assay: **5 µL supernatant + 10 µL crude-glycerol medium**.
6. Add **15 µL assay reagent mix**; incubate **1 h protected from light**.
7. Add **15 µL 1 M acetic acid** (explicitly reconfirmed); read absorbance.

There is no additional washing, filtration, extraction, deproteinisation or
cleanup. Do not retrospectively add mixing, bubble removal, pH adjustment,
weighing, growth measurements or validation experiments from a planning SOP.
Source dilution with glycerol-free screen medium is recorded in the design
recipes; verify run-specific execution rather than assuming every written detail.

The routine assay is 45 µL after stopping. The original culture-derived sample
fraction is 1/3 before reagents and 1/9 after stop. The added fresh medium brings
its own crude-glycerol matrix. Unwashed inocula can carry nutrient-medium solutes
and metabolites; assess the magnitude from dilution/transfer volumes rather than
declaring carryover the explanation without a quantitative check.

The SOP/plots use **A490**, but the user's recent recollection was "400 nm or sum".
October CSV headers contain method names, not a numerical wavelength. Thus 490 nm
is a recorded protocol assumption, not independently verified instrument metadata.
`SPR` is a historical method label whose settings are unknown. Do not guess them.

## Read these first

Paths are relative to the project root:

- `protocols/SOP — Community Assembly and Crude-Glycerol Fermentation Screen.md`
- `protocols/SOP — Colorimetric L-Lactate Assay for Bacterial Supernatants.md`
- `protocols/archive/20261023_before_operator_correction/`: preserved full media
  and reagent recipes. Its old wash/dose/layout instructions are not execution.
- `protocols/Strain Identity and Lab Mapping.md`: lab IDs to candidate identities.
- `protocols/echo/experiment_setup.ipynb` and `adaptive_run.py`: dated input blocks.

Recorded screen-medium formulation: modified M9-K, 100 mM phosphate, target 3 g/L
crude-glycerol equivalent. Recorded preculture recipe: 10 g/L yeast extract and
10 g/L soytone. Actual batch chemistry, crude-glycerol composition, assay stock
activity, reaction pH and oxygen conditions are not established by these recipes.
The 48-hour endpoint does not establish stationary phase or its production route.

## Data chronology

| Assembly / readout | Primary records | Important context |
|---|---|---|
| Aug 28 / Aug 30 | `protocols/echo/readout/Karl_20260830*`; top-level Aug 28 design/picklist CSVs in `experiment_setup_output` | Initial diversity screen; 13 species; triplicate culture wells; two reader runs after bottom cleaning. Edge/media and standard placement must follow operator notes. |
| Sep 2 / Sep 4 | `readout/04-09/`; `experiment_setup_output/20260902_adaptive_12strains/` | Four acquisition arms, 12 species. No standards reported. Ordinary and SPR reader methods. Major logged transfer failures; use reconstructed inocula. |
| Sep 9 / Sep 11 | `readout/11-09/`; `experiment_setup_output/20260909_paths/` | Four arms, transfer failures, conflicting report timestamps. Operator said assembly occurred after midnight; this does not resolve every internal report timestamp. |
| Sep 16 / Sep 19 | `readout/19-09-20h-delay-read/`; `experiment_setup_output/20260916_corrected_paths/20260916T155952_260216Z/` | Four source rows. Operator reports no warnings. Approximately 20 h from acid stop to read. A1–F1 crude-medium standards; A24–F24 water standards. Do not include these as edge blanks. |
| Sep 30 / Oct 2, exported Oct 3 | `readout/02-10/`; `experiment_setup_output/20260930_analytical_one_plate/20260930T193441_052059Z/` | Fixed diagnostic panel, 96 cultures, 294 mapped assay wells, two reads of the same assay plate. Detailed map below. |
| Oct 21 / results not yet supplied | `experiment_setup_output/20261021_assay_troubleshooting/20261021T132729_361675Z/` | New small diagnostic design, 30 culture wells, 104 planned assay wells. Do not invent results or treat blank templates as measurements. |

All shortened `readout/...` and `experiment_setup_output/...` paths above are
under `protocols/echo/`. Filenames and acquisition dates alone do not identify a
physical plate. Join using reader IDs, picklist naming and transfer fingerprints.

## Echo assembly reconstruction

Raw reports: `protocols/echo/transfer_logs/usb_import_20260914/`.

Start with `work/echo_assembly_audit_20260914/full_report_output/ASSEMBLY_AUDIT.md`
and its tables; verify against `full_reports.py` and the raw reports. The older
`output/` directory is exception-only and superseded. It must not be merged as
independent execution evidence.

Previous full-report analysis found two August run records: one entirely
"Not processed" and one with all 1,971 planned transfers recorded at requested
volume. The successful report, not the exception-only record, supports the
August intended assembly. September reports contain explicit zero actual volumes.
Reconstruct member identity **and delivered dose** per well. A positive reported
volume is not proof of viable growth, and an exception list alone is not proof
that nothing transferred. Avoid counting exceptions and full reports twice.

The prior audit reports poor repeatability even among fully transferred August
triplicates and some composition-matched September groups. Therefore assembly
failure may explain some disagreement but must not be assumed to explain all.

## October diagnostic map and actual deviation

The authoritative intended map is
`protocols/echo/experiment_setup_output/20260930_analytical_one_plate/20260930T193441_052059Z/plate-1_assay_design_20260930.csv`.

- Twelve monocultures: 256, 75, 269, 294, 275, 41, 158, 145, 200, 196, 291, 264.
- C01: **41 + 200 + 196**.
- C02: **256 + 269 + 294 + 41 + 145 + 200**.
- C03: **294 + 275**.
- Six culture wells per condition plus six incubated medium controls = 96.
- Every culture row B/D/F/H/J/L contains all conditions in rotated order,
  columns 2–17. Rows are NOT separated into monocultures and communities.
- Analytical duplicates: B/C, D/E, F/G, H/I, J/K, L/M at the same column.
- Extra spikes/backgrounds cover culture replicate 1 only. They do not directly
  measure background in the high replicate-5/6 cultures.
- Two CSVs: header times 03 Oct 01:10:19 and 01:16:35 (UTC+2). Second ID says
  `plate_1_2ndmeasure_postshake`. These are repeated reads, not independent plates.
- Operator confirmed planned standard concentrations/positions with a short
  "yes" response. Water and crude-medium series were separate. Avoid inferring
  answers to additional parts of compound questions from that single reply.

Actual October 2 spike: **7.5 µL 12 mM lactate + 5 µL medium + 2.5 µL undiluted
supernatant**, plus normal assay mix and stop. That is +6 mM before reagents and
+2 mM in the final 45 µL well. It has half the ordinary sample volume. Stock
matrix assignment per spike well is unresolved. This is NOT matched recovery.

The current clarification confirms the routine 5+10 µL unspiked preparation and
acid stop. Historical analytical outputs still describe these as less certain;
keep their numerical data, but use this later clarification when interpreting.

Existing October outputs and scripts:

- `protocols/echo/readout_analysis/20261002_diagnostic/analysis.json`
- `protocols/echo/readout_analysis/20261002_diagnostic/report.html`
- `protocols/echo/readout_analysis/20261002_diagnostic/visualizations/`
- `protocols/echo/analyze_readout_20261002.py`
- `protocols/echo/visualize_readout_20261002.py`

Recompute rather than anchor on their narrative. Notable raw postshake pairs
to verify: J8/K8 = 0.349/0.359 (C01), J3/K3 = 0.497/0.297 (C02),
L11/M11 = 0.515/0.287 (C01), L6/M6 = 0.090/0.354 (C02).
These are absorbance values, not lactate concentrations. Flat standard series
and different sample fractions block straightforward concentration/recovery claims.

No populated endpoint OD record was supplied for that analysis. Visible pellets
were reported in some wells, not quantified growth. Later OD templates are not
evidence that growth was measured.

## Latest planned troubleshooting run

October 21: monocultures 41, 200, 196 and C01, six cultures each plus six blanks.
196 was restarted from original stock because the personal stock grew poorly.
296 (a different strain) had previously failed to grow repeatedly and was excluded.
OD readings: 41=0.77, 200=0.25, 196=0.10, interpreted using the established 1:10
measurement convention; that dilution was not separately reconfirmed in the reply.

The design retains 0.100 OD600·µL per member with 25 nL rounding and 50 µL
medium plus inoculum. A three-member community is NOT total-biomass-matched to a
monoculture. Do not use apparent superiority as proof of biomass-independent synergy.
New matched spikes, pH fields and OD records are plans until populated.

## Independent measurement option

KU Leuven replied that a free feasibility pilot is possible: lactate standard
first, then a few samples; detection/quantification in this crude matrix is not
yet demonstrated. Their Echo-MS method would measure **total D+L lactate**, not
distinguish enantiomers by ordinary MS/MS. They request at least 50 µL per submitted
well, homogeneous supernatant without bacteria/precipitate, and advised -20 °C
storage/shipping. They request the full medium composition and experimental details.

They quoted €150/day or €500/5 consecutive days, €16/plate, €0.40/seal,
€5/day carrier solvent and €40/hour staff time. These are the supplied email
terms, not an independent current price verification. The exact email date and
sender signature were not included in the pasted reply. No samples were sent.

Agreement with the current L-specific assay is not guaranteed even if both
methods work. Preserve untreated supernatant for a pilot, not completed assay
mixture. Existing 50-µL cultures may not yield 50 µL recoverable sample after
sampling/evaporation; ask about volumes before pooling independent cultures.

## Reproducibility and scope

`PROJECT_FILE_MANIFEST.json` records SHA-256 hashes of the included non-cache
project files. Frozen experiment snapshot manifests are retained and should be
verified. Paths inside old JSON manifests may be absolute; resolve them relative
to the included project, without changing the original records.

The requested repository commit includes available experimental data, genomic
inputs, model files, saved designs/analyses, results and supporting code.
Use the connected repository for Pro; the operator requested no ZIP archive.
External files not present locally cannot be supplied; the USB volume is not
currently mounted. Do not interpret missing raw evidence as a recorded negative.

Never run the entire notebook or `adaptive_run.py` blindly: historical entry
points can regenerate old designs. Run analyses in a copy or new output folder,
retain raw data and do not train new models until label validity is assessed.
