# SOP — Community Assembly and Crude-Glycerol Fermentation Screen

384-well fermentation screen. 14-member panel, presence/absence, one fixed medium.
Community assembly by Echo, all bulk liquid by VIAFLO.

Produces cleared supernatant for `SOP — Colorimetric L-Lactate Assay for Bacterial
Supernatants.md`. Design rationale: `results/planning/EXPERIMENT_PLAN_V5.md`. Panel:
`pool_v5_panel.yaml`. Physical strain retrieval: `Strain Identity and Lab Mapping.md`.

---

## 1. Instrument assignment

| Step | Instrument |
|---|---|
| member combinations per well | Echo |
| medium into culture plate | VIAFLO |
| supernatant 384 → 384 | 384-compatible multichannel |
| assay buffer, stop solution | VIAFLO |

Crude glycerol is not Echo-transferable (viscosity, opacity). It is a constant across the
plate, so it goes in the bulk medium. Every well receives **48 µL** bulk screen medium plus its
variable cell inoculum (§5.2), with a hard maximum of **50 µL total per 384-well well**.

---

## 2. Media

### 2.1 Preculture — Nutrient Media

```
10 g/L  yeast extract
10 g/L  soytone
        water to 980 mL
        adjust pH 7.0
        fill to 1 L
        autoclave
```

### 2.2 Screen medium — M9-K raised to 100 mM phosphate + crude glycerol

For this fermentation screen, modify that M9-K further by raising phosphate to
**100 mM total PO4**, then add crude glycerol. M9-K supplies N as NH4Cl.

Per 1 L final screen medium:

```
 2.5 g   NaCl                               (~43 mM)
 4.0 g   KCl                                (~54 mM)
 5.0 g   NH4Cl                              (~93 mM)
100 mmol PO4 total, as sodium dihydrogen phosphate
           = 10 mmol/L already in M9-K + 90 mmol/L screen supplement
 60 mL   50 g/L crude glycerol stock  (see 2.3) -> 3 g/L glycerol
  2 mL   1 M MgSO4
  1 mL   1000x Teknova trace element mix
100 µL   1 M CaCl2                          -- add last, prevents precipitation
          de-ionized water to 1 L
```

Dissolve the M9-K salts plus the additional 90 mmol/L sodium dihydrogen phosphate in less than
the final volume, adjust to **pH 7.0 with NaOH**, and autoclave. After cooling, add the sterile
glycerol, trace element, MgSO4 and CaCl2 stocks aseptically; bring to 1 L. Add CaCl2 last.
Measure and record the final pH.

This is **not unmodified 10 mM-phosphate M9-K and not standard M9**: it uses the M9-K salt
formulation with phosphate deliberately raised from 10 to 100 mM for the screen.

| Screen-medium exclusions | |
|---|---|
| yeast extract, soytone | complex nitrogen sources are excluded from the screen |
| NiSO4 | required only by *Ca*, not in this panel. 4 mg/L is ~10³× the trace-mix background and is the LarA racemase cofactor — no panel member carries `lar`. **Exception: add 4 mg/L to the `P1_B3`/`P1_B8` purity wells**, or nickel starvation returns a false-clean chiral result |
| MnCl2 | required only by *Lp*, not in this panel |

Make one batch for the whole campaign. Aliquot into five sterile bottles, one per culture-plate run.
Store 4 °C.

### 2.3 Glycerol concentration and pH

Use **3 g/L crude-glycerol equivalent**, corresponding to a maximum theoretical lactate titre
of 32.6 mM at one mole lactate per mole glycerol. Start at pH 7.0. The screen-specific phosphate
increase to 100 mM is intended to prevent lactate production from causing an early pH crash;
record final pH rather than assuming that the starting pH is maintained.

### 2.4 Buffer capacity in a 50 µL culture

At 100 mM total phosphate, a 48–50 µL culture contains **4.8–5.0 µmol total phosphate**. At
pH 7.0, approximately **1.83–1.90 µmol** is present as HPO4²⁻ and available as the immediate
proton-accepting reserve. Assuming one acid equivalent is generated per lactate accumulated,
the phosphate-only estimate for a 50 µL well is:

| Approximate pH | Lactate concentration | Lactate in 50 µL |
|---:|---:|---:|
| 6.5 | 21.8 mM | 1.09 µmol |
| 6.0 | 32.3 mM | 1.62 µmol |
| ~5 / outside useful phosphate-buffer range | ~37.5 mM | ~1.88 µmol |

The 48 µL screen-medium dose contains 0.144 mg glycerol (1.56 µmol), so its **maximum theoretical
lactate is 1.56 µmol per well (32.6 mM)**. In the actual 48–50 µL culture, the 100 mM phosphate
has approximately **1.55–1.62 µmol acid capacity before pH 6** and approximately 1.80–1.88 µmol
before pH ~5. Therefore the buffer is sized to support essentially complete theoretical
glycerol-to-lactate conversion without an earlier phosphate-buffer crash. This calculation
assumes one acid equivalent per lactate; final pH still must be measured because ammonium uptake
and pathway proton balance affect net acidification.

---

## 3. Crude glycerol stock

One lot for the entire campaign. Record the lot number.

Glycerol content (% w/w) from the supplier CoA. If unavailable, run the stock once on the
enzymatic glycerol assay (§9). No other characterisation.

For 200 mL of 50 g/L glycerol-equivalent stock:

1. Warm the jar to 50–60 °C until it pours.
2. Weigh — do not pipette:
   `mass = (50 g/L × 0.2 L) ÷ glycerol fraction` → **12.5 g** at 80% w/w.
3. Add water to 200 mL. Stir warm until homogeneous.
4. If pH > 8, adjust to 7.0 with HCl. Record volume used.
5. Autoclave 121 °C, 20 min. Dilute enough not to foam or darken.
6. Store 4 °C, dark. Shake before use.

---

## 4. Precultures

All at 30 °C in Nutrient Media (§2.1).

| Members | Notes |
|---|---|
| `P1_B11`, `P1_F7`, `P3_F11`, `P2_D9` | 12–18 h |
| `P3_A9`, `P3_A8`, `P3_C5`, `P3_G9` | 18–24 h |
| `P2_C2`, `P3_D2`, `P2_C10` | 24–48 h |
| `P3_G8` | 24–48 h; clumps, see §5.3 |
| `P2_B6` | 12–18 h |
| `P3_C3` | + 5 g/L glucose, static, no shaking |
| `P1_B3`, `P1_B8` | purity wells only |

Harvest all strains in **early stationary**, at a fixed hours-after-inoculation per strain set
in the pilot (§11) and not changed thereafter.

Staggered schedule, back from harvest at T0:

```
T−48 h   inoculate P2_C2, P3_D2, P2_C10, P3_G8
T−24 h   inoculate P3_A9, P3_A8, P3_C5, P3_G9
T−16 h   inoculate P1_B11, P1_F7, P3_F11, P2_D9, P2_B6, P3_C3
T0       harvest all
```

Record harvest OD600 per strain.

---

## 5. Inoculum

### 5.1 Wash and measure

Carrier = screen medium (§2.2) **without crude glycerol**.

Per strain, from 5–10 mL preculture:

1. Pellet 4,000 × g, 10 min.
2. Discard supernatant completely. Resuspend in 1 volume carrier.
3. Repeat once (**two washes total**).
4. Resuspend the final pellet in a convenient volume of carrier and measure OD600 (dilute an
   aliquot into the reader's linear range if necessary). Record the undiluted suspension OD.

The Echo transfer volume is calculated from each suspension's measured OD (5.2).
Only dilute or concentrate a suspension if its calculated transfers fall outside the validated
Echo range.

Two washes are required because unwashed `P3_C3` would introduce preculture-derived lactate
correlated with a design variable. Confirm adequate removal with the pilot's t = 0 samples
(§11.1).

### 5.2 Dose

Normalise **OD-equivalent bacterial input, not inoculum volume**. Use a total biomass budget of
**2.00 OD600·µL per culture well**, equivalent to a starting OD600 contribution of approximately
**0.040 in a 50 µL well**. There is no single inoculum volume: it depends on the measured OD of
the source suspensions.

For a community containing *k* members, give every member an equal share of that biomass
budget:

```
transfer per member = 2000 / (k × measured source OD600)    nL
```

If all source suspensions had the same OD, the total inoculum would be 2.00 µL at OD600 1.0,
1.00 µL at OD600 2.0, or 0.67 µL at OD600 3.0. With unequal source ODs, calculate every member
separately and sum the executed transfers to obtain that well's actual total inoculum volume.

The summed inoculum must be **≤2.00 µL per well**, because every well already contains 48 µL
screen medium and the 384-well culture limit is 50 µL. If a calculated community exceeds this
limit, concentrate the relevant source suspension(s) and regenerate the picklist; do not reduce
the biomass target and do not overfill the well.

Calculate the transfer independently for every strain from its washed-suspension OD. Round to
the nearest valid droplet multiple (nominally 2.5 nL with the aqueous class; confirm on the
instrument) and record both the calculated and executed volumes. This keeps total starting
biomass approximately constant across communities without requiring equal inoculum volumes
across wells. OD600 is only an approximate cross-strain biomass proxy, especially for the
aggregating strains in §5.3; consistency of measurement and mixing is more important than
forcing identical suspension concentrations.

If a calculated transfer is below the Echo minimum, dilute that source suspension. If it is
above the validated maximum, concentrate it or split the dose into supported transfers.

### 5.3 Aggregating strains

`P3_G8` clumps and floats. `P2_C2`, `P3_D2`, `P2_C10` may run filamentous.

1. Vortex hard immediately before loading the source plate.
2. If clumps persist: bath-sonicate 30 s. Not probe.
3. No surfactant.
4. Re-mix source plate between destination plates.

---

## 6. Echo source plate

- Use the lab-designated, instrument-compatible Echo source plate; this SOP does not specify its
  format.
- Load each strain according to the validated operating volume for that source plate and the
  picklist demand plus dead volume.
- Spin 200 × g, 1 min before loading — surface bubbles fail the survey.
- Standard aqueous fluid class.
- Load fresh for every run. Do not store cells in a source plate.

Picklist generated from the design matrix. Save per round: design matrix, well map,
randomisation seed, picklist, Echo transfer report.

> `build_picklist.py` does not exist yet. Required before round 1.

---

## 7. Plate assembly

Use a 384-well culture destination plate. Destination plate is held inverted by the Echo. Do
not Echo into a dry plate and do not exceed **50 µL total volume** in any well.

| Step | Instrument | Volume |
|---|---|---|
| 1 | VIAFLO | **5 µL** screen medium into every used 384-well destination well |
| 2 | Echo | OD-based variable inoculum per well, per picklist (§5.2) |
| 3 | VIAFLO | **43 µL** screen medium |
| 4 | — | seal, optically clear non-breathable film |

Steps 2 → 4 within 30 min.

Do not backfill low-inoculum wells with carrier to equalise volume. Every well receives **48 µL
screen medium plus its variable inoculum (≤2 µL), for ≤50 µL total**. Biomass equivalence is the
controlled quantity; exact total-volume equivalence is not.

### 7.1 Layout — design wells + 8 controls per culture plate

| Wells | |
|---:|---|
| as specified in §7.2 | design conditions, randomised positions |
| 2 | medium blank, no inoculum |
| 2 | glycerol-free control — reference community, carrier medium |
| 4 | reference community, dispersed positions |

Reference community: `P2_C2` + `P3_C3` + `P1_B11` + `P2_B6`. Same on every plate of the
campaign. Fix before round 1.

Do not reserve the edge ring as a moat; disperse controls into edge positions.

### 7.2 Plate budget

| Round | Plates | Conditions | n | Design wells |
|---|---:|---:|---:|---:|
| 0 — pilot / interaction | 1 | 36 | 2 | 72 |
| 1 — random | 1 | 88 | 2 | 176 |
| 2 — EI | 1 | 44 | 2 | 88 |
| 3 — EI | 1 | 44 | 2 | 88 |
| 4 — confirm | 1 | 22 | 4 | 88 |
| | **5** | | | **512** + 40 control = **552** |

`P1_B3` / `P1_B8` purity monocultures occupy 2 wells of the round-1 control block, + Ni (§2.2).

---

## 8. Incubation and reads

- 30 °C, **48 µL screen medium plus variable inoculum, ≤50 µL total per well**, flat-bottom
  384-well culture plate, sealed.
- Static, or one fixed shaking speed/orbit for the entire campaign. Record it.
- Incubate for a fixed **48 h**. Read OD600 through the seal at **T0 and T48 only**.
- **Weigh the plate before sealing and at harvest.**
  `corrected conc = measured × (V_final / V_initial)`

Oxygen exposure is fixed by the chosen 384-well plate, seal, fill volume and shaking condition;
record those specifications and do not vary them across the campaign. Do not assign an oxygen
class from nominal well volume alone.

The pilot (§11) must confirm that 48 h is a suitable endpoint. Once fixed, apply the same
incubation duration to every campaign plate.

---

## 9. Harvest

1. Centrifuge the sealed plate 4,000 × g, 15 min, swing-out rotor.
2. Unseal. With a 384-compatible multichannel, transfer **35 µL clear supernatant** to a fresh
   384-well plate. Set aspiration height explicitly; leave ≥10 µL over the pellet.
3. Assay same day, or seal and freeze −20 °C.

Keep the plate map with the supernatant plate. Verify orientation at every transfer.

### 9.1 Analytes

| Analyte | Coverage | Status |
|---|---|---|
| L-lactate | every well | ready — colorimetric SOP |
| OD600 | every well, T0 and endpoint | ready |
| **glycerol consumed** | every well | **required — yield denominator** |
| **D-lactate** | wells above a titre floor | **required — purity constraint** |
| final pH | 12-well subset, or all if the pilot shows spread | |
| acetoin | every well | tier 2 |
| HPLC spectrum | ~15% subset | tier 3 |
| composition | top/bottom ~20 | tier 4 |

Objective is mol lactate / mol glycerol consumed. Optical purity is a hard constraint in
acquisition, not a term in the objective.

### 9.2 Assay dilution and layout

Expected range 0–32.6 mM against a 0.375–12 mM calibration. **Run every sample at 1:3**
(1 volume supernatant + 2 volumes uninoculated screen medium from the same batch). The
theoretical maximum then enters the assay at 10.9 mM. Re-run below-range wells undiluted.

One 384-well assay plate holds at most **120 culture samples in triplicate**. Each assay well
has a **45 µL final volume** (15 µL sample + 15 µL assay buffer + 15 µL stop solution):

| Wells | |
|---:|---|
| 360 | 120 samples × 3 technical replicates, 1:3 |
| 21 | standards, 7 levels × 3 |
| 3 | medium blank |

Use two assay plates when a culture plate contains more than 120 used wells. Put undiluted
re-reads and other repeats on a separate assay plate.

Everything else per the assay SOP unchanged.

---

## 10. Records per plate

Plate barcode, round, date. Glycerol lot and stock batch. Medium batch. Per-strain harvest OD,
washed-source OD and hours. Calculated and executed inoculum volumes and realised OD-equivalent
input. Echo picklist, transfer report, flagged wells. Design matrix, well map,
randomisation seed. Plate mass before and after. Seal lot. Incubator settings. T0 and endpoint
OD files. Harvest time. Assay plate map and raw A490.

---

## 11. Pilot — must pass before the campaign

### 11.1 Pilot plate

Reference community + three spanning communities, n = 4. Read OD at T0, T24 and T48; do not
run 30-minute kinetics. Read t = 0 lactate on four wells.

Run every pilot culture only in the **screen medium: M9-K raised to 100 mM total phosphate,
pH 7.0, with 3 g/L crude glycerol (§2.2)**. There is no alternative-medium control arm.

Confirms:

- growth on crude glycerol in M9-K. If nothing grows, stop and diagnose the crude-glycerol lot,
  salt load, inoculum recovery and seal.
- whether growth between T24 and T48 is small enough to use the fixed 48 h endpoint; if not,
  add a T72 pilot before fixing the campaign endpoint
- seal integrity (plate mass)
- final pH on 12 wells — this verifies the 100 mM screen buffer under the actual fermentation
- **t = 0 lactate below 0.375 mM in every well** — confirms §5.1 carryover

### 11.2 Assay in this matrix

Per the assay SOP Part J, in screen medium:

- standard curve linearity
- spike–recovery on three pilot supernatants
- **A490 of cleared supernatants with no assay reagent, at 1:3.** If well-to-well spread
  exceeds the A490 increment of the 0.375 mM standard, add per-sample background wells
  (15 µL sample + 15 µL plain 0.2 M Tris pH 8.2 + 15 µL acetic acid). Pair every measured
  sample with its background well and reduce the number of samples per assay plate accordingly.
  Otherwise use §9.2 unchanged.

### 11.3 Echo dose response

For each strain, use its measured washed-source OD to prepare four twofold-spaced
OD-equivalent inputs spanning the routine §5.2 dose range. Transfer into **48 µL Nutrient Media
plus ≤2 µL inoculum (≤50 µL total)** in a 384-well plate, triplicate, and read OD600 at T0, T12
and T24. Confirm recovery and a dose-dependent shift in growth. Flag any strain that does not
show it.

---

## Appendix — one plate

| | |
|---|---|
| T−72 h | streak slow strains |
| T−48 h | inoculate `P2_C2`, `P3_D2`, `P2_C10`, `P3_G8` |
| T−24 h | inoculate Paenibacillus |
| T−16 h | inoculate remaining strains |
| T−2 h | harvest, wash ×2, measure source OD, calculate picklist, load source plate |
| T−1 h | VIAFLO 5 µL into used culture wells |
| T−0.5 h | Echo run |
| T−0.25 h | VIAFLO 43 µL; weigh; seal |
| T0 | read OD600 once; incubate at 30 °C |
| T+48 h | weigh; centrifuge; supernatant transfer; assay or freeze |
