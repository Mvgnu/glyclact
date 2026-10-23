# Lactate assay workflow

Actual routine clarified by Magnus on 23 October 2026. Reagent preparation is
retained in the [original recipe reference](archive/20261023_before_operator_correction/lactate_assay_planned.md).
Proposed validation steps in that reference are not evidence they were performed.

## Procedure

1. Spin down the bacterial cultures with the centrifuge set to **0 °C**.
2. Transfer **5 µL supernatant + 10 µL crude-glycerol screen medium** into each
   ordinary assay well. No filtration or additional cleanup.
3. Add **15 µL lactate assay reagent mix** (the "screen stuff").
4. Incubate **1 hour protected from light**.
5. Add **15 µL 1 M acetic acid** to stop the reaction. Magnus explicitly
   reconfirmed this step on October 23.
6. Measure absorbance using the reader method used for that run.

Final assay volume is **45 µL**. The 15-µL sample aliquot contains one-third
supernatant; the final stopped well contains one-ninth supernatant. The added
10 µL is crude-glycerol medium, **not water**. It dilutes culture-derived changes
but also introduces fresh medium into the assay.

## Reader wavelength

The protocol and existing analysis use **490 nm (A490)**. Magnus's latest
recollection was "400 nm or sum". The October 2/3 CSV headers identify
`MAGNUS_ABS`, one chromatic channel and the read time, but do **not** state a
numerical wavelength. Verify the actual instrument method/export before treating
490 nm as independently confirmed; do not change it to 400 nm from that remark.
`SPR` is a separate historical reader-method label whose exact setting still
needs verification. Do not infer its meaning from the label alone.

## Assay reagent recipe on record

The high-LDH working recipe listed previously is:

| Component | Stock | Volume |
|---|---|---:|
| Tris | 0.2 M, pH 8.2 | 5,440 µL |
| β-NAD⁺ | 22 mg/mL | 640 µL |
| INT | 10 mg/mL in 67% methanol | 320 µL |
| L1254 LDH | 2.75 U/µL | 7.68 µL |
| Methoxy-PMS | 25 mg/mL in DMSO | 1.41 µL |

Those additions sum to **6,409.09 µL**, described previously as approximately
6.4 mL. Keep this recipe available for auditing stock concentrations, enzyme
activity and preparation arithmetic. It is not a measurement of each historical
reagent batch's composition or activity. Stock identities and original preparation
instructions remain in the archived reference, Parts A–G. The published method
reference is https://doi.org/10.1371/journal.pone.0271818.

## Standards and spikes vary by run

Use the saved assay map and operator corrections, not a universal plate layout.
Not every earlier plate had standards. Plate edges are not automatically blanks.

For the October 2 assay, the reported spike recipe was **7.5 µL of 12 mM lactate
stock + 5 µL medium + 2.5 µL undiluted supernatant**, followed by the same reagent
and acid volumes. It adds **6 mM before reagents / 2 mM in the final well**, but
contains half the supernatant of ordinary wells. Raw spiked-minus-unspiked
absorbance is therefore not a matched recovery measurement. The reply about
water/crude stock preparations did not assign a stock matrix to every spike well.

The October 21 diagnostic design proposes matched spikes with 5 µL supernatant,
3.75 µL of 12 mM lactate in matched medium, and 6.25 µL matched medium blank.
That is a **planned correction**, not a retrospective description of October 2.

Some runs had prolonged handling/read delays. September 19 was reported as
approximately **20 hours after acid stop**. A later comment described a couple
of hours on the bench during pipetting without specifying which preparation or
reaction stage. Do not turn either statement into a measured reaction duration.

## Keep measurement and interpretation separate

Save the raw values and the actual read/stop times when available. Re-reading
one physical plate is not a new culture replicate. Analytical duplicates are
two aliquots of one culture; culture replicates are separate wells.

A rise in a spike well does not by itself validate the standard curve. Where
standards are flat, do not report fitted lactate concentrations or percent
recovery as if calibration worked. No washing, sample extraction, filtration,
deproteinisation, pH adjustment, bubble removal or extra clarification step should
be inserted into the performed procedure without an operator record.
