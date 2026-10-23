# Excluded reconstructions

Kept rather than deleted so the CarveMe compute is not lost, but moved out of
`models/` because `pipeline/comets_simulation/convert_models.py` globs `*.xml` and
would otherwise pull them into the harmonised community metabolite space.

- `P3_A3.xml` — highest producer score in the pool (43.8) but 99.78% Staphylococcus,
  nearest reference *S. aureus* USA300. Not an industrially usable host and a
  biosafety review case.
- `P3_H7.xml` — carves to a reaction set identical to `P3_C11.xml` (2484 reactions,
  zero symmetric difference). `P3_C9` replaces it as the independent Gluconobacter.
- `P3_G9.xml` — cannot use glycerol. `GLYCtex` delivers glycerol to the periplasm and
  `GLYK` phosphorylates it in the cytosol, but no inner-membrane transport connects
  them (`GLYCtpp` and `GLYCt` both absent), and the genome carries no `glpF` or any
  MIP-family aquaglyceroporin. Grows on glucose (1.347) and glycerol-3-phosphate
  (0.703), zero on glycerol. Adding `GLYCtpp` would be an assumption rather than
  curation from genomic evidence.
- `P1_G5.xml` — ANI identifies this as *Brucella pseudogrignonensis* (99.56%, 94%
  aligned, strain-level), an opportunistic human pathogen. The screen had cleared it
  because its genus read as *Ochrobactrum*, which was folded into *Brucella* in 2020
  and so was absent from the risk-genus list. Replaced by `P1_A3`
  (*Lysinibacillus xylanilyticus*) in the same lactate-sink role.

---

## Pool v3 exclusions

The panel was widened from 9 to 19 evaluated strains and cut back to 12. Seven moved
here. All have finished models, species calls and gate rows, so re-adding any is free.

**Correction to the `P3_A3` entry above.** It read "99.78% Staphylococcus, nearest
reference *S. aureus* USA300". That was the artefact of a genus-restricted search: with
the reference set escalated to Staphylococcaceae, ANI places it at **99.50% to
*Mammaliicoccus sciuri*** with 96% aligned — no match anywhere among 68 *Staphylococcus*
species. *M. sciuri* left the genus in 2020, the third post-2020 reclassification to
matter here. "Nearest reference *S. aureus*" simply meant it was the least-bad hit in a
reference set that did not contain the right organism. It is **not** an *S. aureus*.

- `P3_A3.xml` — *Mammaliicoccus sciuri*. Re-evaluated and re-excluded for two reasons
  that are not the original one. (1) Its model **fails the redox gate**: electron
  recovery 1.0909, and its 1.00 mol/mol lactate yield collapses to **0.676** once oxygen
  is held to the stoichiometric 0.5 per glycerol, so the headline number was bought with
  electrons the model does not have. (2) *M. sciuri* is a documented `mecA` reservoir and
  needs a real AMR screen (AMRFinderPlus/CARD) rather than eggNOG name matching, which
  found only collisions. **Worth revisiting**: the strain's genome is genuinely clean —
  one fermentative `ldh`, no `lldD`/`dld`, and no ethanol or acetate exchange in the
  model at all — so if the AMR screen clears it and the leak is repaired it is a strong
  bench candidate.
- `P1_A3.xml` — *Lysinibacillus xylanilyticus*. Pure lactate-sink negative control,
  deprioritised under the wider brief. Also **fails the redox gate** at 1.0244; that leak
  is localisable to `SUCD2`, a succinate dehydrogenase written to reduce *two* electron
  carriers (`fad` and `mqn7`) while oxidising one substrate. Blocking `SUCD2` restores
  recovery to exactly 1.0000.
- `P2_C2.xml` — *Oerskovia merdavium*. **Fails the redox gate** at 1.0714, and unlike
  `P1_A3` the leak does not localise: the per-reaction imbalances cancel, so it runs
  through a metabolite whose formula cannot be scored. Was also dominated as a supporter
  by `P3_G8` (6 amino-acid exporters versus 1, catalase 3 versus 1). Its one advantage is
  that it consumes no lactate — re-add if a non-consuming Actinobacterium is wanted.
- `P1_F12.xml` — *Leuconostoc mesenteroides*. **Cannot grow on glycerol**: µ 0.0000 both
  aerobically and anaerobically against 0.3102 on glucose. No `glpK`, no `glpD`, no
  `gldA`, and no lactate exchange either. Viable only as an obligate cross-feeder, which
  makes a null result likely.
- `P3_C9.xml` — *Gluconobacter cerinus*. Near-duplicate of `P3_C11`, catalase-negative,
  and not a genuine producer: its `LDH_L` is built from `lldD` (EC 1.1.2.3), the
  respiratory L-lactate *oxidase*.
- `P1_G3.xml` — unplaced Erwiniaceae (best hit *Pantoea coffeiphila*, 84.6%, below the
  genus threshold). Functional duplicate of `P2_D5`. Its 1.00 lactate yield is also an
  artifact — `LDH_L` from an unnamed FMN-dependent dehydrogenase with no EC number.
- `P3_A10.xml` — *Paenibacillus taichungensis*. The unannotated-export test case; lower
  priority now that the pool has four confirmed producers with working export.

---

## Reinstated for the pure-D fork

- `P1_G3.xml` — **moved back into `models/`.** Two of its three exclusion reasons were
  specific to the L-lactate pool and do not survive the fork. The "1.00 lactate yield is
  an artifact — `LDH_L` from an unnamed FMN-dependent dehydrogenase" was a spurious *L*-LDH,
  and the isomer-aware screen agrees there is none: it qualifies on `K03778` / EC 1.1.1.28
  with no L-LDH and no racemase. "Functional duplicate of `P2_D5`" also fails on the axis
  that now separates members — `P1_G3` makes acids only, `P2_D5` diverts to 2,3-butanediol.
  The identity problem stands and is unresolved: ANI put it at 84.6% to *Pantoea
  coffeiphila*, below the genus threshold, while MEGAN reads *Enterobacter*. It is one of
  six bench strains with no competing reductive product.
- `MYb535.xml` — **fails the redox gate.** Electron recovery 1.153 against a 1.015
  threshold, driven by an ethanol yield of 1.259 mol per mol glycerol. The ceiling is
  1.167, from the degree of reduction: glycerol carries 14 available electrons and ethanol
  12, so 14/12 = 1.167 mol/mol is the most that can be made even with perfect conservation.
  A model above it is manufacturing electrons, and every yield it reports is unsafe.
  Excluded after the D-lactate export curation, which cleared the only other flagged
  strain (`P1_C12`, 1.167 -> 1.000) rather than masking it.
