# Lactate export: why this cannot be a selection criterion

Pool v4 marks most producers `unannotated_gap` for lactate transport, and "carrying an
exporter" was proposed as a pool requirement. This screen was run to fill that gap. It
does not fill it, and the reason is worth stating plainly because the alternative is a
filter that confidently selects the wrong strains.

## There is no characterized bacterial lactate exporter

Swiss-Prot holds **13** manually reviewed bacterial lactate transporters in total. Every
one of them, read from its own UniProt function annotation:

| Gene | Organism | Annotated function |
|---|---|---|
| `lldP`/`lctP` | *E. coli*, *Salmonella* | "**Uptake** of L-lactate across the membrane" |
| `lutP` | *B. subtilis* | "the principal permease for the **uptake** of L-lactate" |
| `lctP` | *B. subtilis*, *S. iniae*, *Halalkalibacterium* | "role in L-lactate transport" / "L-lactate **utilization**" |
| `Sfum_3364` | *Syntrophobacter* | proton-coupled L-lactate transporter |
| `larD` | *L. plantarum* | channel, facilitates **diffusion** of D/L-lactic acid — part of the **lactate racemization operon** |
| `yhjX` | *E. coli* | uncharacterized MFS, pyruvate-sensing regulon |

Not one is a dedicated efflux system. The permease family is an uptake family: `lldP`
exists so a cell can eat lactate for respiration. `larD` is a bidirectional channel
whose characterized role is feeding substrate to a racemase — the opposite of what is
wanted here.

Lactate efflux in lactic acid bacteria proceeds by undissociated lactic acid diffusing
across the membrane, plus an electrogenic lactate/H⁺ symport that has never been
assigned to a gene in *L. lactis*. There is nothing to screen for.

## What the screen actually found

1,687,322 proteins searched against the 13 characterized transporters, plus an HMM search
for PF02652 `Lactate_perm`.

- **241 of 384** genomes carry a PF02652 permease domain.
- 1,532 `lldP`-family hits at 65% median identity — real permeases, and by annotation
  these import.
- 412 `larD` hits at **31% median identity**. `LarD` sits in the formate/nitrite
  transporter family (PF01226) alongside `FocA` (formate channel) and `NirC` (nitrite
  channel), which are near-universal. At 31% identity these are FocA/NirC, not lactate
  channels. The single 100% hit is `P2_A5`, the *L. plantarum* that genuinely has one.
- `yhjX` hits at 27.6% median — MFS background.

## The filter would invert the answer

Applying "must carry an annotated lactate transporter" to the producer bench:

- It **excludes `P3_B12`** (*L. lactis*, zero PF02652 domains) — the most industrially
  characterized homolactic organism there is, which exports lactate at commercial rates.
- It **includes the Bacillaceae** (`P1_A1`, `P1_A9`, `P1_B4`, `P3_B5`, `P3_C8` … each
  with 3 PF02652 domains) whose permeases are, on every characterized homolog, uptake
  systems — i.e. the genomic signal most likely indicates a strain that can *re-consume*
  the product.

So the criterion is not merely uninformative; it is anti-correlated with what is wanted.

## The filter had already been applied, inside the models

This is not hypothetical. The transporter annotation was gating a model curation step,
and it had silently removed the lactate node from a third of the bench.

Checking every candidate reconstruction for `EX_lac__L_e`:

**12 of 36 candidate models carry `LDH_L` and have no L-lactate exchange at all.** They
can make lactate and can never excrete it, so at steady state the reaction is blocked and
the model produces none. Affected: `P2_D9`, `P3_A9`, `P3_A10`, `P3_C5`, `P3_D2`, and
**all seven *Lactococcus*** (`P3_A12`, `P3_A7`, `P3_B1`, `P3_B12`, `P3_C1`, `P3_C3`,
`P3_H11`).

That is why `P3_B12`'s lactate yield columns were blank in every model gate, and why the
provision matrix reports it secreting **0.00** lactate while the other producers report
8.5–9.1 mmol/gDW/h.

`modeling/curate_models.py` exists to fix exactly this — it adds `L_LACt3`, the electrogenic
lactate/H⁺ symport, plus the exchange. But it was gated on:

```python
if evidence["l_ldh"] < 1 or evidence["lactate_transporter"] < 1:
    return changed          # no annotated lctP/lldP -> no export
```

So the amendment was withheld from *L. lactis*, which has zero annotated lactate
permeases and exports lactate at commercial rates, and granted to Bacillaceae carrying
three uptake permeases each. The gate inverted the very thing it was protecting.

**The guard has been removed.** Export is now gated on a confirmed fermentative L-LDH,
taken from the sequence-level calls rather than the annotation. The justification is
mechanistic rather than convenient: an organism that reduces pyruvate to lactate must
move it out, and `L_LACt3` is the mechanism the LAB literature describes.

This also retires the `P3_A10` "unannotated export" test case, which asked whether a
missing transporter annotation predicts a real export gap. The answer is no — the
annotation tracks uptake, so its absence carries no information about export.

## Recommendation

Drop lactate export as a selection criterion. Keep the PF02652 count in the bench table
as a **re-uptake risk** annotation rather than a capability, and read it alongside the
`lldD` count, which is the enzyme that makes re-uptake metabolically worthwhile.

Export is a phenotype. It is measured as titre in the supernatant, and every strain that
goes into culture reports it for free.
