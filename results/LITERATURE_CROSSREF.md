# Literature cross-reference against the ANI-identified species

Once species were pinned by ANI, each was checked against published physiology and
the result compared to what our genome annotation and GEM say. Three of the checks
changed a conclusion. Sources are linked inline.

The purpose is the one the CarveMe/COMETS limitations forced: **amend model defaults
using organism-level knowledge instead of working around annotation gaps.**

## 1. *Lactococcus lactis* (`P3_B12`) — the model has the wrong glycerol route

Documented physiology: glycerol-3-phosphate oxidase (GlpO) uses **O₂ as electron
sink and produces H₂O₂**; *L. lactis* ATCC 19435 grows aerobically but **H₂O₂
accumulation caps its growth rate**, and **no NADH peroxidase or other enzymatic
H₂O₂-reducing activity was detected**
([AEM 68:4350](https://journals.asm.org/doi/full/10.1128/aem.68.9.4350-4356.2002);
[Microbiology 158:2661](https://www.microbiologyresearch.org/docserver/fulltext/micro/158/10/2661_mic061663.pdf?expires=1772496210&id=id&accname=guest&checksum=D5A3CB0C05451C96A09238A6F5246B44)).
In oral streptococci, conserved `glpKOF` supports only low glycerol growth and
glycerol catabolism drives H₂O₂ production
([J Bacteriol](https://journals.asm.org/doi/10.1128/jb.00227-24)).

What our data say:

- The annotation **cannot resolve the route**: `P3_B12`'s `glpD` hit carries **both**
  `K00105` (glpO, H₂O₂-forming) and `K00111` (quinone-linked), described as
  "C-terminal domain of alpha-glycerophosphate oxidase". `P1_B3` and `P3_F11` carry
  the identical dual assignment.
- The **model picked the peroxide-free option**: G3P is consumed by `G3PD2` (NADP)
  and `G3PD4` (menaquinone). No O₂-dependent oxidase, so no H₂O₂ from glycerol.
- `h2o2_c` exists in the model with seven consuming reactions, and the strain has no
  catalase — `K11065` maps to `tpx` (thiol peroxidase), not a true NADH peroxidase,
  matching the literature.

**Consequence.** The GEM overstates `P3_B12`'s aerobic glycerol performance and
structurally cannot represent the peroxide burden that is the actual reason a
catalase-positive partner would help it. Adding a GlpO reaction
(`glyc3p + o2 -> dhap + h2o2`) is the single amendment that makes the pool's
catalase-rescue hypothesis testable.

## 2. *Priestia aryabhattai* (`P1_B3`) — B12 route is the wrong one for an aerobic process

*Priestia megaterium* (its closest relative here at 95.1% ANI) is a long-standing
industrial B12 producer, but it uses the **anaerobic pathway**, and production
requires an anaerobic or microaerobic stage — published processes use two- or
three-stage cultivation with an anaerobic phase to form the precursor
([Microb Cell Fact 13:102](https://pmc.ncbi.nlm.nih.gov/articles/PMC4105875/);
[review](https://pubmed.ncbi.nlm.nih.gov/34263356/)).

Our genomes agree, and split cleanly by route:

| Strain | Species | `cob*` (aerobic) | `cbi*` (anaerobic) | Route |
|---|---|---:|---:|---|
| `P1_B3` | *Priestia aryabhattai* | 8 | **12** | anaerobic |
| `P3_C11` | *Gluconobacter frateurii* | **17** | 2 | aerobic |
| `P3_C9` | *Gluconobacter cerinus* | **17** | 2 | aerobic |
| others | — | 1–8 | 0–1 | incomplete |

**Consequence.** This settles the open question of who supplies B12. `P1_B3` carries
the anaerobic route, which will not run under aeration, so **the Gluconobacters
remain the only plausible B12 source** — and now for a mechanistic reason rather than
a gene count. Note the GEMs cannot test it regardless: no model has an `EX_cbl1_e`.

## 3. *Gluconobacter frateurii* (`P3_C11`) — DHA is a by-product, glyceric acid is the main one

Direct species-level support for the donor role: in *G. frateurii*, a **ΔsldA mutant
produced 89.1 g/L glyceric acid over 4 days with no DHA accumulation**
([Appl Microbiol Biotechnol](https://pubmed.ncbi.nlm.nih.gov/21071844/)), so SldA is
the enzyme responsible for DHA from glycerol — exactly the marker our screen used to
assign the donor role. *G. frateurii* NBRC 103465 is characterised primarily as a
**high glyceric-acid producer** with several membrane-bound dehydrogenases
([Genome Announc](https://ncbi.nlm.nih.gov/pmc/articles/PMC3735062)). Engineered
*G. oxydans* reaches 96 g/L DHA from 100 g/L glycerol
([Bioresour Technol](https://www.sciencedirect.com/science/article/abs/pii/S0960852410009235)).

Our genomes: `sldA`/`sldB` present in both, plus a membrane-bound ADH subunit and six
`glycerate`-annotated proteins in each.

**Consequence, and it is a caution.** The donor assignment is confirmed at species
level, but **wild-type *G. frateurii* directs most glycerol to glyceric acid, with
DHA as a by-product.** Glyceric acid is a carbon sink competing directly with the
shuttle, and it is absent from our models and from the planned product labels. The
DHA shuttle should be expected to leak carbon to glyceric acid in reality even if the
curated model shows clean DHA transfer.

## 4. Species that raised no conflict

- ***Bacillus velezensis*** (`P1_B11`) — well-characterised plant-associated and
  biocontrol organism, widely used industrially; no conflict with its role as the
  low-aerobic-drain producer.
- ***Bacillus altitudinis*** (`P3_F11`), ***Paenibacillus taichungensis*** (`P3_A10`),
  ***Rahnella inusitata*** (`P2_D12`), ***Lysinibacillus xylanilyticus*** (`P1_A3`) —
  environmental organisms, no documented physiology contradicting their assigned roles.

## What this pass changed

1. **A concrete GEM amendment for `P3_B12`** — add GlpO, so glycerol catabolism
   produces H₂O₂ and the catalase-rescue mechanism becomes representable.
2. **The B12 question is settled on route, not gene count** — Gluconobacters only.
3. **A new competing product to track**: glyceric acid from the Gluconobacters.
4. **Confirmation that `sldA` is the right marker** for the DHA donor role, from a
   knockout in the exact species we have.

The general lesson is the one that motivated the pass: for `P3_B12` the annotation
carried both KOs and could not choose, and only organism-level knowledge resolves it.
Gene names alone were also wrong about `noxE` in the same strain, while the species
literature independently supports an NADH oxidase phenotype. **Annotation sets the
hypothesis; the species literature adjudicates it.**
