# Genomic modification opportunities — a running reference

Not a plan. A place to record, as they surface, the strain edits the analysis keeps
pointing at — so that if strain engineering is ever revisited alongside or after
community design, the reasoning does not have to be reconstructed.

Everything here is derived from annotation plus FBA, so each entry states what the
evidence actually is and what would have to be measured before acting. Ranked by expected
effect on lactate yield.

---

## 1. `alsSD` — acetoin / 2,3-butanediol. **Best-evidenced single edit.**

**What.** Delete `alsS` (acetolactate synthase) or `alsD` (acetolactate decarboxylase).

**Why.** Acetolactate synthase consumes **two** pyruvate and regenerates **no** NAD⁺, so
it competes for pyruvate without doing the redox job lactate does. In *Bacillus* the
dominant aerobic overflow route is acetate, acetoin and 2,3-butanediol — not lactate.

**Evidence.** Inactivating `alsS` in *B. subtilis* raised L-lactate yield by **15%**, and
that strain reached **95% of theoretical yield** at 105 g/L
([Romero et al., Microb Cell Fact 8:23](https://pmc.ncbi.nlm.nih.gov/articles/PMC2680810/)).
In our pool, `alsS` **and** `alsD` are annotated in `P1_B3`, `P1_B11` and `P3_F11`.

**Caution — the in-silico test is misleading here.** Blocking `ACLS` is *lethal* for
`P1_B11` and `P3_F11` in the models, because CarveMe merged the catabolic AlsS with the
anabolic IlvBH acetolactate synthase, so deleting it also deletes branched-chain amino
acid biosynthesis. Real *B. subtilis* `alsS` mutants are viable. **Target the catabolic
copy specifically**, or delete `alsD` instead, which has no biosynthetic role.

**Measure first.** Acetoin and 2,3-BD in the baseline runs. If the strains make little of
either on glycerol, this edit buys nothing.

---

## 2. Ethanol route — `adhE` / alcohol dehydrogenase

**What.** Delete the acetaldehyde/alcohol dehydrogenase.

**Why.** Ethanol is the largest competing sink in the models: **~17 mmol/gDW/h of
capacity** sits alongside the lactate in `P1_B3`, `P1_B11` and `P3_C11`, against a
lactate capacity of ~17. And it is genuinely *forced* under under-aeration — ethanol
minimum rises above zero once O₂/glycerol falls below ~0.4, so it is the specific route
by which carbon leaks when aeration is short.

**Strain note.** `P3_F11` already has by far the lowest ethanol capacity in the pool
(**4.98** against ~17 for the others). If ethanol proves to be the main loss at the
bench, `P3_F11` is the strain that needs this edit least — which is an argument for
choosing it rather than engineering the others.

**Caution.** `adhE` is not annotated in `P1_B3`, `P1_B11` or `P3_F11`; the ethanol
capacity runs through other aldehyde/alcohol dehydrogenases (`P1_B3` has eleven
aldehyde-dehydrogenase family proteins). So this is not a clean single deletion in these
strains, unlike in *E. coli*.

---

## 3. Acetate route — `pta` / `ackA`, and `poxB`/`spxB`

**What.** Delete phosphotransacetylase / acetate kinase, and/or pyruvate oxidase.

**Why.** Two separate drains. `pta`–`ackA` takes acetyl-CoA to acetate; pyruvate oxidase
takes pyruvate *directly* to acetate + CO₂, bypassing the lactate node entirely and, in
the `spxB` form, producing H₂O₂ as well.

**Annotated in the pool.** `pta` and `ackA` in `P1_B3`, `P1_B11`, `P3_F11`, `P2_A5`.
Pyruvate oxidase in `P1_B3` (2 copies), `P1_B11`, `P3_F11`, and **`P2_A5` (5 copies)**.

**Precedent.** The 85%-of-theoretical *E. coli* glycerol→D-lactate strain deleted `pta`
alongside `pflB`, `adhE`, `frdA` and `dld`
([Mazumdar et al., AEM 2010](https://journals.asm.org/doi/10.1128/aem.00664-10)).

**`P2_A5` is the standout case**: five pyruvate-oxidase copies is the largest pyruvate
drain in the pool, in the strain that is otherwise the most interesting producer.

---

## 4. `nox` — H₂O-forming NADH oxidase. **Specific to the two LAB.**

**What.** Delete the water-forming NADH oxidase.

**Why.** It reoxidises NADH straight to water, doing LDH's job without making product.
Under aeration it is a direct competitor for the NADH that lactate production depends on.

**Where.** `nox` (K00359/K17869) is **absent from every Firmicute producer in the pool** —
a favourable finding recorded earlier — and present only in `P2_A5` and `P1_F12`, the two
lactic acid bacteria. So this is a targeted edit for `P2_A5`, not a general one.

**Precedent.** Controlled NADH-oxidase expression is the classic cofactor-engineering
lever in *L. lactis*, used in the opposite direction to *reduce* lactate
([J Bacteriol 180:3804](https://journals.asm.org/doi/10.1128/jb.180.15.3804-3808.1998)).
Running it backwards — deletion to *raise* lactate — is the obvious inverse.

---

## 5. `ldhA` + `lldD` in the four enterobacterial supporters. **Two deletions each, and it recovers the whole economic case for a community.**

**Why this is the highest-value entry here.** Optical purity forced the four best
amino-acid providers out of the pool (`POOL_V4.md`): total supporter AA-exporter count fell
from **36 to 11** and the peak from 9 to 6, because the D-LDH risk sits precisely in the
Enterobacteriaceae, which is exactly where the provisioning lives. That cost the
precursor-supply hypothesis (H5) its strain-based test — the one worth up to 38% of
production cost as a yeast-extract substitute. **Two deletions per strain reverses it.**

| Supporter | D-LDH to delete | L-lactate oxidase to delete | AA exporters recovered |
|---|---|---|---:|
| `P2_G9` (unplaced *Leclercia*-like) | `ldhA` KEHDBFFF_02799 | `lldD` KEHDBFFF_01360 | **9** |
| `P3_F8` (*Ewingella americana*) | `ldhA` PALJNNCD_04400 | `lldD` PALJNNCD_03016 | **7** |
| `P2_D5` (*Pantoea pleuroti*) | `ldhA` NAIGPMMC_01317 | `lldD` NAIGPMMC_04245 | 6 |
| `P2_D12` (*Rahnella inusitata*) | `ldhA` JOLBCOOK_00601 | `lldD` JOLBCOOK_02571 | 6 |

**Single copies, no redundancy, both dispensable.** Fermentative D-lactate production and
respiratory lactate oxidation are non-essential. No `dld` in any of the four, so D-lactate
consumption is not a separate target. Removing `ldhA` makes them purity-safe; removing
`lldD` stops them eating the product — so the pair converts each from "purity risk that
also consumes" to clean in the strict sense.

**And these are the most tractable hosts in the pool.** λ-Red recombineering or CRISPR-Cas9
double deletion in an enterobacterium is routine, where editing *Priestia*, *Paenibacillus*
or *Chishuiella* is not. The strains that need engineering are the ones easiest to engineer,
which inverts the usual difficulty ordering.

**Two caveats.** `P2_G9` and `P3_F8` carry a genuine L-LDH as well (real `Ldh_1_N`/`Ldh_1_C`
domains, not an MDH mis-call), so post-deletion they would still make **L**-lactate — on-spec,
and arguably a bonus rather than a problem. And `P3_F8` still needs the *Ewingella*
biosafety conversation regardless of any edit (`LITERATURE_POOL_V3.md`).

**Precedent.** `dld` deletion was part of the 85%-of-theoretical *E. coli* glycerol→lactate
strain, for exactly this reason
([Mazumdar et al.](https://journals.asm.org/doi/10.1128/aem.00664-10)).

---

## 5b. Background: L- and D-LDH are different protein families, and the glycerol route does not set the isomer

Two things worth having written down, because both bear on how much the annotation-based
isomer calls can be trusted and on why the published *E. coli* work produced the wrong
isomer.

**The two lactate dehydrogenases are non-homologous.** They are convergent enzymes from
entirely different folds, so an L/D call rests on domain architecture rather than on a gene
name:

| | EC | KO | PFAM domains | Superfamily |
|---|---|---|---|---|
| **L-LDH** | 1.1.1.27 | K00016 | `Ldh_1_N`, `Ldh_1_C` | LDH/MDH, Rossmann fold |
| **D-LDH** | 1.1.1.28 | K03778 | `2-Hacid_dh`, `2-Hacid_dh_C` | D-isomer specific 2-hydroxyacid dehydrogenase |

Verified across all three producers and all four D-carrying supporters. So the purity calls
in `results/purity_scan.tsv` are about as robust as annotation gets.

**The confusable pair is not L versus D — it is L-LDH versus malate dehydrogenase.** Both
sit in the LDH/MDH superfamily, both Rossmann fold, EC 1.1.1.27 against 1.1.1.37. That is
the ambiguity that made `P3_G8` and `P2_C10` unreliable. The three producers are clean on
it: `Ldh_1_C,Ldh_1_N`, EC 1.1.1.27 **only**, and the descriptor "LDH **family**" rather
than merely "superfamily".

**The glycerol dissimilation route does not determine the isomer.** `GlpK`-`GlpD` versus
`GldA`-`DhaK` decides where glycerol's surplus electrons land — quinone pool versus NADH —
and therefore the oxygen requirement (`OVERFLOW_TRIGGERS.md` H1). The isomer is set purely
by which LDH reduces pyruvate at the end.

That resolves an apparent contradiction with the literature. *E. coli*'s only fermentative
LDH is **`ldhA`, EC 1.1.1.28, D-specific** — it has no L-LDH at all. So Mazumdar et al. got
D-lactate because that is the only enzyme the host has, and the companion L-lactate paper
had to introduce a heterologous L-LDH. **Our producers run the same favourable upper route
(`GLYK` → `G3PD4`, menaquinone-linked, carrying all the flux in the model solution) but
terminate in a genuine L-LDH — so they get *E. coli*'s redox advantage with the L isomer,
with no engineering.**

---

## 6. Upregulation targets, rather than deletions

- **`glpK` / `glpD`.** Overexpressing the *respiratory* glycerol route gave a **~5-fold**
  lactate improvement in *E. coli*, because GlpD puts glycerol's surplus electrons on the
  quinone pool and leaves the NAD pool free for LDH. **Our Firmicutes already use this
  route natively** (`GLYK` → `G3PD4`, menaquinone-linked, carrying all the flux in the
  model solution), so the *qualitative* benefit is already banked — but glycerol uptake
  rate is the classic bottleneck on this substrate, so raising `glpK`/`glpF` expression
  is a rate lever even where the routing is already right.
- **`ldh` copy number / promoter strength.** `P2_A5` carries six annotated L-LDH copies
  against one in `P1_B3` — a natural experiment in whether copy number matters here.
- **Lactate export.** `lctP` is annotated in `P1_B3` (2), `P1_B11`, `P3_F11` and `P2_A5`,
  and export has been the binding constraint in **eight** separate findings in this
  project. Export is also electrogenic H⁺ symport whose driving force collapses at high
  titre and low pH, so overexpression is a plausible titre lever late in a fed-batch.

---

## 7. A note on what *not* to touch

Pyruvate dehydrogenase and the TCA cycle are the obvious "competing" routes, and both are
effectively off-limits: PDH deletion is lethal or near-lethal in most aerobes, and the TCA
cycle supplies biosynthetic precursors. The literature route to suppressing PDH is
**nutritional, not genetic** — withhold lipoate, PDH's essential cofactor. That is already
in the pool design as a medium variable (see `POOL_V3.md`), and it is reversible, which a
deletion is not.

---

## 8. PHB storage in `P1_B3` — a sink no model in the pool can see

**What.** `phaC` / the polyhydroxybutyrate polymerase, if the strain proves to accumulate it.

**Why.** *Priestia* species are documented PHB producers, with accumulation strongly
dependent on carbon source, and *P. aryabhattai* specifically has served as the
PHB-producing member of a synthetic consortium. PHB is intracellular carbon and energy
storage — a direct competitor for the carbon we want in lactate.

**The reason it is invisible.** Checked directly: **no PHB polymer species and no
3-hydroxybutyryl-CoA sink exists in any of the twelve models.** So `P1_B3`'s 0.850 carbon
yield — the highest in the pool — is computed in a model structurally blind to a drain its
species is known for.

**Measure before considering any edit.** If bench yield undershoots the model for `P1_B3`
specifically, PHB is the first thing to look for. See `LITERATURE_POOL_V3.md`.

---

## Ranking, revised after the pool v3 literature pass

0. **`ldhA` + `lldD` in `P2_G9`, `P3_F8`, `P2_D5`, `P2_D12`** — new top entry. Two routine
   deletions per strain, in the pool's most tractable hosts, recovering the four best
   amino-acid providers (9, 7, 6, 6 exporters) that the purity constraint had forced out.
   It is the only edit here that restores an *economic* argument rather than a yield
   increment, and unlike the rest it does not depend on a hypothesis being right.
1. **`poxB`/`poxF` in `P2_A5`** — promoted to the top. *L. plantarum* **consumes its own
   lactate aerobically** via reverse LDH → pyruvate oxidase → acetate kinase, and "acetate
   is the major fermentation end product... produced at the expense of lactate" under
   aeration. In a pyruvate-oxidase-null mutant, **lactate remains the major end product**
   ([J Bacteriol 186:3749](https://jb.asm.org/content/186/12/3749)). Five annotated copies,
   the pool's largest pyruvate drain, in the strain that is otherwise the most interesting
   producer. **Non-genetic alternative first**: pyruvate oxidase is repressed by glucose,
   so the small sugar co-feed proposed under hypothesis H3 suppresses this route *and*
   raises FBP — one intervention, two mechanisms, no strain engineering.
2. `alsD` (acetoin) — best evidence for the *Bacillus* producers, +15% precedent, no
   biosynthetic role
3. `lldD`/`dld` in a chosen supporter — creates a partner that feeds without eating
4. `nox` in `P2_A5` — targeted, strong precedent in the reverse direction
5. `pta`/`ackA` — precedent from the *E. coli* strain, but a smaller drain here
6. Ethanol route — largest capacity, but not a clean single deletion in these strains
7. `phaC` in `P1_B3` — only if PHB is measured and turns out to matter

**The gate on all of it stays the same**: measure the baseline product spectrum first.
Every entry above is a *capacity* argument, and capacity is not a prediction that the
carbon actually goes there.
