# The pure-D fork: a producer-only guild, screened on genomes rather than genus

Three scripts, run in order:

```
screen_d_guild.py     trait screen over all 384 annotated genomes
virulence_screen.py   DIAMOND vs VFDB core + NCBI AMRFinderPlus, guild members only
d_guild_bench.py      joins the two, ranks the usable bench
```

## What changed in the design

The producer/supporter split is gone. The community landscape falsified it twice over:
supporters never secrete their provisions at maximum growth, and they suppress producer
biomass ~80-fold (0.257 alone to 0.0033 with four supporters) because Enterobacteriaceae
run at mu ~0.55-0.60 on glycerol against 0.117-0.129 for the Bacillus producers. A partner
that eats the feedstock and returns nothing measurable is a cost, not a support.

So every member of this pool is a producer, following LS30 (PMC6387908): Enterobacter
57.97%, Escherichia/Shigella 39.25%, Klebsiella 0.67%, Citrobacter 0.41% -- four glycerol
fermenters, a redundant guild rather than a division of labour. Redundancy also solves the
kinetic problem for free, since a guild drawn from one physiology is growth-rate matched
by construction.

## Attrition

| Stage | Genomes |
|---|---:|
| annotated genomes scanned | 384 |
| carry D-LDH (K03778 / EC 1.1.1.28) | 203 |
| and no L-LDH | 116 |
| and ferment pyruvate (pflB or adhE) | 52 |
| and glycerol-capable, no racemase | 52 |
| excluding GCF reference assemblies | 47 |
| carrying a lactate transporter | 43 |

The bench is 40 lab-plate isolates and 7 from the worm collection. 81 of the 384 genomes
are downloaded reference assemblies rather than freezer strains; they are scanned and
reported but cannot be pooled.

## The fermentation gate is the one that matters

Carrying D-LDH does not make a strain a producer. Lactate is overflow from a fermentative
pyruvate node, and an obligate respirer routes pyruvate through PDH to CO2 instead.

Before this gate the bench was headed by 14 *Pseudomonas*, plus *Macromonas*,
*Chishuiella*, *Agrobacterium*, *Mixta* and *Hydrogenophilus*. **Not one of them carries
pflB or adhE.** Every Enterobacteriaceae genus present carries both:

| genus | n | pflB | adhE |
|---|---:|---:|---:|
| Pseudomonas | 14 | 0 | 0 |
| Citrobacter | 13 | 12 | 11 |
| Enterobacter | 9 | 8 | 8 |
| Cedecea | 5 | 5 | 5 |
| Leclercia | 4 | 4 | 4 |
| Macromonas | 3 | 0 | 0 |
| Buttiauxella | 3 | 3 | 3 |

Applying it reproduces the LS30 membership without being told to: Citrobacter 12,
Enterobacter 8, Cedecea 5, Leclercia 4, Buttiauxella 3, Arsenophonus 3, plus single
Escherichia, Klebsiella, Pantoea and Serratia.

## Virulence, by homology rather than by genus

The previous round rejected 48 of 61 candidates for sitting in "opportunistic pathogen
genera". That is a statement about a genus, and the genus is not what gets cultured.

The ANI-analogue for this question is the standard one: search the predicted proteome
against curated references with DIAMOND and call a gene present above an identity and
coverage threshold, which is what ABRicate does. Two references were used -- **VFDB core
(setA)**, 4,759 experimentally verified virulence proteins, and **NCBI AMRFinderPlus
AMRProt**, 9,998 resistance/stress/virulence proteins -- at >=80% identity over >=80% of
the reference protein.

**45 of 47 carry nothing decisive.** Two are flagged:

| strain | tier | evidence |
|---|---|---|
| P2_E3 (*Escherichia*) | review required | intact K1 capsule locus: kpsC/D/E/F/M/S/U, 7/7 genes at 96-99% identity, full coverage |
| P2_H1 | elevated | `ant(3'')-Ic`, a genuinely acquired aminoglycoside nucleotidyltransferase |

Absent across the whole bench: **no verified exotoxin, no T3SS injectisome, no
carbapenemase, no ESBL, no `mcr` colistin resistance, no `vanA`/`vanB`.**

What the bench does carry, and why none of it counts:

- **T6SS** in 34 strains (T6SS-1, T6SS-2, HSI-1/2/3). An antibacterial weapon carried by
  most environmental Gram-negatives, not a determinant of human disease.
- **Chromosomal class C AmpC matching each genus**: `blaCMY` (22, the *Citrobacter
  freundii* enzyme), `blaACT` (7, the *Enterobacter* one), `blaEC` (*E. coli*), `blaBUT`
  (*Buttiauxella*). Finding these is a statement about species identity, not resistance.
  Likewise `qnrB46` (14) and `fosA2` (14), both native to these genera.
- **Environmental metal and efflux systems**: arsenic, copper/silver, cadmium, chromate.
- Adherence, motility, siderophores and ompA-type proteins, which are ubiquitous.

### Three false positives this screen had to be corrected for

Each was found by checking the evidence behind a call, not by inspection of the code.

1. **`ompA` is not an invasion determinant.** Treating VFDB's whole "Invasion" class as
   decisive flagged 43 of 85 strains, driven by ompA -- the major outer membrane protein of
   every enterobacterium, catalogued by VFDB because it crosses the blood-brain barrier in
   *E. coli* K1. `ibeB`, `ibeC` and `aslA` behave the same way. Fixed by requiring
   locus-level evidence: >=3 distinct genes of a named high-consequence system.
2. **AMRProt doubles as the point-mutation reference set**, so a plain homology search hits
   `gyrA`, `parC`, `rpsL`, `folP` and `murA` in every genome. And AMRFinder's `_rel` suffix
   marks proteins *related to but not conferring* resistance -- an unfiltered run called
   `qseB_rel`/`qseC_rel` a colistin-resistant strain.
3. **Gene count alone over-calls a locus.** Three strains hit "K1 capsule": P2_E3 with all
   seven kps genes at 96-99%, and two others carrying only the kpsD/E/M export apparatus at
   83-87%. That machinery is shared across capsule types, so a >=90% identity floor was
   added for high-consequence calls.

This does not replace ANI. ANI says what the organism is, which is what maps onto a
risk-group listing (TRBA 466); this says what the isolate carries. Both are needed, and
species-level ANI for the 47 is the outstanding step.

## What actually limits yield: the reductive branch

`dld` is present in 46 of 47, so gating on it would leave a bench of one. It is treated as
a **process parameter and a knockout target, not a selection criterion**.

The reason is *not* that it can be ignored. An earlier draft of this document claimed
*E. coli* carries `dld` and still reaches 85% of theoretical D-lactate from glycerol. That
is wrong and the error mattered, because it was the argument for treating `dld` as benign.
Mazumdar, Clomburg and Gonzalez (AEM 2010) **deleted it**: both their strains carry
Δ*dld*, described there as "a mutation that blocked the aerobic d-lactate dehydrogenase".
The best-documented process for this exact conversion removed the gene rather than
tolerating it, which is evidence in the opposite direction to the one originally cited.

What survives is narrower and independent of that paper: `dld` is quinone-linked
(EC 1.1.5.12), so it needs a terminal acceptor to turn over and is thermodynamically idle
in a strictly fermentative regime. Under the microaerobic conditions LS30 used, or with
nitrate present, it is active and will re-consume product. Treat it as a gene to delete,
or as a reason to control aeration, and not as a strain-selection filter.

The full precedent is worth stating, because it sets the expectation for this project:

> LA02Δ*dld* (Δ*frdA* Δ*pta* Δ*adhE* Δ*dld*) with GlpK-GlpD overexpressed produced
> **32 g/L D-lactate from 40 g/L glycerol, 85% of the theoretical maximum, >99.9% chiral
> purity.** The sister strain LA01 was Δ*pflB* Δ*frdA* Δ*dld*.

Four deletions and one overexpression. **Homofermentative D-lactate from glycerol was an
engineering result, not a strain-selection result**, and no wild isolate in this screen
should be expected to reproduce it off the shelf. What the screen contributes is the
starting point: which strains need the fewest of those deletions. That is exactly what the
`competing_product` column measures.

The axis that does separate strains is what competes for the same NADH:

| competing product | no dld | dld present |
|---|---:|---:|
| acids only | 1 | 5 |
| 2,3-butanediol | 0 | 24 |
| 1,3-propanediol | 0 | 15 |

**39 of 45 clean strains divert carbon to 1,3-PDO or 2,3-BD.** This is the mechanism behind
LS30's own numbers: 27.77 g/L 1,3-PDO against 14.68 g/L lactate. On glycerol the reductive
branch, not acetate, is the main rival for reducing power.

Six strains make neither, and carry a lactate transporter:

| strain | genus | glycerol route | dld knockout needed | grows anaerobically in silico |
|---|---|---|---|---|
| P2_A11 | Citrobacter | both | no | yes |
| P1_G3 | Enterobacter | both | yes | yes |
| MYb174 | unassigned | glpK only | yes | yes |
| MYb176 | unassigned | glpK only | yes | yes |
| P3_E8 | Buttiauxella | both | yes | **no** |
| P3_E10 | Buttiauxella | both | yes | **no** |

P2_A11 is the only strain on the bench that is biosafety-clean, makes no competing
reductive product, exports lactate, and needs no knockout.

The last column is a correction from the modelling below: `P3_E8` and `P3_E10` reach the
D-lactate yield ceiling in an FBA probe but cannot grow anaerobically in their models, and
an FBA yield probe does not require growth. They are not equivalent to the other four.

## Models

45 CarveMe reconstructions, from the 45-strain panel (`d_guild_panel.yaml`, generated by
`d_guild_bench.py --emit-panel --exclude-flagged`). The two biosafety-flagged strains were
never built. Genomes and models live alongside the L pool; separation is carried by
`candidate_genomes/manifest_d_guild.tsv`, not by directory.

### Every model was structurally unable to secrete D-lactate

The first audit returned a D-lactate yield of exactly zero in all 45 — from strains whose
genomes carry D-LDH, in a conversion where *E. coli* reaches 85% of theoretical. The cause
was one bound:

```
D_LACt2pp   h_p + lac__D_p --> h_c + lac__D_c    bounds (0, 1000)
```

CarveMe wrote the lactate proton-symporter inward-only in all 45. Nothing else moves
`lac__D_c` outward, so the models could consume D-lactate and never produce it. The rest of
the route was already intact: `D_LACtex` (periplasm to extracellular) is reversible in all
45 and `EX_lac__D_e` is present.

`curate_d_lactate_export` reverses that single bound where the genome carries both a
D-LDH confirmed by EC 1.1.1.28 and a lactate permease (K03303). No reaction is added.
The justification is that these are proton symporters, whose direction is set by the
proton-motive force and the lactate gradient rather than by the enzyme, so
unidirectionality is the artifact.

`D_LACt2` is deliberately left alone. It is BiGG's two-compartment form
(`h_e + lac__D_e <=> h_c + lac__D_c`), which collapses both membranes; in a model with a
periplasm it bypasses it, and reversing it alongside the pp/tex pair would open a proton
path from periplasm to extracellular space. 17 models carry it inward-only, harmlessly.

### Result

| | strains |
|---|---:|
| curated (D-LDH + permease) | 41 of 45 |
| secrete D-lactate at **1.000 mol/mol**, the theoretical maximum | 41 |
| ...and can also grow anaerobically | **37** |
| no permease, export gap left in place | 4 |

The 4 uncurated are `BIGb0393`, `MYb121`, `P1_A7`, `P2_F2`.

Two signs the amendment is a correction rather than a convenience. It lands on exactly
1.000 mol/mol — the theoretical carbon maximum for a C3 to C3 conversion, not an arbitrary
improvement. And it made the redox gate *stricter on itself*: `P1_C12` had been flagged and
now reads electron recovery 1.000, because a one-way proton symporter is itself a
free-energy source and making it reversible closed that. With every uptake closed the
objective stays at 0.0000, so no free lunch was opened.

### One model excluded

`MYb535` fails the redox gate: electron recovery 1.153 against a 1.015 threshold, from an
ethanol yield of 1.259 mol/mol against a ceiling of 1.167. Moved to `models/excluded/`.

### What the models cannot answer

`yield_13ppd` is nonzero in 1 of 45, while the genomic screen found the glycerol
dehydratase route in 15. CarveMe did not reconstruct it, most likely the B12-dependent
`dhaB`. **So the central yield question — whether the reductive branch steals carbon from
lactate — cannot be tested in silico with these models.** Adding the pathway would be an
assumption rather than curation. It stays an experimental question.

Related: yields for different products are independent maxima from separate probes, not
co-products. A strain reporting lactate 1.000 and formate 2.11 does not make both at once,
and a 1.000 anaerobic lactate yield still leaves +1 NADH per glycerol to be paid for.

## Identity, by ANI

`identify_strains.py --manifest candidate_genomes/manifest_d_guild.tsv`.

34 of 45 resolved to species.

| species | n |
|---|---:|
| *Citrobacter braakii* | 14 |
| *Enterobacter ludwigii* | 7 |
| *Lelliottia nimipressuralis* | 7 |
| *Pantoea pleuroti* / *P. rwandensis* / *P. coffeiphila* | 2 / 2 / 1 |
| *Serratia plymuthica*, *Citrobacter freundii*, *Kluyvera ascorbata*, *Rouxiella chamberiensis* | 1 each |

**28 of 45 are three species.** The pool is more redundant than the genus view showed.

### MEGAN was wrong about several strains this project had built arguments on

| strain | MEGAN said | ANI says | |
|---|---|---|---|
| P1_C12 | *Pantoea* | ***Citrobacter braakii*** | 99.0%, different family |
| P2_A11 | *Citrobacter* | ***Kluyvera ascorbata*** | 96.9% |
| P2_F2 | *Cedecea* | ***Rouxiella chamberiensis*** | 98.8% |
| P3_E8 / P3_E10 | *Buttiauxella* | novel *Pantoea* | 89.4%, below threshold |

`P1_G3` remains unplaced at *Pantoea coffeiphila* **84.6%**, aligned 0.72 — the identical
figure to its original exclusion note. Reinstating it for the D fork did not resolve the
one exclusion reason that survived the fork.

## Second-level signals: what separates members of a redundant guild

`guild_diversity.py`. 45 strains is 2^45 subsets, so the question is not which strain is
best but what the smallest panel is that spans the diversity.

### Gene-content clustering collapses the pool 45 -> 13

Jaccard distance on KEGG ortholog profiles, single-linkage at <=0.10. **The clusters are
species boundaries**, recovered independently -- the clustering never saw the ANI:

| cluster | n | species |
|---|---:|---|
| 0 | 15 | *C. braakii* x14, *C. freundii* |
| 1 | 14 | *L. nimipressuralis* x7, *E. ludwigii* x7 |
| 4, 5 | 2 each | *P. pleuroti* x2 / *P. rwandensis* x2 |
| 9-12 | 1 each | *Kluyvera*, *Rouxiella*, *P. coffeiphila*, *Serratia* |
| 2, 3, 6, 7 | 7 total | worm-collection isolates, identity pending |

29 of 45 sit in two clusters. A 13-strain panel spans the same functional range.

### Distinctness is confounded by assembly quality

A short or fragmented assembly loses genes, which inflates gene-content distance and
manufactures auxotrophies. Five strains are caveated, **including the top two ranks**:
`P2_A11` (3.85 Mb against a 4.93 Mb pool median), `P2_F2` (516 contigs), `P1_D10` (1021),
`P1_F1` (663), `P2_A4` (short). `guild_diversity.tsv` carries
`distinctness_trustworthy` so a rank is never read without it.

`P2_A11` is the instructive case. Its distinctness is **partly real** -- ANI makes it
*Kluyvera ascorbata*, the only member of that genus in the pool, which is exactly why it
sits alone in cluster 11. But roughly a megabase is missing against a typical *Kluyvera*,
so its six apparent amino-acid auxotrophies are not safe to believe. Real outlier,
incomplete assembly, both true at once.

### Biosynthesis: self-sufficient, and B12 is not a lever

**43 of 45 are complete for all 20 amino acids.** Biotin, pantothenate and B6 are 45/45;
folate and niacin 44; riboflavin 43; thiamine 42. The pool is biosynthetically
self-sufficient, which is consistent with having dropped the provisioner concept rather
than in tension with it.

**15 of 45 make B12 de novo** by the anaerobic `cbi` route, all in cluster 0. The aerobic
`cobN` route is absent pool-wide. But **no strain needs it**: `metE`, the B12-independent
methionine synthase, is present in 45/45. B12 is a lineage marker here, not a cross-feeding
currency, and designing an experiment around it would be designing around nothing.

### Accessory traits: which vary, and which of those actually mean something

Catalase, superoxide dismutase and cytochrome *bd* are 45/45 and carry no information at
all. Four traits do vary across the pool, but varying is not the same as mattering, and
only two have a mechanism that touches this process:

| trait | n | why it matters here |
|---|---:|---|
| nitrate respiration | 41/45 | **Supplies the terminal acceptor that switches `dld` on.** A strain with `narGHI` will re-oxidise its own D-lactate as soon as nitrate is available. This is the same reason nitrate was removed from the defined medium after it let models respire around every oxygen cap. |
| hydrogenase | 15/45 | **A route for disposing of surplus reducing equivalents.** Glycerol to lactate leaves +1 NADH per glycerol, which has to go somewhere; H2 evolution is one of the few sinks that costs no carbon, unlike 1,3-PDO or ethanol. Worth noting the KO set used here (`K06281`/`K06282`/`K18016`/`K18017`) is the uptake hydrogenases, not the formate-hydrogen-lyase complex, so this is suggestive rather than a direct FHL call. |
| urease | 32/45 | Nitrogen flexibility. Crude glycerol carries no urea, so no established relevance -- reported for completeness. |
| cellulose biofilm | 29/45 | Aggregation and cell retention, which would matter in a continuous reactor but not in the batch work planned. Also reported for completeness. |

The first two are levers on the redox problem that runs through this whole document. The
last two are variance without a demonstrated consequence, and are labelled as such rather
than presented as findings.

### Three marker sets that had to be corrected

The first run reported His and Tyr incomplete in 45/45 and B12 complete in 44/45. Both are
impossible for Enterobacteriaceae, which is how they were caught -- an implausible number
is a better error detector than reading the code.

* **His** -- the library carries the bifunctional `hisIE` (**K11755**, 44/45); scoring only
  the split K01523/K01496 forms called every strain a histidine auxotroph.
* **Tyr** -- `tyrA` is **K04517** or the bifunctional **K14187**. The plausible-looking
  K00220/K24018/K15226 are absent from all 45.
* **B12** -- the first attempt used K02302, the `cysG`-family methyltransferase shared with
  siroheme biosynthesis and present in all 45, plus `btuR`, which is a *salvage* enzyme
  that adenosylates imported cobinamide. Neither implies de novo synthesis. Switching to
  the committed steps `cbiK`/`cbiX`, `cbiA` and `cobS` took the count from 44 to 15. The
  script's own docstring had warned about this exact trap.

### The worm-collection isolates, and a scope error

Seven strains carry no MEGAN taxonomy at all, so neither the genus nor the family branch
could fire and all seven fell through as `no_reference_set` -- never searched, but
indistinguishable in the output from searched-and-not-found. `--fallback-family` was added
to give them a reference set, and the reporting now separates `no_ani_hit` from
`no_reference_set`.

Against Enterobacteriaceae, two resolved:

| strain | best hit | ANI | aligned | call |
|---|---|---:|---:|---|
| MYb174 | *Scandinavium goeteborgense* | 94.85 | 0.867 | genus match, species unresolved |
| MYb176 | *Scandinavium goeteborgense* | 94.83 | 0.865 | genus match, species unresolved |

*Scandinavium* was described in 2021, which is why no earlier method placed them.

The other five best-hit at **79-84% ANI with alignment as low as 0.24**. That is not a weak
match; it says the reference set is the wrong one. The cause is the **2016 reclassification**
that split the old broad Enterobacteriaceae into Yersiniaceae, Morganellaceae, Hafniaceae,
Pectobacteriaceae, Budviciaceae and Erwiniaceae -- a fermentative enterobacterium can sit in
any of those and be invisible to an Enterobacteriaceae-only search.

`--fallback-family` now accepts ORDER as well as FAMILY. Enterobacterales carries **568
species against Enterobacteriaceae's 203**, so the first search had 64% of the plausible
clade out of scope. Against the order, all five resolved -- and every one is Erwiniaceae,
exactly the family the split removed:

| strain | species | ANI | aligned |
|---|---|---:|---:|
| MYb416 | *Erwinia rhapontici* | 98.55 | 0.892 |
| MYb121 | *Erwinia billingiae* | 97.68 | 0.925 |
| MYb375 | *Erwinia rhapontici* | 96.58 | 0.919 |
| MYb535 | *Erwinia rhapontici* | 96.49 | 0.908 |
| BIGb0393 | ***Pantoea nemavictus*** | **100.0** | 0.999 |

`BIGb0393` is a strain-level identity: it *is* the deposited *P. nemavictus*, a species
described from the *C. elegans* microbiome, which is where the BIGb collection comes from.

### Reclassification is a standing hazard here, not a run of bad luck

`P1_A10` then failed for a related but distinct reason: its MEGAN genus reads *Raoultella*,
which **merged into *Klebsiella* in 2021** and no longer exists at genus rank, so the tax id
lookup returned nothing and the strain was skipped rather than searched. The fallback was
extended to fire on an unresolvable genus as well as a missing one. It then identified it as
***Klebsiella terrigena*** at 99.0% -- precisely the post-merge name for *Raoultella
terrigena*, so the fix confirmed its own diagnosis.

That is five reclassifications that have changed a conclusion in this project:
Brucella/Ochrobactrum (`P1_G5`), Mammaliicoccus/Staphylococcus (`P3_A3`),
Methylobacterium/Methylorubrum, the 2016 Enterobacterales split, and Raoultella/Klebsiella.
**A name absent from a reference set is not evidence that the organism is unusual.**

### Final identity: 45/45 searched, 40 to species

| | |
|---|---:|
| strain-level match | 1 |
| species match | 39 |
| genus match, species unresolved | 2 |
| no species match | 3 |

*C. braakii* 14, *E. ludwigii* 7, *L. nimipressuralis* 7, *E. rhapontici* 3, *Scandinavium
goeteborgense* 2, *P. pleuroti* 2, *P. rwandensis* 2, and one each of *P. nemavictus*,
*E. billingiae*, *K. terrigena*, *S. plymuthica*, *C. freundii*, *P. coffeiphila*,
*K. ascorbata*, *R. chamberiensis*.

Still unresolved, and each is a candidate novel species rather than a failure: `MYb174` and
`MYb176` at 94.85% to *Scandinavium goeteborgense* (a genus only described in 2021),
`P3_E8`/`P3_E10` at 89.4% to *Pantoea rwandensis*, and `P1_G3` at 84.6% to *P. coffeiphila*.

### The bench is Enterobacterales, not Enterobacteriaceae

Worth correcting an earlier statement. The fermentation gate was described above as
reproducing the LS30 membership, and at genus level it looked that way. With species names
attached it is broader: LS30 was Enterobacter, Escherichia, Klebsiella and Citrobacter --
all Enterobacteriaceae -- whereas this bench also carries seven Erwiniaceae (*Erwinia*,
*Pantoea*) and two *Scandinavium*. The physiology that the gate selected for is shared
across the order, so the guild is wider than its precedent.

## Open

- `P2_A11` needs its assembly resolved before anything is built on it: 3.85 Mb against
  ~4.9 Mb typical for *Kluyvera ascorbata*.
- Four candidate novel species to confirm or dismiss.
- Whether the 2,3-BD and 1,3-PDO carriers are worth keeping for impurity tolerance despite
  the carbon diversion is a fermentation question the models cannot settle.
- Nothing in the 384-genome library oxidises methanol (zero methanol dehydrogenases; 26
  genomes carry hxlA/hxlB and can detoxify formaldehyde without being able to make it). If
  the crude glycerol carries methanol, that role has to be bought in.
