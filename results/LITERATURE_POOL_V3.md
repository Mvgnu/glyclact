# Literature pass over the pool v3 strains, and the multi-producer question

Two things this document answers. First, what is published about each newly identified
strain — monoculture, since that is nearly all that exists, but in community context
wherever it does. Second, the question that has been sitting under the whole design:
**does adding more lactate-producing strains actually raise lactate?**

Three findings change project decisions. They are marked **⚠**.

---

# Part 1 — Does a second LDH help? Mostly no, and the exceptions say why

## The direct evidence

[Two-species community design of lactic acid bacteria for optimal production of lactate](https://pmc.ncbi.nlm.nih.gov/articles/PMC8605394/)
is the closest published study to this project's design question. It screened pairwise LAB
communities in silico and found genuine winners — up to a **10-fold increase in product
yield** over expectation. But the mechanism is not what "more LDH" would predict:

- **Substrate division did the work.** The best pairs put a glucose fermenter with a
  xylose-capable species. Where *both* members used *both* substrates, yield was **lower**
  than the glucose specialist alone.
- **Lactate itself was cross-fed in only 35% of communities** — the members were
  complementary, not additively productive.
- **Metabolic distance above 0.7 (70% dissimilarity) predicted success.** Similar strains
  did not help each other.
- ***L. plantarum* WCFS1 was the keystone organism, present in every top-performing pair.**

Supporting cases follow the same pattern: *L. brevis* + *L. plantarum* on glucose **and**
xylose reached 0.8 g/g, beating either monoculture — again two substrates, two specialists
([review](https://pmc.ncbi.nlm.nih.gov/articles/PMC7500094/)).

## What that means for this pool

**Our members share one substrate.** Glycerol is the only carbon input, so the mechanism
that made those co-cultures work is not available by default. That is consistent with —
and now explains — the negative result in `CAPACITY_RANKING.md`: total lactate is bounded
above by the best member's yield, so adding a lower-yield producer can only dilute.

**The pool does contain exactly one genuine substrate-division axis**, and it is the DHA
branch: `P3_C11` oxidises glycerol periplasmically to DHA, which `P1_B3` and `P3_F11` can
take up via DHA kinase. That splits one feedstock into two niches, which is structurally
the same move the successful LAB pairs made. It is the one multi-producer configuration
with published logic behind it, and it deserves priority in the presence/absence screen
over any pairing of two direct glycerol producers.

**⚠ A direct warning against our own medium recommendation.** The same paper found that
**parasitism dominated under minimal nutrient conditions, and hindered lactate yields**.
`CAPACITY_RANKING.md` recommends glycerol minimal medium because it is the most
carbon-efficient in FBA. Those two statements are in tension: FBA sees carbon accounting,
not interaction mode. The resolution is to treat medium richness as a *community* variable
rather than a per-strain optimisation — minimal may be right for a monoculture and wrong
for a consortium, and that is testable directly.

## The best published precedent for this project

[Bioconversion of raw glycerol from waste-cooking-oil biodiesel to 1,3-propanediol and
lactate by a microbial consortium](https://pmc.ncbi.nlm.nih.gov/articles/PMC6387908/) is
the closest thing in the literature to what we are attempting:

| | |
|---|---|
| Consortium | LS30, ~58% *Enterobacter* + ~39% *Escherichia*, enriched from activated sludge over 30 serial transfers |
| Conditions | **microaerobic**, **non-sterile**, pH 7.0, 37 °C |
| Products | 27.77 g/L 1,3-PDO **and 14.68 g/L lactate** in fed-batch |
| Versus monoculture | Four single isolates showed **long lag phases and much lower product** |

**The consortium's advantage was tolerance of crude glycerol's toxic impurities, not
yield.** That is precisely the benchmark this project's README already sets — a community
must earn its place through productivity, robustness, medium cost or crude-glycerol
tolerance, *not* yield. It is worth noting that the winning consortium was **enriched, not
designed**, and that non-sterile operation was a deliberate feature.

Two design points follow: **microaerobic** operation appears again (a third independent
literature source after the *E. coli* glycerol→D-lactate work and the *B. coagulans*
aeration data), and 1,3-PDO is a competing product on glycerol that our pool should be
watched for.

---

# Part 2 — The strains

## `P2_A5` — *Lactiplantibacillus plantarum*

The best-documented organism in the pool by a wide margin, and the literature is
double-edged.

**In its favour.** Keystone organism of every top-performing pair in the two-species study
above. Food-grade, GRAS, decades of process data. Its LDH is the documented
FBP-**independent** one, which is the whole basis of hypothesis H3.

**⚠ Against it, and this is serious for an aerobic process.** *L. plantarum* **consumes its
own lactate aerobically**. The route is well characterised: NAD-dependent LDH runs
*backwards* (lactate → pyruvate), pyruvate oxidase decarboxylates to acetyl-phosphate, and
acetate kinase releases acetate
([J Bacteriol 186:6661](https://journals.asm.org/doi/10.1128/jb.186.19.6661-6666.2004);
[AEM 72](https://journals.asm.org/doi/10.1128/aem.00659-06)).
**"Acetate is the major fermentation end product when *L. plantarum* is cultivated under
aerobic conditions and sugar limitation, and is produced at the expense of lactate."**

Our own data agree from three independent directions, which is reassuring for the pipeline
and bad news for the strain:

- the annotation carries **five pyruvate-oxidase copies**, the most in the pool, plus `nox`;
- the screen scored it `lactate_consumption = 1`;
- its capacity **peaks at O₂/glycerol 0.6 and falls to 4.00 at ratio 2.0** — uniquely among
  the producers, it is harmed by excess oxygen.

**Two levers, both already in our notes and both now literature-backed.** `poxB`/`poxF`
deletion: in a pyruvate-oxidase-null mutant "lactate remain[s] the major end product"
([J Bacteriol 186:3749](https://jb.asm.org/content/186/12/3749)) — this moves to the top of
`GENOMIC_MODIFICATION_OPPORTUNITIES.md`. And **pyruvate oxidase is repressed by glucose**,
so the small sugar co-feed proposed under H3 to raise FBP would *also* suppress the acetate
route. One intervention, two independent mechanisms — that raises the co-feed experiment's
value considerably.

Also relevant: lactate production is growth-associated with **product inhibition above
~15 g/L** ([Fermentation](https://pmc.ncbi.nlm.nih.gov/articles/PMC11051871/)), which
matches the export-thermodynamics argument in `OVERFLOW_TRIGGERS.md` H6.

## `P1_B3` — *Priestia aryabhattai*

Plant-growth-promoting and endophytic in essentially all the literature; no lactate
fermentation record. It does produce **butanoic acid**
([PGPR study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10909845/)), and one strain was
isolated from **Vietnamese pickled mustard greens** with probiotic characterisation
([Sci Rep](https://www.nature.com/articles/s41598-025-26744-0)) — a genuine fermentation
context, and the only such record for our top-yield producer.

**⚠ The important finding is a carbon sink our models cannot see.** *Priestia* species are
documented **polyhydroxybutyrate producers**, with PHB accumulation strongly dependent on
carbon source ([Priestia PHB](https://www.thepharmajournal.com/archives/2023/vol12issue12/PartAY/12-12-147-136.pdf)),
and *P. aryabhattai* specifically has been used as the PHB-producing partner in a synthetic
consortium ([ACS ES&T Eng](https://pubs.acs.org/doi/10.1021/acsestengg.5c01078)).

PHB is intracellular carbon and energy storage — a direct competitor for the carbon we want
in lactate. I checked whether the models can represent it: **no PHB polymer species and no
3-hydroxybutyryl-CoA sink exists in any of them.** So `P1_B3`'s 0.850 carbon yield, the
highest in the pool, is computed in a model that is structurally blind to a carbon drain
its species is known for. That does not invalidate the ranking, but it means **PHB should
be measured** — and it is a plausible explanation in advance if bench yield undershoots.

## `P3_G8` — *Rhodococcus erythropolis*

Valued precisely for what it was added for. "Extraordinary metabolic versatility and unique
environmental persistence and robustness"
([Appl Microbiol Biotechnol](https://link.springer.com/article/10.1007/s00253-020-10861-z)),
and documented consortium behaviour: complementary degradation roles with other
*Rhodococcus*
([carbendazim](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11767360/)), and co-cultivation
with *E. coli* or *B. subtilis* **increasing the biological activity of its surface-active
substances** ([Enzyme Microb Technol](https://www.sciencedirect.com/science/article/abs/pii/S0141022920301708)).

That last point is a caution as well as a feature: it produces surfactants, and it modifies
what it makes in the presence of partners. Not a lactate producer — our ambiguous-LDH call
is consistent with the literature.

## `P3_F8` — *Ewingella americana*

**⚠ A biosafety note that should be raised, not buried.** Our top provisioner (23 pathways,
7 amino-acid exporters, best oxidative defence) was **described in 1983 from clinical
specimens**, and is a **plant/fungal pathogen — the causal agent of internal stipe necrosis
in cultivated mushrooms**
([Phytopathology](https://apsjournals.apsnet.org/doi/10.1094/PHYTO-08-22-0299-R)). A
comparative genomics paper exists specifically on the genus's **resistome and virulome**
([PMC7281767](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7281767/)).

Our screen gave it zero virulence name matches and no risk-genus flag, which is why it
passed. That was the screen working as designed — and it is also exactly the pattern that
missed *Brucella* in `P1_G5`. It is a mushroom pathogen rather than a human one, so this is
a lower-stakes call than `P3_A3`, but it should be a conscious decision rather than a
default. It produces chitinase and behaves as a mycoparasite, which is worth knowing if the
lab handles fungal cultures.

## `P2_D5` — *Pantoea pleuroti*

Plant-associated; the genus is used both as a biocontrol agent and, industrially, as a host
for **L-amino acid production**
([review](https://www.sciencedirect.com/science/article/pii/S1049964424000033)) — which is
direct support for the provisioning role it was picked for. Also documented that
**co-culture induces cryptic biosynthetic pathways** in *Pantoea*
([PMC7072716](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7072716/)): it makes different
things in company than alone, which is an argument for testing it in community rather than
inferring from monoculture.

## `P2_E7` — *Chryseobacterium pennae*

Little strain-level literature, but the genus supports the role it was chosen for.
*Chryseobacterium proteolyticum* is the source organism of **protein-glutaminase**, which
deamidates glutamine residues in proteins and peptides
([Eur J Biochem](https://www.researchgate.net/publication/12101049_Protein-glutaminase_from_Chryseobacterium_proteolyticum_an_enzyme_that_deamidates_glutaminyl_residues_in_proteins_Purification_characterization_and_gene_cloning)).
A proteolytic Bacteroidete that liberates peptides and amino acids is a plausible
in-situ substitute for the yeast extract that hypothesis H5 prices at up to 38% of
lactic-acid production cost — and it is the one supporter in the pool that **consumes no
lactate**, so it can provision without competing for product.

## `P3_C11` — *Gluconobacter frateurii*

The DHA route is industrially established: high-yield DHA from **crude** glycerol in
fed-batch *G. oxydans*
([Molecules 29:2932](https://www.mdpi.com/1420-3049/29/12/2932)). And the genus is the
textbook case of a **two-stage consortium** — the classical two-step vitamin C fermentation
([Sci Rep](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4727326/)). So the DHA-donor role
has both the monoculture chemistry and the consortium precedent behind it.

## Cut strains, briefly

*Leuconostoc mesenteroides* (`P1_F12`) is heterofermentative and food-grade, and
*Leuconostoc* was the **xylose-capable partner** in the successful LAB pairs above — its
value there was substrate division, which our single-substrate design cannot offer, and our
gate showed it cannot grow on glycerol at all. Consistent with the cut.
*Mammaliicoccus sciuri* (`P3_A3`) is unchanged: AMR screen first.

---

# What I would change

1. **Prioritise the DHA-division communities** (`P3_C11` + `P1_B3`/`P3_F11`) over
   two-direct-producer pairings in the presence/absence screen. It is the pool's only
   substrate-division axis, and substrate division is what made every published
   lactate co-culture work.
2. **Treat medium richness as a community variable, not a per-strain optimum.** Minimal
   medium is most carbon-efficient in FBA but the LAB study found minimal nutrients drive
   parasitic interactions that *reduce* lactate. Test minimal versus modestly supplemented
   at community level.
3. **Move `poxB`/`poxF` to the top of the engineering list** for `P2_A5`, and raise the
   priority of the glucose co-feed experiment — it suppresses pyruvate oxidase *and*
   raises FBP, two independent routes to the same outcome.
4. **Add PHB to the measurement list** for `P1_B3`. Its species is a documented producer
   and no model in the pool can represent the drain.
5. **Put `P3_F8` in front of the supervisor** alongside `P3_A3` — a mushroom pathogen
   originally described from clinical specimens is a conscious-decision strain, even
   though it is a lower-stakes call.
6. **Expect the win to be robustness, not yield.** The one published crude-glycerol
   consortium beat its own isolates on lag phase and impurity tolerance while producing
   14.68 g/L lactate — and it was enriched from sludge, not designed. That is the
   benchmark to measure this pool against.
