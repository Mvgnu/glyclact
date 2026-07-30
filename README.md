# Glycerol-to-lactate community assembly

## Objective

Construct minimal communities that convert glycerol into optically defined
lactate in inexpensive mineral medium. Community value must come from a
measurable function such as nutrient rescue, DHA carbon transfer, crude-glycerol
tolerance, redox balancing, or process stability.

The benchmark is a single strain that performs the complete conversion. A
community is useful only if it improves yield, productivity, medium cost,
robustness, or the need for strain engineering.

## Candidate architectures

1. **Direct producer:** one strain carries glycerol entry, lower glycolysis, and
   a productive NAD-dependent LDH.
2. **Direct producer plus supporter:** the direct producer is nutritionally
   demanding and a glycerol-grown prototroph supplies required metabolites.
3. **DHA shuttle:** one strain oxidizes glycerol and releases dihydroxyacetone
   (DHA); another imports/phosphorylates DHA and converts it to lactate.

The annotation screen is triage, not validation. A pathway in a genome does not
prove flux, extracellular metabolite release, growth on glycerol, or safe use.

## Process regime

Operation is **aerobic**. That decision drives the screen: under aeration the
yield-limiting competition is pyruvate dehydrogenase, pyruvate oxidase and NADH
oxidase draining pyruvate and NADH away from lactate, so those are scored
explicitly. Anaerobic infeasibility in a reconstruction is a gap to note, not a
blocker.

Candidate choice is driven by **medium simplicity**: which strains grow on
glycerol plus mineral salts, what each one requires, and which partners can
supply the difference. The GEMs evaluate those choices; DHA, glycerol uptake and
lactate production are layered on afterwards in COMETS.

## Directory layout

- `screen_genomes.py` — eggNOG annotation screen (v2). Marker sets live in
  `markers.py` so they can be reviewed independently of the scoring.
- `identify_strains.py` — species/strain assignment by ANI against one RefSeq
  genome per species, escalating to family when no species matches.
- `audit_models.py` — model gate: mass balance, free-lunch loops, growth on
  glycerol/glucose/DHA, and product ceilings per unit glycerol.
- `results/` — score table, marker hits, ranked summary, `MODEL_GATE.md`,
  `SCREEN_V2_FINDINGS.md`, `strain_identity.tsv`.
- `candidate_genomes/`, `models/`, `references/`, `comets_runs/` — gitignored.
- `docs/` — research decisions, validation criteria, and experimental plan.

## Workflow

```text
384 annotated local genomes
        ↓
annotation screen: glycerol entry, lactate route AND export, aerobic pyruvate
competition, per-pathway requirements and provisions, assembly QC, biosafety
        ↓
ANI species assignment for shortlisted strains
        ↓
shortlist producers, DHA donors/receivers, supporters, and controls
        ↓
CarveMe reconstructions
        ↓
model gate: mass balance, no free lunch, glycerol growth, product ceilings
        ↓
GEM minimal-medium and dropout analysis — the authoritative requirement list
        ↓
match each requirement to a partner that can supply it
        ↓
reality gate: measured single-strain substrate use and products
        ↓
COMETS community × medium screen with product labels
        ↓
minimal lab-ready communities and defined media
```

Annotation is an approximation and a community's behaviour is not the sum of its
parts, so every stage above is a prior for the next, and the wet-lab gates in
`docs/PROJECT_BRIEF.md` remain the arbiter.

## Primary validation measurements

- glycerol disappearance;
- L- and D-lactate separately;
- DHA;
- acetate, ethanol, succinate, pyruvate, and 1,3-propanediol;
- biomass and viable counts;
- final community composition;
- ideally `13C3`-glycerol incorporation into M+3 lactate.

The eventual optimization target should reward lactate yield, productivity,
optical purity, glycerol conversion, and medium simplicity while penalizing
competing products, residual glycerol, and community instability.

