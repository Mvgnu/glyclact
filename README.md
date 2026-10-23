# Glycerol → L-lactate community

This repository records the literature research, genome mining, community designs,
experimental readouts and assay troubleshooting for a glycerol-to-lactate project.

**For the October 23 data audit, start with [the handoff](docs/handoff/START_HERE.md)
and [the analysis prompt](docs/handoff/GPT_PRO_PROMPT.md).** The concise actual
workflows are in `protocols/`. They supersede unperformed wash/cleanup steps in
historical planning documents.

Start with the [research and decision trail](docs/RESEARCH_TRAIL.md): it connects the literature, screening evidence, changing candidate pools and reasons for retaining, excluding or holding strains for validation.

## Scope of this repository

On October 23, 2026, Magnus requested a complete **project commit and push**, including
data. The snapshot now includes available genome/model inputs, wet-lab measurements,
strain mappings, Echo records, acquisition runs and analysis outputs. Runtime
caches, installed dependencies, credentials and OS metadata are excluded.
Inspect source-data permissions before any further redistribution.
These records do not establish validated lactate production by themselves.

Historical documents retain the hypotheses and wording used at each stage, including conclusions later revised. A historical claim about growth, secretion, optical purity or community benefit must not be read as a measured phenotype. The current panel rationale supersedes the earlier pools. This is a curated research record, not a complete transcript of every analysis or an untouched time-stamped pre-bench snapshot.

The source genome library and some reference databases are external inputs. The code and decision records are inspectable here, but reproducing every result requires those inputs and their applicable permissions. No rights to redistribute the source strain library are implied.

## Authoritative state

[panel/pool_v5_panel.yaml](panel/pool_v5_panel.yaml) is the only authoritative pool definition: **14 producers plus 2 strains held for chiral-HPLC resolution of a `lar`-operon risk**. There is no non-producing supporter.

The panel is intentionally manual. LDH, racemase, identity, physiology, and virulence screens produce evidence to inspect; they do not regenerate the YAML. In particular, virulence homology calls require locus-by-locus interpretation and are not safe as an automated exclusion rule. There is also no “lab plate only” eligibility rule.

The staged-input manifest, `candidate_genomes/manifest_pool_v5.tsv`, is derived
from the YAML. Available local genome files are included in this local snapshot;
see [input provenance](candidate_genomes/README.md). The genomic candidate panel is
not the executed species list for every culture run; use the dated Echo maps.

## Directory map

| Directory | Purpose |
|---|---|
| [panel/](panel/) | current manually curated panel |
| [screening/](screening/) | marker, LDH/fold, racemase, and virulence evidence |
| [identity/](identity/) | ANI hit generation and conservative taxonomic calls |
| [modeling/](modeling/) | input staging, CarveMe reconstruction, curation, and diagnostics |
| [physiology/](physiology/) | product capacity/consumption and gene-content diagnostics |
| [results/](results/) | generated and interpreted evidence, with cohort scopes indexed |
| [results/planning/](results/planning/) | proposed validation experiments; not completed experimental results |
| `references/` (external, untracked) | LDH, VFDB, and AMR reference material |
| [comets/](comets/) | historical COMETS inputs from the retired community enumerator |
| [archive/](archive/) | superseded panel logic and historical work |
| [protocols/](protocols/) | actual workflows, historical plans, raw readouts, Echo logs and designs |
| [work/](work/) | assembly reconstruction, September readout analyses and presentation work |
| [docs/handoff/](docs/handoff/) | experiment context, data inventory and independent-analysis prompt |

Each topical code directory has its own README describing what to run and how to interpret it.

## Decision boundary

The current panel combines:

1. Sequence-adjudicated L-LDH/D-LDH evidence and `lar`-system review.
2. Conservative identity evidence. A closest ANI hit below the species threshold remains a closest reference, not a species call.
3. Manual virulence/resistance review of the underlying homology hits.
4. Deliberate functional and taxonomic spread, using gene-content diversity as an aid rather than an automatic selector.
5. Model outputs as diagnostics only. CarveMe/FBA results are useful for finding transport, exchange, gap-fill, and redox problems; they are not organism growth-rate or medium-validation measurements.

## Process caveat: oxygen and alternative acceptors

The existing “anaerobic” lactate-consumption table means **oxygen-free with no nitrate, fumarate, or other alternative respiratory acceptor**. Its 35/36 non-consumer result must not be generalized to nitrate- or fumarate-containing production medium.

Quinone/cytochrome-linked L-lactate oxidation can feed anaerobic respiration when a compatible terminal acceptor is available. Nitrate or fumarate may help redox balance and growth while also reopening product consumption. Therefore acceptor choice is an experimental axis: measure L-lactate retention in the exact proposed medium rather than treating “low oxygen” or “nitrate added” as a solved condition.

## Typical workflow

From this directory:

```bash
python screening/screen_genomes.py \
  --genome-root /path/to/genome-library \
  --output-dir results/screening
python screening/assign_folds.py
python screening/adjudicate_ldh.py
python identity/identify_strains.py

# Requires the strain-library volume to refresh staged inputs.
python modeling/prepare_candidates.py
python modeling/reconstruct_models.py
python modeling/curate_models.py
python modeling/audit_models.py

python physiology/guild_diversity.py
python physiology/lactate_capacity.py
python physiology/lactate_consumption.py
```

Read [results/pool_v5/POOL_V5.md](results/pool_v5/POOL_V5.md) for the current selection rationale and [results/planning/EXPERIMENT_PLAN_V5.md](results/planning/EXPERIMENT_PLAN_V5.md) for experimental validation.
