#!/usr/bin/env python3
"""Does a producer eat its own product?

Gene counts answer this badly. `lldD` is the wrong question asked three ways: the
FMN-dependent alpha-hydroxyacid oxidase family holds L-lactate dehydrogenase (`lldD`),
L-lactate oxidase (`lctO`), glycolate oxidase, mandelate dehydrogenase and long-chain
hydroxyacid oxidase in one fold, so a homology hit names the family. And even a real
`lldD` only turns over when it has somewhere to put the electrons.

The model answers a narrower structural question, so this asks two questions per strain:

1. **Can it grow on L-lactate as its only carbon source?** Aerobic and anaerobic. A yes
   means the strain has a complete, connected route from extracellular lactate to biomass
   -- transport, oxidation, and a downstream sink for the electrons. This is the binary
   the pool decision needs.

2. **How fast, relative to glycerol?** A strain that grows on lactate at 5% of its
   glycerol rate is a different proposition from one that prefers it.

The medium is the one `audit_models.py` uses, so the numbers are comparable to the gate:
same mineral base, same 10 mmol/gDW/h carbon uptake, same 20 mmol/gDW/h oxygen.

The oxygen-free condition contains no alternative terminal electron acceptor. It therefore
does not predict consumption in the proposed nitrate-, fumarate-, or acetate-containing
production media; those conditions must be tested explicitly.

**Transport direction has to be neutralised or the answer is circular.** `curate_models.py`
adds `L_LACt3` with bounds (0, 1000) -- efflux only, which is the right default for a
producer but makes lactate uptake structurally impossible. The first run of this script
found exactly the twelve curated strains to be non-consumers, which measured the curation
rather than the organism. So every membrane reaction carrying `lac__L` is opened
reversibly for the duration of the test. That asks the question that can actually be
answered from a genome-scale model: *given that lactate reaches the cytoplasm, is there a
route from it to biomass?* Whether the transporter runs inward in a real vessel is a
kinetic question about gradients that FBA cannot answer either way.

Both directions are reported: `eats_L_lactate` uses reversible transport, and
`eats_L_lactate_asbuilt` uses the model as stored, so the gap between them is visible.

Writes results/pool_v5/lactate_consumption.tsv.
"""

from __future__ import annotations

import argparse
import csv
import logging
import warnings
from pathlib import Path

logging.disable(logging.CRITICAL)
warnings.filterwarnings("ignore")

import cobra  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]

# Kept identical to audit_models.py so the two tables can be read side by side.
MINERAL_BASE = [
    "h2o", "pi", "so4", "nh4", "k", "na1", "cl", "mg2", "ca2", "fe2", "fe3",
    "mn2", "zn2", "cu2", "cobalt2", "mobd", "ni2", "sel", "slnt", "tungs", "cd2", "hg2",
]
CARBON_UPTAKE = 10.0
OXYGEN_UPTAKE = 20.0


def apply_medium(model, carbon_ids, aerobic):
    for reaction in model.exchanges:
        reaction.lower_bound = 0.0
    for salt in MINERAL_BASE:
        rid = f"EX_{salt}_e"
        if rid in model.reactions:
            model.reactions.get_by_id(rid).lower_bound = -1000.0
    if aerobic and "EX_o2_e" in model.reactions:
        model.reactions.get_by_id("EX_o2_e").lower_bound = -OXYGEN_UPTAKE
    ok = True
    for cid in carbon_ids:
        rid = f"EX_{cid}_e"
        if rid not in model.reactions:
            ok = False
            continue
        model.reactions.get_by_id(rid).lower_bound = -CARBON_UPTAKE
    return ok


def solve(model) -> float:
    v = model.slim_optimize()
    return 0.0 if v is None or v != v else float(v)


def open_lactate_transport(model) -> list[str]:
    """Make every lac__L membrane reaction reversible, and say which ones they were.

    A transport reaction is one that has lac__L on both sides in different compartments.
    Exchanges are excluded -- their direction is set by the medium, not by the cell.
    """
    opened = []
    for reaction in model.reactions:
        if reaction.boundary:
            continue
        comps = {m.compartment for m in reaction.metabolites if m.id.startswith("lac__L")}
        if len(comps) > 1:
            reaction.lower_bound = -1000.0
            reaction.upper_bound = 1000.0
            opened.append(reaction.id)
    return opened


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--models-dir", type=Path, default=ROOT / "models")
    ap.add_argument("--manifest", type=Path,
                    default=ROOT / "candidate_genomes" / "manifest_pool_v5.tsv")
    ap.add_argument("--genomes", type=Path,
                    help="optional file with one genome id per line; overrides --manifest")
    ap.add_argument("--output", type=Path,
                    default=ROOT / "results" / "pool_v5" / "lactate_consumption.tsv")
    args = ap.parse_args()

    if args.genomes:
        ids = [ln.strip() for ln in args.genomes.read_text().splitlines() if ln.strip()]
    else:
        with args.manifest.open(encoding="utf-8") as handle:
            ids = [row["genome"] for row in csv.DictReader(handle, delimiter="\t")]
    rows = []
    for genome in ids:
        path = args.models_dir / f"{genome}.xml"
        if not path.exists():
            print(f"{genome}: no model")
            continue
        model = cobra.io.read_sbml_model(str(path))

        has_lac = "EX_lac__L_e" in model.reactions
        if not has_lac:
            # No exchange means the question cannot be asked of this model, which is a
            # different answer from "no". curate_models.py should have fixed this.
            rows.append({"genome": genome, "has_lac_exchange": "0"})
            print(f"{genome}: no EX_lac__L_e")
            continue

        # As stored, so the effect of the curated efflux direction stays visible.
        with model:
            apply_medium(model, ["lac__L"], aerobic=True)
            lac_asbuilt = solve(model)

        with model:
            apply_medium(model, ["glyc"], aerobic=True)
            gly_aer = solve(model)

        # Transport neutralised: the question is whether the metabolism can use lactate.
        with model:
            transporters = open_lactate_transport(model)
            apply_medium(model, ["lac__L"], aerobic=True)
            lac_aer = solve(model)
        with model:
            open_lactate_transport(model)
            apply_medium(model, ["lac__L"], aerobic=False)
            lac_ana = solve(model)

        rows.append({
            "genome": genome,
            "has_lac_exchange": "1",
            "mu_lactate_aerobic": f"{lac_aer:.4f}",
            "mu_lactate_anaerobic": f"{lac_ana:.4f}",
            "mu_glycerol_aerobic": f"{gly_aer:.4f}",
            "lactate_vs_glycerol": f"{(lac_aer / gly_aer):.2f}" if gly_aer > 1e-6 else "",
            "eats_L_lactate": "1" if lac_aer > 1e-4 else "0",
            "eats_L_lactate_asbuilt": "1" if lac_asbuilt > 1e-4 else "0",
            "lac_transporters": ";".join(transporters) or "none",
        })
        print(f"{genome}: lac_aer={lac_aer:.4f} (as-built {lac_asbuilt:.4f}) "
              f"lac_ana={lac_ana:.4f} gly={gly_aer:.4f}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = ["genome", "has_lac_exchange", "mu_lactate_aerobic", "mu_lactate_anaerobic",
              "mu_glycerol_aerobic", "lactate_vs_glycerol", "eats_L_lactate",
              "eats_L_lactate_asbuilt", "lac_transporters"]
    with args.output.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})
    eaters = sum(1 for r in rows if r.get("eats_L_lactate") == "1")
    print(f"\nwrote {args.output}")
    print(f"{eaters}/{len(rows)} grow on L-lactate as sole carbon source")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
