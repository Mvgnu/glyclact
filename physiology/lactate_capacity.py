#!/usr/bin/env python3
"""Per-strain lactate capacity across medium conditions.

Titre from dFBA is not a prediction for these models: at maximum growth the minimum
lactate is zero and the maximum is ~8 mmol/gDW/h in every producer, so which value a
simulation reports is the solver's arbitrary choice of vertex
(see results/models/PRODUCT_PREDICTION_LIMITS.md).

Capacity is the well-posed question instead. Holding growth at a fixed fraction of its
maximum and then maximising lactate has a unique answer, and it means something
concrete: **the best this strain could do on this medium while still growing.** That is
the right quantity for ranking community x medium combinations to take to the bench,
because the bench is what will decide where in the feasible range the strain actually
sits.

Two numbers are reported per strain and condition:

* ``max_growth`` — whether the strain grows at all here, which gates everything.
* ``lactate_capacity`` — max lactate flux at ``--growth-fraction`` of max growth,
  and the same divided by glycerol uptake to give a yield.

The competing-product columns are independent capacity envelopes: each product is
maximised in its own solution at the same growth constraint. They do not describe
co-production in the L-lactate-maximising solution.

Writes results/lactate_biology/capacity_table.tsv.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import logging
import warnings
from pathlib import Path

logging.disable(logging.CRITICAL)
warnings.filterwarnings("ignore")

import cobra  # noqa: E402

MINERAL_BASE = [
    "nh4", "pi", "so4", "h2o", "k", "na1", "mg2", "ca2", "fe2", "fe3", "cl",
    "cobalt2", "cu2", "mn2", "zn2", "mobd", "ni2", "sel", "slnt", "tungs",
    "cd2", "hg2",
]

# Nitrate is deliberately *not* in the mineral base. It was, and it silently broke every
# oxygen conclusion: nitrate is a respiratory electron acceptor, Bacillus respires it, and
# every model here carries a nitrite efflux exchange. With nitrate open at -1000 the cell
# simply switched acceptor, so growth was flat against oxygen and oxygen looked like no
# lever at all. Ammonium already supplies nitrogen, so nothing is lost by closing it.
# These are closed explicitly rather than merely omitted, in case a model's own bounds
# leave them open.
ALT_ELECTRON_ACCEPTORS = [
    "EX_no3_e", "EX_no2_e", "EX_dmso_e", "EX_tmao_e", "EX_fum_e", "EX_so3_e", "EX_tsul_e",
]

# Sulfide and free protons are absent from the base above for the same reason, one class
# further on: both act as electron *donors* here. Sulfide oxidation pushed ethanol to 1.50
# mol/mol glycerol against a redox ceiling of 1.167, and unbounded H+ uptake let the
# models' several hundred charge-unbalanced reactions behave as a free electron source
# (P2_E7 took up 57 H+ per 10 glycerol). See audit_models.py for the full note.

# Supplements the annotation flags as requirements somewhere in the pool. The GEM
# cannot see biotin, pantothenate, lipoate, heme or folate because the CarveMe biomass
# equation omits them, so supplying them changes nothing in silico — they are swept
# anyway so the table stays valid if the biomass equation is later amended.
SUPPLEMENT_SETS = {
    "none": [],
    "vitamins": ["btn", "pnto__R", "thm", "ribflv", "fol", "pydxn", "nac", "lipoate"],
    "vitamins_aa": ["btn", "pnto__R", "thm", "ribflv", "fol", "pydxn", "nac", "lipoate",
                    "met__L", "cys__L", "ile__L", "leu__L", "val__L", "his__L"],
}

# Acetoin and 2,3-butanediol are here because they were missing and they are the textbook
# Bacillus overflow sinks -- two pyruvate into acetolactate, no NADH consumed. Leaving them
# out meant the competition for pyruvate was never scored against lactate.
PRODUCTS = ["EX_lac__L_e", "EX_lac__D_e", "EX_etoh_e", "EX_for_e", "EX_ac_e", "EX_succ_e",
            "EX_dha_e", "EX_actn__R_e", "EX_btd_RR_e"]
SUPPLEMENT_UPTAKE = 1.0

# Carbon per molecule, for an honest yield. Dividing lactate by *glycerol* alone makes
# supplement-rich media look best: fourteen amino acids at 1 mmol each carry more carbon
# than 20 mmol of glycerol, so the strain is partly making lactate out of the
# supplements. Yield is therefore reported per carbon actually taken up.
CARBON_NUMBER = {
    "glyc": 3, "btn": 10, "pnto__R": 9, "thm": 12, "ribflv": 17, "fol": 19,
    "pydxn": 8, "nac": 6, "lipoate": 8, "met__L": 5, "cys__L": 3, "ile__L": 6,
    "leu__L": 6, "val__L": 5, "his__L": 6,
}


def carbon_uptake(model: cobra.Model, solution, supplements: list[str]) -> tuple[float, float]:
    """(glycerol carbon, supplement carbon) actually consumed at this solution."""
    def uptake(metabolite_id: str) -> float:
        reaction_id = f"EX_{metabolite_id}_e"
        if reaction_id not in model.reactions:
            return 0.0
        return max(-float(solution.fluxes.get(reaction_id, 0.0)), 0.0)

    glycerol_carbon = uptake("glyc") * CARBON_NUMBER["glyc"]
    supplement_carbon = sum(uptake(s) * CARBON_NUMBER.get(s, 0) for s in supplements)
    return glycerol_carbon, supplement_carbon


def apply_condition(
    model: cobra.Model,
    glycerol: float,
    oxygen: float,
    supplements: list[str],
) -> None:
    """Close all uptake, then open the mineral base, glycerol and oxygen at given rates.

    Oxygen is capped as an uptake *rate*, not a supplied amount. That distinction is
    what makes it a real lever: limiting the amount merely starves the culture once it
    is used up, whereas limiting the rate makes respiratory capacity the binding
    constraint and pushes carbon into overflow products.

    The cap only binds once every *other* electron acceptor is closed -- see
    ALT_ELECTRON_ACCEPTORS.
    """
    for reaction in model.exchanges:
        reaction.lower_bound = 0.0
    for metabolite_id in MINERAL_BASE:
        reaction_id = f"EX_{metabolite_id}_e"
        if reaction_id in model.reactions:
            model.reactions.get_by_id(reaction_id).lower_bound = -1000.0
    for metabolite_id in supplements:
        reaction_id = f"EX_{metabolite_id}_e"
        if reaction_id in model.reactions:
            model.reactions.get_by_id(reaction_id).lower_bound = -SUPPLEMENT_UPTAKE
    for reaction_id in ALT_ELECTRON_ACCEPTORS:
        if reaction_id in model.reactions:
            model.reactions.get_by_id(reaction_id).lower_bound = 0.0
    if "EX_o2_e" in model.reactions:
        model.reactions.EX_o2_e.lower_bound = -oxygen
    if "EX_glyc_e" in model.reactions:
        model.reactions.EX_glyc_e.lower_bound = -glycerol


def solve(model: cobra.Model) -> float:
    value = model.slim_optimize()
    return 0.0 if value is None or value != value else float(value)


def capacity_at(model: cobra.Model, target: str, growth_fraction: float, max_growth: float):
    """Max flux through ``target`` while holding growth at ``growth_fraction`` of max.

    Returns (flux, solution) so the caller can read the uptake side of the same
    operating point rather than assuming it.
    """
    if target not in model.reactions or max_growth <= 1e-9:
        return 0.0, None
    with model:
        biomass = next(r for r in model.reactions if r.objective_coefficient != 0)
        biomass.lower_bound = growth_fraction * max_growth
        model.objective = target
        solution = model.optimize()
        if solution.status != "optimal":
            return 0.0, None
        return float(solution.objective_value or 0.0), solution


def evaluate(
    sbml_path: Path,
    glycerol_levels: list[float],
    oxygen_ratios: list[float],
    growth_fraction: float,
) -> list[dict[str, object]]:
    genome = sbml_path.stem
    base = cobra.io.read_sbml_model(str(sbml_path))
    rows: list[dict[str, object]] = []
    # Oxygen is swept as a ratio to glycerol, not as an absolute rate, because the ratio
    # is what the stoichiometry fixes: glycerol -> lactate leaves two surplus electrons
    # per glycerol, so 0.5 mol O2 per mol glycerol is the exact respiratory demand and
    # the same absolute O2 cap means something different at each glycerol level.
    for glycerol, oxygen_ratio, (supplement_name, supplements) in itertools.product(
        glycerol_levels, oxygen_ratios, SUPPLEMENT_SETS.items()
    ):
        oxygen = oxygen_ratio * glycerol
        model = base.copy()
        apply_condition(model, glycerol, oxygen, supplements)
        max_growth = solve(model)
        row: dict[str, object] = {
            "genome": genome,
            "glycerol_uptake": glycerol,
            "oxygen_per_glycerol": oxygen_ratio,
            "oxygen_uptake": oxygen,
            "supplements": supplement_name,
            "max_growth": round(max_growth, 5),
            "grows": int(max_growth > 1e-6),
        }
        if max_growth <= 1e-6:
            for target in PRODUCTS:
                short = target.replace("EX_", "").replace("_e", "")
                row[f"{short}_capacity"] = 0.0
            row["glycerol_carbon_uptake"] = 0.0
            row["supplement_carbon_uptake"] = 0.0
            row["lactate_capacity_yield"] = 0.0
            row["lactate_carbon_yield"] = 0.0
            row["supplement_carbon_share"] = 0.0
            row["growth_at_operating_point"] = 0.0
            rows.append(row)
            continue

        lactate_solution = None
        for target in PRODUCTS:
            short = target.replace("EX_", "").replace("_e", "")
            flux, solution = capacity_at(model, target, growth_fraction, max_growth)
            row[f"{short}_capacity"] = round(flux, 4)
            if target == "EX_lac__L_e":
                lactate_solution = solution
        # This project targets optically pure L-lactate. D-lactate is reported as an
        # independent diagnostic envelope and must not be added to the target yield.
        lactate = float(row["lac__L_capacity"])

        glycerol_carbon, supplement_carbon = (
            carbon_uptake(model, lactate_solution, supplements) if lactate_solution is not None else (0.0, 0.0)
        )
        total_carbon = glycerol_carbon + supplement_carbon
        row["glycerol_carbon_uptake"] = round(glycerol_carbon, 3)
        row["supplement_carbon_uptake"] = round(supplement_carbon, 3)
        # Per glycerol, kept for comparability with the classic yield figure.
        row["lactate_capacity_yield"] = round(lactate / glycerol, 4) if glycerol > 0 else 0.0
        # Per carbon actually consumed — the number that survives adding supplements.
        row["lactate_carbon_yield"] = (
            round(lactate * 3 / total_carbon, 4) if total_carbon > 1e-9 else 0.0
        )
        row["supplement_carbon_share"] = (
            round(supplement_carbon / total_carbon, 4) if total_carbon > 1e-9 else 0.0
        )
        row["growth_at_operating_point"] = round(growth_fraction * max_growth, 5)
        rows.append(row)
    return rows


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--models-dir", type=Path, default=root / "models")
    parser.add_argument("--manifest", type=Path,
                        default=root / "candidate_genomes" / "manifest_pool_v5.tsv")
    parser.add_argument("--output", type=Path, default=root / "results" / "lactate_biology" / "capacity_table.tsv")
    parser.add_argument("--glycerol", type=float, nargs="+", default=[5.0, 10.0, 20.0])
    parser.add_argument("--oxygen-ratio", type=float, nargs="+",
                        default=[2.0, 1.0, 0.6, 0.5, 0.4, 0.25, 0.1],
                        help="O2 uptake cap as a multiple of the glycerol uptake rate. "
                             "0.5 is the stoichiometric demand of glycerol -> lactate.")
    parser.add_argument("--growth-fraction", type=float, default=0.9,
                        help="Fraction of maximum growth held while maximising product.")
    options = parser.parse_args()

    with options.manifest.open(encoding="utf-8") as handle:
        genomes = {row["genome"] for row in csv.DictReader(handle, delimiter="\t")}
    rows: list[dict[str, object]] = []
    for path in sorted(p for p in options.models_dir.glob("*.xml") if p.stem in genomes):
        print(f"evaluating {path.stem}", flush=True)
        rows.extend(evaluate(path, options.glycerol, options.oxygen_ratio, options.growth_fraction))

    options.output.parent.mkdir(parents=True, exist_ok=True)
    with options.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nwrote {options.output} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
