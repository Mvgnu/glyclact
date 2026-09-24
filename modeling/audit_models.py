#!/usr/bin/env python3
"""Model diagnostics for the glycerol-to-lactate candidate panel.

Runs the checks that decide whether a CarveMe reconstruction is fit to enter a
COMETS community screen at all:

* mass/charge balance (element-level imbalance is the disqualifying kind);
* no growth on a strictly carbon-free mineral base (catches free carbon);
* no product formation with zero carbon uptake (catches free-lunch loops);
* growth on glycerol, glucose and DHA, aerobic and anaerobic;
* maximum secretion of lactate, DHA and the competing products, per glycerol.

The mineral base is deliberately carbon-free. Opening the usual "trace"
exchanges (nucleobases, vitamins, folate) feeds carbon and silently invalidates
both the no-carbon control and every yield number.

Writes results/models/model_gate.tsv and prints a per-model summary.
"""

from __future__ import annotations

import argparse
import csv
import logging
import re
import warnings
from pathlib import Path

logging.disable(logging.CRITICAL)
warnings.filterwarnings("ignore")

import cobra  # noqa: E402
import pandas as pd  # noqa: E402

# Strictly carbon-free mineral salts. co2/hco3 are deliberately excluded: leaving
# them open lets a model pull inorganic carbon and push apparent product yields
# just past the carbon ceiling (seen as 1.05 mol lactate/mol glycerol in P3_A3).
MINERAL_BASE = [
    "nh4", "pi", "so4", "h2o", "k", "na1", "mg2", "ca2", "fe2", "fe3", "cl",
    # Nitrate is deliberately absent: ammonium already covers nitrogen, and nitrate
    # doubles as a respiratory electron acceptor, which silently rescues anaerobic
    # growth and makes any oxygen-limited test meaningless.
    "cobalt2", "cu2", "mn2", "zn2", "mobd", "ni2", "sel", "slnt", "tungs",
    "cd2", "hg2",
]
# Sulfide is out for the same reason as nitrate, one step further on: it is an electron
# *donor*, and models that can oxidise it pushed ethanol to 1.50 mol/mol glycerol against
# a redox ceiling of 1.167 (glycerol degree of reduction 14, ethanol 12). Closing it puts
# P3_F8 and P1_G3 back on the ceiling exactly. Sulfate still supplies sulfur.
#
# Free protons are out for a subtler version of the same problem. These reconstructions
# carry 200-500 charge-unbalanced reactions each, and an unbounded H+ *uptake* lets that
# imbalance act as a free electron source: P2_E7 reached 1.50 ethanol/glycerol by taking
# up 57 H+ per 10 glycerol, and six models sat at an identical, artificial mu of 0.7031.
# Closing proton uptake drops them to 0.55-0.60 and makes them distinguishable. This is a
# workaround for the charge imbalance, not a claim that cells do not import protons --
# fixing it properly would mean recharging several hundred reactions per model.
# Omitting "h" here leaves proton *export* open (only lower bounds are closed), so
# fermentation products can still be secreted with their protons.

CARBON_SOURCES = {"glycerol": "glyc", "glucose": "glc__D", "DHA": "dha", "glycerol3p": "glyc3p"}

# Products that matter for the economics: the target, the optical-purity partner,
# the DHA intermediate, and the competing sinks named in the brief.
PRODUCTS = [
    "EX_lac__L_e", "EX_lac__D_e", "EX_dha_e", "EX_ac_e", "EX_etoh_e",
    "EX_succ_e", "EX_pyr_e", "EX_for_e", "EX_13ppd_e",
]

GLYCEROL_UPTAKE = 10.0
OXYGEN_UPTAKE = 20.0


def apply_medium(model: cobra.Model, carbon: str | None, aerobic: bool) -> bool:
    """Close every uptake, then open the mineral base plus one carbon source.

    Returns False when the carbon source has no exchange reaction, which is a
    finding in its own right rather than a zero-growth result.
    """
    for reaction in model.exchanges:
        reaction.lower_bound = 0.0

    def open_uptake(metabolite_id: str, lower_bound: float) -> bool:
        reaction_id = f"EX_{metabolite_id}_e"
        if reaction_id not in model.reactions:
            return False
        model.reactions.get_by_id(reaction_id).lower_bound = lower_bound
        return True

    for salt in MINERAL_BASE:
        open_uptake(salt, -1000.0)
    if aerobic:
        open_uptake("o2", -OXYGEN_UPTAKE)
    if carbon is None:
        return True
    return open_uptake(carbon, -GLYCEROL_UPTAKE)


def solve(model: cobra.Model) -> float:
    value = model.slim_optimize()
    return 0.0 if value is None or value != value else float(value)


def element_imbalances(model: cobra.Model) -> list[str]:
    """Reaction ids with a non-charge, non-biomass mass imbalance."""
    offenders = []
    for reaction in model.reactions:
        if reaction.boundary or reaction.id == "Growth":
            continue
        try:
            balance = reaction.check_mass_balance()
        except Exception:
            continue
        if any(abs(value) > 1e-6 for key, value in balance.items() if key != "charge"):
            offenders.append(reaction.id)
    return offenders


def charge_imbalance_count(model: cobra.Model) -> int:
    count = 0
    for reaction in model.reactions:
        if reaction.boundary:
            continue
        try:
            balance = reaction.check_mass_balance()
        except Exception:
            continue
        if abs(balance.get("charge", 0.0)) > 1e-6:
            count += 1
    return count


ELECTRON_WEIGHTS = {"C": 4, "H": 1, "O": -2, "N": -3, "S": 6, "P": 5}

# Reduced products to drive the redox test with, in preference order. Ethanol is the most
# reduced and so the most sensitive probe, but five models in this panel have no ethanol
# exchange at all -- testing only on ethanol left them silently untested, which is exactly
# the kind of coverage gap that makes a gate worse than useless.
REDOX_PROBES = ("EX_etoh_e", "EX_lac__L_e", "EX_lac__D_e", "EX_ac_e", "EX_succ_e")


def available_electrons(metabolite: cobra.Metabolite) -> float | None:
    """Available electrons of a metabolite, or None when the formula cannot be used.

    For CcHhOoNnSsPp carrying net charge z, referenced to CO2/H2O/NH3/SO4/H3PO4:

        gamma = 4c + h - 2o - 3n + 6s + 5p - z

    This is a conservation law, not a fitted constant. It reproduces every textbook
    value: glycerol 14, ethanol 12, lactate 12, acetate 8, CO2 0, H2O 0, O2 -4, NH4+ 0,
    sulfate 0 -- and, decisively for this project, **H+ = 1 - 1 = 0**. A proton carries no
    available electrons, so in a model whose reactions all conserve charge, proton uptake
    cannot supply reducing power.

    Metabolites containing metals or R-groups return None: they sit outside this
    reference system, so counting them would be worse than skipping them. Skipping is
    safe for the *rejection* test below, because a skipped uptake can only make the
    measured recovery too low, never too high.
    """
    if not metabolite.formula or metabolite.charge is None:
        return None
    counts: dict[str, int] = {}
    for element, number in re.findall(r"([A-Z][a-z]?)(\d*)", metabolite.formula):
        if element:
            counts[element] = counts.get(element, 0) + (int(number) if number else 1)
    if set(counts) - set(ELECTRON_WEIGHTS):
        return None
    return sum(ELECTRON_WEIGHTS[e] * n for e, n in counts.items()) - metabolite.charge


def electron_recovery(model: cobra.Model, products: tuple[str, ...] = REDOX_PROBES) -> float | None:
    """Electrons leaving the cell divided by electrons entering it, maximising ``product``.

    This is the redox counterpart of the free-lunch test above, and it catches a whole
    class of defect that a mass balance cannot see. The mass test asks "can carbon appear
    from nothing"; every model in this panel passed it while three of them were quietly
    manufacturing *electrons* -- from sulfide left in the mineral base, from unbounded
    proton uptake against 200-500 charge-unbalanced reactions, and in one case from an
    internal loop that survives closing both.

    A value above 1.0 is impossible and disqualifying. Below 1.0 is normal: the remainder
    went into biomass. Reported at maximum reduced-product formation because that is the
    condition under which an electron leak is easiest to see.
    """
    if "EX_glyc_e" not in model.reactions:
        return None
    worst: float | None = None
    for product in products:
        if product not in model.reactions:
            continue
        probe = model.copy()
        if not apply_medium(probe, "glyc", aerobic=True):
            return None
        probe.objective = product
        solution = probe.optimize()
        if solution.status != "optimal":
            continue
        entering = leaving = 0.0
        for reaction in probe.exchanges:
            flux = solution.fluxes[reaction.id]
            if abs(flux) < 1e-9:
                continue
            electrons = available_electrons(next(iter(reaction.metabolites)))
            if electrons is None:
                continue
            if flux < 0:
                entering += -flux * electrons
            else:
                leaving += flux * electrons
        if entering <= 1e-9:
            continue
        ratio = leaving / entering
        # Keep the worst case: one leaking route is enough to disqualify a model, and
        # reporting the average would let a clean probe mask a broken one.
        if worst is None or ratio > worst:
            worst = ratio
    return None if worst is None else round(worst, 4)


def free_lunch(model: cobra.Model) -> dict[str, float]:
    """Products still reachable with zero carbon uptake — a hard rejection."""
    leaks = {}
    for target in [*PRODUCTS, "Growth"]:
        if target not in model.reactions:
            continue
        probe = model.copy()
        apply_medium(probe, None, aerobic=True)
        probe.objective = target
        value = solve(probe)
        if abs(value) > 1e-6:
            leaks[target] = value
    return leaks


def audit_model(sbml_path: Path) -> dict[str, object]:
    model = cobra.io.read_sbml_model(str(sbml_path))
    exchange_ids = {reaction.id for reaction in model.exchanges}
    row: dict[str, object] = {
        "genome": sbml_path.stem,
        "reactions": len(model.reactions),
        "metabolites": len(model.metabolites),
        "genes": len(model.genes),
        "exchanges": len(exchange_ids),
    }

    element_bad = element_imbalances(model)
    row["element_unbalanced_reactions"] = len(element_bad)
    row["element_unbalanced_examples"] = ";".join(element_bad[:5])
    row["charge_unbalanced_reactions"] = charge_imbalance_count(model)

    leaks = free_lunch(model)
    row["free_lunch"] = ";".join(f"{k}={v:.3g}" for k, v in leaks.items())

    recovery = electron_recovery(model)
    row["electron_recovery"] = recovery if recovery is not None else float("nan")
    # 1.5% tolerance absorbs rounding and the handful of metabolites whose formulas cannot
    # be scored; anything past it is manufacturing reducing power.
    row["redox_free_lunch"] = int(recovery is not None and recovery > 1.015)

    for label, metabolite_id in CARBON_SOURCES.items():
        for aerobic in (True, False):
            probe = model.copy()
            present = apply_medium(probe, metabolite_id, aerobic)
            suffix = "aerobic" if aerobic else "anaerobic"
            row[f"growth_{label}_{suffix}"] = solve(probe) if present else float("nan")
        row[f"has_EX_{metabolite_id}_e"] = int(f"EX_{metabolite_id}_e" in exchange_ids)

    for aerobic in (True, False):
        probe = model.copy()
        apply_medium(probe, None, aerobic)
        row[f"growth_no_carbon_{'aerobic' if aerobic else 'anaerobic'}"] = solve(probe)

    for target in PRODUCTS:
        short = target.replace("EX_", "").replace("_e", "")
        if target not in model.reactions:
            row[f"yield_{short}_aerobic"] = float("nan")
            row[f"yield_{short}_anaerobic"] = float("nan")
            continue
        for aerobic in (True, False):
            probe = model.copy()
            apply_medium(probe, "glyc", aerobic)
            probe.objective = target
            suffix = "aerobic" if aerobic else "anaerobic"
            row[f"yield_{short}_{suffix}"] = round(solve(probe) / GLYCEROL_UPTAKE, 4)
    return row


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--models-dir", type=Path, default=root / "models")
    parser.add_argument("--manifest", type=Path,
                        default=root / "candidate_genomes" / "manifest_pool_v5.tsv")
    parser.add_argument("--output", type=Path, default=root / "results" / "models" / "model_gate.tsv")
    options = parser.parse_args()

    with options.manifest.open(encoding="utf-8") as handle:
        genomes = {row["genome"] for row in csv.DictReader(handle, delimiter="\t")}
    rows = [audit_model(path) for path in sorted(options.models_dir.glob("*.xml"))
            if path.stem in genomes]
    frame = pd.DataFrame(rows)
    options.output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(options.output, sep="\t", index=False)

    for row in rows:
        print(f"=== {row['genome']}")
        print(f"    element-unbalanced {row['element_unbalanced_reactions']}, "
              f"charge-unbalanced {row['charge_unbalanced_reactions']}, "
              f"free-lunch {row['free_lunch'] or 'none'}")
        print(f"    growth glycerol aer/ana {row['growth_glycerol_aerobic']:.4f}/"
              f"{row['growth_glycerol_anaerobic']:.4f}   "
              f"no-carbon aer {row['growth_no_carbon_aerobic']:.4f}")
        missing = [
            key.replace("has_EX_", "EX_").replace("_e", "_e")
            for key, value in row.items()
            if key.startswith("has_EX_") and value == 0
        ]
        absent_products = [
            key.replace("yield_", "").replace("_aerobic", "")
            for key, value in row.items()
            if key.startswith("yield_") and key.endswith("_aerobic") and value != value
        ]
        print(f"    missing carbon-source exchanges: {', '.join(missing) or 'none'}")
        print(f"    absent product exchanges: {', '.join(absent_products) or 'none'}")
        print(f"    max lactate_L per glycerol aer/ana: {row['yield_lac__L_aerobic']}/"
              f"{row['yield_lac__L_anaerobic']}")
    print(f"\nwrote {options.output}")


if __name__ == "__main__":
    main()
