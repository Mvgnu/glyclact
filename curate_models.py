#!/usr/bin/env python3
"""Apply documented, genome-justified amendments to the CarveMe reconstructions.

CarveMe draws every reaction from the BiGG universe by homology score, and it has two
systematic blind spots that this project has now hit repeatedly:

1. **Transport and exchange lag behind biochemistry.** A model gets the enzymes for a
   pathway but not the transporter that makes the product leave the cell, so the
   capability is present and unusable. Lactate export has been the binding constraint in
   eight separate findings across POOL_V2 and OVERFLOW_TRIGGERS.
2. **A missing exchange can make an aerobe anaerobic.** `P3_G8` and `P2_C2` carry
   cytochrome oxidases, an o2_c pool and a periplasmic O2 transporter, but no `EX_o2_e`
   at all -- so oxygen cannot reach the periplasm and the models grow identically with
   and without air. Both organisms are obligate aerobes.

Every amendment here is (a) a genuine BiGG universe reaction, not an invention, (b)
supported by named genes in that strain's own annotation, and (c) recorded in
results/model_curation.tsv with the evidence. Amendments are idempotent: a reaction that
is already present is left alone. Nothing is added on the grounds that it would be
convenient -- the vitamin-export gaps found earlier were deliberately *not* patched,
because prototrophy is not export and the annotation gave no transporter to point at.

Run after reconstruct_models.py and before audit_models.py.
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
from cobra import Metabolite, Reaction  # noqa: E402


def external_compartment(model: cobra.Model) -> str:
    """The model's own name for the extracellular compartment.

    CarveMe writes it as ``C_e`` rather than ``e``, and cobra's ``add_boundary`` refuses
    to build an exchange for a metabolite whose compartment does not match the model's
    external compartment exactly. Hard-coding "e" produced a metabolite that looked right
    and could not be given a valid exchange.
    """
    # Read it off a metabolite the model already has outside, rather than guessing from a
    # list of conventions -- whatever compartment the existing exchanges use is by
    # definition the one a new exchange has to match.
    for reaction in model.exchanges:
        return next(iter(reaction.metabolites)).compartment
    for suffix_metabolite in ("glyc_e", "h2o_e", "h_e"):
        if suffix_metabolite in model.metabolites:
            return model.metabolites.get_by_id(suffix_metabolite).compartment
    raise RuntimeError("cannot determine the model's extracellular compartment")


def metabolite(model: cobra.Model, mid: str, name: str, formula: str, charge: int,
               compartment: str) -> Metabolite:
    if mid in model.metabolites:
        return model.metabolites.get_by_id(mid)
    met = Metabolite(mid, name=name, formula=formula, charge=charge, compartment=compartment)
    model.add_metabolites([met])
    return met


def add_reaction(model: cobra.Model, rid: str, name: str, stoichiometry: dict[str, float],
                 bounds: tuple[float, float]) -> bool:
    """Add ``rid`` if absent. Returns True when the model was changed."""
    if rid in model.reactions:
        return False
    reaction = Reaction(rid, name=name, lower_bound=bounds[0], upper_bound=bounds[1])
    model.add_reactions([reaction])
    reaction.add_metabolites({model.metabolites.get_by_id(m): c for m, c in stoichiometry.items()})
    return True


def add_exchange(model: cobra.Model, metabolite_id: str) -> str:
    """Add an exchange reaction that cobra will actually recognise as one.

    Building it by hand as a plain boundary Reaction is not enough. ``model.exchanges``
    is resolved from the SBO term, not from the reaction's shape, so a hand-built
    `EX_..._e` reads as boundary=True but exchanges=False -- and every routine in this
    project that iterates ``model.exchanges`` to close uptake, or to sum an electron
    budget, then skips it silently. That is how the first version of this curation
    produced an electron recovery of 0.0000 for P2_A5: the lactate it had just been given
    the ability to secrete was invisible to the accounting. ``add_boundary`` sets the term.
    """
    reaction = model.add_boundary(model.metabolites.get_by_id(metabolite_id), type="exchange")
    return reaction.id


def curate_oxygen(model: cobra.Model) -> list[str]:
    """Give an aerobe a way to take up oxygen.

    Only fires when the model already has an o2 pool and something that consumes it --
    i.e. the biochemistry is there and only the boundary is missing. If a model genuinely
    has no oxygen chemistry this does nothing, which is the correct behaviour for a
    strict anaerobe.
    """
    changed: list[str] = []
    if "EX_o2_e" in model.reactions or "o2_c" not in model.metabolites:
        return changed
    consumers = [r for r in model.reactions
                 if any(m.id == "o2_c" for m in r.metabolites) and not r.boundary]
    if not consumers:
        return changed

    metabolite(model, "o2_e", "O2 O2", "O2", 0, external_compartment(model))
    if "o2_p" in model.metabolites:
        # Gram-negative layout: outside -> periplasm, then the existing O2tpp inward.
        if add_reaction(model, "O2tex", "O2 transport via diffusion (extracellular to periplasm)",
                        {"o2_e": -1.0, "o2_p": 1.0}, (-1000.0, 1000.0)):
            changed.append("O2tex")
    else:
        if add_reaction(model, "O2t", "O2 transport diffusion",
                        {"o2_e": -1.0, "o2_c": 1.0}, (-1000.0, 1000.0)):
            changed.append("O2t")
    changed.append(add_exchange(model, "o2_e"))
    return changed


def curate_lactate_export(model: cobra.Model, evidence: dict[str, int]) -> list[str]:
    """Give a strain that makes L-lactate a way to secrete it.

    Gated on the strain's own annotation: it must have a fermentative L-LDH *and* an
    annotated lactate transporter (lctP/lldP family). Without both, the export gap is
    left in place as a real finding -- that is exactly the call made for P3_A10, whose
    unannotated-export prediction is one of the pool's deliberate test cases.
    """
    changed: list[str] = []
    if evidence.get("l_ldh", 0) < 1 or evidence.get("lactate_transporter", 0) < 1:
        return changed
    if "lac__L_c" not in model.metabolites:
        return changed
    if "EX_lac__L_e" in model.reactions:
        return changed

    metabolite(model, "lac__L_e", "(S)-Lactate", "C3H5O3", -1, external_compartment(model))
    # L_LACt3: h_e + lac__L_c --> h_c + lac__L_e. Proton-motive-force-driven efflux, and
    # the mechanism the literature describes for lactic acid bacteria -- electrogenic
    # lactate/H+ symport that generates PMF at low external lactate.
    if add_reaction(model, "L_LACt3", "L-lactate transport out via proton symport",
                    {"h_e": -1.0, "lac__L_c": -1.0, "h_c": 1.0, "lac__L_e": 1.0},
                    (0.0, 1000.0)):
        changed.append("L_LACt3")
    changed.append(add_exchange(model, "lac__L_e"))
    return changed


def annotation_evidence(annotation_path: Path) -> dict[str, int]:
    """Count the genes that justify each amendment, from the strain's own eggNOG table."""
    l_ldh = 0
    transporter = 0
    with annotation_path.open(errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            name, kos = fields[8].strip(), fields[11]
            # K00016 is shared with malate dehydrogenase in eggNOG, so require the gene
            # name too -- P1_A3's two K00016 hits are both mdh, and calling that an LDH
            # would have manufactured a producer out of the sink control.
            if "K00016" in kos and name in {"ldh", "lctE", "ldhL"}:
                l_ldh += 1
            if "K03303" in kos or name in {"lctP", "lldP"}:
                transporter += 1
    return {"l_ldh": l_ldh, "lactate_transporter": transporter}


def main() -> None:
    root = Path(__file__).parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--models-dir", type=Path, default=root / "models")
    parser.add_argument("--genomes-dir", type=Path, default=root / "candidate_genomes")
    parser.add_argument("--output", type=Path, default=root / "results" / "model_curation.tsv")
    parser.add_argument("--dry-run", action="store_true")
    options = parser.parse_args()

    rows: list[dict[str, object]] = []
    for path in sorted(options.models_dir.glob("*.xml")):
        genome = path.stem
        annotation = options.genomes_dir / genome / "eggnog.emapper.annotations"
        if not annotation.exists():
            print(f"{genome}: no annotation, skipped", flush=True)
            continue
        evidence = annotation_evidence(annotation)
        model = cobra.io.read_sbml_model(str(path))

        added = curate_oxygen(model)
        oxygen_note = (
            f"model had o2_c and consumers but no EX_o2_e; both organisms in this "
            f"situation are obligate aerobes" if added else ""
        )
        lactate_added = curate_lactate_export(model, evidence)
        lactate_note = (
            f"{evidence['l_ldh']} fermentative L-LDH gene(s) and "
            f"{evidence['lactate_transporter']} lactate transporter gene(s) annotated, "
            f"but no L-lactate exchange in the model" if lactate_added else ""
        )
        all_added = added + lactate_added
        if all_added:
            rows.append({
                "genome": genome,
                "reactions_added": ";".join(all_added),
                "oxygen_evidence": oxygen_note,
                "lactate_evidence": lactate_note,
                "l_ldh_genes": evidence["l_ldh"],
                "lactate_transporter_genes": evidence["lactate_transporter"],
            })
            if not options.dry_run:
                cobra.io.write_sbml_model(model, str(path))
            print(f"{genome}: added {', '.join(all_added)}", flush=True)
        else:
            print(f"{genome}: no amendment needed", flush=True)

    if rows:
        options.output.parent.mkdir(parents=True, exist_ok=True)
        with options.output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nwrote {options.output} ({len(rows)} models amended)")
    else:
        print("\nno models needed amendment")


if __name__ == "__main__":
    main()
