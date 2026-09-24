#!/usr/bin/env python3
"""GEM medium and exchange diagnostics for the candidate pool.

This is a reconstruction diagnostic, not an organismal minimal-medium prediction.
CarveMe was gapfilled on M9 and its biomass reaction omits several vitamins, so the
results describe what the models can represent in COMETS.

Three questions, in the order they matter for medium design:

1. Does the reconstruction grow on the defined base used elsewhere in this project?
2. Which supplied components are individually essential in that reconstruction?
3. Which exchange reactions can carry secretion while the model still grows?

Single-component dropout is not a true medium minimisation: alternative nutrients can
be jointly necessary without either being individually essential. Provision columns are
capacity envelopes, not evidence that a culture actually secretes the metabolite.

A fourth output matters as much as the biology: **which cross-feeding mechanisms are
representable at all**. A metabolite with no exchange reaction in a model cannot be
traded no matter what the genome says, so any community mechanism resting on it is
untestable until the model is curated.

Writes results/models/minimal_medium.tsv, results/models/provision_matrix.tsv
and results/models/exchangeable_mechanisms.tsv.
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

# Carbon-free mineral salts, identical to audit_models.py so the two gates agree.
MINERAL_BASE = [
    "nh4", "pi", "so4", "h2o", "k", "na1", "mg2", "ca2", "fe2", "fe3", "cl",
    # Nitrate omitted: ammonium already covers nitrogen, and nitrate doubles as a
    # respiratory electron acceptor, which distorts anything oxygen-related.
    "cobalt2", "cu2", "mn2", "zn2", "mobd", "ni2", "sel", "slnt", "tungs",
    "cd2", "hg2",
]
# Sulfide and free protons are also absent: both act as electron donors in these
# reconstructions and inflate reduced-product yields past the redox ceiling. Full
# reasoning in audit_models.py.

AMINO_ACIDS = [
    "ala__L", "arg__L", "asn__L", "asp__L", "cys__L", "gln__L", "glu__L", "gly",
    "his__L", "ile__L", "leu__L", "lys__L", "met__L", "phe__L", "pro__L", "ser__L",
    "thr__L", "trp__L", "tyr__L", "val__L",
]
COFACTORS = [
    "thm", "ribflv", "pnto__R", "pydxn", "pydx", "pydam", "btn", "fol", "nac",
    "lipoate", "pheme", "sheme", "4abz",
]
NUCLEOBASES = ["ade", "gua", "ura", "thym", "csn", "ins", "hxan", "xan"]
OTHER = ["spmd", "ptrc", "26dap__M", "glyc3p", "dha", "ac", "pyr", "succ", "mal__L", "cit"]
# Products are scanned alongside supplements. Omitting lac__L made the mechanism table
# report "lactate: no_secretor" when the production envelope in fact reaches
# 0.81-0.92 mol/mol at half-maximal growth.
PRODUCTS = ["lac__L", "lac__D", "etoh", "for", "13ppd"]
SUPPLEMENT_PANEL = [*AMINO_ACIDS, *COFACTORS, *NUCLEOBASES, *OTHER, *PRODUCTS]

# End metabolites whose de novo synthesis the annotation flagged as missing in at
# least one pool member, used to reconcile model prototrophy against annotation.
KEY_COFACTORS = {
    "biotin": "btn_c", "heme": "pheme_c", "thiamine": "thmpp_c",
    "pyridoxal": "pydx5p_c", "pantothenate": "pnto__R_c", "methionine": "met__L_c",
    "riboflavin": "ribflv_c", "folate": "fol_c", "lipoate": "lipoate_c",
}

# Named community mechanisms, each with the exchanged metabolite, how the receiver
# would benefit, and which biomass components it feeds.
#
# The benefit route matters for judging testability. A *cofactor* cross-feed can only
# show up if the biomass equation contains it or something downstream of it — and the
# demand must be checked against the right species: biomass holds pydx5p, not pyridoxine,
# and thf, not folate. A *carbon* cross-feed needs no biomass entry at all; the receiver
# benefits by growing on it, so the test is whether a consumer grows. A *product* needs
# only a secretor.
MECHANISMS = {
    "biotin rescue":                       ("btn",     "cofactor", ["btn_c"]),
    "heme rescue (catalase, respiration)": ("pheme",   "cofactor", ["pheme_c"]),
    "lipoate supply (PDH cofactor)":       ("lipoate", "cofactor", ["lipoate_c"]),
    "methionine supply to DHA donors":     ("met__L",  "cofactor", ["met__L_c"]),
    "pantothenate supply":                 ("pnto__R", "cofactor", ["pnto__R_c", "coa_c"]),
    "pyridoxal B6 supply":                 ("pydxn",   "cofactor", ["pydxn_c", "pydx5p_c"]),
    "thiamine B1 supply":                  ("thm",     "cofactor", ["thm_c", "thmpp_c"]),
    "folate B9 supply":                    ("fol",     "cofactor", ["fol_c", "thf_c"]),
    "riboflavin B2 supply":                ("ribflv",  "cofactor", ["ribflv_c", "fad_c"]),
    "cobalamin B12 provision":             ("cbl1",    "cofactor", ["cbl1_c", "adocbl_c"]),
    "respiratory quinone":                 ("mqn8",    "cofactor", ["mqn8_c", "q8_c"]),
    "DHA transfer":                        ("dha",     "carbon",   []),
    "glycerol-3-phosphate transfer":       ("glyc3p",  "carbon",   []),
    "lactate (product)":                   ("lac__L",  "product",  []),
}

GLYCEROL_UPTAKE = 10.0
OXYGEN_UPTAKE = 20.0
SUPPLEMENT_UPTAKE = 1.0
GROWTH_EPSILON = 1.0e-6
SECRETION_GROWTH_FRACTION = 0.5


def open_medium(model: cobra.Model, supplements: list[str], carbon: str = "glyc") -> None:
    """Close all uptake, then open mineral salts, oxygen, the carbon source and supplements.

    Supplements are capped at ``SUPPLEMENT_UPTAKE`` rather than opened wide. Opening
    twenty amino acids at 1000 each lets the model eat them as carbon — growth then
    reads 12-82/h against 0.1-0.7 on glycerol, which measures nothing about
    supplementation. A trace cap keeps them in the role a supplement actually plays.
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
    if "EX_o2_e" in model.reactions:
        model.reactions.EX_o2_e.lower_bound = -OXYGEN_UPTAKE
    carbon_id = f"EX_{carbon}_e"
    if carbon_id in model.reactions:
        model.reactions.get_by_id(carbon_id).lower_bound = -GLYCEROL_UPTAKE


def solve(model: cobra.Model) -> float:
    value = model.slim_optimize()
    return 0.0 if value is None or value != value else float(value)


def essential_components(model: cobra.Model, components: list[str], baseline: float) -> list[str]:
    """Components whose removal costs more than 99% of growth."""
    essential = []
    for metabolite_id in components:
        reaction_id = f"EX_{metabolite_id}_e"
        if reaction_id not in model.reactions:
            continue
        reaction = model.reactions.get_by_id(reaction_id)
        original = reaction.lower_bound
        reaction.lower_bound = 0.0
        if solve(model) < 0.01 * baseline:
            essential.append(metabolite_id)
        reaction.lower_bound = original
    return essential


def secretion_capacity(model: cobra.Model, baseline: float) -> dict[str, float]:
    """Max secretion of each panel metabolite while holding growth at >= half of max.

    Requiring growth is what separates a real provision from a stoichiometric
    curiosity: a strain that can only make a metabolite by not growing will not feed
    a partner in a community.
    """
    provisions: dict[str, float] = {}
    biomass = next((r for r in model.reactions if r.objective_coefficient != 0), None)
    if biomass is None or baseline <= GROWTH_EPSILON:
        return provisions
    original_bounds = (biomass.lower_bound, biomass.upper_bound)
    biomass.lower_bound = SECRETION_GROWTH_FRACTION * baseline
    try:
        for metabolite_id in SUPPLEMENT_PANEL:
            reaction_id = f"EX_{metabolite_id}_e"
            if reaction_id not in model.reactions:
                continue
            with model:
                model.objective = reaction_id
                value = solve(model)
            if value > 1.0e-6:
                provisions[metabolite_id] = round(value, 4)
    finally:
        biomass.lower_bound, biomass.upper_bound = original_bounds
    return provisions


def biosynthesis_evidence(model: cobra.Model) -> dict[str, str]:
    """For each key cofactor, whether de novo synthesis has gene support.

    This reconciles the model against the annotation. The models grow on glycerol
    minimal with no supplements, which contradicts the annotation's auxotrophy calls,
    so the question is whether the model's route is backed by genes.

    Two distinctions make the difference between signal and noise:

    * **Transport is not synthesis.** A reaction producing ``btn_c`` from ``btn_e`` is
      biotin *uptake*. In P3_B12 the only gene-backed producers of biotin are
      ``BTNt2i`` and ``BTNt``, both transporters — which agrees with the annotation
      calling it an auxotroph rather than contradicting it. Reactions spanning more
      than one compartment are classified as transport and excluded.
    * **Peptidases are not synthesis either**, but they are single-compartment, so a
      "gene-supported" verdict here still means "can obtain", not "can synthesise de
      novo from central metabolism". A full pathway-level audit is the stronger test.

    Verdicts: ``gene_supported`` (>=1 gene-backed cytosolic producer),
    ``evidence_free`` (producers exist but none gene-backed — the model's prototrophy
    rests on reactions with no genomic support, so the annotation's requirement call
    is corroborated), ``transport_only``, or ``absent``.
    """
    verdicts: dict[str, str] = {}
    for label, metabolite_id in KEY_COFACTORS.items():
        if metabolite_id not in model.metabolites:
            verdicts[label] = "absent"
            continue
        metabolite = model.metabolites.get_by_id(metabolite_id)
        producers = [
            r for r in metabolite.reactions
            if not r.boundary and r.metabolites.get(metabolite, 0) > 0
        ]
        internal = [r for r in producers if len({m.compartment for m in r.metabolites}) == 1]
        transport = [r for r in producers if r not in internal]
        gene_backed = [r for r in internal if r.gene_reaction_rule.strip()]
        if gene_backed:
            verdicts[label] = "gene_supported"
        elif internal:
            verdicts[label] = "evidence_free"
        elif transport:
            verdicts[label] = "transport_only"
        else:
            verdicts[label] = "absent"
    return verdicts


def analyse(sbml_path: Path) -> dict[str, object]:
    genome = sbml_path.stem
    model = cobra.io.read_sbml_model(str(sbml_path))
    exchange_ids = {r.id for r in model.exchanges}

    # 1. Glycerol minimal medium: salts + glycerol only, no supplements.
    probe = model.copy()
    open_medium(probe, [])
    growth_minimal = solve(probe)

    # 2. Same medium plus the full supplement panel, to see whether supplements help
    #    at all and to give the dropout a baseline when minimal growth is zero.
    rich = model.copy()
    open_medium(rich, SUPPLEMENT_PANEL)
    growth_rich = solve(rich)

    # 3. Requirements. If the strain grows without supplements it has none by
    #    definition, and the informative question becomes which salts are essential.
    if growth_minimal > GROWTH_EPSILON:
        probe2 = model.copy()
        open_medium(probe2, [])
        required_supplements: list[str] = []
        essential_salts = essential_components(probe2, MINERAL_BASE, growth_minimal)
    else:
        rich2 = model.copy()
        open_medium(rich2, SUPPLEMENT_PANEL)
        required_supplements = essential_components(rich2, SUPPLEMENT_PANEL, growth_rich)
        essential_salts = essential_components(rich2, MINERAL_BASE, growth_rich)

    # 4. Provisions, measured on the medium the strain will actually be grown in.
    provider = model.copy()
    open_medium(provider, [])
    baseline = growth_minimal if growth_minimal > GROWTH_EPSILON else 0.0
    provisions = secretion_capacity(provider, baseline) if baseline else {}

    # 5. Reconcile against the annotation: is the prototrophy gene-supported?
    # Which panel metabolites the biomass equation actually demands. CarveMe's biomass
    # reaction (57 metabolites, identical across the pool) omits biotin, pantothenate,
    # lipoate, heme and folate, so those auxotrophies cannot be predicted by these models
    # at all — which is the root reason every strain reads prototrophic against an
    # annotation that calls several of them auxotrophs.
    biomass_reaction = next((r for r in model.reactions if r.objective_coefficient != 0), None)
    biomass_metabolites = {m.id for m in biomass_reaction.metabolites} if biomass_reaction else set()
    biomass_demand = {m for m in SUPPLEMENT_PANEL if f"{m}_c" in biomass_metabolites}

    # Growth on each carbon-route mechanism metabolite, since a carbon cross-feed is
    # judged by whether a partner can grow on it rather than by biomass composition.
    grows_on: dict[str, float] = {}
    for _metabolite, _route, _forms in MECHANISMS.values():
        if _route != "carbon" or f"EX_{_metabolite}_e" not in model.reactions:
            continue
        _probe = model.copy()
        open_medium(_probe, [], carbon=_metabolite)
        grows_on[_metabolite] = round(solve(_probe), 4)

    evidence = biosynthesis_evidence(model)
    evidence_free = sorted(k for k, v in evidence.items() if v == "evidence_free")
    transport_only = sorted(k for k, v in evidence.items() if v == "transport_only")

    return {
        "genome": genome,
        "growth_glycerol_minimal": round(growth_minimal, 4),
        "synthesis_evidence_free_count": len(evidence_free),
        "synthesis_evidence_free": ";".join(evidence_free) or "none",
        "synthesis_transport_only": ";".join(transport_only) or "none",
        "growth_glycerol_plus_supplements": round(growth_rich, 4),
        "supplements_improve_growth": int(growth_rich > growth_minimal * 1.05),
        "grows_on_glycerol_minimal": int(growth_minimal > GROWTH_EPSILON),
        "required_supplements_count": len(required_supplements),
        "required_supplements": ";".join(required_supplements) or "none",
        "essential_salts_count": len(essential_salts),
        "essential_salts": ";".join(essential_salts) or "none",
        "can_secrete_count": len(provisions),
        "can_secrete": ";".join(f"{k}:{v:g}" for k, v in sorted(provisions.items())) or "none",
        "biomass_demands_count": len(biomass_demand),
        "biomass_demands": ";".join(sorted(biomass_demand)) or "none",
        "_provisions": provisions,
        "_biomass_demand": biomass_demand,
        "_biomass_forms": biomass_metabolites,
        "_grows_on": grows_on,
        "_exchanges": exchange_ids,
    }


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--models-dir", type=Path, default=root / "models")
    parser.add_argument("--manifest", type=Path,
                        default=root / "candidate_genomes" / "manifest_pool_v5.tsv")
    parser.add_argument("--output-dir", type=Path, default=root / "results" / "models")
    options = parser.parse_args()

    with options.manifest.open(encoding="utf-8") as handle:
        genomes = {row["genome"] for row in csv.DictReader(handle, delimiter="\t")}
    paths = sorted(p for p in options.models_dir.glob("*.xml") if p.stem in genomes)
    rows = []
    for path in paths:
        print(f"analysing {path.stem}", flush=True)
        rows.append(analyse(path))
    options.output_dir.mkdir(parents=True, exist_ok=True)

    fields = [k for k in rows[0] if not k.startswith("_")]
    with (options.output_dir / "minimal_medium.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    # Provision matrix: strain x metabolite, only metabolites someone can secrete.
    secretable = sorted({m for row in rows for m in row["_provisions"]})
    with (options.output_dir / "provision_matrix.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metabolite", *(row["genome"] for row in rows)])
        for metabolite in secretable:
            writer.writerow([metabolite, *(f"{row['_provisions'].get(metabolite, 0):g}" for row in rows)])

    # Which named mechanisms are representable at all. A cross-feed needs BOTH a
    # supplier and a consumer with demand: if the biomass equation does not contain the
    # metabolite, no partner can benefit from receiving it and the mechanism cannot show
    # up as a growth effect however well the supply side works.
    with (options.output_dir / "exchangeable_mechanisms.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "mechanism", "metabolite", "benefit_route", "models_with_exchange",
            "models_that_can_secrete", "receivers_that_benefit", "testable", "blocked_by",
        ])
        for mechanism, (metabolite, route, biomass_forms) in MECHANISMS.items():
            with_exchange = [r["genome"] for r in rows if f"EX_{metabolite}_e" in r["_exchanges"]]
            can_secrete = [r["genome"] for r in rows if r["_provisions"].get(metabolite, 0) > 1e-6]
            if route == "cofactor":
                receivers = [r["genome"] for r in rows
                             if any(f in r["_biomass_forms"] for f in biomass_forms)]
                missing_demand = "biomass_demands_no_form_of_it"
            elif route == "carbon":
                receivers = [r["genome"] for r in rows if r["_grows_on"].get(metabolite, 0) > 1e-6]
                missing_demand = "no_consumer_grows_on_it"
            else:
                receivers = ["(product, no receiver needed)"]
                missing_demand = ""

            if not with_exchange:
                testable, blocked = "no", "no_exchange_reaction"
            elif not can_secrete:
                testable, blocked = "no", "no_secretor"
            elif not receivers:
                testable, blocked = "no", missing_demand
            else:
                testable, blocked = "yes", ""
            writer.writerow([
                mechanism, metabolite, route, f"{len(with_exchange)}/{len(rows)}",
                ";".join(can_secrete) or "none",
                ";".join(receivers) if route != "cofactor" else f"{len(receivers)}/{len(rows)}",
                testable, blocked,
            ])

    print("\n=== minimal medium ===")
    for row in rows:
        note = "" if row["grows_on_glycerol_minimal"] else "  <-- NO GROWTH on glycerol minimal"
        print(f"  {row['genome']:8} growth {row['growth_glycerol_minimal']:>7} "
              f"(+supplements {row['growth_glycerol_plus_supplements']:>7})  "
              f"requires: {row['required_supplements']}{note}")
    print(f"\nwrote {options.output_dir}/minimal_medium.tsv, provision_matrix.tsv, exchangeable_mechanisms.tsv")


if __name__ == "__main__":
    main()
