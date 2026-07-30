#!/usr/bin/env python3
"""Rank local annotated genomes for glycerol-to-lactate community roles.

Version 2. The process is aerobic, so the screen is organised around the three
questions that actually decide a candidate:

* **Can it make lactate from glycerol under aeration?** Which needs glycerol
  entry, lower glycolysis, an NAD-dependent LDH, a way to get lactate *out*, and
  as little pyruvate dehydrogenase / pyruvate oxidase / NADH oxidase competition
  as possible. v1 checked none of the last two.
* **What does it require?** Per-pathway auxotrophy priors, reported individually
  rather than averaged, because averaging hides the single missing pathway that
  defines a strain's medium.
* **What can it provide?** The same pathways read in the other direction, plus
  catalase and heme, which gate a partner's aerobic survival.

Requirements and provisions are reported separately and never collapsed into one
number: v1's ``9*aa_fraction + 9*vitamin_fraction`` gave a strain complete in six
of seven pathways and empty in the seventh a high supporter score, when it is in
fact an auxotroph that needs rescue.

These are annotation-based priors. Absence of a named gene is weak evidence of
absent function, so pathways carried by promiscuous enzymes are excluded from
auxotrophy calls, and the authoritative requirement list comes from the GEMs.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import markers as M

PROVIDES_THRESHOLD = 0.80
REQUIRES_THRESHOLD = 0.40
MIN_LOWER_GLYCOLYSIS = 5
MIN_CONTIG_N50 = 20_000
MAX_CONTIGS = 500


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def split_values(value: str) -> set[str]:
    return {item.removeprefix("ko:") for item in re.split(r"[,;]", value) if item and item != "-"}


GENUS_FAMILY = {
    "Aneurinibacillus": "Aneurinibacillaceae",
    "Bacillus": "Bacillaceae",
    "Priestia": "Bacillaceae",
    "Paenibacillus": "Paenibacillaceae",
    "Lactococcus": "Streptococcaceae",
    "Staphylococcus": "Staphylococcaceae",
    "Gluconobacter": "Acetobacteraceae",
    "Rahnella": "Yersiniaceae",
    "Ochrobactrum": "Brucellaceae",
    "Brucella": "Brucellaceae",
    "Ureibacillus": "Planococcaceae",
    "Pseudomonas": "Pseudomonadaceae",
}


def taxonomy(genome_dir: Path) -> tuple[str, str, str]:
    """Best genus and its support percentage, checked against its parent family.

    v1 returned the pipe-joined genus list and the summary table printed only
    that, so a call with 99.78% support for one genus and three genera below 0.2%
    read as "ambiguous". The percentages were already in the TSV; they just were
    not used.

    The percentages are normalised *within each rank row*, so a genus row holding a
    single sparse assignment reads as 100% however few proteins support it. P1_B3
    reported "Aneurinibacillus 100%" while its own family row said Bacillaceae
    99.98% — and Aneurinibacillus is not in Bacillaceae. ANI later placed it in
    *Priestia* (Bacillaceae) at 99.5%. So a genus is only trusted when the dominant
    family is consistent with it; otherwise the call is downgraded to the family and
    flagged, because the genus percentage is measuring the wrong denominator.
    """
    for base in (genome_dir / "diamond_megan", genome_dir):
        taxa_path, pct_path = base / "unique_taxa.txt", base / "unique_taxa_percentage.txt"
        if not taxa_path.exists():
            continue
        taxa = [x.strip().strip("|") for x in taxa_path.read_text(errors="replace").splitlines() if x.strip()]
        if not taxa:
            continue
        pcts = ([x.strip().strip("|") for x in pct_path.read_text(errors="replace").splitlines() if x.strip()]
                if pct_path.exists() else [])
        genera = taxa[-1].split("|")
        values: list[float] = []
        if pcts:
            try:
                values = [float(x) for x in pcts[-1].split("|") if x]
            except ValueError:
                values = []
        if values and len(values) == len(genera):
            ranked = sorted(zip(genera, values), key=lambda pair: pair[1], reverse=True)
            best_genus, best_pct = ranked[0]
            runner_up = ranked[1][1] if len(ranked) > 1 else 0.0
            flag = "clean" if best_pct >= 90.0 and runner_up < 10.0 else "mixed"

            # Consistency check against the dominant family, one rank up.
            dominant_family = ""
            if len(taxa) >= 2:
                family_names = taxa[-2].split("|")
                family_values: list[float] = []
                if len(pcts) >= 2:
                    try:
                        family_values = [float(x) for x in pcts[-2].split("|") if x]
                    except ValueError:
                        family_values = []
                if family_values and len(family_values) == len(family_names):
                    dominant_family = max(zip(family_names, family_values),
                                          key=lambda pair: pair[1])[0]
            expected_family = GENUS_FAMILY.get(best_genus)
            if dominant_family and expected_family and dominant_family != expected_family:
                # The genus is kept, not overruled. A mismatch has two causes and this
                # test cannot separate them: a genuinely spurious genus (P1_B3, where
                # "Aneurinibacillus" sat under Bacillaceae and ANI found *Priestia*), or
                # simply an older reference taxonomy (MEGAN places Rahnella in
                # Enterobacteriaceae and Ochrobactrum in Rhizobiaceae, both since
                # reclassified — those genus calls were right). Overruling the genus
                # would have discarded two correct calls, so the flag only says
                # "resolve this by ANI", which was the correct action in all three cases.
                return best_genus, f"{best_pct:.2f}", "verify_by_ani_family_mismatch"
            return best_genus, f"{best_pct:.2f}", flag
        return genera[0], "", "no_percentages"
    return "", "", "missing"


def assembly_quality(genome_dir: Path) -> dict[str, object]:
    """QUAST metrics, which sit unused next to every genome on the drive.

    Total length and GC% independently corroborate a taxonomy call, and a
    fragmented assembly degrades both annotation and reconstruction.
    """
    report = genome_dir / "quast_output" / "report.tsv"
    wanted = {"Total length": "assembly_bp", "GC (%)": "assembly_gc_percent",
              "# contigs": "contigs", "N50": "n50"}
    values: dict[str, object] = {name: "" for name in wanted.values()}
    values["assembly_quality"] = "no_quast"
    if not report.exists():
        return values
    for line in report.read_text(errors="replace").splitlines():
        parts = line.split("\t")
        if len(parts) >= 2 and parts[0] in wanted:
            values[wanted[parts[0]]] = parts[1]
    try:
        n50, contigs = int(values["n50"]), int(values["contigs"])
        if n50 >= 100_000 and contigs <= 100:
            values["assembly_quality"] = "good"
        elif n50 >= MIN_CONTIG_N50 and contigs <= MAX_CONTIGS:
            values["assembly_quality"] = "acceptable"
        else:
            values["assembly_quality"] = "fragmented"
    except (ValueError, TypeError):
        values["assembly_quality"] = "unparsed"
    return values


def pathway_status(names: set[str], definitions: dict[str, dict[str, object]]) -> dict[str, object]:
    """Per-pathway completeness, plus provides/requires lists.

    Pathways flagged ``unreliable`` keep their fraction but never enter the
    requires list, because their absence reflects naming, not biology.
    """
    fractions: dict[str, float] = {}
    provides: list[str] = []
    requires: list[str] = []
    for pathway, definition in definitions.items():
        # A function served by several unrelated routes counts as present if any
        # one route is complete, so the best route sets the fraction.
        routes = [set(definition.get("genes", set())), *definition.get("alternatives", [])]
        fraction = max(
            (len(names & route) / len(route) if route else 0.0) for route in routes
        )
        fractions[pathway] = round(fraction, 3)
        if fraction >= PROVIDES_THRESHOLD:
            provides.append(pathway)
        elif fraction <= REQUIRES_THRESHOLD and not definition.get("unreliable"):
            salvage = set(definition.get("salvage", set()))
            requires.append(f"{pathway}(salvage_only)" if salvage and names & salvage else pathway)
    return {"fractions": fractions, "provides": sorted(provides), "requires": sorted(requires)}


def scan(annotation_path: Path) -> tuple[dict[str, object], list[dict[str, str]]]:
    genome_dir = annotation_path.parent
    names: set[str] = set()
    hits: list[dict[str, str]] = []
    ldh_l: set[str] = set()
    ldh_d: set[str] = set()
    lactate_oxidases: set[str] = set()
    nadh_oxidase_confirmed: set[str] = set()
    descriptions_dha_kinase = 0
    descriptions_quinoprotein = 0
    header: list[str] | None = None
    index: dict[str, int] = {}
    proteins = 0
    marker_names = M.all_marker_names()

    with annotation_path.open(errors="replace") as handle:
        for line in handle:
            if line.startswith("#query"):
                header = line[1:].rstrip("\n").split("\t")
                index = {name: position for position, name in enumerate(header)}
                continue
            if line.startswith("#") or not line.strip() or header is None:
                continue
            fields = line.rstrip("\n").split("\t")
            fields.extend(["-"] * (len(header) - len(fields)))
            proteins += 1
            query = fields[index["query"]]
            preferred = fields[index["Preferred_name"]].strip()
            description = fields[index["Description"]].strip()
            lowered = description.lower()
            kegg = fields[index["KEGG_ko"]].strip()
            kos = split_values(kegg)
            categories: set[str] = set()

            if preferred not in {"", "-"}:
                names.add(preferred)
                if preferred in marker_names:
                    categories.add("marker")

            if M.L_LDH_KOS & kos and M.MALATE_DH_KO not in kos:
                ldh_l.add(query)
                categories.add("l_ldh")
            if M.D_LDH_KOS & kos:
                ldh_d.add(query)
                categories.add("d_ldh")
            if M.LACTATE_OXIDATION_KOS & kos:
                lactate_oxidases.add(query)
                categories.add("lactate_oxidation")
            if "nad-dependent l-lactate dehydrogenase" in lowered:
                ldh_l.add(query)
                categories.add("l_ldh")
            if "nad-dependent d-lactate dehydrogenase" in lowered:
                ldh_d.add(query)
                categories.add("d_ldh")
            if "dihydroxyacetone kinase" in lowered:
                descriptions_dha_kinase += 1
                categories.add("dha_kinase_description")
            if M.NADH_OXIDASE_KOS & kos:
                nadh_oxidase_confirmed.add(query)
                categories.add("nadh_oxidase_confirmed")
            if "quinoprotein" in lowered or "pqq" in lowered:
                descriptions_quinoprotein += 1
                categories.add("quinoprotein_description")
            if "lactate permease" in lowered or "lactate transport" in lowered:
                categories.add("lactate_transport_description")
                names.add("lctP")

            if categories:
                hits.append({
                    "genome": genome_dir.name, "protein": query, "gene_name": preferred,
                    "description": description, "kegg_ko": kegg,
                    "categories": ",".join(sorted(categories)),
                })

    # ---- glycerol entry
    glp_abc = len(names & {"glpA", "glpB", "glpC"})
    g3p_route = "glpK" in names and ("glpD" in names or glp_abc >= 2)
    facilitator = bool(names & M.GLYCEROL_FACILITATOR)
    cytoplasmic_dh = len(names & M.GLYCEROL_DH_CYTOPLASMIC)
    periplasmic_dh = len(names & M.GLYCEROL_DH_PERIPLASMIC)
    pqq = len(names & M.PQQ_BIOSYNTHESIS)

    # ---- DHA handling. Two kinase families; either is receiver capability.
    pts_subunits = len(names & M.DHA_KINASE_PTS)
    atp_kinase = bool(names & M.DHA_KINASE_ATP)
    # A named gene is stronger evidence than a description match. Only the
    # confident call penalises a donor: description-only hits are near-universal
    # (all six panel genomes have one), so scoring them would have made the donor
    # role unassignable to anything.
    dha_kinase_confident = pts_subunits >= 2 or atp_kinase
    dha_uptake = dha_kinase_confident or descriptions_dha_kinase > 0
    if periplasmic_dh >= 1 and pqq >= 1:
        # Periplasmic oxidation puts DHA outside the cell by construction.
        dha_export = "periplasmic_oxidation"
    elif cytoplasmic_dh >= 1 and facilitator:
        dha_export = "cytoplasmic_via_facilitator"
    elif cytoplasmic_dh >= 1:
        dha_export = "cytoplasmic_no_exporter"
    else:
        dha_export = "none"

    # ---- lactate
    lower_steps = len(names & M.LOWER_GLYCOLYSIS)
    productive_ldh = bool(ldh_l or ldh_d)
    lactate_transport = len(names & M.LACTATE_TRANSPORT)
    consumption = len(names & M.LACTATE_CONSUMPTION) + len(lactate_oxidases)
    # Which isomer, if any. Optical purity is the economic differentiator, and
    # L-lactate is the polymer-grade one, so a D-only strain is a different product.
    if ldh_l and ldh_d:
        isomer = "L_and_D"
    elif ldh_l:
        isomer = "L"
    elif ldh_d:
        isomer = "D"
    else:
        isomer = "none"

    # LctP/LldP-family permeases are catabolic lactate *uptake* carriers. Crediting
    # them as export while only weakly penalising consumption ranked three lactate
    # catabolisers (no L-LDH, 1-4 consumption genes) among the top producers. A
    # strain with consumption machinery and no NAD-dependent L-LDH is a consumer.
    if consumption >= 2 and not ldh_l:
        lactate_direction = "consumer_favoured"
    elif ldh_l and consumption <= 1:
        lactate_direction = "producer_favoured"
    else:
        lactate_direction = "ambiguous"

    if lactate_transport:
        lactate_export = "annotated"
    elif productive_ldh:
        # LDH with no annotated carrier: the strain may still export lactate, but
        # nothing in the genome supports it and a reconstruction will strand
        # lactate in the cytosol as a dead-end metabolite.
        lactate_export = "unannotated_gap"
    else:
        lactate_export = "no_ldh"

    # ---- aerobic pyruvate competition
    pdh = len(names & M.PYRUVATE_DEHYDROGENASE)
    pyruvate_oxidase = len(names & M.PYRUVATE_OXIDASE)
    # Name-only hits are recorded separately from corroborated ones; only corroborated
    # hits carry the aerobic-drain penalty.
    nadh_oxidase_named = len(names & M.NADH_OXIDASE)
    nadh_oxidase = len(nadh_oxidase_confirmed)
    competing = len(names & M.COMPETING_SINKS)
    pdo = len(names & M.PDO_SINK)
    methylglyoxal = len(names & M.METHYLGLYOXAL)

    # ---- oxygen handling
    catalase = len(names & M.CATALASE)
    peroxidase = len(names & M.PEROXIDASE_OTHER)
    sod = len(names & M.SUPEROXIDE_DISMUTASE)
    oxidase = len(names & M.TERMINAL_OXIDASE)

    amino = pathway_status(names, M.AMINO_ACID_PATHWAYS)
    vitamin = pathway_status(names, M.VITAMIN_PATHWAYS)

    # Cobalamin needs the conditional treatment described in markers.py: a missing
    # synthesis pathway only becomes a requirement when something in the cell
    # actually depends on B12 and has no B12-independent substitute.
    b12_dependencies = sorted(
        enzyme for enzyme, substitutes in M.B12_DEPENDENT_ENZYMES.items()
        if enzyme in names and not (names & substitutes)
    )
    b12_nonessential = sorted(names & M.B12_DEPENDENT_NONESSENTIAL)
    if not b12_dependencies:
        vitamin["requires"] = [item for item in vitamin["requires"] if not item.startswith("cobalamin_B12")]

    requires = [*amino["requires"], *vitamin["requires"]]
    provides = [*amino["provides"], *vitamin["provides"]]
    heme_dependent = any(item.startswith("heme") for item in vitamin["requires"])
    lipoate_salvage_only = any(item.startswith("lipoate") for item in vitamin["requires"])
    # H2O2 is generated by pyruvate oxidase and NADH oxidase; without catalase the
    # strain depends on a partner to remove it under aeration. Heme auxotrophy is
    # the usual underlying cause, since catalase is a haemoprotein.
    peroxide_exposed = int((pyruvate_oxidase or nadh_oxidase) and not catalase)

    genus, genus_pct, tax_flag = taxonomy(genome_dir)
    virulence = sorted(names & M.VIRULENCE_REVIEW_NAMES)

    # ---- scores
    glycerol_entry = g3p_route or (cytoplasmic_dh and dha_uptake) or periplasmic_dh >= 1
    # Aerobic penalties are what v1 lacked entirely: PDH and NADH oxidase are the
    # reason an aerated culture makes CO2 and biomass instead of lactate.
    aerobic_drain = 1.2 * min(pdh, 4) + 2.5 * min(pyruvate_oxidase, 2) + 3.0 * min(nadh_oxidase, 2)
    producer = (
        16 * int(bool(glycerol_entry))
        + 1.6 * lower_steps
        + 11 * int(productive_ldh)
        + 8 * int(lactate_export == "annotated")
        + 2 * min(len(ldh_l) + len(ldh_d), 2)
        - 5 * min(consumption, 4)
        - aerobic_drain
        - 2.5 * min(pdo, 2)
        - 0.35 * competing
    )
    dha_donor = (
        20 * int(dha_export == "periplasmic_oxidation")
        + 8 * int(dha_export == "cytoplasmic_via_facilitator")
        + min(5, pqq)
        + 4 * int(periplasmic_dh >= 2)
        - 3 * int(dha_kinase_confident)  # a donor that phosphorylates DHA competes with its own partner
        - 2 * min(pdo, 2)
    )
    dha_receiver = (
        18 * int(bool(dha_uptake))
        + 1.6 * lower_steps
        + 11 * int(productive_ldh)
        + 8 * int(lactate_export == "annotated")
        - 4 * min(consumption, 3)
        - aerobic_drain
        - 1.5 * min(pdo, 2)
    )
    supporter = (
        1.6 * len(provides)
        - 1.2 * len(requires)
        + 5 * int("cobalamin_B12" in provides)
        + 4 * int(catalase > 0)
        + 3 * int("heme" in provides)
        + 1.0 * len(names & M.AMINO_ACID_EXPORTERS)
        + 4 * int(bool(glycerol_entry))
    )

    roles = []
    if (producer >= 34 and productive_ldh and glycerol_entry
            and lower_steps >= MIN_LOWER_GLYCOLYSIS
            and lactate_direction != "consumer_favoured"):
        roles.append("direct_producer")
    if dha_donor >= 20:
        roles.append("dha_donor")
    if dha_receiver >= 34 and dha_uptake and productive_ldh and lower_steps >= MIN_LOWER_GLYCOLYSIS:
        roles.append("dha_receiver")
    if supporter >= 20:
        roles.append("supporter")
    if catalase and not peroxide_exposed:
        roles.append("peroxide_scavenger")

    row: dict[str, object] = {
        "genome": genome_dir.name,
        "genus": genus,
        "genus_percent": genus_pct,
        "taxonomy_flag": tax_flag,
        **assembly_quality(genome_dir),
        "glycerol_kinase_route": int(g3p_route),
        "glycerol_facilitator": int(facilitator),
        "glycerol_dh_cytoplasmic": cytoplasmic_dh,
        "glycerol_dh_periplasmic": periplasmic_dh,
        "pqq_genes": pqq,
        "quinoprotein_descriptions": descriptions_quinoprotein,
        "dha_kinase_pts_subunits": pts_subunits,
        "dha_kinase_atp": int(atp_kinase),
        "dha_kinase_confident": int(dha_kinase_confident),
        "dha_uptake_possible": int(bool(dha_uptake)),
        "dha_export_route": dha_export,
        "lower_glycolysis_steps": lower_steps,
        "l_ldh_count": len(ldh_l),
        "d_ldh_count": len(ldh_d),
        "lactate_isomer": isomer,
        "lactate_transport": lactate_transport,
        "lactate_export_call": lactate_export,
        "lactate_consumption": consumption,
        "lactate_direction": lactate_direction,
        "pyruvate_dehydrogenase": pdh,
        "pyruvate_oxidase": pyruvate_oxidase,
        "nadh_oxidase": nadh_oxidase,
        "nadh_oxidase_name_only": nadh_oxidase_named,
        "competing_sinks": competing,
        "propanediol_sinks": pdo,
        "methylglyoxal_genes": methylglyoxal,
        "catalase": catalase,
        "other_peroxidase": peroxidase,
        "superoxide_dismutase": sod,
        "terminal_oxidase_genes": oxidase,
        "peroxide_exposed_no_catalase": peroxide_exposed,
        "requires_count": len(requires),
        "requires": ";".join(requires) or "none",
        "provides_count": len(provides),
        "provides": ";".join(provides) or "none",
        "heme_dependent": int(heme_dependent),
        "lipoate_salvage_only": int(lipoate_salvage_only),
        "b12_required_because": ";".join(b12_dependencies) or "not_required",
        "b12_dependent_nonessential": ";".join(b12_nonessential) or "none",
        "amino_acid_exporters": len(names & M.AMINO_ACID_EXPORTERS),
        "nitrogen_assimilation": len(names & M.NITROGEN_ASSIMILATION),
        "risk_genus": int(genus in M.RISK_GENERA),
        "biosafety_review_required": int(bool(virulence) or genus in M.RISK_GENERA),
        "virulence_name_matches": ";".join(virulence) or "none",
        "producer_score": round(producer, 2),
        "dha_donor_score": round(dha_donor, 2),
        "dha_receiver_score": round(dha_receiver, 2),
        "supporter_score": round(supporter, 2),
        "predicted_roles": ",".join(roles) or "none",
        "annotated_proteins": proteins,
        "annotation_path": str(annotation_path),
        "nucleotide_fasta": next((str(p) for p in sorted(genome_dir.glob("prokka/*.fna"))), ""),
        "protein_fasta": next((str(p) for p in sorted(genome_dir.glob("prokka/*.faa"))), ""),
    }
    for pathway, fraction in {**amino["fractions"], **vitamin["fractions"]}.items():
        row[f"pathway_{pathway}"] = fraction
    return row, hits


def summary(rows: list[dict[str, object]], path: Path) -> None:
    sections = [
        ("Direct glycerol-to-lactate producers (aerobic)", "producer_score"),
        ("Glycerol-to-DHA donors", "dha_donor_score"),
        ("DHA-to-lactate receivers", "dha_receiver_score"),
        ("Nutrient supporters", "supporter_score"),
    ]
    lines = [
        "# Annotation-based candidate ranking (screen v2)",
        "",
        "Priors for choosing models and experiments. They do not establish growth,",
        "secretion, flux, yield or safety. Requirements are annotation-based and must",
        "be confirmed against the GEMs; absence of a named gene is weak evidence.",
        "",
    ]
    for title, metric in sections:
        lines += [
            f"## {title}", "",
            "| Rank | Genome | Genus | % | Score | Isomer | Lactate direction | Export | Drain (PDH/POX/NOX) | Assembly | Requires | Roles |",
            "|---:|---|---|---:|---:|---|---|---|---|---|---|---|",
        ]
        ranked = sorted(rows, key=lambda row: (-float(row[metric]), str(row["genome"])))[:25]
        for rank, row in enumerate(ranked, 1):
            drain = f"{row['pyruvate_dehydrogenase']}/{row['pyruvate_oxidase']}/{row['nadh_oxidase']}"
            requires = str(row["requires"])
            if len(requires) > 55:
                requires = requires[:52] + "..."
            lines.append(
                f"| {rank} | `{row['genome']}` | {row['genus'] or 'unresolved'} | "
                f"{row['genus_percent'] or '-'} | {row[metric]} | {row['lactate_isomer']} | "
                f"{row['lactate_direction']} | {row['lactate_export_call']} | "
                f"{drain} | {row['assembly_quality']} | {requires} | {row['predicted_roles']} |"
            )
        lines.append("")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    options = parse_args()
    annotations = sorted(options.genome_root.glob("*/eggnog.emapper.annotations"))
    if not annotations:
        raise SystemExit(f"No eggNOG annotation files found under {options.genome_root}")
    options.output_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    hits: list[dict[str, str]] = []
    for position, annotation in enumerate(annotations, 1):
        row, marker_hits = scan(annotation)
        rows.append(row)
        hits.extend(marker_hits)
        if position % 25 == 0 or position == len(annotations):
            print(f"Scanned {position}/{len(annotations)} genomes", flush=True)

    fields = list(rows[0].keys())
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with (options.output_dir / "genome_scores.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", restval="")
        writer.writeheader()
        writer.writerows(rows)
    with (options.output_dir / "marker_hits.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["genome", "protein", "gene_name", "description", "kegg_ko", "categories"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(hits)
    summary(rows, options.output_dir / "RANKED_CANDIDATES.md")
    print(f"wrote {options.output_dir}/genome_scores.tsv, marker_hits.tsv, RANKED_CANDIDATES.md")


if __name__ == "__main__":
    main()
