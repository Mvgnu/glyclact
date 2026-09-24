#!/usr/bin/env python3
"""Summarise crude-glycerol detoxification and competition markers.

Why this exists rather than another column on the old screen.

The existing supporter score is provisioning: completed biosynthetic pathways, amino-acid
exporters, catalase. Two results make that the wrong target.

1. **The provisioning never happens in silico.** Every supporter carries all 20 amino-acid
   exchanges and could secrete several mmol/gDW/h, but at maximum growth secretes exactly
   zero -- export is not growth-coupled, so FBA sits at the minimum. Same degeneracy as
   lactate.
2. **Supporters competitively exclude the producer.** In the community landscape, producer
   biomass falls from 0.257 alone to 0.0033 with four supporters, an ~80x suppression,
   because the Enterobacteriaceae grow at mu ~0.55-0.60 on glycerol against 0.117-0.129
   for the Bacillus producers. On one shared substrate the fast grower simply takes it.

The literature says the same thing from the other side: every documented case where a
consortium helps lactic acid production is a *heterogeneity* case -- mixed sugars, whey,
molasses, raw glycerol -- and for crude glycerol specifically the proposal is to "partition
the detoxification and production tasks". Crude glycerol carries residual methanol, soap,
free long-chain fatty acids, 50-200 g/L Na/K salts and heavy metals.

So the criterion inverts. A useful partner here should:

* **grow on the impurities** -- methanol and LCFA are carbon sources, so a detoxifier can
  fund itself from what is poisoning the producer;
* **not grow on glycerol** -- which was previously a rejection reason and is now a feature,
  because it removes the competition that suppressed the producer 80-fold;
* tolerate salt, and not consume lactate.

KO-first throughout. This project has repeatedly been bitten by gene-name collisions
(`cat` matching cation efflux, `mecA` matching a ClpC adaptor, `noxE` matching an Old
Yellow Enzyme), so names are only ever a fallback where KEGG is genuinely ambiguous.

Writes results/crude_glycerol/detox_screen.tsv.

This is an annotation-evidence report. It does not reinstate a supporter or select panel
members; the current panel deliberately has no non-producing partner.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOCAL_GENOME_ROOT = ROOT / "candidate_genomes"

# Each entry: label -> (KO set, fallback gene names)
TRAITS: dict[str, dict[str, tuple[set[str], set[str]]]] = {
    # Methanol is a carbon source for the partner and a poison for the producer. The
    # dehydrogenase alone is not enough -- it makes formaldehyde, which is worse -- so
    # detoxification capacity is scored separately below.
    "methanol_oxidation": {
        "mdh_pqq": ({"K14028", "K14029", "K23995"}, {"mxaF", "mxaI", "xoxF"}),
        "mdh_nad": ({"K17066"}, {"mdo", "mdh2"}),
    },
    # The step that actually determines survival on methanol-containing feedstock.
    "formaldehyde_detox": {
        "gsh_dependent": ({"K00121", "K01070"}, {"frmA", "frmB", "adhC", "fghA"}),
        "rump_pathway": ({"K08093", "K08094"}, {"hxlA", "hxlB"}),
        "formaldehyde_dh": ({"K00148"}, set()),
        "formate_dh": ({"K00122", "K00123", "K00126"}, {"fdhA", "fdhF"}),
    },
    # Soap and free long-chain fatty acids: both a toxicity and a carbon source.
    "lcfa_utilisation": {
        "transport": ({"K06076"}, {"fadL"}),
        "activation": ({"K01897"}, {"fadD"}),
        "beta_oxidation": ({"K06445", "K00632", "K01782", "K01825"}, {"fadE", "fadA", "fadB"}),
        "lipase": ({"K01046", "K07001"}, set()),
    },
    # 50-200 g/L Na/K is the largest single impurity load in crude glycerol.
    "halotolerance": {
        "betaine": ({"K00108", "K00130", "K02168"}, {"betA", "betB", "betT"}),
        "betaine_transport": ({"K02000", "K02001", "K02002"}, {"proV", "proW", "proX"}),
        "ectoine": ({"K06718", "K00836", "K06720"}, {"ectA", "ectB", "ectC"}),
        "trehalose": ({"K00697", "K01087"}, {"otsA", "otsB"}),
        "na_antiport": ({"K03313", "K03314", "K05565"}, {"nhaA", "nhaB", "mrpA"}),
    },
    "metal_resistance": {
        "efflux": ({"K15726", "K17686", "K01534"}, {"czcA", "copA", "zntA"}),
        "arsenic_mercury": ({"K03325", "K00520"}, {"arsB", "merA"}),
    },
    # Now scored as a LIABILITY for a detoxifier: competing for glycerol is what suppressed
    # the producer 80-fold. Reported, not summed into the score.
    "glycerol_competition": {
        "uptake": ({"K00864", "K02440"}, {"glpK", "glpF"}),
        "dissimilation": ({"K00111", "K00005"}, {"glpD", "gldA"}),
    },
    # Still disqualifying: a partner that eats the product.
    "lactate_consumption": {
        "oxidase": ({"K00101", "K00102"}, {"lldD", "dld"}),
    },
    "lactate_purity_risk": {
        "d_ldh": ({"K03778"}, set()),
    },
}


def count_traits(annotation_path: Path) -> dict[str, int]:
    kos: dict[str, int] = {}
    names: dict[str, int] = {}
    with annotation_path.open(errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            name = fields[8].strip()
            for ko in fields[11].replace("ko:", "").split(","):
                ko = ko.strip()
                if ko and ko != "-":
                    kos[ko] = kos.get(ko, 0) + 1
            if name and name != "-":
                names[name] = names.get(name, 0) + 1

    counts: dict[str, int] = {}
    for trait, components in TRAITS.items():
        for component, (ko_set, name_set) in components.items():
            hits = sum(kos.get(k, 0) for k in ko_set)
            if hits == 0:
                hits = sum(names.get(n, 0) for n in name_set)
            counts[f"{trait}__{component}"] = hits
        counts[trait] = sum(counts[f"{trait}__{c}"] for c in components)
    return counts


def summarise(counts: dict[str, int]) -> dict[str, object]:
    """Evidence categories rather than one additive score.

    The first version of this summed weighted trait counts, and the top four hits were the
    three producers -- already in the pool, and scoring well only on generic stress genes.
    Three things were wrong with it, and they are the reason this returns a classification
    instead of a number:

    * **Additive scoring let breadth beat fitness for the role.** A strain with many
      moderate counts outranked one with the single capability that matters.
    * **Halotolerance is near-universal**, so scoring it inflated everything equally and
      discriminated nothing. It belongs as a threshold.
    * **Tolerating an impurity is not removing it.** Only utilisation helps a partner, so
      that is the role evidence; tolerance only indicates whether it may survive to do it.

    Glycerol competition is a *tier*, not a penalty, because the landscape result is closer
    to disqualifying than to costly: a glycerol-capable partner suppressed producer biomass
    ~80-fold. A detoxifier that cannot touch the feedstock is the point, not a compromise.
    """
    glycerol = counts.get("glycerol_competition__uptake", 0) + counts.get(
        "glycerol_competition__dissimilation", 0
    )
    if glycerol == 0:
        competition = "none"            # ideal: eats impurities, leaves the substrate
    elif counts.get("glycerol_competition__uptake", 0) == 0:
        competition = "partial"         # dissimilation genes but no uptake route
    else:
        competition = "competitor"

    methanol = counts.get("methanol_oxidation", 0)
    lcfa = counts.get("lcfa_utilisation__activation", 0) + counts.get(
        "lcfa_utilisation__beta_oxidation", 0
    )
    # Methanol oxidation without formaldehyde detox turns one toxin into a worse one, so
    # the methylotroph role requires both.
    methylotroph = methanol > 0 and counts.get("formaldehyde_detox", 0) > 0
    lipolytic = lcfa >= 2 and counts.get("lcfa_utilisation__lipase", 0) > 0

    roles = []
    if methylotroph:
        roles.append("methanol")
    if lipolytic:
        roles.append("lcfa_soap")
    return {
        "competition_tier": competition,
        "detox_roles": ",".join(roles) or "none",
        "halotolerance_ok": int(counts.get("halotolerance", 0) >= 4),
        "methanol_oxidation": methanol,
        "formaldehyde_detox": counts.get("formaldehyde_detox", 0),
        "lcfa_utilisation": counts.get("lcfa_utilisation", 0),
        "halotolerance": counts.get("halotolerance", 0),
        "metal_resistance": counts.get("metal_resistance", 0),
        "glycerol_genes": glycerol,
    }


def main() -> None:
    root = ROOT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome-root", type=Path, default=LOCAL_GENOME_ROOT)
    parser.add_argument("--manifest", type=Path,
                        default=root / "candidate_genomes" / "manifest_pool_v5.tsv")
    parser.add_argument("--genomes", type=Path,
                        help="optional file with one genome id per line; overrides --manifest")
    parser.add_argument("--output", type=Path, default=root / "results" / "crude_glycerol" / "detox_screen.tsv")
    parser.add_argument("--top", type=int, default=20)
    options = parser.parse_args()

    genome_root = options.genome_root
    if not genome_root.is_dir():
        genome_root = LOCAL_GENOME_ROOT
        print(f"drive unavailable; scanning the local panel copies in {genome_root}")

    if options.genomes:
        genome_ids = [line.strip() for line in options.genomes.read_text().splitlines()
                      if line.strip()]
    else:
        with options.manifest.open(encoding="utf-8") as handle:
            genome_ids = [row["genome"] for row in csv.DictReader(handle, delimiter="\t")]

    rows: list[dict[str, object]] = []
    for genome in genome_ids:
        directory = genome_root / genome
        annotation = directory / "eggnog.emapper.annotations"
        if not annotation.exists():
            continue
        counts = count_traits(annotation)
        rows.append({
            "genome": directory.name,
            **summarise(counts),
            "eats_lactate": int(counts.get("lactate_consumption", 0) > 0),
            "d_ldh": counts.get("lactate_purity_risk", 0),
        })

    if not rows:
        raise SystemExit(f"no annotations found under {genome_root}")
    options.output.parent.mkdir(parents=True, exist_ok=True)
    tier_order = {"none": 0, "partial": 1, "competitor": 2}
    rows.sort(key=lambda r: (r["detox_roles"] == "none",
                             tier_order[str(r["competition_tier"])],
                             -int(r["lcfa_utilisation"])))
    with options.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"scanned {len(rows)} genomes -> {options.output}")

    print(f"\n{'':<24}{'no glycerol':>12}{'partial':>9}{'competitor':>12}")
    for role in ("methanol", "lcfa_soap", "none"):
        line = f"  {role:<22}"
        for tier in ("none", "partial", "competitor"):
            n = len([r for r in rows if role in str(r["detox_roles"]).split(",")
                     and r["competition_tier"] == tier])
            line += f"{n:>12}" if tier != "partial" else f"{n:>9}"
        print(line)

    clean = [r for r in rows if not r["eats_lactate"] and r["d_ldh"] == 0
             and r["detox_roles"] != "none" and r["halotolerance_ok"]]
    print(f"\n{len(clean)} genomes match the displayed marker profile: a detox role, "
          f"halotolerant, no annotated lactate oxidase, no annotated D-LDH")
    print(f"{'genome':<10}{'roles':<14}{'competition':<13}{'MeOH':>5}{'HCHO':>6}"
          f"{'LCFA':>6}{'salt':>6}{'metal':>6}")
    for row in clean[: options.top]:
        print(f"{row['genome']:<10}{row['detox_roles']:<14}{row['competition_tier']:<13}"
              f"{row['methanol_oxidation']:>5}{row['formaldehyde_detox']:>6}"
              f"{row['lcfa_utilisation']:>6}{row['halotolerance']:>6}"
              f"{row['metal_resistance']:>6}")


if __name__ == "__main__":
    main()
