#!/usr/bin/env python3
"""Scan a declared genome list for overflow and redox marker genes.

This is annotation evidence only. It reports which routes are represented; it does not
infer flux or select pool members. The current manifest is the default; another manifest
or explicit genome list must be supplied to inspect a historical cohort.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# KO-first, names only as a fallback where KEGG is ambiguous.
PROBES = [
    ("glpF  glycerol facilitator",        {"K02440"}, {"glpF"}),
    ("glpK  glycerol kinase",             {"K00864"}, {"glpK"}),
    ("glpD  G3P-DH (quinone, aerobic)",   {"K00111"}, {"glpD", "glpA"}),
    ("glpB/C anaerobic G3P-DH",           {"K00112", "K00113"}, {"glpB", "glpC"}),
    ("gpsA  G3P-DH (NAD(P), biosynth)",   {"K00057"}, {"gpsA"}),
    ("gldA  glycerol DH (NAD)",           {"K00005"}, {"gldA", "dhaD"}),
    ("dhaK  DHA kinase",                  {"K00863", "K05878", "K05879", "K05881"}, {"dhaK"}),
    ("ldh   L-LDH (fermentative)",        {"K00016"}, {"ldh", "lctE", "ldhL"}),
    ("ldhA  D-LDH (fermentative)",        {"K03778"}, {"ldhA"}),
    ("lldD  L-lactate ox (respiratory)",  {"K00101"}, {"lldD", "lctO"}),
    ("dld   D-lactate ox (respiratory)",  {"K00102"}, {"dld"}),
    ("lctP  lactate permease",            {"K03303"}, {"lctP", "lldP"}),
    ("rex   NADH/NAD+ repressor",         {"K03710"}, {"rex", "ydiH"}),
    ("cydAB cytochrome bd (high aff.)",   {"K00425", "K00426"}, {"cydA", "cydB"}),
    ("cyo/qox low-affinity oxidase",      {"K02297", "K02298", "K02299", "K02300"},
                                             {"qoxA", "qoxB", "cyoA", "cyoB"}),
    ("pdh   pyruvate dehydrogenase",      {"K00161", "K00162", "K00627", "K00382"}, set()),
    ("alsS  acetolactate synthase",       {"K01652"}, {"alsS", "budB"}),
    ("alsD  acetolactate decarboxylase",  {"K01575"}, {"alsD", "budA"}),
    ("budC  2,3-BD dehydrogenase",        {"K00004"}, {"budC", "bdhA"}),
    ("pta   phosphotransacetylase",       {"K00625", "K13788"}, {"pta"}),
    ("ackA  acetate kinase",              {"K00925"}, {"ackA"}),
    ("poxB/spxB pyruvate oxidase",        {"K00156", "K00158"}, {"poxB", "spxB", "poxL"}),
    ("pflB  pyruvate formate-lyase",      {"K00656"}, {"pflB"}),
    ("adhE  acetaldehyde/alcohol DH",     {"K04072"}, {"adhE"}),
    ("pfkA  phosphofructokinase",         {"K00850", "K21071"}, {"pfkA"}),
    ("fbp   fructose-1,6-bisphosphatase", {"K03841", "K02446", "K11532"}, {"fbp", "glpX"}),
    ("fbaA  FBP aldolase",                {"K01624", "K01623"}, {"fba", "fbaA"}),
    ("pyk   pyruvate kinase",             {"K00873"}, {"pyk", "pykA", "pykF"}),
    ("ccpA  catabolite control A",        set(), {"ccpA"}),
    ("nox   NADH oxidase (H2O)",          {"K00359", "K17869"}, set()),
]


def scan_annotation(path: Path) -> tuple[dict[str, int], dict[str, int]]:
    kos: dict[str, int] = defaultdict(int)
    names: dict[str, int] = defaultdict(int)
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            name, ko_field = fields[8].strip(), fields[11].strip()
            for ko in ko_field.replace("ko:", "").split(","):
                if ko and ko != "-":
                    kos[ko] += 1
            if name and name != "-":
                names[name] += 1
    return kos, names


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path,
                        default=ROOT / "candidate_genomes" / "manifest_pool_v5.tsv")
    parser.add_argument("--genomes", type=Path,
                        help="optional file with one genome id per line; overrides --manifest")
    parser.add_argument("--genomes-dir", type=Path, default=ROOT / "candidate_genomes")
    parser.add_argument("--output", type=Path,
                        default=ROOT / "results" / "pool_v5" / "overflow_gene_scan.tsv")
    args = parser.parse_args()

    if args.genomes:
        genomes = [line.strip() for line in args.genomes.read_text().splitlines() if line.strip()]
    else:
        with args.manifest.open(encoding="utf-8") as handle:
            genomes = [row["genome"] for row in csv.DictReader(handle, delimiter="\t")]
    counts: dict[str, dict[str, int]] = defaultdict(dict)
    for genome in genomes:
        annotation = args.genomes_dir / genome / "eggnog.emapper.annotations"
        if not annotation.exists():
            raise FileNotFoundError(f"missing annotation for {genome}: {annotation}")
        kos, names = scan_annotation(annotation)
        for label, ko_set, name_set in PROBES:
            hits = sum(kos.get(ko, 0) for ko in ko_set)
            counts[label][genome] = hits or sum(names.get(name, 0) for name in name_set)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["gene_probe", *genomes])
        for label, _, _ in PROBES:
            writer.writerow([label, *(counts[label][genome] for genome in genomes)])
    print(f"wrote {args.output} ({len(genomes)} genomes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
