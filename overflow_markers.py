#!/usr/bin/env python3
"""Scan the 9 candidates for the genes that decide *where* overflow carbon goes.

Driven by the literature pass: whether a strain can spill pyruvate to lactate
depends on (a) where glycerol's electrons land (quinone via GlpD vs NADH via GldA),
(b) whether ldh is present and how it is regulated (Rex), and (c) which competing
overflow sink exists (acetoin via AlsSD, acetate via Pta-AckA).
"""
import csv
from collections import defaultdict
from pathlib import Path

ROOT = Path("/Users/magnusohle/cursorprojects/ratzke/glycerol_lactate_community")

# KO-first, names only as a fallback where KEGG is ambiguous.
PROBES = [
    # glycerol entry: which route, and therefore where the electrons go
    ("glpF  glycerol facilitator",        {"K02440"}, {"glpF"}),
    ("glpK  glycerol kinase",             {"K00864"}, {"glpK"}),
    ("glpD  G3P-DH (quinone, aerobic)",   {"K00111"}, {"glpD", "glpA"}),
    ("glpB/C anaerobic G3P-DH",           {"K00112", "K00113"}, {"glpB", "glpC"}),
    ("gpsA  G3P-DH (NAD(P), biosynth)",   {"K00057"}, {"gpsA"}),
    ("gldA  glycerol DH (NAD)",           {"K00005"}, {"gldA", "dhaD"}),
    ("dhaK  DHA kinase",                  {"K00863", "K05878", "K05879", "K05881"}, {"dhaK"}),
    # the lactate node
    ("ldh   L-LDH (fermentative)",        {"K00016"}, {"ldh", "lctE", "ldhL"}),
    ("ldhA  D-LDH (fermentative)",        {"K03778"}, {"ldhA"}),
    ("lldD  L-lactate ox (respiratory)",  {"K00101"}, {"lldD", "lctO"}),
    ("dld   D-lactate ox (respiratory)",  {"K00102"}, {"dld"}),
    ("lctP  lactate permease",            {"K03303"}, {"lctP", "lldP"}),
    # the regulator that gates ldh on redox state
    ("rex   NADH/NAD+ repressor",         {"K03710"}, {"rex", "ydiH"}),
    # terminal oxidases: which one is used sets the P/O and the O2 affinity
    ("cydAB cytochrome bd (high aff.)",   {"K00425", "K00426"}, {"cydA", "cydB"}),
    ("cyo/qox low-affinity oxidase",      {"K02297", "K02298", "K02299", "K02300"}, {"qoxA", "qoxB", "cyoA", "cyoB"}),
    # competing pyruvate sinks
    ("pdh   pyruvate dehydrogenase",      {"K00161", "K00162", "K00627", "K00382"}, set()),
    ("alsS  acetolactate synthase",       {"K01652"}, {"alsS", "budB"}),
    ("alsD  acetolactate decarboxylase",  {"K01575"}, {"alsD", "budA"}),
    ("budC  2,3-BD dehydrogenase",        {"K00004"}, {"budC", "bdhA"}),
    ("pta   phosphotransacetylase",       {"K00625", "K13788"}, {"pta"}),
    ("ackA  acetate kinase",              {"K00925"}, {"ackA"}),
    ("poxB/spxB pyruvate oxidase",        {"K00156", "K00158"}, {"poxB", "spxB", "poxL"}),
    ("pflB  pyruvate formate-lyase",      {"K00656"}, {"pflB"}),
    ("adhE  acetaldehyde/alcohol DH",     {"K04072"}, {"adhE"}),
    # FBP node: the allosteric activator of Firmicute L-LDH and pyruvate kinase
    ("pfkA  phosphofructokinase",         {"K00850", "K21071"}, {"pfkA"}),
    ("fbp   fructose-1,6-bisphosphatase", {"K03841", "K02446", "K11532"}, {"fbp", "glpX"}),
    ("fbaA  FBP aldolase",                {"K01624", "K01623"}, {"fba", "fbaA"}),
    ("pyk   pyruvate kinase",             {"K00873"}, {"pyk", "pykA", "pykF"}),
    ("ccpA  catabolite control A",        set(), {"ccpA"}),
    ("nox   NADH oxidase (H2O)",          {"K00359", "K17869"}, set()),
]

def pool_genomes() -> list[str]:
    """Read the pool from the manifest rather than hard-coding it, so the scan tracks
    the panel as it grows."""
    manifest = ROOT / "candidate_genomes" / "manifest.tsv"
    with manifest.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
        index = header.index("genome")
        return [line.rstrip("\n").split("\t")[index] for line in handle if line.strip()]


GENOMES = pool_genomes()

counts = defaultdict(dict)
for genome in GENOMES:
    path = ROOT / "candidate_genomes" / genome / "eggnog.emapper.annotations"
    kos, names = defaultdict(int), defaultdict(int)
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            name, ko_field = f[8].strip(), f[11].strip()
            for ko in ko_field.replace("ko:", "").split(","):
                if ko and ko != "-":
                    kos[ko] += 1
            if name and name != "-":
                names[name] += 1
    for label, ko_set, name_set in PROBES:
        hits = sum(kos.get(k, 0) for k in ko_set)
        if hits == 0:
            hits = sum(names.get(n, 0) for n in name_set)
        counts[label][genome] = hits

width = max(len(label) for label, _, _ in PROBES)
print(" " * width + "  " + " ".join(f"{g:>7}" for g in GENOMES))
for label, _, _ in PROBES:
    row = counts[label]
    cells = " ".join(f"{row[g]:>7}" if row[g] else f"{'.':>7}" for g in GENOMES)
    print(f"{label:<{width}}  {cells}")

out = ROOT / "results" / "overflow_gene_scan.tsv"
with out.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["gene_probe"] + GENOMES)
    for label, _, _ in PROBES:
        writer.writerow([label] + [counts[label][g] for g in GENOMES])
print(f"\nwrote {out}")
